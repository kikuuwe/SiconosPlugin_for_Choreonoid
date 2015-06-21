/* 
 This file is a part of SiconosPlugin.
 
 Author: Shin'ichiro Nakaoka
 Author: Ryo Kikuuwe
 
 Copyright (c) 2007-2015 Shin'ichiro Nakaoka
 Copyright (c) 2014-2015 Ryo Kikuuwe
 Copyright (c) 2007-2015 National Institute of Advanced Industrial
                         Science and Technology (AIST)
 Copyright (c) 2014-2015 Kyushu University

 SiconosPlugin is a plugin for Choreonoid to provide access to Siconos-Numerics.
 
 SiconosPlugin is a free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 SiconosPlugin is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with SiconosPlugin; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 
 Contact: Ryo Kikuuwe, kikuuwe@ieee.org
*/



#include "SPSimulatorItem.h"
#include "SPConstraintForceSolver.h"
#include <cnoid/ForwardDynamicsCBM>
#include <cnoid/DyWorld>

#include <cnoid/ItemManager>
#include <cnoid/Archive>
#include <cnoid/EigenArchive>
#include <cnoid/DyBody>
#include <cnoid/LeggedBodyHelper>
#include <cnoid/ControllerItem>
#include <cnoid/BodyMotionItem>

#include <cnoid/FloatingNumberString>
#include <cnoid/EigenUtil>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include "gettext.h"

using namespace std;
using namespace cnoid;

// for Windows
#undef min
#undef max

namespace {

const bool TRACE_FUNCTIONS = false;
const bool ENABLE_DEBUG_OUTPUT = false;
const double DEFAULT_GRAVITY_ACCELERATION = 9.80665;

class HighGainControllerItem : public ControllerItem
{
    BodyPtr body;
    MultiValueSeqPtr qseqRef;
    int currentFrame;
    int lastFrame;
    int numJoints;

public:
    HighGainControllerItem(BodyItem* bodyItem, BodyMotionItem* bodyMotionItem) {
        qseqRef = bodyMotionItem->jointPosSeq();
        setName(str(fmt(_("HighGain Controller with %1%")) % bodyMotionItem->name()));
    }

    virtual bool start(Target* target) {
        body = target->body();
        currentFrame = 0;
        lastFrame = std::max(0, qseqRef->numFrames() - 1);
        numJoints = std::min(body->numJoints(), qseqRef->numParts());
        if(qseqRef->numFrames() == 0){
            putMessage(_("Reference motion is empty()."));
            return false;
        }
        if(fabs(qseqRef->frameRate() - (1.0 / target->worldTimeStep())) > 1.0e-6){
            putMessage(_("The frame rate of the reference motion is different from the world frame rate."));
            return false;
        }
        control();
        return true;
    }

    virtual double timeStep() const {
        return qseqRef->getTimeStep();
    }
        
    virtual void input() { }

    virtual bool control() {

        if(++currentFrame > lastFrame){
            currentFrame = lastFrame;
            return false;
        }
        return true;
    }
        
    virtual void output() {

        int prevFrame = std::max(currentFrame - 1, 0);
        int nextFrame = std::min(currentFrame + 1, lastFrame);
            
        MultiValueSeq::Frame q0 = qseqRef->frame(prevFrame);
        MultiValueSeq::Frame q1 = qseqRef->frame(currentFrame);
        MultiValueSeq::Frame q2 = qseqRef->frame(nextFrame);

        double dt = qseqRef->getTimeStep();
        double dt2 = dt * dt;

        for(int i=0; i < numJoints; ++i){
            Link* joint = body->joint(i);
            joint->q() = q1[i];
            joint->dq() = (q2[i] - q1[i]) / dt;
            joint->ddq() = (q2[i] - 2.0 * q1[i] + q0[i]) / dt2;
        }
    }
        
    virtual void stop() { }

};
}


namespace cnoid {
  
class SPSimulatorItemImpl
{
public:
    SPSimulatorItem* self;

    World<SPConstraintForceSolver> world;
        
    Selection dynamicsMode;
    Selection integrationMode;
    Vector3 gravity;
    double staticFriction;
    double slipFriction;
    FloatingNumberString contactCullingDistance;
    FloatingNumberString contactCullingDepth;
    FloatingNumberString errorCriterion;
    int maxNumIterations;
    FloatingNumberString contactCorrectionDepth;
    FloatingNumberString contactCorrectionVelocityRatio;
    double epsilon;
    bool is2Dmode;
    bool isKinematicWalkingEnabled;

    typedef std::map<Body*, int> BodyIndexMap;
    BodyIndexMap bodyIndexMap;

    SPSimulatorItemImpl(SPSimulatorItem* self);
    SPSimulatorItemImpl(SPSimulatorItem* self, const SPSimulatorItemImpl& org);
    bool initializeSimulation(const std::vector<SimulationBody*>& simBodies);
    void addBody(SimulationBody* simBody);
    void doPutProperties(PutPropertyFunction& putProperty);
    bool store(Archive& archive);
    bool restore(const Archive& archive);

    // for debug
    ofstream os;
};

}


void SPSimulatorItem::initializeClass(ExtensionManager* ext)
{
    ext->itemManager().registerClass<SPSimulatorItem>(N_("SiconosSimulatorItem"));
    ext->itemManager().addCreationPanel<SPSimulatorItem>();
}


SPSimulatorItem::SPSimulatorItem()
{
    impl = new SPSimulatorItemImpl(this);
}
 

SPSimulatorItemImpl::SPSimulatorItemImpl(SPSimulatorItem* self)
    : self(self),
      dynamicsMode(SPSimulatorItem::N_DYNAMICS_MODES, CNOID_GETTEXT_DOMAIN_NAME),
      integrationMode(SPSimulatorItem::N_INTEGRATION_MODES, CNOID_GETTEXT_DOMAIN_NAME)
{
    dynamicsMode.setSymbol(SPSimulatorItem::FORWARD_DYNAMICS,  N_("Forward dynamics"));
    dynamicsMode.setSymbol(SPSimulatorItem::HG_DYNAMICS,       N_("High-gain dynamics"));
    dynamicsMode.setSymbol(SPSimulatorItem::KINEMATICS,        N_("Kinematics"));

    integrationMode.setSymbol(SPSimulatorItem::EULER_INTEGRATION,  N_("Euler"));
//    integrationMode.setSymbol(SPSimulatorItem::RUNGE_KUTTA_INTEGRATION,  N_("Runge Kutta"));
    integrationMode.select(SPSimulatorItem::EULER_INTEGRATION);
    
    gravity << 0.0, 0.0, -DEFAULT_GRAVITY_ACCELERATION;

    SPConstraintForceSolver& cfs = world.constraintForceSolver;
    staticFriction = cfs.staticFriction();
    slipFriction = cfs.slipFriction();
    contactCullingDistance = cfs.contactCullingDistance();
    contactCullingDepth = cfs.contactCullingDepth();
    epsilon = cfs.coefficientOfRestitution();
    
    errorCriterion = cfs.gaussSeidelErrorCriterion();
    maxNumIterations = cfs.gaussSeidelMaxNumIterations();
    contactCorrectionDepth = cfs.contactCorrectionDepth();
    contactCorrectionVelocityRatio = cfs.contactCorrectionVelocityRatio();

    isKinematicWalkingEnabled = false;
    is2Dmode = false;
}


SPSimulatorItem::SPSimulatorItem(const SPSimulatorItem& org)
    : SimulatorItem(org),
      impl(new SPSimulatorItemImpl(this, *org.impl))
{

}


SPSimulatorItemImpl::SPSimulatorItemImpl(SPSimulatorItem* self, const SPSimulatorItemImpl& org)
    : self(self),
      dynamicsMode(org.dynamicsMode),
      integrationMode(org.integrationMode)
{
    gravity = org.gravity;
    staticFriction = org.staticFriction;
    slipFriction = org.slipFriction;
    contactCullingDistance = org.contactCullingDistance;
    contactCullingDepth = org.contactCullingDepth;
    errorCriterion = org.errorCriterion;
    maxNumIterations = org.maxNumIterations;
    contactCorrectionDepth = org.contactCorrectionDepth;
    contactCorrectionVelocityRatio = org.contactCorrectionVelocityRatio;
    epsilon = org.epsilon;
    isKinematicWalkingEnabled = org.isKinematicWalkingEnabled;
    is2Dmode = org.is2Dmode;
}


SPSimulatorItem::~SPSimulatorItem()
{
    delete impl;
}


void SPSimulatorItem::setDynamicsMode(int mode)
{
    impl->dynamicsMode.select(mode);
}


void SPSimulatorItem::setIntegrationMode(int mode)
{
    impl->integrationMode.select(mode);
}


void SPSimulatorItem::setGravity(const Vector3& gravity)
{
    impl->gravity = gravity;
}


void SPSimulatorItem::setStaticFriction(double value)
{
    impl->staticFriction = value; 
}


void SPSimulatorItem::setSlipFriction(double value)
{
    impl->slipFriction = value;
}


void SPSimulatorItem::setContactCullingDistance(double value)    
{
    impl->contactCullingDistance = value;
}


void SPSimulatorItem::setContactCullingDepth(double value)    
{
    impl->contactCullingDepth = value;
}

    
void SPSimulatorItem::setErrorCriterion(double value)    
{
    impl->errorCriterion = value;
}

    
void SPSimulatorItem::setMaxNumIterations(int value)
{
    impl->maxNumIterations = value;   
}


void SPSimulatorItem::setContactCorrectionDepth(double value)
{
    impl->contactCorrectionDepth = value;
}


void SPSimulatorItem::setContactCorrectionVelocityRatio(double value)
{
    impl->contactCorrectionVelocityRatio = value;
}


void SPSimulatorItem::setEpsilon(double epsilon)
{
    impl->epsilon = epsilon;
}


void SPSimulatorItem::set2Dmode(bool on)
{
    impl->is2Dmode = on;
}


void SPSimulatorItem::setKinematicWalkingEnabled(bool on)
{
    impl->isKinematicWalkingEnabled = on;
}


ItemPtr SPSimulatorItem::doDuplicate() const
{
    return new SPSimulatorItem(*this);
}


SimulationBodyPtr SPSimulatorItem::createSimulationBody(BodyPtr orgBody)
{
    return new SimulationBody(new DyBody(*orgBody));
}


ControllerItem* SPSimulatorItem::createBodyMotionController(BodyItem* bodyItem, BodyMotionItem* bodyMotionItem)
{
    return new HighGainControllerItem(bodyItem, bodyMotionItem);
}


bool SPSimulatorItem::initializeSimulation(const std::vector<SimulationBody*>& simBodies)
{
    return impl->initializeSimulation(simBodies);
}


bool SPSimulatorItemImpl::initializeSimulation(const std::vector<SimulationBody*>& simBodies)
{
    if(ENABLE_DEBUG_OUTPUT){
        static int ntest = 0;
        os.open((string("test-log-") + boost::lexical_cast<string>(ntest++) + ".log").c_str());
        os << setprecision(30);
    }

    // if(integrationMode.is(SPSimulatorItem::EULER_INTEGRATION)){
        world.setEulerMethod();
    //} else if(integrationMode.is(SPSimulatorItem::RUNGE_KUTTA_INTEGRATION)){
    //    world.setRungeKuttaMethod();
    //}
    world.setGravityAcceleration(gravity);
    world.enableSensors(true);
    world.setTimeStep(self->worldTimeStep());
    world.setCurrentTime(0.0);

    SPConstraintForceSolver& cfs = world.constraintForceSolver;

    cfs.setGaussSeidelErrorCriterion(errorCriterion.value());
    cfs.setGaussSeidelMaxNumIterations(maxNumIterations);
    cfs.setContactDepthCorrection(
        contactCorrectionDepth.value(), contactCorrectionVelocityRatio.value());

    world.clearBodies();
    bodyIndexMap.clear();
    for(size_t i=0; i < simBodies.size(); ++i){
        addBody(simBodies[i]);
    }

    cfs.setFriction(staticFriction, slipFriction);
    cfs.setContactCullingDistance(contactCullingDistance.value());
    cfs.setContactCullingDepth(contactCullingDepth.value());
    cfs.setCoefficientOfRestitution(epsilon);
    cfs.setCollisionDetector(self->collisionDetector());
    
    if(is2Dmode){
        cfs.set2Dmode(true);
    }

    world.initialize();

    return true;
}


void SPSimulatorItemImpl::addBody(SimulationBody* simBody)
{
    DyBody* body = static_cast<DyBody*>(simBody->body());

    DyLink* rootLink = body->rootLink();
    rootLink->v().setZero();
    rootLink->dv().setZero();
    rootLink->w().setZero();
    rootLink->dw().setZero();
    rootLink->vo().setZero();
    rootLink->dvo().setZero();

    bool isHighGainMode = dynamicsMode.is(SPSimulatorItem::HG_DYNAMICS);
    if(dynamic_cast<HighGainControllerItem*>(simBody->controller())){
        isHighGainMode = true;
    }

    for(int i=0; i < body->numLinks(); ++i){
        Link* link = body->link(i);
        link->u() = 0.0;
        link->dq() = 0.0;
        link->ddq() = 0.0;
    }
    
    body->clearExternalForces();
    body->calcForwardKinematics(true, true);

    if(isHighGainMode){
        ForwardDynamicsCBMPtr cbm = make_shared_aligned<ForwardDynamicsCBM>(body);
        cbm->setHighGainModeForAllJoints();
        bodyIndexMap[body] = world.addBody(body, cbm);
    } else {
        bodyIndexMap[body] = world.addBody(body);
    }
}


bool SPSimulatorItem::stepSimulation(const std::vector<SimulationBody*>& activeSimBodies)
{
    impl->world.constraintForceSolver.clearExternalForces();

    if(!impl->dynamicsMode.is(KINEMATICS)){
        impl->world.calcNextState();
        return true;
    }

    // Kinematics mode
    if(!impl->isKinematicWalkingEnabled){
        for(size_t i=0; i < activeSimBodies.size(); ++i){
            activeSimBodies[i]->body()->calcForwardKinematics(true, true);
        }
    } else {
        for(size_t i=0; i < activeSimBodies.size(); ++i){
            Body* body = activeSimBodies[i]->body();
            LeggedBodyHelper* legged = getLeggedBodyHelper(body);
            if(!legged->isValid()){
                body->calcForwardKinematics(true, true);
            } else {
                Link* supportFoot = 0;
                const int n = legged->numFeet();
                for(int i=0; i < n; ++i){
                    Link* foot = legged->footLink(i);
                    if(!supportFoot || foot->p().z() < supportFoot->p().z()){
                        supportFoot = foot;
                    }
                }
                LinkTraverse traverse(supportFoot, true, true);
                traverse.calcForwardKinematics(true, true);
            }
        }
    }
    return true;
}


void SPSimulatorItem::finalizeSimulation()
{
    if(ENABLE_DEBUG_OUTPUT){
        impl->os.close();
    }
}


void SPSimulatorItem::doPutProperties(PutPropertyFunction& putProperty)
{
    SimulatorItem::doPutProperties(putProperty);
    impl->doPutProperties(putProperty);
}


void SPSimulatorItemImpl::doPutProperties(PutPropertyFunction& putProperty)
{
    putProperty(_("Dynamics mode"), dynamicsMode,
                boost::bind((bool(Selection::*)(int))&Selection::select, &dynamicsMode, _1));
    putProperty(_("Integration mode"), integrationMode,
                boost::bind((bool(Selection::*)(int))&Selection::select, &integrationMode, _1));
    putProperty(_("Gravity"), str(gravity), boost::bind(toVector3, _1, boost::ref(gravity)));
    putProperty.decimals(3).min(0.0);
    putProperty(_("Static friction"), staticFriction, changeProperty(staticFriction));
    putProperty(_("Slip friction"), slipFriction, changeProperty(slipFriction));
    putProperty(_("Contact culling distance"), contactCullingDistance,
                (boost::bind(&FloatingNumberString::setNonNegativeValue, boost::ref(contactCullingDistance), _1)));
    putProperty(_("Contact culling depth"), contactCullingDepth,
                (boost::bind(&FloatingNumberString::setNonNegativeValue, boost::ref(contactCullingDepth), _1)));
    putProperty(_("Error criterion"), errorCriterion,
                boost::bind(&FloatingNumberString::setPositiveValue, boost::ref(errorCriterion), _1));
    putProperty.min(1.0)(_("Max iterations"), maxNumIterations, changeProperty(maxNumIterations));
    putProperty(_("CC depth"), contactCorrectionDepth,
                boost::bind(&FloatingNumberString::setNonNegativeValue, boost::ref(contactCorrectionDepth), _1));
    putProperty(_("CC v-ratio"), contactCorrectionVelocityRatio,
                boost::bind(&FloatingNumberString::setNonNegativeValue, boost::ref(contactCorrectionVelocityRatio), _1));
    putProperty(_("Kinematic walking"), isKinematicWalkingEnabled,
                changeProperty(isKinematicWalkingEnabled));
    putProperty(_("2D mode"), is2Dmode, changeProperty(is2Dmode));
}


bool SPSimulatorItem::store(Archive& archive)
{
    SimulatorItem::store(archive);
    return impl->store(archive);
}


bool SPSimulatorItemImpl::store(Archive& archive)
{
    archive.write("dynamicsMode", dynamicsMode.selectedSymbol());
    archive.write("integrationMode", integrationMode.selectedSymbol());
    write(archive, "gravity", gravity);
    archive.write("staticFriction", staticFriction);
    archive.write("slipFriction", slipFriction);
    archive.write("cullingThresh", contactCullingDistance);
    archive.write("contactCullingDepth", contactCullingDepth);
    archive.write("errorCriterion", errorCriterion);
    archive.write("maxNumIterations", maxNumIterations);
    archive.write("contactCorrectionDepth", contactCorrectionDepth);
    archive.write("contactCorrectionVelocityRatio", contactCorrectionVelocityRatio);
    archive.write("kinematicWalking", isKinematicWalkingEnabled);
    archive.write("2Dmode", is2Dmode);
    return true;
}


bool SPSimulatorItem::restore(const Archive& archive)
{
    SimulatorItem::restore(archive);
    return impl->restore(archive);
}


bool SPSimulatorItemImpl::restore(const Archive& archive)
{
    string symbol;
    if(archive.read("dynamicsMode", symbol)){
        dynamicsMode.select(symbol);
    }
    if(archive.read("integrationMode", symbol)){
        integrationMode.select(symbol);
    }
    read(archive, "gravity", gravity);
    archive.read("staticFriction", staticFriction);
    archive.read("slipFriction", slipFriction);
    contactCullingDistance = archive.get("cullingThresh", contactCullingDistance.string());
    contactCullingDepth = archive.get("contactCullingDepth", contactCullingDepth.string());
    errorCriterion = archive.get("errorCriterion", errorCriterion.string());
    archive.read("maxNumIterations", maxNumIterations);
    contactCorrectionDepth = archive.get("contactCorrectionDepth", contactCorrectionDepth.string());
    contactCorrectionVelocityRatio = archive.get("contactCorrectionVelocityRatio", contactCorrectionVelocityRatio.string());
    archive.read("kinematicWalking", isKinematicWalkingEnabled);
    archive.read("2Dmode", is2Dmode);
    return true;
}
