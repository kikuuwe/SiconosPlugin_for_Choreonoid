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



#ifndef CNOID_SICONOSPLUGIN_SICONOS_SIMULATOR_ITEM_H_INCLUDED
#define CNOID_SICONOSPLUGIN_SICONOS_SIMULATOR_ITEM_H_INCLUDED

#include <cnoid/SimulatorItem>
#include <cnoid/EigenTypes>

namespace cnoid {

class SPSimulatorItemImpl;
        
class SPSimulatorItem : public SimulatorItem
{
public:
    static void initializeClass(ExtensionManager* ext);
        
    SPSimulatorItem();
    SPSimulatorItem(const SPSimulatorItem& org);
    virtual ~SPSimulatorItem();

    enum DynamicsMode { FORWARD_DYNAMICS = 0, HG_DYNAMICS, KINEMATICS, N_DYNAMICS_MODES };
    enum IntegrationMode { EULER_INTEGRATION = 0, /* RUNGE_KUTTA_INTEGRATION,*/ N_INTEGRATION_MODES };

    void setDynamicsMode(int mode);
    void setIntegrationMode(int mode);
    void setGravity(const Vector3& gravity);
    void setStaticFriction(double value);
    void setSlipFriction(double value);
    void setContactCullingDistance(double value);        
    void setContactCullingDepth(double value);        
    void setErrorCriterion(double value);        
    void setMaxNumIterations(int value);
    void setContactCorrectionDepth(double value);
    void setContactCorrectionVelocityRatio(double value);
    void setEpsilon(double epsilon);
    void set2Dmode(bool on);
    void setKinematicWalkingEnabled(bool on); 

protected:
    virtual SimulationBodyPtr createSimulationBody(BodyPtr orgBody);
    virtual ControllerItem* createBodyMotionController(BodyItem* bodyItem, BodyMotionItem* bodyMotionItem);
    virtual bool initializeSimulation(const std::vector<SimulationBody*>& simBodies);
    virtual bool stepSimulation(const std::vector<SimulationBody*>& activeSimBodies);
    virtual void finalizeSimulation();
        
    virtual ItemPtr doDuplicate() const;
    virtual void doPutProperties(PutPropertyFunction& putProperty);
    virtual bool store(Archive& archive);
    virtual bool restore(const Archive& archive);

private:
    SPSimulatorItemImpl* impl;
    friend class SPSimulatorItemImpl;
};

typedef ref_ptr<SPSimulatorItem> SPSimulatorItemPtr;
}

#endif
