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



#include "SPForwardDynamicsCBM.h"
#include <cnoid/DyBody>
#include <cnoid/LinkTraverse>
#include <cnoid/EigenUtil>
#include <boost/bind.hpp>
#include <iostream>

using namespace std;
using namespace boost;
using namespace cnoid;

static const bool CALC_ALL_JOINT_TORQUES = false;
static const bool ROOT_ATT_NORMALIZATION_ENABLED = false;

static const bool debugMode = false;

#include <boost/format.hpp>

template<class TMatrix>
static void putMatrix(TMatrix& M, char* name)
{
    if(M.cols() == 1){
        std::cout << "Vector " << name << M << std::endl;
    } else {
        std::cout << "Matrix " << name << ": \n";
        for(size_t i=0; i < M.rows(); i++){
            for(size_t j=0; j < M.cols(); j++){
                std::cout << boost::format(" %6.3f ") % M(i, j);
            }
            std::cout << std::endl;
        }
    }
}

template<class TVector>
static void putVector(TVector& M, char* name)
{
    std::cout << "Vector " << name << M << std::endl;
}


SPForwardDynamicsCBM::SPForwardDynamicsCBM(const DyBodyPtr& body) :
    SPForwardDynamics(body)
{
    highGainModeLinkFlag.resize(body->numLinks(), false);
}


SPForwardDynamicsCBM::~SPForwardDynamicsCBM()
{

}


void SPForwardDynamicsCBM::setHighGainModeForAllJoints()
{
    highGainModeLinkFlag.set();
    highGainModeLinkFlag.reset(0);
}


void SPForwardDynamicsCBM::setHighGainMode(int linkIndex, bool on)
{
    highGainModeLinkFlag[linkIndex] = on;
}


void SPForwardDynamicsCBM::initialize()
{
    DyLink* root = body->rootLink();
    unknown_rootDof = (root->isFreeJoint() && !highGainModeLinkFlag[0]) ? 6 : 0;
    given_rootDof = (root->isFreeJoint() && highGainModeLinkFlag[0]) ? 6 : 0;

    torqueModeJoints.clear();
    highGainModeJoints.clear();

    const int numLinks = body->numLinks();

    for(int i=1; i < numLinks; ++i){
        DyLink* link = body->link(i);
        if(link->isRotationalJoint() || link->isSlideJoint()){
            if(highGainModeLinkFlag[i]){
                highGainModeJoints.push_back(link);
            } else {
                torqueModeJoints.push_back(link);
            }
        }
    }

    const int n = unknown_rootDof + torqueModeJoints.size();
    const int m = highGainModeJoints.size();
    
    isNoUnknownAccelMode = (n == 0);

    M11.resize(n, n);
    M12.resize(n, given_rootDof + m);
    b1. resize(n, 1);
    c1. resize(n);
    d1. resize(n);

    qGiven.  resize(m);
    dqGiven. resize(m);
    ddqGiven.resize(given_rootDof + m);

    qGivenPrev. resize(m);
    dqGivenPrev.resize(m);

    ddqorg.resize(numLinks);
    uorg.  resize(numLinks);

    calcPositionAndVelocityFK();

    if(!isNoUnknownAccelMode){
        calcMassMatrix();
    }

    ddqGivenCopied = false;

    initializeSensors();

    if(integrationMode == RUNGEKUTTA_METHOD){

        q0. resize(numLinks);
        dq0.resize(numLinks);
        dq. resize(numLinks);
        ddq.resize(numLinks);

        preserveHighGainModeJointState();
    }
}


void SPForwardDynamicsCBM::solveUnknownAccels()
{
    if(!isNoUnknownAccelMode){
        initializeAccelSolver();
        solveUnknownAccels(fextTotal, tauextTotal);
    }
}


inline void SPForwardDynamicsCBM::calcAccelFKandForceSensorValues()
{
    Vector3 f, tau;
    calcAccelFKandForceSensorValues(body->rootLink(), f, tau);
}


void SPForwardDynamicsCBM::calcNextState()
{
    if(isNoUnknownAccelMode && !sensorHelper.isActive()){

        calcPositionAndVelocityFK();

    } else {

        switch(integrationMode){

        case EULER_METHOD:
            calcMotionWithEulerMethod();
            break;
			
        case RUNGEKUTTA_METHOD:
            calcMotionWithRungeKuttaMethod();
            break;
        }
		
        if(ROOT_ATT_NORMALIZATION_ENABLED && unknown_rootDof){
            normalizeRotation(body->rootLink()->T());
        }
    }
		
    if(sensorHelper.isActive()){
        updateForceSensors();
        sensorHelper.updateGyroAndAccelSensors();
    }
	
    ddqGivenCopied = false;
}


void SPForwardDynamicsCBM::calcMotionWithEulerMethod()
{
    sumExternalForces();
    solveUnknownAccels();
    calcAccelFKandForceSensorValues();

    DyLink* root = body->rootLink();

    if(unknown_rootDof){
        Position T;
        SE3exp(T, root->T(), root->w(), root->vo(), timeStep);
        root->T() = T;
        root->vo() += root->dvo() * timeStep;
        root->w()  += root->dw()  * timeStep;
    }

    const int n = torqueModeJoints.size();
    for(int i=0; i < n; ++i){
        DyLink* link = torqueModeJoints[i];
        link->q()  += link->dq()  * timeStep;
        link->dq() += link->ddq() * timeStep;
    }

    calcPositionAndVelocityFK();
    calcMassMatrix();
}


void SPForwardDynamicsCBM::calcMotionWithRungeKuttaMethod()
{
    const int numHighGainJoints = highGainModeJoints.size();
    DyLink* root = body->rootLink();

    if(given_rootDof){
        pGiven = root->p();
        RGiven = root->R();
        voGiven = root->vo();
        wGiven = root->w();
        root->p() = pGivenPrev;
        root->R() = RGivenPrev;
        root->vo() = voGivenPrev;
        root->w() = wGivenPrev;
    }

    for(int i=0; i < numHighGainJoints; ++i){
        DyLink* link = highGainModeJoints[i];
        qGiven [i] = link->q();
        dqGiven[i] = link->dq();
        link->q()  = qGivenPrev[i];
        link->dq() = dqGivenPrev[i];
    }
  
    if(unknown_rootDof || given_rootDof){
        T0  = root->T();
        vo0 = root->vo();
        w0  = root->w();
    }

    vo.setZero();
    w.setZero();
    dvo.setZero();
    dw.setZero();

    const int numLinks = body->numLinks();
    for(int i=1; i < numLinks; ++i){
        const DyLink* link = body->link(i);
        q0 [i] = link->q();
        dq0[i] = link->dq();
        dq [i] = 0.0;
        ddq[i] = 0.0;
    }

    sumExternalForces();

    solveUnknownAccels();
    calcAccelFKandForceSensorValues();

    integrateRungeKuttaOneStep(1.0 / 6.0, timeStep / 2.0);

    calcPositionAndVelocityFK();
    calcMassMatrix();
    solveUnknownAccels();

    integrateRungeKuttaOneStep(2.0 / 6.0, timeStep / 2.0);

    calcPositionAndVelocityFK();
    calcMassMatrix();
    solveUnknownAccels();

    integrateRungeKuttaOneStep(2.0 / 6.0, timeStep);

    calcPositionAndVelocityFK();
    calcMassMatrix();
    solveUnknownAccels();

    if(unknown_rootDof){
        SE3exp(root->T(), T0, w0, vo0, timeStep);
        root->vo() = vo0 + (dvo + root->dvo() / 6.0) * timeStep;
        root->w()  = w0  + (dw  + root->dw()  / 6.0) * timeStep;
    }
    if(given_rootDof){
        root->p() = pGiven;
        root->R() = RGiven;
        root->vo() = voGiven;
        root->w() = wGiven;
    }

    for(size_t i=0; i < torqueModeJoints.size(); ++i){
        DyLink* link = torqueModeJoints[i];
        const int index = link->index();
        link->q()  = q0 [index] + (dq [index] + link->dq()  / 6.0) * timeStep;
        link->dq() = dq0[index] + (ddq[index] + link->ddq() / 6.0) * timeStep;
    }
    for(size_t i=0; i < highGainModeJoints.size(); ++i){
        DyLink* link = highGainModeJoints[i];
        link->q()  = qGiven [i];
        link->dq() = dqGiven[i];
    }

    calcPositionAndVelocityFK();
    calcMassMatrix();

    preserveHighGainModeJointState();
}


void SPForwardDynamicsCBM::integrateRungeKuttaOneStep(double r, double dt)
{
    DyLink* root = body->rootLink();

    if(unknown_rootDof || given_rootDof){
        SE3exp(root->T(), T0, root->w(), root->vo(), dt);
        root->vo() = vo0 + root->dvo() * dt;
        root->w()  = w0  + root->dw()  * dt;

        vo  += r * root->vo();
        w   += r * root->w();
        dvo += r * root->dvo();
        dw  += r * root->dw();
    }

    const int n = body->numLinks();
    for(int i=1; i < n; ++i){
        DyLink* link = body->link(i);
        link->q()  = q0 [i] + dt * link->dq();
        link->dq() = dq0[i] + dt * link->ddq();
        dq [i] += r * link->dq();
        ddq[i] += r * link->ddq();
    }
}


void SPForwardDynamicsCBM::preserveHighGainModeJointState()
{
    if(given_rootDof){
        const DyLink* root = body->rootLink();
        pGivenPrev = root->p();
        RGivenPrev = root->R();
        voGivenPrev = root->vo();
        wGivenPrev = root->w();
    }
    for(size_t i=0; i < highGainModeJoints.size(); ++i){
        const DyLink* link = highGainModeJoints[i];
        qGivenPrev [i] = link->q();
        dqGivenPrev[i] = link->dq();
    }
}


void SPForwardDynamicsCBM::calcPositionAndVelocityFK()
{
    DyLink* root = body->rootLink();
    root_w_x_v = root->w().cross(root->vo() + root->w().cross(root->p()));
    if(given_rootDof){
        root->vo() = root->v() - root->w().cross(root->p());
    }

    const LinkTraverse& traverse = body->linkTraverse();
    const int n = traverse.numLinks();

    for(int i=0; i < n; ++i){
        DyLink* link = static_cast<DyLink*>(traverse[i]);
        const DyLink* parent = link->parent();

        if(parent){

            switch(link->jointType()){

            case Link::SLIDE_JOINT:
                link->p().noalias() = parent->R() * (link->b() + link->q() * link->d()) + parent->p();
                link->R() = parent->R();
                link->sw().setZero();
                link->sv().noalias() = parent->R() * link->d();
                link->w() = parent->w();
                break;

            case Link::ROTATIONAL_JOINT:
                link->R().noalias() = parent->R() * AngleAxisd(link->q(), link->a());
                link->p().noalias() = parent->R() * link->b() + parent->p();
                link->sw().noalias() = parent->R() * link->a();
                link->sv() = link->p().cross(link->sw());
                link->w() = link->dq() * link->sw() + parent->w();
                break;

            case Link::FIXED_JOINT:
            default:
                link->p().noalias() = parent->R() * link->b() + parent->p();
                link->R() = parent->R();
                link->w() = parent->w();
                link->vo() = parent->vo();
                link->sw().setZero();
                link->sv().setZero();
                link->cv().setZero();
                link->cw().setZero();
                goto COMMON_CALCS_FOR_ALL_JOINT_TYPES;
            }

            link->vo() = link->dq() * link->sv() + parent->vo();

            const Vector3 dsv = parent->w().cross(link->sv()) + parent->vo().cross(link->sw());
            const Vector3 dsw = parent->w().cross(link->sw());
            link->cv() = link->dq() * dsv;
            link->cw() = link->dq() * dsw;
        }

COMMON_CALCS_FOR_ALL_JOINT_TYPES:

        /// \todo remove this  equation
        link->v() = link->vo() + link->w().cross(link->p());

        link->wc().noalias() = link->R() * link->c() + link->p();
        const Matrix3 Iw = link->R() * link->I() * link->R().transpose();
        const Matrix3 c_hat = hat(link->wc());
        link->Iww().noalias() = link->m() * c_hat * c_hat.transpose() + Iw;
        link->Iwv() = link->m() * c_hat;

        const Vector3 P = link->m() * (link->vo() + link->w().cross(link->wc()));
        const Vector3 L = link->Iww() * link->w() + link->m() * link->wc().cross(link->vo());

        link->pf()   = link->w().cross(P);
        link->ptau() = link->vo().cross(P) + link->w().cross(L);
    }
}


/**
   calculate the mass matrix using the unit vector method
   \todo replace the unit vector method here with
   a more efficient method that only requires O(n) computation time
*/
void SPForwardDynamicsCBM::calcMassMatrix()
{
    DyLink* root = body->rootLink();
    const int numLinks = body->numLinks();

    // preserve and clear the joint accelerations
    for(int i=1; i < numLinks; ++i){
        DyLink* link = body->link(i);
        ddqorg[i] = link->ddq();
        uorg  [i] = link->u();
        link->ddq() = 0.0;
    }

    // preserve and clear the root link acceleration
    dvoorg = root->dvo();
    dworg  = root->dw();
    root->dvo() = -g - root_w_x_v;   // dv = g, dw = 0
    root->dw().setZero();
	
    setColumnOfMassMatrix(b1, 0);

    if(unknown_rootDof){
        for(int i=0; i < 3; ++i){
            root->dvo()[i] += 1.0;
            setColumnOfMassMatrix(M11, i);
            root->dvo()[i] -= 1.0;
        }
        for(int i=0; i < 3; ++i){
            root->dw()[i] = 1.0;
            Vector3 dw_x_p = root->dw().cross(root->p());
            root->dvo() -= dw_x_p;
            setColumnOfMassMatrix(M11, i + 3);
            root->dvo() += dw_x_p;
            root->dw()[i] = 0.0;
        }
    }
    if(given_rootDof){
        for(int i=0; i < 3; ++i){
            root->dvo()[i] += 1.0;
            setColumnOfMassMatrix(M12, i);
            root->dvo()[i] -= 1.0;
        }
        for(int i=0; i < 3; ++i){
            root->dw()[i] = 1.0;
            Vector3 dw_x_p = root->dw().cross(root->p());
            root->dvo() -= dw_x_p;
            setColumnOfMassMatrix(M12, i + 3);
            root->dvo() += dw_x_p;
            root->dw()[i] = 0.0;
        }
    }

    for(size_t i=0; i < torqueModeJoints.size(); ++i){
        DyLink* link = torqueModeJoints[i];
        link->ddq() = 1.0;
        int j = i + unknown_rootDof;
        setColumnOfMassMatrix(M11, j);
        M11(j, j) += link->Jm2(); // motor inertia
        link->ddq() = 0.0;
    }
    for(size_t i=0; i < highGainModeJoints.size(); ++i){
        DyLink* link = highGainModeJoints[i];
        link->ddq() = 1.0;
        int j = i + given_rootDof;
        setColumnOfMassMatrix(M12, j);
        link->ddq() = 0.0;
    }

    // subtract the constant term
    for(int i=0; i < M11.cols(); ++i){
        M11.col(i) -= b1;
    }
    for(int i=0; i < M12.cols(); ++i){
        M12.col(i) -= b1;
    }

    for(int i=1; i < numLinks; ++i){
        DyLink* link = body->link(i);
        link->ddq() = ddqorg[i];
        link->u()   = uorg  [i];
    }
    root->dvo() = dvoorg;
    root->dw()  = dworg;

    accelSolverInitialized = false;
}


void SPForwardDynamicsCBM::setColumnOfMassMatrix(MatrixXd& M, int column)
{
    Vector3 f;
    Vector3 tau;
    DyLink* root = body->rootLink();
    calcInverseDynamics(root, f, tau);

    MatrixXd::ColXpr col = M.col(column);

    if(unknown_rootDof){
        tau -= root->p().cross(f);
        col << f, tau;
        /*
          col.head(3) = f;
          col.segment(3, 3) = tau;
        */
    }

    for(size_t i=0; i < torqueModeJoints.size(); ++i){
        col(i + unknown_rootDof) = torqueModeJoints[i]->u();
    }
}


void SPForwardDynamicsCBM::calcInverseDynamics(DyLink* link, Vector3& out_f, Vector3& out_tau)
{
    const DyLink* parent = link->parent();
    if(parent){
        link->dvo() = parent->dvo() + link->cv() + link->sv() * link->ddq();
        link->dw()  = parent->dw()  + link->cw() + link->sw() * link->ddq();
    }

    out_f = link->pf();
    out_tau = link->ptau();

    if(link->child()){
        Vector3 f_c;
        Vector3 tau_c;
        calcInverseDynamics(link->child(), f_c, tau_c);
        out_f += f_c;
        out_tau += tau_c;
    }

    out_f  .noalias() += link->m()   * link->dvo() + link->Iwv().transpose() * link->dw();
    out_tau.noalias() += link->Iwv() * link->dvo() + link->Iww()             * link->dw();

    link->u() = link->sv().dot(out_f) + link->sw().dot(out_tau);

    if(link->sibling()){
        Vector3 f_s;
        Vector3 tau_s;
        calcInverseDynamics(link->sibling(), f_s, tau_s);
        out_f += f_s;
        out_tau += tau_s;
    }
}


void SPForwardDynamicsCBM::sumExternalForces()
{
    fextTotal.setZero();
    tauextTotal.setZero();

    const int n = body->numLinks();
    for(int i=0; i < n; ++i){
        const DyLink* link = body->link(i);
        fextTotal   += link->f_ext();
        tauextTotal += link->tau_ext();
    }

    tauextTotal -= body->rootLink()->p().cross(fextTotal);
}

void SPForwardDynamicsCBM::calcd1(DyLink* link, Vector3& out_f, Vector3& out_tau)
{
    out_f = -link->f_ext();
    out_tau = -link->tau_ext();

    if(link->child()){
        Vector3 f_c;
        Vector3 tau_c;
        calcd1(link->child(), f_c, tau_c);
        out_f += f_c;
        out_tau += tau_c;
    }

    link->u() = link->sv().dot(out_f) + link->sw().dot(out_tau);

    if(link->sibling()){
        Vector3 f_s;
        Vector3 tau_s;
        calcd1(link->sibling(), f_s, tau_s);
        out_f += f_s;
        out_tau += tau_s;
    }
}

void SPForwardDynamicsCBM::initializeAccelSolver()
{
    if(!accelSolverInitialized){

        if(!ddqGivenCopied){
            if(given_rootDof){
                DyLink* root = body->rootLink();
                root->dvo() = root->dv() - root->dw().cross(root->p()) - root->w().cross(root->v());
                ddqGiven.head(3) = root->dvo();
                ddqGiven.segment(3, 3) = root->dw();
            }
            for(size_t i=0; i < highGainModeJoints.size(); ++i){
                ddqGiven(given_rootDof + i) = highGainModeJoints[i]->ddq();
            }
            ddqGivenCopied = true;
        }

        b1.noalias() += M12 * ddqGiven;
        
        for(int i=1; i < body->numLinks(); ++i){
            const DyLink* link = body->link(i);
            uorg  [i] = link->u();
        }
        Vector3 f, tau;
        DyLink* root = body->rootLink();
        calcd1(root, f, tau);
        for(int i=0; i < unknown_rootDof; i++){
            d1(i, 0) = 0;
        }
        for(size_t i=0; i < torqueModeJoints.size(); ++i){
            d1(i + unknown_rootDof, 0) = torqueModeJoints[i]->u();
        }
        for(int i=1; i < body->numLinks(); ++i){
            DyLink* link = body->link(i);
            link->u() = uorg  [i];
        }

        accelSolverInitialized = true;
    }
}

//for ConstraintForceSolver 
void SPForwardDynamicsCBM::solveUnknownAccels
(DyLink* link, const Vector3& fext, const Vector3& tauext, const Vector3& rootfext, const Vector3& roottauext)
{
    const Vector3 fextorg = link->f_ext();
    const Vector3 tauextorg = link->tau_ext();
    link->f_ext() = fext;
    link->tau_ext() = tauext;
    for(int i=1; i < body->numLinks(); ++i){
        const DyLink* link = body->link(i);
        uorg[i] = link->u();
    }
    Vector3 f, tau;
    DyLink* root = body->rootLink();
    calcd1(root, f, tau);
    for(int i=0; i < unknown_rootDof; i++){
        d1(i, 0) = 0;
    }
    for(size_t i=0; i < torqueModeJoints.size(); ++i){
        d1(i + unknown_rootDof, 0) = torqueModeJoints[i]->u();
    }
    for(int i=1; i < body->numLinks(); ++i){
        DyLink* link = body->link(i);
        link->u() = uorg[i];
    }
    link->f_ext() = fextorg;
    link->tau_ext() = tauextorg;

    solveUnknownAccels(rootfext, roottauext);
}

void SPForwardDynamicsCBM::solveUnknownAccels(const Vector3& fext, const Vector3& tauext)
{
    if(unknown_rootDof){
        c1.head(3) = fext;
        c1.segment(3, 3) = tauext;
    }

    for(size_t i=0; i < torqueModeJoints.size(); ++i){
        c1(i + unknown_rootDof) = torqueModeJoints[i]->u();
    }

    c1 -= d1;
    c1 -= b1.col(0);

    const VectorXd a(M11.colPivHouseholderQr().solve(c1));
    
    if(unknown_rootDof){
        DyLink* root = body->rootLink();
        root->dw() = a.segment(3, 3);
        const Vector3 dv = a.head(3);
        root->dvo() = dv - root->dw().cross(root->p()) - root_w_x_v;
    }

    for(size_t i=0; i < torqueModeJoints.size(); ++i){
        DyLink* link = torqueModeJoints[i];
        link->ddq() = a(i + unknown_rootDof);
    }
}


void SPForwardDynamicsCBM::calcAccelFKandForceSensorValues(DyLink* link, Vector3& out_f, Vector3& out_tau)
{
    const DyLink* parent = link->parent();
    if(parent){
        link->dvo() = parent->dvo() + link->cv() + link->sv() * link->ddq();
        link->dw()  = parent->dw()  + link->cw() + link->sw() * link->ddq();
    }

    out_f   = link->pf();
    out_tau = link->ptau();

    for(DyLink* child = link->child(); child; child = child->sibling()){
        Vector3 f, tau;
        calcAccelFKandForceSensorValues(child, f, tau);
        out_f   += f;
        out_tau += tau;
    }

    ForceSensorInfo& info = forceSensorInfo[link->index()];

    if(CALC_ALL_JOINT_TORQUES || info.hasSensorsAbove){

        Vector3 fg(link->m() * g);
        Vector3 tg(link->wc().cross(fg));

        out_f   -= fg;
        out_tau -= tg;

        out_f   -= link->f_ext();
        out_tau -= link->tau_ext();

        out_f  .noalias() += link->m()   * link->dvo() + link->Iwv().transpose() * link->dw();
        out_tau.noalias() += link->Iwv() * link->dvo() + link->Iww()             * link->dw();

        if(CALC_ALL_JOINT_TORQUES && highGainModeLinkFlag[link->index()]){
            link->u() = link->sv().dot(out_f) + link->sw().dot(out_tau);
        }

        if(info.hasSensor){
            info.f = -out_f;
            info.tau = -out_tau;
        }
    }
}

void SPForwardDynamicsCBM::initializeSensors()
{
    SPForwardDynamics::initializeSensors();

    const int n = body->numLinks();

    forceSensorInfo.resize(n);

    if(sensorsEnabled){
        const DeviceList<ForceSensor>& forceSensors = sensorHelper.forceSensors();
        for(size_t i=0; i < forceSensors.size(); ++i){
            const ForceSensor* sensor = forceSensors.get(i);
            forceSensorInfo[sensor->link()->index()].hasSensor = true;
        }
        updateForceSensorInfo(body->rootLink(), false);
    }
}


void SPForwardDynamicsCBM::updateForceSensorInfo(DyLink* link, bool hasSensorsAbove)
{
    ForceSensorInfo& info = forceSensorInfo[link->index()];
    hasSensorsAbove |= info.hasSensor;
    info.hasSensorsAbove = hasSensorsAbove;

    for(DyLink* child = link->child(); child; child = child->sibling()){
        updateForceSensorInfo(child, hasSensorsAbove);
    }
}


void SPForwardDynamicsCBM::updateForceSensors()
{
    const DeviceList<ForceSensor>& sensors = sensorHelper.forceSensors();

    for(size_t i=0; i < sensors.size(); ++i){
        
        ForceSensor* sensor = sensors.get(i);
        const Link* link = sensor->link();

        const ForceSensorInfo& info = forceSensorInfo[sensor->link()->index()];
        const Matrix3 R = link->R() * sensor->R_local();
        const Vector3 p = link->p() + link->R() * sensor->p_local();
        
        sensor->f()   = R.transpose() * info.f;
        sensor->tau() = R.transpose() * (info.tau - p.cross(info.f));
        sensor->notifyStateChange();
    }
}
