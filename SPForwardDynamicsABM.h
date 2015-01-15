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



#ifndef CNOID_SICONOSPLUGIN_SICONOS_FORWARD_DYNAMICS_ABM_H_INCLUDED
#define CNOID_SICONOSPLUGIN_SICONOS_FORWARD_DYNAMICS_ABM_H_INCLUDED

#include "SPForwardDynamics.h"
 
namespace cnoid
{
/**
   Forward dynamics calculation using Featherstone's Articulated Body Method (ABM)
*/
class SPForwardDynamicsABM : public SPForwardDynamics
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        
    SPForwardDynamicsABM(DyBodyPtr body);
    ~SPForwardDynamicsABM();
        
    virtual void initialize();
    virtual void calcNextState();

private:
        
    void calcMotionWithEulerMethod();
    void integrateRungeKuttaOneStep(double r, double dt);
    void calcMotionWithRungeKuttaMethod();

    /**
       compute position/orientation/velocity
    */
    void calcABMPhase1();

    /**
       compute articulated inertia
    */
    void calcABMPhase2();
    void calcABMPhase2Part1();
    void calcABMPhase2Part2();

    /**
       compute joint acceleration/spatial acceleration
    */
    void calcABMPhase3();

    inline void calcABMFirstHalf();
    inline void calcABMLastHalf();

    void updateForceSensors();

    // Buffers for the Runge Kutta Method
    Position T0;
    Vector3 vo0;
    Vector3 w0;
    std::vector<double> q0;
    std::vector<double> dq0;
		
    Vector3 vo;
    Vector3 w;
    Vector3 dvo;
    Vector3 dw;
    std::vector<double> dq;
    std::vector<double> ddq;
};
	
};

#endif
