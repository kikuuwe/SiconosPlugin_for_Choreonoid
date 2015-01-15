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



#ifndef CNOID_CICONOSPLUGIN_SICONOS_FORWARD_DYNAMICS_H
#define CNOID_CICONOSPLUGIN_SICONOS_FORWARD_DYNAMICS_H

#include <cnoid/BasicSensorSimulationHelper>
#include <cnoid/Link>
#include <boost/shared_ptr.hpp>


namespace cnoid
{
class DyBody;
typedef ref_ptr<DyBody> DyBodyPtr;

/**
   This class calculates the forward dynamics of a Body object
   by using the Featherstone's articulated body algorithm.
   The class also integrates motion using the Euler method or RungeKutta method.
*/
class SPForwardDynamics {

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        
    SPForwardDynamics(const DyBodyPtr& body);
    virtual ~SPForwardDynamics();
        
    void setGravityAcceleration(const Vector3& g);
    void setEulerMethod();
    void setRungeKuttaMethod();
    void setTimeStep(double timeStep);
    void enableSensors(bool on);

    virtual void initialize() = 0;
    virtual void calcNextState() = 0;

protected:

    virtual void initializeSensors();

    /**
       @brief update position/orientation using spatial velocity
       @param out_p p(t+dt)
       @param out_R R(t+dt)
       @param p0 p(t)
       @param R0 R(t)
       @param w angular velocity
       @param v0 spatial velocity
       @param dt time step[s]
    */
    static void SE3exp(Position& out_T, const Position& T0, const Vector3& w, const Vector3& vo, double dt);
		
    DyBodyPtr body;
    Vector3 g;
    double timeStep;
    bool sensorsEnabled;
    BasicSensorSimulationHelper sensorHelper;

    enum { EULER_METHOD, RUNGEKUTTA_METHOD } integrationMode;
};

typedef boost::shared_ptr<SPForwardDynamics> SPForwardDynamicsPtr;
};

#endif
