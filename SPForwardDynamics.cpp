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

#include "SPForwardDynamics.h"
#include <cnoid/DyBody>

using namespace cnoid;

SPForwardDynamics::SPForwardDynamics(const DyBodyPtr& body)
    : body(body)
{
    g.setZero();
    timeStep = 0.005;

    integrationMode = RUNGEKUTTA_METHOD;
    sensorsEnabled = false;
}


SPForwardDynamics::~SPForwardDynamics()
{

}


void SPForwardDynamics::setTimeStep(double ts)
{
    timeStep = ts;
}


void SPForwardDynamics::setGravityAcceleration(const Vector3& g)
{
    this->g = g;
}


void SPForwardDynamics::setEulerMethod()
{
    integrationMode = EULER_METHOD;
}


void SPForwardDynamics::setRungeKuttaMethod()
{
    integrationMode = RUNGEKUTTA_METHOD;
}


void SPForwardDynamics::enableSensors(bool on)
{
    sensorsEnabled = on;
}


/// function from Murray, Li and Sastry p.42
void SPForwardDynamics::SE3exp(Position& out_T, const Position& T0, const Vector3& w, const Vector3& vo, double dt)
{
    double norm_w = w.norm();
	
    if(norm_w < std::numeric_limits<double>::epsilon()) {
        out_T.linear() = T0.linear();
        out_T.translation() = T0.translation() + vo * dt;
    } else {
        double th = norm_w * dt;
        Vector3 w_n = w / norm_w;
        Vector3 vo_n = vo / norm_w;
        const Matrix3 R(AngleAxisd(th, w_n));
        out_T.translation() =
            R * T0.translation() + (Matrix3::Identity() - R) * w_n.cross(vo_n) + (w_n * w_n.transpose()) * vo_n * th;
        out_T.linear() = R * T0.linear();
    }
}


void SPForwardDynamics::initializeSensors()
{
    body->initializeDeviceStates();

    if(sensorsEnabled){
        sensorHelper.initialize(body, timeStep, g);
    }
}
