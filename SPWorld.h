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

#ifndef CNOID_SICONOSPLUGIN_WORLD_H_INCLUDED
#define CNOID_SICONOSPLUGIN_WORLD_H_INCLUDED

#include "SPForwardDynamics.h"
#include <map>

namespace cnoid {

class DyLink;
class DyBody;
typedef ref_ptr<DyBody> DyBodyPtr;

class SPWorldBase
{
public:
    SPWorldBase();
    virtual ~SPWorldBase();

    /**
       @brief get the number of bodies in this world
       @return the number of bodies
    */ 
    int numBodies() const { return bodyInfoArray.size(); }

    /**
       @brief get body by index
       @param index of the body
       @return body
    */
    const DyBodyPtr& body(int index) const;

    /**
       @brief get body by name
       @param name of the body
       @return body
    */
    const DyBodyPtr& body(const std::string& name) const;

    /**
       @brief get forward dynamics computation method for body
       @param index index of the body
       @return forward dynamics computation method
    */
    const SPForwardDynamicsPtr& forwardDynamics(int index) const {
        return bodyInfoArray[index].forwardDynamics;
    }

    /**
       @brief get index of body by name
       @param name of the body
       @return index of the body
    */
    int bodyIndex(const std::string& name) const;

    /**
       @brief add body to this world
       @param body
       @return index of the body
       @note This must be called before initialize() is called.
    */
    int addBody(const DyBodyPtr& body);

    /**
       Use this method instead of addBody(const DyBodyPtr& body) when you want to specify
       a forward dynamics calculater.
    */
    int addBody(const DyBodyPtr& body, const SPForwardDynamicsPtr& forwardDynamics);
        
    /**
       @brief clear bodies in this world
    */
    void clearBodies();

    /**
       @brief clear collision pairs
    */
    void clearCollisionPairs();

    /**
       @brief set time step
       @param dt time step[s]
    */
    void setTimeStep(double dt);

    /**
       @brief get time step
       @return time step[s]
    */
    double timeStep(void) const { return timeStep_; }
	
    /**
       @brief set current time
       @param tm current time[s]
    */
    void setCurrentTime(double tm);

    /**
       @brief get current time
       @return current time[s]
    */
    double currentTime(void) const { return currentTime_; }
	
    /**
       @brief set gravity acceleration
       @param g gravity acceleration[m/s^2]
    */
    void setGravityAcceleration(const Vector3& g);

    /**
       @brief get gravity acceleration
       @return gravity accleration
    */
    inline const Vector3& gravityAcceleration() const { return g; }

        
    /**
       @brief enable/disable sensor simulation
       @param on true to enable, false to disable
       @note This must be called before initialize() is called.
    */
    void enableSensors(bool on);

    /**
       @brief choose euler method for integration
    */
    void setEulerMethod();

    /**
       @brief choose runge-kutta method for integration
    */
    void setRungeKuttaMethod();

    /**
       @brief initialize this world. This must be called after all bodies are registered.
    */
    virtual void initialize();

    void setVirtualJointForces();
        
    /**
       @brief compute forward dynamics and update current state
    */
    virtual void calcNextState();

    /**
       @brief get index of link pairs
       @param link1 link1
       @param link2 link2
       @return pair of index and flag. The flag is true if the pair was already registered, false othewise.
    */
    std::pair<int,bool> getIndexOfLinkPairs(DyLink* link1, DyLink* link2);

protected:
 
    double currentTime_;
    double timeStep_;

    struct BodyInfo {
        DyBodyPtr body;
        SPForwardDynamicsPtr forwardDynamics;
        bool hasVirtualJointForces;
    };
    std::vector<BodyInfo> bodyInfoArray;

    bool sensorsAreEnabled;

private:
    typedef std::map<std::string, int> NameToIndexMap;
    NameToIndexMap nameToBodyIndexMap;

    typedef std::map<DyBodyPtr, int> BodyToIndexMap;
    BodyToIndexMap bodyToIndexMap;

    Vector3 g;

    bool isEulerMethod; // Euler or Runge Kutta ?

    struct LinkPairKey {
        DyLink* link1;
        DyLink* link2;
        bool operator<(const LinkPairKey& pair2) const;
    };
    typedef std::map<LinkPairKey, int> LinkPairKeyToIndexMap;
    LinkPairKeyToIndexMap linkPairKeyToIndexMap;

    int numRegisteredLinkPairs;
		
};

template <class TConstraintForceSolver> class SPWorld : public SPWorldBase
{
public:
    TConstraintForceSolver constraintForceSolver;

    SPWorld() : constraintForceSolver(*this) { }

    virtual void initialize() {
        SPWorldBase::initialize();
        constraintForceSolver.initialize();
    }

    virtual void calcNextState(){
        SPWorldBase::setVirtualJointForces();
        constraintForceSolver.solve();
        SPWorldBase::calcNextState();
    }
};

};

#endif
