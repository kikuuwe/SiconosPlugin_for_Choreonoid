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



#include "SPWorld.h"
#include "SPForwardDynamicsABM.h"
#include "SPForwardDynamicsCBM.h"
#include <cnoid/DyBody>
#include <cnoid/EigenUtil>
#include <string>
#include <iostream>

using namespace std;
using namespace boost;
using namespace cnoid;

static const double DEFAULT_GRAVITY_ACCELERATION = 9.80665;

static const bool debugMode = false;


SPWorldBase::SPWorldBase()
{
    currentTime_ = 0.0;
    timeStep_ = 0.005;

    g << 0.0, 0.0, -DEFAULT_GRAVITY_ACCELERATION;

    isEulerMethod =false;
    sensorsAreEnabled = false;
    numRegisteredLinkPairs = 0;
}


SPWorldBase::~SPWorldBase()
{

}


int SPWorldBase::bodyIndex(const std::string& name) const
{
    NameToIndexMap::const_iterator p = nameToBodyIndexMap.find(name);
    return (p != nameToBodyIndexMap.end()) ? p->second : -1;
}


const DyBodyPtr& SPWorldBase::body(int index) const
{
    static const DyBodyPtr null;
    
    if(index < 0 || (int)bodyInfoArray.size() <= index){
        return null;
    }
    return bodyInfoArray[index].body; 
}


const DyBodyPtr& SPWorldBase::body(const std::string& name) const
{
    static const DyBodyPtr null;

    int idx = bodyIndex(name);
    if(idx < 0 || (int)bodyInfoArray.size() <= idx){
        return null;
    }
    return bodyInfoArray[idx].body;
}


void SPWorldBase::setTimeStep(double ts)
{
    timeStep_ = ts;
}


void SPWorldBase::setCurrentTime(double time)
{
    currentTime_ = time;
}


void SPWorldBase::setGravityAcceleration(const Vector3& g)
{
    this->g = g;
}


void SPWorldBase::enableSensors(bool on)
{
    sensorsAreEnabled = on;
}


void SPWorldBase::initialize()
{
    const int n = bodyInfoArray.size();

    for(int i=0; i < n; ++i){

        BodyInfo& info = bodyInfoArray[i];
        DyBodyPtr& body = info.body;

        if(!info.forwardDynamics){
            info.forwardDynamics = make_shared_aligned<SPForwardDynamicsABM>(body);
        }
        
        if(isEulerMethod){ 
            info.forwardDynamics->setEulerMethod();
        } else {
            info.forwardDynamics->setRungeKuttaMethod();
        }
        info.forwardDynamics->setGravityAcceleration(g);
        info.forwardDynamics->setTimeStep(timeStep_);
        info.forwardDynamics->enableSensors(sensorsAreEnabled);
        info.forwardDynamics->initialize();
    }
}


void SPWorldBase::setVirtualJointForces()
{
    for(size_t i=0; i < bodyInfoArray.size(); ++i){
        BodyInfo& info = bodyInfoArray[i];
        if(info.hasVirtualJointForces){
            info.body->setVirtualJointForces();
        }
    }
}


void SPWorldBase::calcNextState()
{
    if(debugMode){
        cout << "World current time = " << currentTime_ << endl;
    }
    const int n = bodyInfoArray.size();

    for(int i=0; i < n; ++i){
        BodyInfo& info = bodyInfoArray[i];
        info.forwardDynamics->calcNextState();
    }
    currentTime_ += timeStep_;
}


int SPWorldBase::addBody(const DyBodyPtr& body)
{
    if(!body->name().empty()){
        nameToBodyIndexMap[body->name()] = bodyInfoArray.size();
    }
    BodyInfo info;
    info.body = body;
    info.hasVirtualJointForces = body->hasVirtualJointForces();
    bodyInfoArray.push_back(info);

    return bodyInfoArray.size() - 1;
}


int SPWorldBase::addBody(const DyBodyPtr& body, const SPForwardDynamicsPtr& forwardDynamics)
{
    int index = addBody(body);
    bodyInfoArray[index].forwardDynamics = forwardDynamics;
    return index;
}


void SPWorldBase::clearBodies()
{
    nameToBodyIndexMap.clear();
    bodyInfoArray.clear();
}


void SPWorldBase::clearCollisionPairs()
{
    linkPairKeyToIndexMap.clear();
    numRegisteredLinkPairs = 0;
}


void SPWorldBase::setEulerMethod()
{
    isEulerMethod = true;
}


void SPWorldBase::setRungeKuttaMethod()
{
    isEulerMethod = false;
}


std::pair<int,bool> SPWorldBase::getIndexOfLinkPairs(DyLink* link1, DyLink* link2)
{
    int index = -1;
    bool isRegistered = false;

    if(link1 != link2){

        LinkPairKey linkPair;
        if(link1 < link2){
            linkPair.link1 = link1;
            linkPair.link2 = link2;
        } else {
            linkPair.link1 = link2;
            linkPair.link2 = link1;
        }

        LinkPairKeyToIndexMap::iterator p = linkPairKeyToIndexMap.find(linkPair);

        if(p != linkPairKeyToIndexMap.end()){
            index = p->second;
            isRegistered = true;
        } else {
            index = numRegisteredLinkPairs++;
            linkPairKeyToIndexMap[linkPair] = index;
        }
    }

    return std::make_pair(index, isRegistered);
}


bool SPWorldBase::LinkPairKey::operator<(const LinkPairKey& pair2) const
{
    if(link1 < pair2.link1){
        return true;
    } else if(link1 == pair2.link1){
        return (link2 < pair2.link2);
    } else {
        return false;
    }
}
