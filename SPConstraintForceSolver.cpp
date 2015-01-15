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


#ifdef __WIN32__
#define NOMINMAX
#endif

#include "SPWorld.h"
#include "SPForwardDynamicsCBM.h"
#include "SPConstraintForceSolver.h"

#include <cnoid/LinkTraverse>
#include <cnoid/DyBody>
#include <cnoid/BodyCollisionDetectorUtil>
#include <cnoid/IdPair>
#include <cnoid/EigenUtil>
#include <cnoid/AISTCollisionDetector>
#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/random.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>
#include <limits>

#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include "SiconosNumerics.h"
#include "FrictionContact3D_compute_error.h"

using namespace std;
using namespace cnoid;
// settings

static const double VEL_THRESH_OF_DYNAMIC_FRICTION = 1.0e-4;

static const bool SKIP_REDUNDANT_ACCEL_CALC = true;

static const bool ASSUME_SYMMETRIC_MATRIX = true;

static const int DEFAULT_MAX_NUM_GAUSS_SEIDEL_ITERATION = 1000;


static const int DEFAULT_NUM_GAUSS_SEIDEL_INITIAL_ITERATION = 0;
static const double DEFAULT_GAUSS_SEIDEL_ERROR_CRITERION = 1.0e-3;

static const double THRESH_TO_SWITCH_REL_ERROR = 1.0e-8;
//static const double THRESH_TO_SWITCH_REL_ERROR = numeric_limits<double>::epsilon();


static const bool ENABLE_CONTACT_DEPTH_CORRECTION = true; // KIKUUWE

// normal setting
static const double DEFAULT_CONTACT_CORRECTION_DEPTH          = 0.0001;
static const double DEFAULT_CONTACT_CORRECTION_VELOCITY_RATIO = 1.0;
static const double DEFAULT_CONTACT_CULLING_DISTANCE          = 0.005;
static const double DEFAULT_CONTACT_CULLING_DEPTH             = 0.05;


// test for mobile robots with wheels
//static const double DEFAULT_CONTACT_CORRECTION_DEPTH = 0.005;
//static const double PENETRATION_A = 500.0;
//static const double PENETRATION_B = 80.0;
//static const double NEGATIVE_VELOCITY_RATIO_FOR_PENETRATION = 10.0;
//static const bool ENABLE_CONTACT_POINT_THINNING = false;

// experimental options
static const bool ENABLE_RANDOM_STATIC_FRICTION_BASE = false;

// CP options

static const bool KIK_USE_FULL_MATRIX = false;

// debug options
static const bool CFS_DEBUG           = false;
static const bool CFS_DEBUG_VERBOSE   = false; 
static const bool CFS_DEBUG_VERBOSE_2 = false;
static const bool CFS_DEBUG_LCPCHECK = false;
static const bool CFS_MCP_DEBUG = false;
static const bool CFS_MCP_DEBUG_SHOW_ITERATION_STOP = false;

static const bool CFS_PUT_NUM_CONTACT_POINTS = false;

static const Vector3 local2dConstraintPoints[3] = {
    Vector3( 1.0, 0.0, (-sqrt(3.0) / 2.0)),
    Vector3(-1.0, 0.0, (-sqrt(3.0) / 2.0)),
    Vector3( 0.0, 0.0, ( sqrt(3.0) / 2.0))
};


/**** KIKUUWE   ***/
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;
typedef VectorXd VectorX;
/**** KIKUUWE   ***/


namespace cnoid
{
class SPCFSImpl
{
public:

    SPCFSImpl(SPWorldBase& world);
    ~SPCFSImpl();

    void initialize(void);
    void solve();
    inline void clearExternalForces();

    SPWorldBase& world;

    bool isConstraintForceOutputMode;
        
    struct ConstraintPoint {
        int globalIndex;
        Vector3 point;
        Vector3 normalTowardInside[2];
        Vector3 defaultAccel[2];
        double normalProjectionOfRelVelocityOn0;
        double depth; // position error in the case of a connection point

        double mu;
        Vector3 relVelocityOn0;
        int globalFrictionIndex;
        int numFrictionVectors;
        Vector3 frictionVector[4][2];
    };
    typedef std::vector<ConstraintPoint> ConstraintPointArray;

    struct LinkData
    {
        Vector3 dvo;
        Vector3 dw;
        Vector3 pf0;
        Vector3 ptau0;
        double uu;
        double uu0;
        double ddq;
        int numberToCheckAccelCalcSkip;
        int parentIndex;
        DyLink* link;
    };
    typedef std::vector<LinkData> LinkDataArray;

    struct BodyData
    {
        DyBodyPtr body;
        bool isStatic;
        bool hasConstrainedLinks;
        bool isTestForceBeingApplied;
        int geometryId;
        LinkDataArray linksData;

        Vector3 dpf;
        Vector3 dptau;

        /**
           If the body includes high-gain mode joints,
           the ForwardDynamisCBM object of the body is set to this pointer.
           The pointer is null when all the joints are torque mode and
           the forward dynamics is calculated by ABM.
        */
        SPForwardDynamicsCBMPtr forwardDynamicsCBM;
    };

    std::vector<BodyData> bodiesData;

    class LinkPair
    {
    public:
        virtual ~LinkPair() { }
        bool isSameBodyPair;
        int bodyIndex[2];
        BodyData* bodyData[2];
        DyLink* link[2];
        LinkData* linkData[2];
        ConstraintPointArray constraintPoints;
        bool isNonContactConstraint;
        double muStatic;
        double muDynamic;
        double contactCullingDistance;
        double contactCullingDepth;
        double epsilon;
    };
    typedef boost::shared_ptr<LinkPair> LinkPairPtr;

    CollisionDetectorPtr collisionDetector;
    vector<int> geometryIdToBodyIndexMap;
    //typedef std::map<IdPair<>, LinkPairPtr> GeometryPairToLinkPairMap;
    typedef std::map<IdPair<>, LinkPair> GeometryPairToLinkPairMap;
    GeometryPairToLinkPairMap geometryPairToLinkPairMap;
        
    double defaultStaticFriction;
    double defaultSlipFriction;
    double defaultContactCullingDistance;
    double defaultContactCullingDepth;
    double defaultCoefficientOfRestitution;

    class ExtraJointLinkPair : public LinkPair
    {
    public:
        Vector3 jointPoint[2];
        Vector3 jointConstraintAxes[3];
    };
    typedef boost::shared_ptr<ExtraJointLinkPair> ExtraJointLinkPairPtr;
    vector<ExtraJointLinkPairPtr> extraJointLinkPairs;

    bool is2Dmode;
    DyBodyPtr bodyFor2dConstraint;
    BodyData bodyDataFor2dConstraint;

    class Constrain2dLinkPair : public LinkPair
    {
    public:
        double globalYpositions[3];
    };
    typedef boost::shared_ptr<Constrain2dLinkPair> Constrain2dLinkPairPtr;
    vector<Constrain2dLinkPairPtr> constrain2dLinkPairs;
        

    std::vector<LinkPair*> constrainedLinkPairs;

    int globalNumConstraintVectors;

    int globalNumContactNormalVectors;
    int globalNumFrictionVectors;

    int prevGlobalNumConstraintVectors;
    int prevGlobalNumFrictionVectors;

    bool areThereImpacts;
    int numUnconverged;

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;
    typedef VectorXd VectorX;
        
    // Mlcp * solution + b   _|_  solution

    MatrixX Mlcp;

    // constant acceleration term when no external force is applied
    VectorX an0;
    VectorX at0;

    // constant vector of LCP
    VectorX b;

    // contact force solution: normal forces at contact points
    VectorX solution;

    // random number generator
    boost::variate_generator<boost::mt19937, boost::uniform_real<> > randomAngle;

    // for special version of gauss sidel iterative solver
    std::vector<int> frictionIndexToContactIndex;
    VectorX contactIndexToMu;
    VectorX mcpHi;

    int  maxNumGaussSeidelIteration;
    int  numGaussSeidelInitialIteration;
    double gaussSeidelErrorCriterion;
    double contactCorrectionDepth;
    double contactCorrectionVelocityRatio;

    int numGaussSeidelTotalLoops;
    int numGaussSeidelTotalCalls;
    int numGaussSeidelTotalLoopsMax;

    void initBody(const DyBodyPtr& body, BodyData& bodyData);
    void initExtraJoints(int bodyIndex);
    void init2Dconstraint(int bodyIndex);
    void setConstraintPoints();
    void extractConstraintPoints(const CollisionPair& collisionPair);
    bool setContactConstraintPoint(LinkPair& linkPair, const Collision& collision);
    void setFrictionVectors(ConstraintPoint& constraintPoint);
    void setExtraJointConstraintPoints(const ExtraJointLinkPairPtr& linkPair);
    void set2dConstraintPoints(const Constrain2dLinkPairPtr& linkPair);
    void putContactPoints();
    void solveImpactConstraints();
    void initMatrices();
    void setAccelCalcSkipInformation();
    void setDefaultAccelerationVector();
    void setAccelerationMatrix();
    void initABMForceElementsWithNoExtForce(BodyData& bodyData);
    void calcABMForceElementsWithTestForce(BodyData& bodyData, DyLink* linkToApplyForce, const Vector3& f, const Vector3& tau);
    void calcAccelsABM(BodyData& bodyData, int constraintIndex);
    void calcAccelsMM(BodyData& bodyData, int constraintIndex);

    void extractRelAccelsOfConstraintPoints
    (Eigen::Block<MatrixX>& Kxn, Eigen::Block<MatrixX>& Kxt, int testForceIndex, int constraintIndex);

    void extractRelAccelsFromLinkPairCase1
    (Eigen::Block<MatrixX>& Kxn, Eigen::Block<MatrixX>& Kxt, LinkPair& linkPair, int testForceIndex, int constraintIndex);
    void extractRelAccelsFromLinkPairCase2
    (Eigen::Block<MatrixX>& Kxn, Eigen::Block<MatrixX>& Kxt, LinkPair& linkPair, int iTestForce, int iDefault, int testForceIndex, int constraintIndex);
    void extractRelAccelsFromLinkPairCase3
    (Eigen::Block<MatrixX>& Kxn, Eigen::Block<MatrixX>& Kxt, LinkPair& linkPair, int testForceIndex, int constraintIndex);

    void copySymmetricElementsOfAccelerationMatrix
    (Eigen::Block<MatrixX>& Knn, Eigen::Block<MatrixX>& Ktn, Eigen::Block<MatrixX>& Knt, Eigen::Block<MatrixX>& Ktt);

    void clearSingularPointConstraintsOfClosedLoopConnections();
		
    void setConstantVectorAndMuBlock();
    void addConstraintForceToLinks();
    void addConstraintForceToLink(LinkPair* linkPair, int ipair);

    void solveMCPByProjectedGaussSeidel
    (const MatrixX& M, const VectorX& b, VectorX& x);
    void solveMCPByProjectedGaussSeidelMainStep
    (const MatrixX& M, const VectorX& b, VectorX& x);
    void solveMCPByProjectedGaussSeidelInitial
    (const MatrixX& M, const VectorX& b, VectorX& x, const int numIteration);

    void checkLCPResult(MatrixX& M, VectorX& b, VectorX& x);
    void checkMCPResult(MatrixX& M, VectorX& b, VectorX& x);

    ofstream os;

    template<class TMatrix>
    void putMatrix(TMatrix& M, const char *name) {
        if(M.cols() == 1){
            os << "Vector " << name << M << std::endl;
        } else {
            os << "Matrix " << name << ": \n";
            for(int i=0; i < M.rows(); i++){
                for(int j=0; j < M.cols(); j++){
                    //os << boost::format(" %6.3f ") % M(i, j);
                    os << boost::format(" %.50g ") % M(i, j);
                }
                os << std::endl;
            }
        }
    }

    template<class TVector>
    void putVector(const TVector& M, const char *name) {
        os << "Vector " << name << M << std::endl;
    }

    template<class TMatrix>
    void debugPutMatrix(const TMatrix& M, const char *name) {
        if(CFS_DEBUG_VERBOSE) putMatrix(M, name);
    }

    template<class TVector>
    void debugPutVector(const TVector& M, const char *name) {
        if(CFS_DEBUG_VERBOSE) putVector(M, name);
    }
    
	bool kik_callSiconosSolver(MatrixX& Mlcp, VectorX& b, VectorX& solution, VectorX& contactIndexToMu);
	void kik_NewBuffer(int aNC) ;
	void kik_DeleteBuffer(int prv);
	FrictionContactProblem kik_prob;
	NumericsOptions  kik_numops;
	double          *kik_reaction  ;
	double          *kik_velocity  ;
	SolverOptions   *kik_solops    ;
 
};
/*
  #ifdef _MSC_VER
  const double SPCFSImpl::PI   = 3.14159265358979323846;
  const double SPCFSImpl::PI_2 = 1.57079632679489661923;
  #endif
*/
};



SPCFSImpl::SPCFSImpl(SPWorldBase& world) :
    world(world),
    randomAngle(boost::mt19937(), boost::uniform_real<>(0.0, 2.0 * PI))
{
  defaultStaticFriction = 1.0;
  defaultSlipFriction = 1.0;
  defaultContactCullingDistance = DEFAULT_CONTACT_CULLING_DISTANCE;
  defaultContactCullingDepth = DEFAULT_CONTACT_CULLING_DEPTH;
  defaultCoefficientOfRestitution = 0.0;
  
  maxNumGaussSeidelIteration = DEFAULT_MAX_NUM_GAUSS_SEIDEL_ITERATION;
  numGaussSeidelInitialIteration = DEFAULT_NUM_GAUSS_SEIDEL_INITIAL_ITERATION;
  gaussSeidelErrorCriterion = DEFAULT_GAUSS_SEIDEL_ERROR_CRITERION;
  contactCorrectionDepth = DEFAULT_CONTACT_CORRECTION_DEPTH;
  contactCorrectionVelocityRatio = DEFAULT_CONTACT_CORRECTION_VELOCITY_RATIO;

  isConstraintForceOutputMode = false;
  is2Dmode = false;
	kik_solops    = new SolverOptions;
	kik_NewBuffer(0);
	kik_prob.dimension    = 3;
	kik_prob.M            =  new NumericsMatrix;
	kik_prob.M->matrix1    = new SparseBlockStructuredMatrix;
	kik_numops.verboseMode = 0;  /*************/
//	SICONOS_FRICTION_3D_projectionOnCylinder;
	frictionContact3D_setDefaultSolverOptions(kik_solops, SICONOS_FRICTION_3D_NSGS   );
	kik_solops->iparam[0] = maxNumGaussSeidelIteration;
//	kik_solops->iparam[1] = 1; // This makes computation slower
	kik_solops->dparam[0] = gaussSeidelErrorCriterion;//1e-10;
	kik_solops->internalSolvers->iparam[0] = maxNumGaussSeidelIteration;
	kik_solops->internalSolvers->dparam[0] = gaussSeidelErrorCriterion;
}


SPCFSImpl::~SPCFSImpl()
{
	kik_DeleteBuffer(Mlcp.rows());
	deleteSolverOptions(kik_solops);
	delete  kik_solops ;
//	freeNumericMatrix(kik_prob.M);
	delete  kik_prob.M ->matrix1;
	delete  kik_prob.M     ;
  if(CFS_DEBUG) os.close();
}


void SPCFSImpl::initBody(const DyBodyPtr& body, BodyData& bodyData)
{
    body->clearExternalForces();
    bodyData.body = body;
    bodyData.linksData.resize(body->numLinks());
    bodyData.hasConstrainedLinks = false;
    bodyData.isTestForceBeingApplied = false;
    bodyData.isStatic = body->isStaticModel();

    LinkDataArray& linksData = bodyData.linksData;
    const int n = body->numLinks();
    for(int i=0; i < n; ++i)
    {
        DyLink* link = body->link(i);
        linksData[link->index()].link = link;
        linksData[link->index()].parentIndex = link->parent() ? link->parent()->index() : -1;
    }
}


// initialize extra joints for making closed links
void SPCFSImpl::initExtraJoints(int bodyIndex)
{
    const DyBodyPtr& body = world.body(bodyIndex);
    BodyData& bodyData = bodiesData[bodyIndex];
    for(int j=0; j < body->numExtraJoints(); ++j){
        Body::ExtraJoint& bodyExtraJoint = body->extraJoint(j);
        ExtraJointLinkPairPtr linkPair;
        if(bodyExtraJoint.type == Body::EJ_PISTON){
            linkPair = boost::make_shared<ExtraJointLinkPair>();
            linkPair->isSameBodyPair = true;
            linkPair->isNonContactConstraint = true;
        
            // generate two vectors orthogonal to the joint axis
            Vector3 u = Vector3::Zero();
            int minElem = 0;
            Vector3& axis = bodyExtraJoint.axis;
            for(int k=1; k < 3; k++){
                if(fabs(axis(k)) < fabs(axis(minElem))){
                    minElem = k;
                }
            }
            u(minElem) = 1.0;
            linkPair->constraintPoints.resize(2);
            const Vector3 t1 = axis.cross(u).normalized();
            linkPair->jointConstraintAxes[0] = t1;
            linkPair->jointConstraintAxes[1] = axis.cross(t1).normalized();
            
        } else if(bodyExtraJoint.type == Body::EJ_BALL){
            
        }
        
        if(linkPair){
            int numConstraints = linkPair->constraintPoints.size();
            for(int k=0; k < numConstraints; ++k){
                ConstraintPoint& constraint = linkPair->constraintPoints[k];
                constraint.numFrictionVectors = 0;
                constraint.globalFrictionIndex = numeric_limits<int>::max();
            }
            for(int k=0; k < 2; ++k){
                linkPair->bodyIndex[k] = bodyIndex;
                linkPair->bodyData[k] = &bodiesData[bodyIndex];
                DyLink* link = static_cast<DyLink*>(bodyExtraJoint.link[k]);
                linkPair->link[k] = link;
                linkPair->linkData[k] = &(bodyData.linksData[link->index()]);
                linkPair->jointPoint[k] = bodyExtraJoint.point[k];
            }
            extraJointLinkPairs.push_back(linkPair);
        }
    }
}    


void SPCFSImpl::init2Dconstraint(int bodyIndex)
{
    if(!bodyFor2dConstraint){
        bodyFor2dConstraint = new DyBody();
        Link* link = bodyFor2dConstraint->createLink();
        bodyFor2dConstraint->setRootLink(link);
        link->p().setZero();
        link->R().setIdentity();
        initBody(bodyFor2dConstraint, bodyDataFor2dConstraint);
        LinkData& linkData0 = bodyDataFor2dConstraint.linksData[0];
        linkData0.dw.setZero();
        linkData0.dvo.setZero();
    }

    DyLink* rootLink = world.body(bodyIndex)->rootLink();

    Constrain2dLinkPairPtr linkPair = boost::make_shared<Constrain2dLinkPair>();
    linkPair->isSameBodyPair = false;
    linkPair->isNonContactConstraint = true;
    
    linkPair->constraintPoints.resize(3);
    for(int i=0; i < 3; ++i){
        ConstraintPoint& constraint = linkPair->constraintPoints[i];
        constraint.numFrictionVectors = 0;
        constraint.globalFrictionIndex = numeric_limits<int>::max();
        linkPair->globalYpositions[i] = (rootLink->R() * local2dConstraintPoints[i] + rootLink->p()).y();
    }
        
    linkPair->bodyIndex[0] = -1;
    linkPair->bodyData[0] = &bodyDataFor2dConstraint;
    linkPair->link[0] = bodyFor2dConstraint->rootLink();
    linkPair->linkData[0] = &bodyDataFor2dConstraint.linksData[0];
    
    linkPair->bodyIndex[1] = bodyIndex;
    linkPair->bodyData[1] = &bodiesData[bodyIndex];
    linkPair->link[1] = rootLink;
    linkPair->linkData[1] = &bodiesData[bodyIndex].linksData[0];
            
    constrain2dLinkPairs.push_back(linkPair);
}
        
    
void SPCFSImpl::initialize(void)
{
    if(CFS_DEBUG || CFS_MCP_DEBUG){
        static int ntest = 0;
        os.close();
        os.open((string("cfs-log-") + boost::lexical_cast<string>(ntest++) + ".log").c_str());
        //os << setprecision(50);
    }

    if(CFS_MCP_DEBUG){
        numGaussSeidelTotalCalls = 0;
        numGaussSeidelTotalLoops = 0;
        numGaussSeidelTotalLoopsMax = 0;
    }

    int numBodies = world.numBodies();

    bodiesData.resize(numBodies);

    if(!collisionDetector){
        collisionDetector = boost::make_shared<AISTCollisionDetector>();
    } else {
        collisionDetector->clearGeometries();
    }
    geometryIdToBodyIndexMap.clear();
    geometryPairToLinkPairMap.clear();

    extraJointLinkPairs.clear();
    constrain2dLinkPairs.clear();

    for(int bodyIndex=0; bodyIndex < numBodies; ++bodyIndex){

        const DyBodyPtr& body = world.body(bodyIndex);
        BodyData& bodyData = bodiesData[bodyIndex];

        initBody(body, bodyData);

        bodyData.forwardDynamicsCBM =
            dynamic_pointer_cast<SPForwardDynamicsCBM>(world.forwardDynamics(bodyIndex));

        if(bodyData.isStatic && !(bodyData.forwardDynamicsCBM)){
            LinkDataArray& linksData = bodyData.linksData;
            for(size_t j=0; j < linksData.size(); ++j){
                LinkData& linkData = linksData[j];
                linkData.dw.setZero();
                linkData.dvo.setZero();
            }
        }

        bodyData.geometryId = addBodyToCollisionDetector(*body, *collisionDetector, false);
        geometryIdToBodyIndexMap.resize(collisionDetector->numGeometries(), bodyIndex);

        initExtraJoints(bodyIndex);

        if(is2Dmode && !body->isStaticModel()){
            init2Dconstraint(bodyIndex);
        }
    }

    collisionDetector->makeReady();

    prevGlobalNumConstraintVectors = 0;
    prevGlobalNumFrictionVectors = 0;
    numUnconverged = 0;

    randomAngle.engine().seed();
}


inline void SPCFSImpl::clearExternalForces()
{
    for(size_t i=0; i < bodiesData.size(); ++i){
        BodyData& bodyData = bodiesData[i];
        if(bodyData.hasConstrainedLinks){
            bodyData.body->clearExternalForces();
        }
    }
}


void SPCFSImpl::solve()
{
    if(CFS_DEBUG){
        os << "Time: " << world.currentTime() << std::endl;
    }

    for(size_t i=0; i < bodiesData.size(); ++i){
        BodyData& data = bodiesData[i];
        data.hasConstrainedLinks = false;
        DyBodyPtr& body = data.body;
        const int n = body->numLinks();
        for(int j=0; j < n; ++j){
            DyLink* link = body->link(j);
            collisionDetector->updatePosition(data.geometryId + j, link->T());
            link->constraintForces().clear();
        }
    }

    globalNumConstraintVectors = 0;
    globalNumFrictionVectors = 0;
    areThereImpacts = false;
    constrainedLinkPairs.clear();

    setConstraintPoints();

    if(CFS_PUT_NUM_CONTACT_POINTS){
        cout << globalNumContactNormalVectors;
    }
    if(globalNumConstraintVectors > 0){

        if(CFS_DEBUG){
            os << "Num Collisions: " << globalNumContactNormalVectors << std::endl;
        }
        if(CFS_DEBUG_VERBOSE) putContactPoints();

        const bool constraintsSizeChanged = ((globalNumFrictionVectors   != prevGlobalNumFrictionVectors) ||
                                             (globalNumConstraintVectors != prevGlobalNumConstraintVectors));

        if(constraintsSizeChanged){
            initMatrices();
        }

        if(areThereImpacts){
            solveImpactConstraints();
        }

        if(SKIP_REDUNDANT_ACCEL_CALC){
            setAccelCalcSkipInformation();
        }

        setDefaultAccelerationVector();
        setAccelerationMatrix();

        clearSingularPointConstraintsOfClosedLoopConnections();
		
        setConstantVectorAndMuBlock();

        if(CFS_DEBUG_VERBOSE){
            debugPutVector(an0 , "an0");
            debugPutVector(at0 , "at0");
            debugPutMatrix(Mlcp, "Mlcp");
            debugPutVector(b.head(globalNumConstraintVectors), "b1");
            debugPutVector(b.segment(globalNumConstraintVectors, globalNumFrictionVectors), "b2");
        }

        bool isConverged = kik_callSiconosSolver(Mlcp, b, solution,contactIndexToMu);
        if(!isConverged)
        {
            ++numUnconverged;
            if(CFS_DEBUG) os << "LCP didn't converge" << numUnconverged << std::endl;
        }
        else
        {
            if(CFS_DEBUG)
                os << "LCP converged" << std::endl;
            if(CFS_DEBUG_LCPCHECK){
                // checkLCPResult(Mlcp, b, solution);
                // checkMCPResult(Mlcp, b, solution);
            }
            addConstraintForceToLinks();
        }
    }
    prevGlobalNumConstraintVectors = globalNumConstraintVectors;
    prevGlobalNumFrictionVectors   = globalNumFrictionVectors;
}


void SPCFSImpl::setConstraintPoints()
{
    collisionDetector->detectCollisions(boost::bind(&SPCFSImpl::extractConstraintPoints, this, _1));

    globalNumContactNormalVectors = globalNumConstraintVectors;

    for(size_t i=0; i < extraJointLinkPairs.size(); ++i){
        setExtraJointConstraintPoints(extraJointLinkPairs[i]);
    }

    for(size_t i=0; i < constrain2dLinkPairs.size(); ++i){
        set2dConstraintPoints(constrain2dLinkPairs[i]);
    }
}


void SPCFSImpl::extractConstraintPoints(const CollisionPair& collisionPair)
{
    LinkPair* pLinkPair;
    
    const IdPair<> idPair(collisionPair.geometryId);
    GeometryPairToLinkPairMap::iterator p = geometryPairToLinkPairMap.find(idPair);
    
    if(p != geometryPairToLinkPairMap.end()){
        pLinkPair = &p->second;
        pLinkPair->constraintPoints.clear();
    } else {
        LinkPair& linkPair = geometryPairToLinkPairMap.insert(make_pair(idPair, LinkPair())).first->second;
        
        for(int i=0; i < 2; ++i){
            const int id = collisionPair.geometryId[i];
            const int bodyIndex = geometryIdToBodyIndexMap[id];
            linkPair.bodyIndex[i] = bodyIndex;
            BodyData& bodyData = bodiesData[bodyIndex];
            linkPair.bodyData[i] = &bodyData;
            const int linkIndex = id - bodyData.geometryId;
            linkPair.link[i] = bodyData.body->link(linkIndex);
            linkPair.linkData[i] = &bodyData.linksData[linkIndex];
        }
        linkPair.isSameBodyPair = (linkPair.bodyIndex[0] == linkPair.bodyIndex[1]);
        linkPair.isNonContactConstraint = false;
        linkPair.muStatic = defaultStaticFriction;
        linkPair.muDynamic = defaultSlipFriction;
        linkPair.contactCullingDistance = defaultContactCullingDistance;
        linkPair.contactCullingDepth = defaultContactCullingDepth;
        linkPair.epsilon = defaultCoefficientOfRestitution;
        pLinkPair = &linkPair;
    }
    pLinkPair->bodyData[0]->hasConstrainedLinks = true;
    pLinkPair->bodyData[1]->hasConstrainedLinks = true;
    
    const vector<Collision>& collisions = collisionPair.collisions;
    for(size_t i=0; i < collisions.size(); ++i){
        setContactConstraintPoint(*pLinkPair, collisions[i]);
    }

    if(!pLinkPair->constraintPoints.empty()){
        constrainedLinkPairs.push_back(pLinkPair);
    }
}


/**
   @retuen true if the point is actually added to the constraints
*/
bool SPCFSImpl::setContactConstraintPoint(LinkPair& linkPair, const Collision& collision)
{
    // skip the contact which has too much depth
    if(collision.depth > linkPair.contactCullingDepth)  return false;
    
    ConstraintPointArray& constraintPoints = linkPair.constraintPoints;
    constraintPoints.push_back(ConstraintPoint());
    ConstraintPoint& contact = constraintPoints.back();

    contact.point = collision.point;

  // dense contact points are eliminated
    int nPrevPoints = constraintPoints.size() - 1;
    for(int i=0; i < nPrevPoints; ++i)
    {
        if((constraintPoints[i].point - contact.point).norm() < linkPair.contactCullingDistance)
        {
            constraintPoints.pop_back();
            return false;
        }
    }

    contact.normalTowardInside[1] = collision.normal;
    contact.normalTowardInside[0] = -contact.normalTowardInside[1];
    contact.depth = collision.depth;
    contact.globalIndex = globalNumConstraintVectors++;

    // check velocities
    Vector3 v[2];
    for(int k=0; k < 2; ++k){
        DyLink* link = linkPair.link[k];
        if(link->isRoot() && link->isFixedJoint()){
            v[k].setZero();
        } else {
            v[k] = link->vo() + link->w().cross(contact.point);
            if (link->jointType() == Link::CRAWLER_JOINT){
                // tentative
                // invalid depths should be fixed
              //  if (contact.depth > contactCorrectionDepth * 2.0){
              //      contact.depth = contactCorrectionDepth * 2.0;
              //  }
                Vector3 axis = link->R() * link->a();
                Vector3 n(collision.normal);
                Vector3 dir = axis.cross(n);
                if (k) dir *= -1.0;
                dir.normalize();
                v[k] += link->u() * dir;
            }
        }
    }
    contact.relVelocityOn0 = v[1] - v[0]; 

    contact.normalProjectionOfRelVelocityOn0 = contact.normalTowardInside[1].dot(contact.relVelocityOn0);

    if(!areThereImpacts){
        if(contact.normalProjectionOfRelVelocityOn0 < -1.0e-6){
            areThereImpacts = true;
        }
    }
    
    Vector3 v_tangent =
        contact.relVelocityOn0 - contact.normalProjectionOfRelVelocityOn0 * contact.normalTowardInside[1];
    
    contact.globalFrictionIndex = globalNumFrictionVectors;
    
    double vt_square = v_tangent.squaredNorm();
    static const double vsqrthresh = VEL_THRESH_OF_DYNAMIC_FRICTION * VEL_THRESH_OF_DYNAMIC_FRICTION;
    bool isSlipping = (vt_square > vsqrthresh);
    contact.mu = isSlipping ? linkPair.muDynamic : linkPair.muStatic;
    
    contact.numFrictionVectors = 2;
    setFrictionVectors(contact);
    globalNumFrictionVectors += contact.numFrictionVectors;

    return true;
}


void SPCFSImpl::setFrictionVectors(ConstraintPoint& contact)
{
    Vector3 u = Vector3::Zero();
    int minAxis = 0;
    Vector3& normal = contact.normalTowardInside[0];

    for(int i=1; i < 3; i++){
        if(fabs(normal(i)) < fabs(normal(minAxis))){
            minAxis = i;
        }
    }
    u(minAxis) = 1.0;

    Vector3 t1 = normal.cross(u).normalized();
    Vector3 t2 = normal.cross(t1).normalized();

    if(ENABLE_RANDOM_STATIC_FRICTION_BASE){
        double theta = randomAngle();
        contact.frictionVector[0][0] = cos(theta) * t1 + sin(theta) * t2;
        theta += PI_2;
        contact.frictionVector[1][0] = cos(theta) * t1 + sin(theta) * t2;
    } else {
        contact.frictionVector[0][0] = t1;
        contact.frictionVector[1][0] = t2;
    }
    contact.frictionVector[0][1] = -contact.frictionVector[0][0];
    contact.frictionVector[1][1] = -contact.frictionVector[1][0];
}


void SPCFSImpl::setExtraJointConstraintPoints(const ExtraJointLinkPairPtr& linkPair)
{
    ConstraintPointArray& constraintPoints = linkPair->constraintPoints;

    DyLink* link0 = linkPair->link[0];
    DyLink* link1 = linkPair->link[1];

    Vector3 point[2];
    point[0].noalias() = link0->p() + link0->R() * linkPair->jointPoint[0];
    point[1].noalias() = link1->p() + link1->R() * linkPair->jointPoint[1];
    Vector3 midPoint = (point[0] + point[1]) / 2.0;
    Vector3 error = midPoint - point[0];

    // check velocities
    Vector3 v[2];
    for(int k=0; k < 2; ++k){
        DyLink* link = linkPair->link[k];
        if(link->isRoot() && link->isFixedJoint()){
            v[k].setZero();
        } else {
            v[k] = link->vo() + link->w().cross(point[k]);
        }
    }
    Vector3 relVelocityOn0 = v[1] - v[0];

    int n = linkPair->constraintPoints.size();
    for(int i=0; i < n; ++i){
        ConstraintPoint& constraint = constraintPoints[i];
        const Vector3 axis = link0->R() * linkPair->jointConstraintAxes[i];
        constraint.point = midPoint;
        constraint.normalTowardInside[0] =  axis;
        constraint.normalTowardInside[1] = -axis;
        constraint.depth = axis.dot(error);
        constraint.globalIndex = globalNumConstraintVectors++;
        constraint.normalProjectionOfRelVelocityOn0 = constraint.normalTowardInside[1].dot(relVelocityOn0);
    }
    linkPair->bodyData[0]->hasConstrainedLinks = true;
    linkPair->bodyData[1]->hasConstrainedLinks = true;

    constrainedLinkPairs.push_back(linkPair.get());
}


void SPCFSImpl::set2dConstraintPoints(const Constrain2dLinkPairPtr& linkPair)
{
    static const Vector3 yAxis(0.0, 1.0, 0.0);
    
    ConstraintPointArray& constraintPoints = linkPair->constraintPoints;
    
    int n = linkPair->constraintPoints.size();
    for(int i=0; i < n; ++i){
        DyLink* link1 = linkPair->link[1];
        Vector3 point1 = link1->p() + link1->R() * local2dConstraintPoints[i];
        Vector3 relVelocityOn0 = link1->vo() + link1->w().cross(point1);
        ConstraintPoint& constraint = constraintPoints[i];
        constraint.point = point1;
        constraint.normalTowardInside[0] =  yAxis;
        constraint.normalTowardInside[1] = -yAxis;
        constraint.depth = point1.y() - linkPair->globalYpositions[i];
        constraint.globalIndex = globalNumConstraintVectors++;
        constraint.normalProjectionOfRelVelocityOn0 = -(link1->vo() + link1->w().cross(point1)).y();
    }
    linkPair->bodyData[0]->hasConstrainedLinks = true;
    linkPair->bodyData[1]->hasConstrainedLinks = true;

    constrainedLinkPairs.push_back(linkPair.get());
}



void SPCFSImpl::putContactPoints()
{
    os << "Contact Points\n";
    for(size_t i=0; i < constrainedLinkPairs.size(); ++i){
        LinkPair* linkPair = constrainedLinkPairs[i];
        ExtraJointLinkPair* ejLinkPair = dynamic_cast<ExtraJointLinkPair*>(linkPair);

        if(!ejLinkPair){

            os << " " << linkPair->link[0]->name() << " of " << linkPair->bodyData[0]->body->modelName();
            os << "<-->";
            os << " " << linkPair->link[1]->name() << " of " << linkPair->bodyData[1]->body->modelName();
            os << "\n";
            os << " culling thresh: " << linkPair->contactCullingDistance << "\n";

            ConstraintPointArray& constraintPoints = linkPair->constraintPoints;
            for(size_t j=0; j < constraintPoints.size(); ++j){
                ConstraintPoint& contact = constraintPoints[j];
                os << " index " << contact.globalIndex;
                os << " point: " << contact.point;
                os << " normal: " << contact.normalTowardInside[1];
                os << " defaultAccel[0]: " << contact.defaultAccel[0];
                os << " defaultAccel[1]: " << contact.defaultAccel[1];
                os << " normal projectionOfRelVelocityOn0" << contact.normalProjectionOfRelVelocityOn0;
                os << " depth" << contact.depth;
                os << " mu" << contact.mu;
                os << " rel velocity: " << contact.relVelocityOn0;
                os << " friction[0][0]: " << contact.frictionVector[0][0];
                os << " friction[0][1]: " << contact.frictionVector[0][1];
                os << " friction[1][0]: " << contact.frictionVector[1][0];
                os << " friction[1][1]: " << contact.frictionVector[1][1];
                os << "\n";
            }
        }
    }
    os << std::endl;
}


void SPCFSImpl::solveImpactConstraints()
{
    if(CFS_DEBUG){
        os << "Impacts !" << std::endl;
    }
}


void SPCFSImpl::initMatrices()
{
  const int n = globalNumConstraintVectors;
  const int m = globalNumFrictionVectors;

  const int dimLCP =  (n + m);

	int prv_sz = Mlcp.rows();
  Mlcp.resize(dimLCP, dimLCP);
  b   .resize(dimLCP);
  solution.resize(dimLCP);

  frictionIndexToContactIndex.resize(m);
  contactIndexToMu.resize(globalNumContactNormalVectors);
  mcpHi.resize(globalNumContactNormalVectors);

  an0.resize(n);
  at0.resize(m);
  
  kik_DeleteBuffer(prv_sz);
	kik_NewBuffer(Mlcp.rows());

}


void SPCFSImpl::setAccelCalcSkipInformation()
{
    // clear skip check numbers
    for(size_t i=0; i < bodiesData.size(); ++i){
        BodyData& bodyData = bodiesData[i];
        if(bodyData.hasConstrainedLinks){
            LinkDataArray& linksData = bodyData.linksData;
            for(size_t j=0; j < linksData.size(); ++j){
                linksData[j].numberToCheckAccelCalcSkip = numeric_limits<int>::max();
            }
        }
    }

    // add the number of contact points to skip check numbers of the links from a contact target to the root
    int numLinkPairs = constrainedLinkPairs.size();
    for(int i=0; i < numLinkPairs; ++i){
        LinkPair* linkPair = constrainedLinkPairs[i];
        int constraintIndex = linkPair->constraintPoints.front().globalIndex;
        for(int j=0; j < 2; ++j){
            LinkDataArray& linksData = linkPair->bodyData[j]->linksData;
            int linkIndex = linkPair->link[j]->index();
            while(linkIndex >= 0){
                LinkData& linkData = linksData[linkIndex];
                if(linkData.numberToCheckAccelCalcSkip < constraintIndex){
                    break;
                }
                linkData.numberToCheckAccelCalcSkip = constraintIndex;
                linkIndex = linkData.parentIndex;
            }
        }
    }
}


void SPCFSImpl::setDefaultAccelerationVector()
{
    // calculate accelerations with no constraint force
    for(size_t i=0; i < bodiesData.size(); ++i){
        BodyData& bodyData = bodiesData[i];
        if(bodyData.hasConstrainedLinks && ! bodyData.isStatic){

            if(bodyData.forwardDynamicsCBM){
                bodyData.forwardDynamicsCBM->sumExternalForces();
                bodyData.forwardDynamicsCBM->solveUnknownAccels();
                calcAccelsMM(bodyData, numeric_limits<int>::max());

            } else {
                initABMForceElementsWithNoExtForce(bodyData);
                calcAccelsABM(bodyData, numeric_limits<int>::max());
            }
        }
    }

    // extract accelerations
    for(size_t i=0; i < constrainedLinkPairs.size(); ++i){

        LinkPair& linkPair = *constrainedLinkPairs[i];
        ConstraintPointArray& constraintPoints = linkPair.constraintPoints;

        for(size_t j=0; j < constraintPoints.size(); ++j){
            ConstraintPoint& constraint = constraintPoints[j];

            for(int k=0; k < 2; ++k){
                if(linkPair.bodyData[k]->isStatic){
                    constraint.defaultAccel[k].setZero();
                } else {
                    DyLink* link = linkPair.link[k];
                    LinkData* linkData = linkPair.linkData[k];
                    constraint.defaultAccel[k] =
                        linkData->dvo - constraint.point.cross(linkData->dw) +
                        link->w().cross(link->vo() + link->w().cross(constraint.point));
                }
            }

            Vector3 relDefaultAccel(constraint.defaultAccel[1] - constraint.defaultAccel[0]);
            an0[constraint.globalIndex] = constraint.normalTowardInside[1].dot(relDefaultAccel);

            for(int k=0; k < constraint.numFrictionVectors; ++k){
                at0[constraint.globalFrictionIndex + k] = constraint.frictionVector[k][1].dot(relDefaultAccel);
            }
        }
    }
}


void SPCFSImpl::setAccelerationMatrix()
{
    const int n = globalNumConstraintVectors;
    const int m = globalNumFrictionVectors;

    Eigen::Block<MatrixX> Knn = Mlcp.block(0, 0, n, n);
    Eigen::Block<MatrixX> Ktn = Mlcp.block(0, n, n, m);
    Eigen::Block<MatrixX> Knt = Mlcp.block(n, 0, m, n);
    Eigen::Block<MatrixX> Ktt = Mlcp.block(n, n, m, m);

    for(size_t i=0; i < constrainedLinkPairs.size(); ++i){

        LinkPair& linkPair = *constrainedLinkPairs[i];
        int numConstraintsInPair = linkPair.constraintPoints.size();

        for(int j=0; j < numConstraintsInPair; ++j){

            ConstraintPoint& constraint = linkPair.constraintPoints[j];
            int constraintIndex = constraint.globalIndex;

            // apply test normal force
            for(int k=0; k < 2; ++k){
                BodyData& bodyData = *linkPair.bodyData[k];
                if(!bodyData.isStatic){

                    bodyData.isTestForceBeingApplied = true;
                    const Vector3& f = constraint.normalTowardInside[k];

                    if(bodyData.forwardDynamicsCBM){
                        //! \todo This code does not work correctly when the links are in the same body. Fix it.
                        Vector3 arm = constraint.point - bodyData.body->rootLink()->p();
                        Vector3 tau = arm.cross(f);
                        Vector3 tauext = constraint.point.cross(f);
                        bodyData.forwardDynamicsCBM->solveUnknownAccels(linkPair.link[k], f, tauext, f, tau);
                        calcAccelsMM(bodyData, constraintIndex);
                    } else {
                        Vector3 tau = constraint.point.cross(f);
                        calcABMForceElementsWithTestForce(bodyData, linkPair.link[k], f, tau);
                        if(!linkPair.isSameBodyPair || (k > 0)){
                            calcAccelsABM(bodyData, constraintIndex);
                        }
                    }
                }
            }
            extractRelAccelsOfConstraintPoints(Knn, Knt, constraintIndex, constraintIndex);

            // apply test friction force
            for(int l=0; l < constraint.numFrictionVectors; ++l){
                for(int k=0; k < 2; ++k){
                    BodyData& bodyData = *linkPair.bodyData[k];
                    if(!bodyData.isStatic){
                        const Vector3& f = constraint.frictionVector[l][k];

                        if(bodyData.forwardDynamicsCBM){
                            //! \todo This code does not work correctly when the links are in the same body. Fix it.
                            Vector3 arm = constraint.point - bodyData.body->rootLink()->p();
                            Vector3 tau = arm.cross(f);
                            Vector3 tauext = constraint.point.cross(f);
                            bodyData.forwardDynamicsCBM->solveUnknownAccels(linkPair.link[k], f, tauext, f, tau);
                            calcAccelsMM(bodyData, constraintIndex);
                        } else {
                            Vector3 tau = constraint.point.cross(f);
                            calcABMForceElementsWithTestForce(bodyData, linkPair.link[k], f, tau);
                            if(!linkPair.isSameBodyPair || (k > 0)){
                                calcAccelsABM(bodyData, constraintIndex);
                            }
                        }
                    }
                }
                extractRelAccelsOfConstraintPoints(Ktn, Ktt, constraint.globalFrictionIndex + l, constraintIndex);
            }

            linkPair.bodyData[0]->isTestForceBeingApplied = false;
            linkPair.bodyData[1]->isTestForceBeingApplied = false;
        }
    }

    if(ASSUME_SYMMETRIC_MATRIX){
        copySymmetricElementsOfAccelerationMatrix(Knn, Ktn, Knt, Ktt);
    }
}


void SPCFSImpl::initABMForceElementsWithNoExtForce(BodyData& bodyData)
{
    bodyData.dpf.setZero();
    bodyData.dptau.setZero();

    std::vector<LinkData>& linksData = bodyData.linksData;
    const LinkTraverse& traverse = bodyData.body->linkTraverse();
    const int n = traverse.numLinks();

    for(int i = n-1; i >= 0; --i){
        DyLink* link = static_cast<DyLink*>(traverse[i]);
        LinkData& data = linksData[i];

        /*
          data.pf0   = link->pf;
          data.ptau0 = link->ptau;
        */
        data.pf0   = link->pf() - link->f_ext();
        data.ptau0 = link->ptau() - link->tau_ext();

        for(DyLink* child = link->child(); child; child = child->sibling()){

            LinkData& childData = linksData[child->index()];

            data.pf0   += childData.pf0;
            data.ptau0 += childData.ptau0;

            if(!child->isFixedJoint()){
                double uu_dd = childData.uu0 / child->dd();
                data.pf0   += uu_dd * child->hhv();
                data.ptau0 += uu_dd * child->hhw();
            }
        }

        if(i > 0){
            if(!link->isFixedJoint()){
                data.uu0  = link->uu() + link->u() - (link->sv().dot(data.pf0) + link->sw().dot(data.ptau0));
                data.uu = data.uu0;
            }
        }
    }
}


void SPCFSImpl::calcABMForceElementsWithTestForce
(BodyData& bodyData, DyLink* linkToApplyForce, const Vector3& f, const Vector3& tau)
{
    std::vector<LinkData>& linksData = bodyData.linksData;

    Vector3 dpf   = -f;
    Vector3 dptau = -tau;

    DyLink* link = linkToApplyForce;
    while(link->parent()){
        if(!link->isFixedJoint()){
            LinkData& data = linksData[link->index()];
            double duu = -(link->sv().dot(dpf) + link->sw().dot(dptau));
            data.uu += duu;
            double duudd = duu / link->dd();
            dpf   += duudd * link->hhv();
            dptau += duudd * link->hhw();
        }
        link = link->parent();
    }

    bodyData.dpf   += dpf;
    bodyData.dptau += dptau;
}


void SPCFSImpl::calcAccelsABM(BodyData& bodyData, int constraintIndex)
{
    std::vector<LinkData>& linksData = bodyData.linksData;
    LinkData& rootData = linksData[0];
    DyLink* rootLink = rootData.link;

    if(rootLink->isFreeJoint()){

        Eigen::Matrix<double, 6, 6> M;
        M << rootLink->Ivv(), rootLink->Iwv().transpose(),
            rootLink->Iwv(), rootLink->Iww();

        Eigen::Matrix<double, 6, 1> f;
        f << (rootData.pf0   + bodyData.dpf),
            (rootData.ptau0 + bodyData.dptau);
        f *= -1.0;

        Eigen::Matrix<double, 6, 1> a(M.colPivHouseholderQr().solve(f));

        rootData.dvo = a.head<3>();
        rootData.dw  = a.tail<3>();

    } else {
        rootData.dw .setZero();
        rootData.dvo.setZero();
    }

    // reset
    bodyData.dpf  .setZero();
    bodyData.dptau.setZero();

    int skipCheckNumber = ASSUME_SYMMETRIC_MATRIX ? constraintIndex : (numeric_limits<int>::max() - 1);
    int n = linksData.size();
    for(int linkIndex = 1; linkIndex < n; ++linkIndex){

        LinkData& linkData = linksData[linkIndex];

        if(!SKIP_REDUNDANT_ACCEL_CALC || linkData.numberToCheckAccelCalcSkip <= skipCheckNumber){

            DyLink* link = linkData.link;
            LinkData& parentData = linksData[linkData.parentIndex];

            if(!link->isFixedJoint()){
                linkData.ddq = (linkData.uu - (link->hhv().dot(parentData.dvo) + link->hhw().dot(parentData.dw))) / link->dd();
                linkData.dvo = parentData.dvo + link->cv() + link->sv() * linkData.ddq;
                linkData.dw  = parentData.dw  + link->cw() + link->sw() * linkData.ddq;
            }else{
                linkData.ddq = 0.0;
                linkData.dvo = parentData.dvo;
                linkData.dw  = parentData.dw;
            }

            // reset
            linkData.uu = linkData.uu0;
        }
    }
}


void SPCFSImpl::calcAccelsMM(BodyData& bodyData, int constraintIndex)
{
    std::vector<LinkData>& linksData = bodyData.linksData;

    LinkData& rootData = linksData[0];
    DyLink* rootLink = rootData.link;
    rootData.dvo = rootLink->dvo();
    rootData.dw  = rootLink->dw();

    const int skipCheckNumber = ASSUME_SYMMETRIC_MATRIX ? constraintIndex : (numeric_limits<int>::max() - 1);
    const int n = linksData.size();

    for(int linkIndex = 1; linkIndex < n; ++linkIndex){

        LinkData& linkData = linksData[linkIndex];

        if(!SKIP_REDUNDANT_ACCEL_CALC || linkData.numberToCheckAccelCalcSkip <= skipCheckNumber){

            DyLink* link = linkData.link;
            LinkData& parentData = linksData[linkData.parentIndex];
            if(!link->isFixedJoint()){
                linkData.dvo = parentData.dvo + link->cv() + link->ddq() * link->sv();
                linkData.dw  = parentData.dw  + link->cw() + link->ddq() * link->sw();
            }else{
                linkData.dvo = parentData.dvo;
                linkData.dw  = parentData.dw;
            }
        }
    }
}


void SPCFSImpl::extractRelAccelsOfConstraintPoints
(Eigen::Block<MatrixX>& Kxn, Eigen::Block<MatrixX>& Kxt, int testForceIndex, int constraintIndex)
{
    int maxConstraintIndexToExtract = ASSUME_SYMMETRIC_MATRIX ? constraintIndex : globalNumConstraintVectors;


    for(size_t i=0; i < constrainedLinkPairs.size(); ++i){

        LinkPair& linkPair = *constrainedLinkPairs[i];

        BodyData& bodyData0 = *linkPair.bodyData[0];
        BodyData& bodyData1 = *linkPair.bodyData[1];

        if(bodyData0.isTestForceBeingApplied){
            if(bodyData1.isTestForceBeingApplied){
                extractRelAccelsFromLinkPairCase1(Kxn, Kxt, linkPair, testForceIndex, maxConstraintIndexToExtract);
            } else {
                extractRelAccelsFromLinkPairCase2(Kxn, Kxt, linkPair, 0, 1, testForceIndex, maxConstraintIndexToExtract);
            }
        } else {
            if(bodyData1.isTestForceBeingApplied){
                extractRelAccelsFromLinkPairCase2(Kxn, Kxt, linkPair, 1, 0, testForceIndex, maxConstraintIndexToExtract);
            } else {
                extractRelAccelsFromLinkPairCase3(Kxn, Kxt, linkPair, testForceIndex, maxConstraintIndexToExtract);
            }
        }
    }
}


void SPCFSImpl::extractRelAccelsFromLinkPairCase1
(Eigen::Block<MatrixX>& Kxn, Eigen::Block<MatrixX>& Kxt,
 LinkPair& linkPair, int testForceIndex, int maxConstraintIndexToExtract)
{
    ConstraintPointArray& constraintPoints = linkPair.constraintPoints;

    for(size_t i=0; i < constraintPoints.size(); ++i){

        ConstraintPoint& constraint = constraintPoints[i];
        int constraintIndex = constraint.globalIndex;

        if(ASSUME_SYMMETRIC_MATRIX && constraintIndex > maxConstraintIndexToExtract){
            break;
        }

        DyLink* link0 = linkPair.link[0];
        DyLink* link1 = linkPair.link[1];
        LinkData* linkData0 = linkPair.linkData[0];
        LinkData* linkData1 = linkPair.linkData[1];

        //! \todo Can the follwoing equations be simplified ?
        Vector3 dv0 =
            linkData0->dvo - constraint.point.cross(linkData0->dw) +
            link0->w().cross(link0->vo() + link0->w().cross(constraint.point));

        Vector3 dv1 =
            linkData1->dvo - constraint.point.cross(linkData1->dw) +
            link1->w().cross(link1->vo() + link1->w().cross(constraint.point));

        Vector3 relAccel = dv1 - dv0;

        Kxn(constraintIndex, testForceIndex) = constraint.normalTowardInside[1].dot(relAccel) - an0(constraintIndex);

        for(int j=0; j < constraint.numFrictionVectors; ++j){
            const int index = constraint.globalFrictionIndex + j;
            Kxt(index, testForceIndex) = constraint.frictionVector[j][1].dot(relAccel) - at0(index);
        }
    }
}


void SPCFSImpl::extractRelAccelsFromLinkPairCase2
(Eigen::Block<MatrixX>& Kxn, Eigen::Block<MatrixX>& Kxt,
 LinkPair& linkPair, int iTestForce, int iDefault, int testForceIndex, int maxConstraintIndexToExtract)
{
    ConstraintPointArray& constraintPoints = linkPair.constraintPoints;

    for(size_t i=0; i < constraintPoints.size(); ++i){

        ConstraintPoint& constraint = constraintPoints[i];
        int constraintIndex = constraint.globalIndex;

        if(ASSUME_SYMMETRIC_MATRIX && constraintIndex > maxConstraintIndexToExtract){
            break;
        }

        DyLink* link = linkPair.link[iTestForce];
        LinkData* linkData = linkPair.linkData[iTestForce];

        Vector3 dv(linkData->dvo - constraint.point.cross(linkData->dw) + link->w().cross(link->vo() + link->w().cross(constraint.point)));

        if(CFS_DEBUG_VERBOSE_2){
            os << "dv " << constraintIndex << " = " << dv << "\n";
        }

        Vector3 relAccel = constraint.defaultAccel[iDefault] - dv;

        Kxn(constraintIndex, testForceIndex) = constraint.normalTowardInside[iDefault].dot(relAccel) - an0(constraintIndex);

        for(int j=0; j < constraint.numFrictionVectors; ++j){
            const int index = constraint.globalFrictionIndex + j;
            Kxt(index, testForceIndex) = constraint.frictionVector[j][iDefault].dot(relAccel) - at0(index);
        }

    }
}


void SPCFSImpl::extractRelAccelsFromLinkPairCase3
(Eigen::Block<MatrixX>& Kxn, Eigen::Block<MatrixX>& Kxt, LinkPair& linkPair, int testForceIndex, int maxConstraintIndexToExtract)
{
    ConstraintPointArray& constraintPoints = linkPair.constraintPoints;

    for(size_t i=0; i < constraintPoints.size(); ++i){

        ConstraintPoint& constraint = constraintPoints[i];
        int constraintIndex = constraint.globalIndex;

        if(ASSUME_SYMMETRIC_MATRIX && constraintIndex > maxConstraintIndexToExtract){
            break;
        }

        Kxn(constraintIndex, testForceIndex) = 0.0;

        for(int j=0; j < constraint.numFrictionVectors; ++j){
            Kxt(constraint.globalFrictionIndex + j, testForceIndex) = 0.0;
        }
    }
}


void SPCFSImpl::copySymmetricElementsOfAccelerationMatrix
(Eigen::Block<MatrixX>& Knn, Eigen::Block<MatrixX>& Ktn, Eigen::Block<MatrixX>& Knt, Eigen::Block<MatrixX>& Ktt)
{
    for(size_t linkPairIndex=0; linkPairIndex < constrainedLinkPairs.size(); ++linkPairIndex){

        ConstraintPointArray& constraintPoints = constrainedLinkPairs[linkPairIndex]->constraintPoints;

        for(size_t localConstraintIndex = 0; localConstraintIndex < constraintPoints.size(); ++localConstraintIndex){

            ConstraintPoint& constraint = constraintPoints[localConstraintIndex];

            int constraintIndex = constraint.globalIndex;
            int nextConstraintIndex = constraintIndex + 1;
            for(int i = nextConstraintIndex; i < globalNumConstraintVectors; ++i){
                Knn(i, constraintIndex) = Knn(constraintIndex, i);
            }
            int frictionTopOfNextConstraint = constraint.globalFrictionIndex + constraint.numFrictionVectors;
            for(int i = frictionTopOfNextConstraint; i < globalNumFrictionVectors; ++i){
                Knt(i, constraintIndex) = Ktn(constraintIndex, i);
            }

            for(int localFrictionIndex=0; localFrictionIndex < constraint.numFrictionVectors; ++localFrictionIndex){

                int frictionIndex = constraint.globalFrictionIndex + localFrictionIndex;

                for(int i = nextConstraintIndex; i < globalNumConstraintVectors; ++i){
                    Ktn(i, frictionIndex) = Knt(frictionIndex, i);
                }
                for(int i = frictionTopOfNextConstraint; i < globalNumFrictionVectors; ++i){
                    Ktt(i, frictionIndex) = Ktt(frictionIndex, i);
                }
            }
        }
    }
}


void SPCFSImpl::clearSingularPointConstraintsOfClosedLoopConnections()
{
    for(int i = 0; i < Mlcp.rows(); ++i){
        if(Mlcp(i, i) < 1.0e-4){
            for(int j=0; j < Mlcp.rows(); ++j){
                Mlcp(j, i) = 0.0;
            }
            Mlcp(i, i) = numeric_limits<double>::max();
        }
    }
}


void SPCFSImpl::setConstantVectorAndMuBlock()
{
    double dtinv = 1.0 / world.timeStep();
    const int block2 = globalNumConstraintVectors;
    const int block3 = globalNumConstraintVectors + globalNumFrictionVectors;

    for(size_t i=0; i < constrainedLinkPairs.size(); ++i){

        LinkPair& linkPair = *constrainedLinkPairs[i];
        int numConstraintsInPair = linkPair.constraintPoints.size();

        for(int j=0; j < numConstraintsInPair; ++j)
        {
            ConstraintPoint& constraint = linkPair.constraintPoints[j];
            int globalIndex = constraint.globalIndex;

            // set constant vector of LCP
            // constraints for normal acceleration

            if(linkPair.isNonContactConstraint)
            {
                // connection constraint
                const double& error = constraint.depth;
                double v;
                if(error >= 0){v = 0.1 * (-1.0 + exp(-error * 20.0));}
                else          {v = 0.1 * ( 1.0 - exp( error * 20.0));}
                b     (globalIndex) = an0(globalIndex) + (constraint.normalProjectionOfRelVelocityOn0 + v) * dtinv;
            }
            else
            {
							// contact constraint
							if(ENABLE_CONTACT_DEPTH_CORRECTION)
							{
								double velOffset;
								const double depth = constraint.depth - contactCorrectionDepth;
								velOffset = contactCorrectionVelocityRatio * min(0.1,depth);
							  //  if(depth <= 0.0){
							  //      velOffset = contactCorrectionVelocityRatio * depth;
							  //  } else {
							  //      velOffset = contactCorrectionVelocityRatio * (-1.0 / (depth + 1.0) + 1.0);
							  //  }
								b(globalIndex) = an0(globalIndex) + (constraint.normalProjectionOfRelVelocityOn0 - velOffset) * dtinv;
							} 
							else
							{
							    b    (globalIndex) = an0(globalIndex) + constraint.normalProjectionOfRelVelocityOn0 * dtinv;/*** SHIORI of KIKUUWE **/
							}
							contactIndexToMu[globalIndex] = constraint.mu;

                int globalFrictionIndex = constraint.globalFrictionIndex;
                for(int k=0; k < constraint.numFrictionVectors; ++k)
                {
                    // constraints for tangent acceleration
              	    double tangentProjectionOfRelVelocity = constraint.frictionVector[k][1].dot(constraint.relVelocityOn0);
                    b    (block2 + globalFrictionIndex) = at0(globalFrictionIndex);
                    b    (block2 + globalFrictionIndex) += tangentProjectionOfRelVelocity * dtinv;
                    frictionIndexToContactIndex[globalFrictionIndex] = globalIndex;
                    ++globalFrictionIndex;
                }
            }
        }
    }
}


void SPCFSImpl::addConstraintForceToLinks()
{
    int n = constrainedLinkPairs.size();
    for(int i=0; i < n; ++i){
        LinkPair* linkPair = constrainedLinkPairs[i];
        for(int j=0; j < 2; ++j){
            // if(!linkPair->link[j]->isRoot() || linkPair->link[j]->jointType != Link::FIXED_JOINT){
            addConstraintForceToLink(linkPair, j);
            // }
        }
    }
}


void SPCFSImpl::addConstraintForceToLink(LinkPair* linkPair, int ipair)
{
    Vector3 f_total   = Vector3::Zero();
    Vector3 tau_total = Vector3::Zero();

    ConstraintPointArray& constraintPoints = linkPair->constraintPoints;
    int numConstraintPoints = constraintPoints.size();
    DyLink* link = linkPair->link[ipair];

    for(int i=0; i < numConstraintPoints; ++i){

        ConstraintPoint& constraint = constraintPoints[i];
        int globalIndex = constraint.globalIndex;

        Vector3 f = solution(globalIndex) * constraint.normalTowardInside[ipair];

        for(int j=0; j < constraint.numFrictionVectors; ++j){
            f += solution(globalNumConstraintVectors + constraint.globalFrictionIndex + j) * constraint.frictionVector[j][ipair];
        }

        f_total   += f;
        tau_total += constraint.point.cross(f);

        if(isConstraintForceOutputMode){
            link->constraintForces().push_back(DyLink::ConstraintForce(constraint.point, f));
        }
    }
    
    link->f_ext()   += f_total;
    link->tau_ext() += tau_total;


    if(CFS_DEBUG){
        os << "Constraint force to " << link->name() << ": f = " << f_total << ", tau = " << tau_total << std::endl;
    }
}




SPConstraintForceSolver::SPConstraintForceSolver(SPWorldBase& world)
{
    impl = new SPCFSImpl(world);
}


SPConstraintForceSolver::~SPConstraintForceSolver()
{
    delete impl;
}


void SPConstraintForceSolver::setCollisionDetector(CollisionDetectorPtr detector)
{
    impl->collisionDetector = detector;
}


CollisionDetectorPtr SPConstraintForceSolver::collisionDetector()
{
    return impl->collisionDetector;
}


void SPConstraintForceSolver::setFriction(double staticFriction, double slipFliction)
{
    impl->defaultStaticFriction = staticFriction;
    impl->defaultSlipFriction = slipFliction;
}


double SPConstraintForceSolver::staticFriction() const
{
    return impl->defaultStaticFriction;
}


double SPConstraintForceSolver::slipFriction() const
{
    return impl->defaultSlipFriction;
}


void SPConstraintForceSolver::setContactCullingDistance(double distance)
{
    impl->defaultContactCullingDistance = distance;
}


double SPConstraintForceSolver::contactCullingDistance() const
{
    return impl->defaultContactCullingDistance;
}


void SPConstraintForceSolver::setContactCullingDepth(double depth)
{
    impl->defaultContactCullingDepth = depth;
}


double SPConstraintForceSolver::contactCullingDepth()
{
    return impl->defaultContactCullingDepth;
}


void SPConstraintForceSolver::setCoefficientOfRestitution(double epsilon)
{
    impl->defaultCoefficientOfRestitution = epsilon;
}


double SPConstraintForceSolver::coefficientOfRestitution() const
{
    return impl->defaultCoefficientOfRestitution;
}


void SPConstraintForceSolver::setGaussSeidelErrorCriterion(double e)
{
    impl->gaussSeidelErrorCriterion = e;
}


double SPConstraintForceSolver::gaussSeidelErrorCriterion()
{
    return impl->gaussSeidelErrorCriterion;
}


void SPConstraintForceSolver::setGaussSeidelMaxNumIterations(int n)
{
    impl->maxNumGaussSeidelIteration = n;
}


int SPConstraintForceSolver::gaussSeidelMaxNumIterations()
{
    return impl->maxNumGaussSeidelIteration;
}


void SPConstraintForceSolver::setContactDepthCorrection(double depth, double velocityRatio)
{
    impl->contactCorrectionDepth = depth;
    impl->contactCorrectionVelocityRatio = velocityRatio;
}


double SPConstraintForceSolver::contactCorrectionDepth()
{
    return impl->contactCorrectionDepth;
}


double SPConstraintForceSolver::contactCorrectionVelocityRatio()
{
    return impl->contactCorrectionVelocityRatio;
}


void SPConstraintForceSolver::enableConstraintForceOutput(bool on)
{
    impl->isConstraintForceOutputMode = on;
}


void SPConstraintForceSolver::set2Dmode(bool on)
{
    impl->is2Dmode = on;
}


void SPConstraintForceSolver::initialize(void)
{
    impl->initialize();
}


void SPConstraintForceSolver::solve()
{
    impl->solve();
}


void SPConstraintForceSolver::clearExternalForces()
{
    impl->clearExternalForces();
}



void SPCFSImpl::kik_NewBuffer(int NC)
{
	if(NC>0)
	{
		int NC3= 3*NC;
		
		if( KIK_USE_FULL_MATRIX )
		{
			kik_prob.q          = new double  [NC3 + NC + NC3 + NC3 + NC3 * NC3];
			kik_prob.mu         = &(kik_prob.q[NC3                  ]);
			kik_reaction        = &(kik_prob.q[NC3 + NC             ]);
			kik_velocity        = &(kik_prob.q[NC3 + NC + NC3       ]);
			kik_prob.M->matrix0 = &(kik_prob.q[NC3 + NC + NC3 + NC3 ]);
		}
		else
		{
			kik_prob.q                       = new double  [NC3 + NC + NC3 + NC3  + NC3 * NC3];
			kik_prob.mu                      = &(kik_prob.q[NC3                  ]);
			kik_reaction                     = &(kik_prob.q[NC3 + NC             ]);
			kik_velocity                     = &(kik_prob.q[NC3 + NC + NC3       ]);
			kik_prob.M->matrix1->block[0]    = &(kik_prob.q[NC3 + NC + NC3 + NC3 ]);
			kik_prob.M->matrix1->block       = new double*[NC*NC  ];
			kik_prob.M->matrix1->index1_data = new size_t [NC+1+NC*NC];
			kik_prob.M->matrix1->blocksize0  = new unsigned int[NC]; 
		}
	}
}
void SPCFSImpl::kik_DeleteBuffer(int prv_sz)
{
	if(prv_sz>0)
	{
		if( KIK_USE_FULL_MATRIX )
		{
			delete [] kik_prob.q          ;
		}
		else
		{
			delete [] kik_prob.q          ;
			delete [] kik_prob.M->matrix1->block;
			delete [] kik_prob.M->matrix1->blocksize0;
			delete [] kik_prob.M->matrix1->index1_data;
	//		freeNumericMatrix(kik_prob.M);
		}
	}
}

int kik_g_check_zero_block(MatrixX& Mlcp, int NC, int ia,int ja)
{
	for(int i = 0; i<3; i++)for(int j = 0; j<3; j++)
	{
		if(fabs(Mlcp(((i==0)?(ia):(2*ia+i+NC-1)),((j==0)?(ja):(2*ja+j+NC-1))) ) >1e-10)
		{
			return 1;
		}
	}
	return 0;
}

int kik_g_construct_sparsity_matrix(int * ibuf, MatrixX& Mlcp, int NC)
{
	int count=0;
	for(int ia = 0; ia<NC; ia++)for(int ja = ia; ja<NC; ja++) 
	{
		if(ia==ja){ibuf[NC*ia+ja]=1; count ++;}
		else
		{
			int tmp = kik_g_check_zero_block(Mlcp,NC,ia,ja)  ;
			ibuf[NC*ia+ja]= ibuf[NC*ja+ia]= tmp;
			if(tmp==1) count += 2;
		}
	}
	return count;
}


void kik_g_check_offdiag(const int * ibuf, int NC, int* pia2, int *pja2)
{
	(*pia2)=(*pja2)=NC;
	for(int k = 0; k<NC; k++)
	{
		                            if(ibuf[NC* k + (NC-k-1) ]==1){(*pia2)=(*pja2)=k; return;}
		for(int ja=NC-k;ja<NC;ja++) if(ibuf[NC* k + ja       ]==1){(*pia2)=(*pja2)=k; return;}
		for(int ia=k-1 ;ia>=0;ia--) if(ibuf[NC*ia + (NC-k-1) ]==1){(*pia2)=(*pja2)=k; return;}
	}
}

void kik_g_sparsify_A(SparseBlockStructuredMatrix* pmat, MatrixX& Mlcp      , int NC, ofstream* pos);

bool SPCFSImpl::kik_callSiconosSolver(MatrixX& Mlcp, VectorX& b, VectorX& solution, VectorX& contactIndexToMu)
{
	int NC3 = Mlcp.rows();
	int NC = NC3/3;
	if(CFS_DEBUG)
	{
		if(NC3%3 != 0           ){ os << "   warning-1 " << std::endl;return false;}
		if(       b.rows()!= NC3){ os << "   warning-2 " << std::endl;return false;}
		if(solution.rows()!= NC3){ os << "   warning-3 " << std::endl;return false;}
	} 
	for(int ia=0;ia<NC;ia++)for(int i=0;i<3;i ++)kik_prob.q [3*ia+i]= b(((i==0)?(ia):(2*ia+i+NC-1)));
	for(int ia=0;ia<NC;ia++)                     kik_prob.mu[  ia  ]= contactIndexToMu[ia];
	kik_prob.numberOfContacts = NC;

	if( KIK_USE_FULL_MATRIX )
	{
		kik_prob.M->storageType = 0;
		kik_prob.M->size0       = NC3;
		kik_prob.M->size1       = NC3;
		for(int ia=0;ia<NC;ia++)for(int i =0;i <3 ;i ++)
		{
			for(int ja=0;ja<NC;ja++)for(int j =0;j <3;j ++) 
			{
				kik_prob.M->matrix0[NC3*(3*ia+i)+(3*ja+j)]=Mlcp(((i==0)?(ia):(2*ia+i+NC-1)),((j==0)?(ja):(2*ja+j+NC-1)));
			}
		}
	}
	else
	{
		kik_prob.M->storageType = 1;
		kik_prob.M->size0       = NC3;
		kik_prob.M->size1       = NC3;
		kik_g_sparsify_A( kik_prob.M->matrix1 , Mlcp , NC , &os);
	}
	
	frictionContact3D_driver(&kik_prob,kik_reaction,kik_velocity,kik_solops, &kik_numops);
	
	for(int ia=0;ia<NC;ia++)for(int i=0;i<3;i++) solution(((i==0)?(ia):(2*ia+i+NC-1))) = kik_reaction[3*ia+i] ;
  if(CFS_DEBUG_VERBOSE)
  {
		os << "=---------------------------------="<< std::endl; 
    os << "| res_error =" << kik_solops->dparam[1] <<  std::endl;
		os << "=---------------------------------="<< std::endl; 
  }
	return true;
}

void kik_g_copy_block(double * bbuf, MatrixX& Mlcp, int NC, int ia,int ja)
{
	for(int i=0;i<3;i++)for(int j=0;j<3;j++) bbuf[3*j+i]= Mlcp(((i==0)?(ia):(2*ia+i+NC-1)),((j==0)?(ja):(2*ja+j+NC-1))) ;
}

void kik_g_sparsify_A(SparseBlockStructuredMatrix* pmat, MatrixX& Mlcp, int NC, ofstream* pos )
{
	SparseBlockStructuredMatrix& mat = *pmat;
	mat.index2_data = mat.index1_data + (NC+1); 
	int NB=0;
	mat.index1_data[0]=0 ;
	for(int ia=0;ia<NC;ia++)
	{
		mat.index1_data[ia+1]=mat.index1_data[ia];
		for(int ja=0;ja<NC;ja++)
		{
			if(kik_g_check_zero_block(Mlcp,NC,ia,ja)==1)
			{
				mat.block[NB] = mat.block[0]+(9*NB);
				kik_g_copy_block(mat.block[NB], Mlcp, NC,ia,ja);
				mat.index1_data[ia+1] ++ ;
				mat.index2_data[NB] = ja;
				NB++;
			}
		}
	}
	mat.nbblocks     = NB;
	mat.blocknumber0 = NC;
	mat.blocknumber1 = NC;
	for(int i=0;i<NC;i++)mat.blocksize0[i]=(i+1)*3;
	mat.blocksize1 = mat.blocksize0;
	mat.filled1 = NC+1;
	mat.filled2 = NB  ;
  if(CFS_DEBUG)
  {
		int* ibuf = new int[NC*NC];
		kik_g_construct_sparsity_matrix(ibuf,Mlcp,NC);
		(*pos) << "=---------------------------------="<< std::endl; 
		(*pos) << "const int NC=" << NC << ";"<< std::endl; 
		(*pos) << "const int NB=" << NB << ";"<< std::endl; 
		for(int i=0;i<NC;i++){for(int j=0;j<NC;j++){(*pos) << ibuf[NC*i+j];} (*pos)<<std::endl;}
		(*pos) << "=---------------------------------="<< std::endl; 
		delete [] ibuf;
	}

}

