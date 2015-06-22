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
#include <cnoid/EigenUtil>
#include <fstream>
#include <iomanip>


#include "SPCore.h"


using namespace cnoid;

SPCore::SPCore(int maxNumGaussSeidelIteration, double gaussSeidelErrorCriterion)
{
    USE_FULL_MATRIX = false;
    kik_prob      = new FrictionContactProblem;
    kik_numops    = new NumericsOptions       ;
    kik_solops    = new SolverOptions         ;
    kik_prob->dimension     = 3;
    kik_prob->M             = new NumericsMatrix;
    kik_prob->M->matrix1    = new SparseBlockStructuredMatrix;
    kik_prob->q             = 0;
    kik_numops->verboseMode = 0;  /*************/
//  SICONOS_FRICTION_3D_projectionOnCylinder;
    frictionContact3D_setDefaultSolverOptions(kik_solops, SICONOS_FRICTION_3D_NSGS   );
    kik_solops->iparam[0] = maxNumGaussSeidelIteration;
//  kik_solops->iparam[1] = 1; // This makes computation slower
    kik_solops->dparam[0] = gaussSeidelErrorCriterion;//1e-10;
    kik_solops->internalSolvers->iparam[0] = maxNumGaussSeidelIteration;
    kik_solops->internalSolvers->dparam[0] = gaussSeidelErrorCriterion;
}

SPCore::~SPCore()
{
    kik_DeleteBuffer();
    deleteSolverOptions(kik_solops);
    delete  kik_solops ;
    delete  kik_prob->M ->matrix1;
    delete  kik_prob->M     ;
    delete  kik_prob     ;
    delete  kik_numops   ;
}

void SPCore::kik_NewBuffer(int NC)
{
  if(NC<=0){return;}
  int NC3= 3*NC;
  
  if( USE_FULL_MATRIX )
  {
    kik_prob->q          = new double   [NC3 + NC + NC3 + NC3 + NC3 * NC3];
    kik_prob->mu         = &(kik_prob->q[NC3                  ]);
    kik_reaction         = &(kik_prob->q[NC3 + NC             ]);
    kik_velocity         = &(kik_prob->q[NC3 + NC + NC3       ]);
    kik_prob->M->matrix0 = &(kik_prob->q[NC3 + NC + NC3 + NC3 ]);
    for(int i=0;i<NC3 + NC + NC3 + NC3 + NC3 * NC3;i++) kik_prob->q[i]=0;
  }
  else
  {
    kik_prob->M->matrix1->block      = new double*[NC*NC  ];
    kik_prob->M->matrix1->index1_data= new size_t [NC+1+NC*NC];
    kik_prob->M->matrix1->blocksize0 = new unsigned int[NC]; 
    kik_prob->q                      = new double  [NC3 + NC + NC3 + NC3  + NC3 * NC3];
    kik_prob->mu                     = &(kik_prob->q[NC3                  ]);
    kik_reaction                     = &(kik_prob->q[NC3 + NC             ]);
    kik_velocity                     = &(kik_prob->q[NC3 + NC + NC3       ]);
    kik_prob->M->matrix1->block[0]    = &(kik_prob->q[NC3 + NC + NC3 + NC3 ]);
    for(int i=0;i<NC3 + NC + NC3 + NC3  + NC3 * NC3;i++) kik_prob->q[i]=0;
  }
}
void SPCore::kik_DeleteBuffer()
{
  if(kik_prob->q==0) return;
  if( USE_FULL_MATRIX )
  {
    delete [] kik_prob->q          ;
  }
  else
  {
    delete [] kik_prob->M->matrix1->block;
    delete [] kik_prob->M->matrix1->blocksize0;
    delete [] kik_prob->M->matrix1->index1_data;
    delete [] kik_prob->q          ;
  }
  kik_prob->q = 0;
}

bool SPCore::callSiconosSolver(MatrixX& Mlcp, VectorX& b, VectorX& solution, VectorX& contactIndexToMu,ofstream& os)
{
  int NC3 = Mlcp.rows();
  int NC = NC3/3;
  int CFS_DEBUG = 0;
  int CFS_DEBUG_VERBOSE = 0;
  if(CFS_DEBUG)
  {
    if(NC3%3 != 0           ){ os << "   warning-1 " << std::endl;return false;}
    if(       b.rows()!= NC3){ os << "   warning-2 " << std::endl;return false;}
    if(solution.rows()!= NC3){ os << "   warning-3 " << std::endl;return false;}
  } 
  for(int ia=0;ia<NC;ia++)for(int i=0;i<3;i++)kik_prob->q [3*ia+i]= b(((i==0)?(ia):(2*ia+i+NC-1)));
  for(int ia=0;ia<NC;ia++)                    kik_prob->mu[  ia  ]= contactIndexToMu[ia];
  kik_prob->numberOfContacts = NC;

  if( USE_FULL_MATRIX )
  {
    kik_prob->M->storageType = 0;
    kik_prob->M->size0       = NC3;
    kik_prob->M->size1       = NC3;
    for(int ia=0;ia<NC;ia++)for(int i =0;i <3 ;i ++)
    {
      for(int ja=0;ja<NC;ja++)for(int j =0;j <3;j ++) 
      {
        kik_prob->M->matrix0[NC3*(3*ia+i)+(3*ja+j)]=Mlcp(((i==0)?(ia):(2*ia+i+NC-1)),((j==0)?(ja):(2*ja+j+NC-1)));
      }
    }
  }
  else
  {
    kik_prob->M->storageType = 1;
    kik_prob->M->size0       = NC3;
    kik_prob->M->size1       = NC3;
    SPGlobal_sparsify_A( kik_prob->M->matrix1 , Mlcp , NC , &os);
  }
  
  frictionContact3D_driver(kik_prob,kik_reaction,kik_velocity,kik_solops, kik_numops);
  
  for(int ia=0;ia<NC;ia++)for(int i=0;i<3;i++) solution(((i==0)?(ia):(2*ia+i+NC-1))) = kik_reaction[3*ia+i] ;
  if(CFS_DEBUG_VERBOSE)
  {
    os << "=---------------------------------="<< std::endl; 
    os << "| res_error =" << kik_solops->dparam[1] <<  std::endl;
    os << "=---------------------------------="<< std::endl; 
  }
  return true;
}




void SPCore::SPGlobal_sparsify_A(SparseBlockStructuredMatrix* pmat, MatrixX& Mlcp, int NC, ofstream* pos )
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
      if(SPGlobal_check_zero_block(Mlcp,NC,ia,ja)==1)
      {
        mat.block[NB] = mat.block[0]+(9*NB);
        SPGlobal_copy_block(mat.block[NB], Mlcp, NC,ia,ja);
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
  if(0)//CFS_DEBUG)
  {
    int* ibuf = new int[NC*NC];
    SPGlobal_construct_sparsity_matrix(ibuf,Mlcp,NC);
    (*pos) << "=---------------------------------="<< std::endl; 
    (*pos) << "const int NC=" << NC << ";"<< std::endl; 
    (*pos) << "const int NB=" << NB << ";"<< std::endl; 
    for(int i=0;i<NC;i++){for(int j=0;j<NC;j++){(*pos) << ibuf[NC*i+j];} (*pos)<<std::endl;}
    (*pos) << "=---------------------------------="<< std::endl; 
    delete [] ibuf;
  }
}

