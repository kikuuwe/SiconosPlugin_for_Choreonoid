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


#ifndef CNOID_SICONOSPLUGIN_CORE_H
#define CNOID_SICONOSPLUGIN_CORE_H

#include "SiconosNumerics.h"
#include "FrictionContact3D_compute_error.h"

using namespace std;

//class FrictionContactProblem ;
//class NumericsOptions        ;
//class SolverOptions          ;
//class SparseBlockStructuredMatrix ;

namespace cnoid
{

class SPCore
{
  public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;
    typedef VectorXd VectorX;
  
    bool USE_FULL_MATRIX ;
    void kik_NewBuffer   (int aNC) ;
    void kik_DeleteBuffer();
    FrictionContactProblem *kik_prob      ;
    NumericsOptions        *kik_numops    ;
    SolverOptions          *kik_solops    ;
    double                 *kik_reaction  ;
    double                 *kik_velocity  ;
  
    SPCore(int maxNumGaussSeidelIteration, double gaussSeidelErrorCriterion);
    ~SPCore();
    bool   callSiconosSolver(MatrixX& Mlcp, VectorX& b, VectorX& solution, VectorX& contactIndexToMu,ofstream& os);    
    static void SPGlobal_sparsify_A(SparseBlockStructuredMatrix* pmat, MatrixX& Mlcp, int NC, ofstream* pos );

    static int SPGlobal_check_zero_block(MatrixX& Mlcp, int NC, int ia,int ja)
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

    static int SPGlobal_construct_sparsity_matrix(int * ibuf, MatrixX& Mlcp, int NC)
    {
      int count=0;
      for(int ia = 0; ia<NC; ia++)for(int ja = ia; ja<NC; ja++) 
      {
        if(ia==ja){ibuf[NC*ia+ja]=1; count ++;}
        else
        {
          int tmp = SPGlobal_check_zero_block(Mlcp,NC,ia,ja)  ;
          ibuf[NC*ia+ja]= ibuf[NC*ja+ia]= tmp;
          if(tmp==1) count += 2;
        }
      }
      return count;
    }


    void SPGlobal_check_offdiag(const int * ibuf, int NC, int* pia2, int *pja2)
    {
      (*pia2)=(*pja2)=NC;
      for(int k = 0; k<NC; k++)
      {
                                    if(ibuf[NC* k + (NC-k-1) ]==1){(*pia2)=(*pja2)=k; return;}
        for(int ja=NC-k;ja<NC;ja++) if(ibuf[NC* k + ja       ]==1){(*pia2)=(*pja2)=k; return;}
        for(int ia=k-1 ;ia>=0;ia--) if(ibuf[NC*ia + (NC-k-1) ]==1){(*pia2)=(*pja2)=k; return;}
      }
    }

    static void SPGlobal_copy_block(double * bbuf, MatrixX& Mlcp, int NC, int ia,int ja)
    {
      for(int i=0;i<3;i++)for(int j=0;j<3;j++) bbuf[3*j+i]= Mlcp(((i==0)?(ia):(2*ia+i+NC-1)),((j==0)?(ja):(2*ja+j+NC-1))) ;
    }




};

};

#endif
