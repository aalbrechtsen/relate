/*
  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.81
  
  HMMld is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Foobar is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with HMMld.  If not, see <http://www.gnu.org/licenses/>.
*/  


/*
  this is the header that contains the class definition
  The program is using a c-version of bfgs.
  This version requires the paramters to be optimized to be a double*.
  So all combinations of different optimize paramters has to be written in a wrapper funktion.
*/


#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#include <vector>
snpMatrix *snp_pair_range(iMatrix *v, int  depth);


class relateHMM{
private:
  static dArray *Sk1, *Sk2,*Sk3, *t, *alim;
  static  double fixA_val,fixK2_val,calcA_val,phi;
  static  int double_recom,fixA, fixK2, calcA,firstRun,localTimesRun,localTimesConv;
  static pars *calculatedValues;
  static dMatrix *convInfo; //list of converged points
 

  static bArray *pruneList;
  static iMatrix *genos_global;
  static dMatrix *mea_global;
  static dArray *maf1_global;
  static dArray *maf2_global;
  static dArray *pos_global;
  static hapStruct *hap_global;
  static bArray *pre_calc_used_list;
  void pre_calc(const functionPars *pars); // will calculate globalstuff
  
public:
  void init_globals(const functionPars *pars);//will simply copy global stuff
  void dinit_globals(int back);
  static double *full_optim(double* var);
  static double full_optim_fun(double *var);
  static double fixK2_calcA_fun(double *var);
  static double *fixK2_calcA_optim(double* var);
  static double fixK2_fun(double *var);
  static double *fixK2_optim(double* var);
  static double fixA_fun(double *var);
  static double *fixA_optim(double* var);
  static double like_optim_one_dim(double* var);
  static double fixK2_fixA_fun(double *var);
  static double *fixK2_fixA_optim(double* var);
  static double calcA_fun(double *var);
  static double *calcA_optim(double* var);

  fullResults *single_pair(const functionPars *pars);
  dArray *run_optimization(int i,int j,double convTol,const dArray *alim);
  int all_pairs(functionPars *pars);
  int fast_all_pairs(functionPars *pars);
  int joblist_fast_all_pairs(functionPars *pars,std::vector<iArray*> pa);
  fullResults *sgl(const functionPars *pars);
};
