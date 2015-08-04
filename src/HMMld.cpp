/*
  Copyright (C) 2009 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.98
  
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

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <queue>
#include <deque>
#include <map>         
#include <string>
#include "relateHMM.h"
#include "asort.h"
#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif

#include "filereader_and_conversions.h"
#include "r_c_conversions.h"
using namespace std;

#include "conf.h"

void setDefaultPars(functionPars *pars){
  iArray *pair = allocIntArray(2);
  dArray *alim = allocDoubleArray(2);
  dArray *par = allocDoubleArray(3);
  
  pars->pair=pair;
  pars->alim = alim;  
  pars->alim->array[0]=0.001;
  pars->alim->array[1]=0.15;

  pars->par = par;
  pars->chr=NULL;
  pars->position=NULL;
  pars->data=NULL;

  pars->doPar=0;
  pars->min=0;
  pars->double_recom=0;
  pars->LD=0;
  pars->ld_adj=0;
  pars->epsilon=0.01;
  pars->back=5;
  pars->doPrune =0;
  
  pars->fixA=0;
  pars->fixK=0;
  //  pars->abs=0;
  pars->calcA =0;
  pars->fixA_val=0;
  pars->fixK_val=0;

  //added in 0.68 version for wrongfull parsing from pars in R
  pars->prune_val=0;
  pars->doAllPairs =0;
  pars->print_results=0;
}

SEXP getListElement(SEXP list,const char *str)     {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}


/*
  This functions eats the data from a native snp.matrix.pedfile object

 */
fromSnpMatrix *buildData(SEXP snpClass){

  SEXP posList=R_NilValue;
  SEXP snpList=R_NilValue;
  SEXP chrList=R_NilValue;
  SEXP v = snpClass;// R_NilValue; //0.98
  //extract the snp.subject file
  PROTECT(v = getListElement(snpClass,"snp.support")); //0.98
  PROTECT(posList = getListElement(v,"position"));
  PROTECT(snpList = getListElement(v,"snp.names"));
  PROTECT(chrList = getListElement(v,"chromosome"));

  printf("Len of posList %i\n",length(posList));
  printf("Len of chrList %i\n",length(chrList));
  dArray *pos = allocDoubleArray(length(posList));
  iArray *chr = allocIntArray(pos->x);

  for(int i=0;i<pos->x;i++){
    pos->array[i] = INTEGER(posList)[i];
    chr->array[i] = INTEGER(chrList)[i];
  }
  
  //now get the snp.data frame
  SEXP s = R_NilValue;
  PROTECT(s = getListElement(snpClass,"snp.data"));
  SEXP dims = R_NilValue;
  PROTECT(dims = getAttrib(s, R_DimSymbol));
  int rows=INTEGER(dims)[0];
  int cols=INTEGER(dims)[1];
  
  
  SEXP dimnames = R_NilValue;
  PROTECT(dimnames = allocVector(VECSXP, cols));
  SET_VECTOR_ELT(dimnames, 1, getAttrib(s, R_NamesSymbol));
  dimnames = coerceVector(dimnames, STRSXP);
  iMatrix *retVal = allocIntMatrix(rows,cols);
  //  printf("\nlen of dimnames:%d\n",length(dimnames));
  for(int j=0;j<cols;j++){
    //   string myStr = string(CHAR(STRING_ELT(snpList, j)));
    //int index_to_input = snpName_vs_index[myStr];
    // printf("ind_to_input=%d\n",index_to_input);
    for (int i=0;i<rows;i++){
      retVal->matrix[i][j] = RAW(s)[j*rows+i];
    }
  }
  //  print_matrix(retVal->matrix,retVal->x,retVal->y);
  fromSnpMatrix *allData = new fromSnpMatrix();
  allData->data=retVal;
  allData->chr=chr;
  allData->pos=pos;
  UNPROTECT(7);

  return allData;
}

extern "C" {
  SEXP readbed(SEXP Bed, SEXP Id, SEXP Snps); 

  //used for commandline args
  SEXP interface(SEXP Rdata,SEXP Rpair,SEXP Rposition,SEXP Rpar,	\
	       SEXP Rmin,SEXP Rdouble_recom,SEXP RLD,SEXP Repsilon,SEXP Rback, \
	       SEXP Ralim,SEXP Roptim,SEXP Rstart,SEXP Rprune,SEXP Rld_adj, \
		 SEXP Rfix_a,SEXP Rfix_k2,SEXP Rchr,SEXP Rmoment,SEXP calcA,SEXP \
		 myphi,SEXP timesRun,SEXP timesConv,SEXP giveCrap,SEXP convTol,SEXP Rmethod,SEXP back2);

  SEXP doPlink(SEXP bedfilename,SEXP dim1,SEXP dim2, SEXP Rpair,SEXP Rposition,SEXP Rpar, \
	       SEXP Rmin,SEXP Rdouble_recom,SEXP RLD,SEXP Repsilon,SEXP Rback, \
	       SEXP Ralim,SEXP Roptim,SEXP Rstart,SEXP Rprune,SEXP Rld_adj, \
	       SEXP Rfix_a,SEXP Rfix_k2,SEXP Rchr,SEXP Rmoment,SEXP calcA,SEXP\
	       myphi,SEXP timesRun,SEXP timesConv,SEXP giveCrap,SEXP convTol,SEXP Rmethod,SEXP back2, SEXP plinkKeepList);


  SEXP getPedfile(SEXP v);
  //used for args by pedfile thats snp.matrix object
  SEXP interface2(SEXP pedfile,SEXP Rpair,SEXP Rpar,			\
		SEXP Rmin,SEXP Rdouble_recom,SEXP RLD,SEXP Repsilon,SEXP Rback,	\
		SEXP Ralim,SEXP Roptim,SEXP Rstart,SEXP Rprune,SEXP Rld_adj, \
		  SEXP Rfix_a,SEXP Rfix_k2,SEXP Rmoment,SEXP calcA,SEXP myphi,SEXP timesRun,SEXP timesConv,SEXP giveCrap,SEXP convTol);
  SEXP snp_pair_ld(SEXP snpdata,SEXP Rback);

}



SEXP getPedfile(SEXP v){
  SEXP ans = R_NilValue;
  if (v==R_NilValue)
    return ans;
  
  //  SEXP retVal = R_NilValue;
  SEXP ans_name = R_NilValue;
  fromSnpMatrix *vars = buildData(v);
  
  PROTECT(ans    = allocVector(VECSXP, 3)); 
  SET_VECTOR_ELT(ans, 0, c_dArray_to_r_array(vars->pos));
  SET_VECTOR_ELT(ans, 1, c_iMatrix_to_r_matrix(vars->data));
  SET_VECTOR_ELT(ans, 2, c_iArray_to_r_array(vars->chr));
  
  PROTECT(ans_name = allocVector(STRSXP, 3));
  SET_STRING_ELT(ans_name, 0, mkChar("position"));
  SET_STRING_ELT(ans_name, 1, mkChar("data"));
  SET_STRING_ELT(ans_name, 2, mkChar("chromosome"));
  setAttrib(ans, R_NamesSymbol, ans_name);
  UNPROTECT(2);
  return ans;
}
  
  //this is for getting the data through command args
SEXP interface(SEXP Rdata,SEXP Rpair,SEXP Rposition,SEXP Rpar,		\
	       SEXP Rmin,SEXP Rdouble_recom,SEXP RLD,SEXP Repsilon,SEXP Rback, \
	       SEXP Ralim,SEXP Roptim,SEXP Rstart,SEXP Rprune,SEXP Rld_adj, \
	       SEXP Rfix_a,SEXP Rfix_k2,SEXP Rchr,SEXP Rmoment,SEXP calcA,SEXP\
	       myphi,SEXP timesRun,SEXP timesConv,SEXP giveCrap,SEXP convTol,SEXP Rmethod,SEXP back2) { 
  //init returnValue
  //    printf("start of interface\n");
  SEXP ans = R_NilValue;
  SEXP ans_name = R_NilValue;
  SEXP hap_names= R_NilValue;
  SEXP hap = R_NilValue;
  SEXP class_name = R_NilValue;
  
  //these are the vars to be inputet
  //SEXP kResult=R_NilValue, kLike = R_NilValue;
  
  //init structure cointain all values needed to run program
  functionPars *pars = new functionPars();
  setDefaultPars(pars);
  
  //now input values
  //these are required
  if(Rdata==R_NilValue||Rpair==R_NilValue||Rposition==R_NilValue){
    Rprintf("Must supply arguments\n");
    return ans;
  }
  
  pars->data = r_matrix_to_c_iMatrix(Rdata);
  pars->position = r_array_to_c_dArray(Rposition);
  pars->pair = r_array_to_c_iArray(Rpair);
  // c is zero-index therefor subtract one from pair
  pars->pair->array[0]--;
  pars->pair->array[1]--;
  
  //added in 0.68 forgot min value
  pars->min = r_real_to_c_double(Rmin);
  
  //these requires only one change to struct
  pars->back = r_int_to_c_int(Rback);
  pars->double_recom = r_bool_to_c_int(Rdouble_recom);
  pars->LD = r_bool_to_c_int(RLD);
  pars->epsilon = r_real_to_c_double(Repsilon);
  //    pars->abs = r_bool_to_c_int(Rabs);
  pars->alim = r_array_to_c_dArray(Ralim);
  pars->ld_adj = r_bool_to_c_int(Rld_adj);
  
  pars->calcA = r_bool_to_c_int(calcA);
  pars->print_results = r_int_to_c_int(giveCrap);
  pars->timesToRun = r_int_to_c_int(timesRun);
  pars->timesToConverge = r_int_to_c_int(timesConv);
  //    printf("calcA in pack.cpp=%d\n",pars->calcA);
  pars->phi = r_real_to_c_double(myphi);
  pars->convTol = r_real_to_c_double(convTol);
  //these requires more than one change to struct
  if(Rpar!=R_NilValue){
    pars->doPar=1;
    pars->par = r_array_to_c_dArray(Rpar);
  }
  
  if (Rchr!=R_NilValue){
    pars->chr = r_array_to_c_iArray(Rchr);
    mysort(pars,pars->print_results);
  }

  if(Rprune!=R_NilValue){
    pars->prune_val = r_real_to_c_double(Rprune);
    pars->doPrune =1;
  }
  if(Rfix_a!=R_NilValue){
    pars->fixA=1;
    pars->fixA_val = r_real_to_c_double(Rfix_a);
  }
  if(Rfix_k2!=R_NilValue){
    pars->fixK=1;
    pars->fixK_val = r_real_to_c_double(Rfix_k2);
    //printf("fixK2_val in pack.cpp=%f\n",pars->fixK_val);
  }
  if(Rpar!=R_NilValue){
    pars->doPar=1;
    pars->par = r_array_to_c_dArray(Rpar);
  }
  //aded in 0.88
  pars->back2 = r_int_to_c_int(back2);
  //now do calculation and optimization
  relateHMM object;
  //    print_functionPars(pars);
  int newVersion = r_int_to_c_int(Rmethod);
  
  fullResults *fromObject;
  int extras = 0;//change in the end
  
  if(newVersion)
    fromObject = object.single_pair(pars);
  else{
    extras = 5;
    object.init_globals(pars);
    fromObject = object.sgl(pars);
 
  }
  //now prepare results for R
  
  //now the data
  
  int verbose = r_bool_to_c_int(giveCrap);
  
  int itemNumber  = 0 ;
  if(verbose==0)
    PROTECT(ans    = allocVector(VECSXP, 20 + extras)); //added 5
  else{
    PROTECT(ans    = allocVector(VECSXP, 26 + extras)); //added 5
    PROTECT(hap = allocVector(VECSXP,4));
  }
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->kResult));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->kLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->kr));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->a));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->uLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_bool(fromObject->LD));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->t));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->snp));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->position));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_bool(fromObject->double_recom));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->alim));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->poLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->back));
  
  //will continue
  
  SET_VECTOR_ELT(ans, itemNumber++, c_iArray_to_r_array(fromObject->chr));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->post));
  
  // conv info
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->timesRun));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->timesConverged));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->convInfo));
  SET_VECTOR_ELT(ans, itemNumber++, c_iArray_to_r_array(fromObject->keepList));
  SET_VECTOR_ELT(ans, itemNumber++ , c_iArray_to_r_array(fromObject->path));      
  
  //extras
  if(!newVersion){
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->mea_global));
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pba_global));
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pBa_global));
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pbA_global));
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pBA_global));
  
  }

  //print_array(fromObject->path->array,fromObject->path->x);
  if(verbose==1){
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->S));      
    SET_VECTOR_ELT(ans,itemNumber++ , c_iArray_to_r_array(fromObject->choose));      
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->hap->mea));
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->S1));
    //will now load hap
    SET_VECTOR_ELT(hap, 0, c_dMatrix_to_r_matrix(fromObject->hap->pBA));
    SET_VECTOR_ELT(hap, 1, c_dMatrix_to_r_matrix(fromObject->hap->pBa));
    SET_VECTOR_ELT(hap, 2, c_dMatrix_to_r_matrix(fromObject->hap->pbA));
    SET_VECTOR_ELT(hap, 3, c_dMatrix_to_r_matrix(fromObject->hap->pba));
    // now input the hap ind the main list
    SET_VECTOR_ELT(ans,itemNumber++ , hap);
    SET_VECTOR_ELT(ans,itemNumber++ , c_dArray_to_r_array(fromObject->maf));
  }
  
  
  //now input the main names 
  
  itemNumber  = 0 ;
  if(verbose ==1)
    PROTECT(ans_name = allocVector(STRSXP, 26+extras));//added 5
  else
    PROTECT(ans_name = allocVector(STRSXP, 20+extras));//added 5
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kResult"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kLike"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kr"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("a"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("uLike"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("LD"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("t"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("snp"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("position"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("double_recom"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("alim"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("poLike"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("back"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("chr"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("post"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("timesRun"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("timesConverged"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("convergenceInfo"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("usedSnps"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("path"));  
  //extras
  if(!newVersion){
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("mea_global"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("pba_global"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("pBa_global"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("pbA_global"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("pBA_global"));    
    object.dinit_globals(pars->back);
  }


  if(verbose==1){
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("S"));      
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("choose"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("mea"));    
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("S1"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("hap"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("maf"));
    //now input the hap names 
    PROTECT(hap_names = allocVector(STRSXP, 4));
    SET_STRING_ELT(hap_names, 0, mkChar("pBA"));
    SET_STRING_ELT(hap_names, 1, mkChar("pBa"));
    SET_STRING_ELT(hap_names, 2, mkChar("pbA"));
    SET_STRING_ELT(hap_names, 3, mkChar("pba"));
    setAttrib(hap, R_NamesSymbol, hap_names);	
  }
    

  //now install the names
  
  setAttrib(ans, R_NamesSymbol, ans_name);
  
  //now set up class name
  PROTECT(class_name = allocVector(STRSXP, 1));
  SET_STRING_ELT(class_name, 0, mkChar("pack.class"));
  classgets(ans, class_name);  
  
  setAttrib(ans, R_NamesSymbol, ans_name);
  if(verbose==0)
    UNPROTECT(3);
  else
    UNPROTECT(5);
  
  //clean up
 
  killFullResults(fromObject);
  killFunctionPars(pars);
  //printf("end of fun\n");
  return(ans) ;
}

// as of december the 19 interface2 was forcked from interface1 this means v. 0.81
SEXP interface2(SEXP pedfile,SEXP Rpair,SEXP Rpar,			\
		SEXP Rmin,SEXP Rdouble_recom,SEXP RLD,SEXP Repsilon,SEXP Rback,	\
		SEXP Ralim,SEXP Roptim,SEXP Rstart,SEXP Rprune,SEXP Rld_adj, \
		SEXP Rfix_a,SEXP Rfix_k2,SEXP Rmoment,SEXP calcA,SEXP myphi,SEXP timesRun,SEXP timesConv,SEXP giveCrap,SEXP convTol) { 
  //init returnValue
  //     printf("start of interface\n");
  SEXP ans = R_NilValue;
  SEXP ans_name = R_NilValue;
  SEXP hap_names= R_NilValue;
  SEXP hap = R_NilValue;
  SEXP class_name = R_NilValue;
  
  //these are the vars to be inputet
  //SEXP kResult=R_NilValue, kLike = R_NilValue;
  
  //init structure cointain all values needed to run program
  functionPars *pars = new functionPars();
  setDefaultPars(pars);
  
  fromSnpMatrix *vars = buildData(pedfile);
  pars->data = vars->data;
  pars->chr = vars->chr;
  pars->position = vars->pos;
  if(pars->chr->x != pars->position->x || pars->chr->x!=pars->data->y){
    printf("Dimension af datastructures aren't equal\n");
    printf("Length of positionvector:%d\n",pars->position->x);
    printf("Length of chromosomevector:%d\n",pars->chr->x);
    printf("Dimension of genotypes: (%d,%d)\n",pars->data->x,pars->data->y);
    return ans;
  }
  
  pars->pair = r_array_to_c_iArray(Rpair);
  // c is zero-index therefor subtract one from pair
  pars->pair->array[0]--;
  pars->pair->array[1]--;
  if(pars->pair->array[0] > pars->data->y-1 || pars->pair->array[1] >pars->data->y-1){
    printf("Pair must be within dim of genotypes");
    return ans;
  }
  
  
  //added in 0.68 forgot min value
  pars->min = r_real_to_c_double(Rmin);
  
  //these requires only one change to struct
  pars->back = r_int_to_c_int(Rback);
  pars->double_recom = r_bool_to_c_int(Rdouble_recom);
  pars->LD = r_bool_to_c_int(RLD);
  pars->epsilon = r_real_to_c_double(Repsilon);
  //    pars->abs = r_bool_to_c_int(Rabs);
  pars->alim = r_array_to_c_dArray(Ralim);
  pars->ld_adj = r_bool_to_c_int(Rld_adj);
  
  pars->calcA = r_bool_to_c_int(calcA);
  pars->print_results = r_int_to_c_int(giveCrap);
  pars->timesToRun = r_int_to_c_int(timesRun);
  pars->timesToConverge = r_int_to_c_int(timesConv);
  //    printf("calcA in pack.cpp=%d\n",pars->calcA);
  pars->phi = r_real_to_c_double(myphi);
  pars->convTol = r_real_to_c_double(convTol);
  //these requires more than one change to struct
  if(Rpar!=R_NilValue){
    pars->doPar=1;
    pars->par = r_array_to_c_dArray(Rpar);
  }
  
  if(Rprune!=R_NilValue){
    pars->prune_val = r_real_to_c_double(Rprune);
    pars->doPrune =1;
  }
  if(Rfix_a!=R_NilValue){
    pars->fixA=1;
    pars->fixA_val = r_real_to_c_double(Rfix_a);
  }
  if(Rfix_k2!=R_NilValue){
    pars->fixK=1;
    pars->fixK_val = r_real_to_c_double(Rfix_k2);
    //printf("fixK2_val in pack.cpp=%f\n",pars->fixK_val);
  }
  if(Rpar!=R_NilValue){
    pars->doPar=1;
    pars->par = r_array_to_c_dArray(Rpar);
  }
  
  for (int i=0;i<pars->position->x;i++)
    pars->position->array[i] = pars->position->array[i]/PLINK_POS_SCALING;


  //now do calculation and optimization
  relateHMM object;
  //    print_functionPars(pars);
  fullResults *fromObject = object.single_pair(pars);
  
  //now prepare results for R
  
  //now the data
  
  int verbose = r_bool_to_c_int(giveCrap);

  int itemNumber  = 0 ;
  if(verbose==0)
    PROTECT(ans    = allocVector(VECSXP, 20)); 
  else{
    PROTECT(ans    = allocVector(VECSXP, 26)); 
    PROTECT(hap = allocVector(VECSXP,4));
  }
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->kResult));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->kLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->kr));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->a));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->uLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_bool(fromObject->LD));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->t));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->snp));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->position));
  //    print_array(fromObject->position->array,fromObject->position->x);
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_bool(fromObject->double_recom));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->alim));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->poLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->back));
  
  //will continue
  
  SET_VECTOR_ELT(ans, itemNumber++, c_iArray_to_r_array(fromObject->chr));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->post));
  
  // conv info
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->timesRun));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->timesConverged));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->convInfo));
  SET_VECTOR_ELT(ans, itemNumber++, c_iArray_to_r_array(fromObject->keepList));
  SET_VECTOR_ELT(ans, itemNumber++ , c_iArray_to_r_array(fromObject->path));      
  //print_array(fromObject->path->array,fromObject->path->x);
  if(verbose==1){
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->S));      
    SET_VECTOR_ELT(ans,itemNumber++ , c_iArray_to_r_array(fromObject->choose));      
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->hap->mea));
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->S1));
    //will now load hap
    SET_VECTOR_ELT(hap, 0, c_dMatrix_to_r_matrix(fromObject->hap->pBA));
    SET_VECTOR_ELT(hap, 1, c_dMatrix_to_r_matrix(fromObject->hap->pBa));
    SET_VECTOR_ELT(hap, 2, c_dMatrix_to_r_matrix(fromObject->hap->pbA));
    SET_VECTOR_ELT(hap, 3, c_dMatrix_to_r_matrix(fromObject->hap->pba));
    // now input the hap ind the main list
    SET_VECTOR_ELT(ans,itemNumber++ , hap);
      SET_VECTOR_ELT(ans,itemNumber++ , c_dArray_to_r_array(fromObject->maf));
  }
  
  
  //now input the main names 
  
  itemNumber  = 0 ;
  if(verbose ==1)
    PROTECT(ans_name = allocVector(STRSXP, 26));
  else
    PROTECT(ans_name = allocVector(STRSXP, 20));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kResult"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kLike"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kr"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("a"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("uLike"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("LD"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("t"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("snp"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("position"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("double_recom"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("alim"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("poLike"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("back"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("chr"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("post"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("timesRun"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("timesConverged"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("convergenceInfo"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("usedSnps"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("path"));  
  if(verbose==1){
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("S"));      
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("choose"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("mea"));    
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("S1"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("hap"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("maf"));
    //now input the hap names 
    PROTECT(hap_names = allocVector(STRSXP, 4));
    SET_STRING_ELT(hap_names, 0, mkChar("pBA"));
    SET_STRING_ELT(hap_names, 1, mkChar("pBa"));
    SET_STRING_ELT(hap_names, 2, mkChar("pbA"));
    SET_STRING_ELT(hap_names, 3, mkChar("pba"));
    setAttrib(hap, R_NamesSymbol, hap_names);	
  }
  

  //now install the names
  
  setAttrib(ans, R_NamesSymbol, ans_name);
  
  //now set up class name
  PROTECT(class_name = allocVector(STRSXP, 1));
  SET_STRING_ELT(class_name, 0, mkChar("pack.class"));
  classgets(ans, class_name);  
  
  setAttrib(ans, R_NamesSymbol, ans_name);
  if(verbose==0)
    UNPROTECT(3);
  else
    UNPROTECT(5);
  
  //clean up
    
  killFullResults(fromObject);
  killFunctionPars(pars);
  //printf("end of fun\n");
  return(ans) ;
}




SEXP snp_pair_ld(SEXP snpdata,SEXP Rback) { 
  //init returnValue
  //  printf("start of snp_pair_ld\n");
  SEXP ans = R_NilValue;
  SEXP ans_name = R_NilValue;
  //  SEXP hap_names= R_NilValue;
  // SEXP hap = R_NilValue;
  SEXP class_name = R_NilValue;
  
  iMatrix *data = r_matrix_to_c_iMatrix(snpdata);
  int back = r_int_to_c_int(Rback);
  //  printf("dims=(%d,%d)\n",data->x,data->y);  
  
  //now do calculation and optimization
  snpMatrix *fromObject = snp_pair_range(data,back);
  
  int itemNumber  = 0 ;
  //prepare
  PROTECT(ans    = allocVector(VECSXP, 8)); 
  //convert and input the data
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->dprime));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pBA));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pBa));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pbA));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pba));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->rmisc));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->lod));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->D));
    
  
  
  //now input the names 
  
  itemNumber  = 0 ;
  PROTECT(ans_name = allocVector(STRSXP, 8));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("dprime"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("pBA"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("pBa"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("pbA"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("pba"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("rmisc"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("lod"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("D"));
    

  //now install the names
  
  setAttrib(ans, R_NamesSymbol, ans_name);
  
  //now set up class name
  PROTECT(class_name = allocVector(STRSXP, 1));
  SET_STRING_ELT(class_name, 0, mkChar("ld.class"));
  classgets(ans, class_name);  
  setAttrib(ans, R_NamesSymbol, ans_name);
  
  //cleanup
  killSnpMatrix(fromObject);
  UNPROTECT(3);
  //printf("end of snp_pair_ld\n");
  return(ans) ;

}



SEXP readbed(SEXP Bed, SEXP Id, SEXP Snps) {
  int nrow = LENGTH(Id);
  int ncol = LENGTH(Snps);

  const char *file = CHAR(STRING_ELT(Bed, 0));
  
  iMatrix *returnMat = bed_to_iMatrix(file,nrow,ncol);

  SEXP ans,ans_name;
  
  PROTECT(ans    = allocVector(VECSXP, 1)); 
  SET_VECTOR_ELT(ans, 0, c_iMatrix_to_r_matrix(returnMat));
  PROTECT(ans_name = allocVector(STRSXP, 1));
  SET_STRING_ELT(ans_name, 0, mkChar("genotypes"));
  setAttrib(ans, R_NamesSymbol, ans_name);

  UNPROTECT(2);
  killMatrix(returnMat);
  return ans;
}



  //this is for getting the data through command args
SEXP doPlink(SEXP bedfilename,SEXP dim1,SEXP dim2, SEXP Rpair,SEXP Rposition,SEXP Rpar, \
	       SEXP Rmin,SEXP Rdouble_recom,SEXP RLD,SEXP Repsilon,SEXP Rback, \
	       SEXP Ralim,SEXP Roptim,SEXP Rstart,SEXP Rprune,SEXP Rld_adj, \
	       SEXP Rfix_a,SEXP Rfix_k2,SEXP Rchr,SEXP Rmoment,SEXP calcA,SEXP\
	     myphi,SEXP timesRun,SEXP timesConv,SEXP giveCrap,SEXP convTol,SEXP Rmethod,SEXP back2, SEXP plinkKeepList) { 
  //init returnValue
  //    printf("start of interface\n");
  SEXP ans = R_NilValue;
  SEXP ans_name = R_NilValue;
  SEXP hap_names= R_NilValue;
  SEXP hap = R_NilValue;
  SEXP class_name = R_NilValue;
  
  //these are the vars to be inputet
  //SEXP kResult=R_NilValue, kLike = R_NilValue;
  
  //init structure cointain all values needed to run program
  functionPars *pars = new functionPars();
  setDefaultPars(pars);
  
  //now input values
  //these are required
  if(bedfilename==R_NilValue||Rpair==R_NilValue||Rposition==R_NilValue){
    Rprintf("Must supply arguments\n");
    return ans;
  }
  

  pars->position = r_array_to_c_dArray(Rposition);
  pars->pair = r_array_to_c_iArray(Rpair);
  // c is zero-index therefor subtract one from pair
  pars->pair->array[0]--;
  pars->pair->array[1]--;
  
  //added in 0.68 forgot min value
  pars->min = r_real_to_c_double(Rmin);
  
  //these requires only one change to struct
  pars->back = r_int_to_c_int(Rback);
  pars->double_recom = r_bool_to_c_int(Rdouble_recom);
  pars->LD = r_bool_to_c_int(RLD);
  pars->epsilon = r_real_to_c_double(Repsilon);
  //    pars->abs = r_bool_to_c_int(Rabs);
  pars->alim = r_array_to_c_dArray(Ralim);
  pars->ld_adj = r_bool_to_c_int(Rld_adj);
  
  pars->calcA = r_bool_to_c_int(calcA);
  pars->print_results = r_int_to_c_int(giveCrap);
  pars->timesToRun = r_int_to_c_int(timesRun);
  pars->timesToConverge = r_int_to_c_int(timesConv);
  //    printf("calcA in pack.cpp=%d\n",pars->calcA);
  pars->phi = r_real_to_c_double(myphi);
  pars->convTol = r_real_to_c_double(convTol);
  //these requires more than one change to struct
  if(Rpar!=R_NilValue){
    pars->doPar=1;
    pars->par = r_array_to_c_dArray(Rpar);
  }
  
  if(Rprune!=R_NilValue){
    pars->prune_val = r_real_to_c_double(Rprune);
    pars->doPrune =1;
  }
  if(Rfix_a!=R_NilValue){
    pars->fixA=1;
    pars->fixA_val = r_real_to_c_double(Rfix_a);
  }
  if(Rfix_k2!=R_NilValue){
    pars->fixK=1;
    pars->fixK_val = r_real_to_c_double(Rfix_k2);
    //printf("fixK2_val in pack.cpp=%f\n",pars->fixK_val);
  }
  if(Rpar!=R_NilValue){
    pars->doPar=1;
    pars->par = r_array_to_c_dArray(Rpar);
  }
  //aded in 0.88
  pars->back2 = r_int_to_c_int(back2);
  //now do calculation and optimization
  

  // we need to get the genotypes, if the length of positions are greater
  // than the number of snps in the genotypes, we should remove the elements given by plinkKeeplist
  const char *file = CHAR(STRING_ELT(bedfilename, 0));
  iMatrix *tmp = bed_to_iMatrix(file,r_int_to_c_int(dim1),r_int_to_c_int(dim2)); 
  

  if(pars->position->x==tmp->y){
    //  printf("yoyoyoyoyo is smae\n");
    pars->data = tmp;
    
  }
  else{
    
    iArray *toKeep = r_array_to_c_iArray(plinkKeepList);
    //just make sure the length of bed is the same as the original bim
    if(toKeep->x!=tmp->y){
      printf("\t-> Length of bim file isn't the same as the bed file will exit\n");
      exit(0);
    }
    iMatrix *lastTmp = allocIntMatrix(tmp->x,pars->position->x);
    int pos = 0;
    for(int i=0 ; i<tmp->y;i++){
      if(toKeep->array[i]!=0){
	for(int j=0; j<tmp->x; j++ )
	  lastTmp->matrix[j][pos] = tmp->matrix[j][i];
	pos++;
      }
    }
    pars->data = lastTmp;
    killArray(toKeep);
    killMatrix(tmp);
    
  }
    
  if (Rchr!=R_NilValue){
    pars->chr = r_array_to_c_iArray(Rchr);
    mysort(pars,pars->print_results);
  }


  relateHMM object;
  //    print_functionPars(pars);
  int newVersion = r_int_to_c_int(Rmethod);
  
  fullResults *fromObject;
  int extras = 0;//change in the end
  
  if(newVersion)
    fromObject = object.single_pair(pars);
  else{
    extras = 5;
    object.init_globals(pars);
    fromObject = object.sgl(pars);
 
  }
  //now prepare results for R
  
  //now the data
  
  int verbose = r_bool_to_c_int(giveCrap);
  
  int itemNumber  = 0 ;
  if(verbose==0)
    PROTECT(ans    = allocVector(VECSXP, 20 + extras)); //added 5
  else{
    PROTECT(ans    = allocVector(VECSXP, 26 + extras)); //added 5
    PROTECT(hap = allocVector(VECSXP,4));
  }
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->kResult));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->kLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->kr));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->a));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->uLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_bool(fromObject->LD));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->t));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->snp));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->position));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_bool(fromObject->double_recom));
  SET_VECTOR_ELT(ans, itemNumber++, c_dArray_to_r_array(fromObject->alim));
  SET_VECTOR_ELT(ans, itemNumber++, c_double_to_r_real(fromObject->poLike));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->back));
  
  //will continue
  
  SET_VECTOR_ELT(ans, itemNumber++, c_iArray_to_r_array(fromObject->chr));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->post));
  
  // conv info
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->timesRun));
  SET_VECTOR_ELT(ans, itemNumber++, c_int_to_r_int(fromObject->timesConverged));
  SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->convInfo));
  SET_VECTOR_ELT(ans, itemNumber++, c_iArray_to_r_array(fromObject->keepList));
  SET_VECTOR_ELT(ans, itemNumber++ , c_iArray_to_r_array(fromObject->path));      
  
  //extras
  if(!newVersion){
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->mea_global));
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pba_global));
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pBa_global));
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pbA_global));
  
    SET_VECTOR_ELT(ans, itemNumber++, c_dMatrix_to_r_matrix(fromObject->pBA_global));
  
  }

  //print_array(fromObject->path->array,fromObject->path->x);
  if(verbose==1){
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->S));      
    SET_VECTOR_ELT(ans,itemNumber++ , c_iArray_to_r_array(fromObject->choose));      
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->hap->mea));
    SET_VECTOR_ELT(ans,itemNumber++ , c_dMatrix_to_r_matrix(fromObject->S1));
    //will now load hap
    SET_VECTOR_ELT(hap, 0, c_dMatrix_to_r_matrix(fromObject->hap->pBA));
    SET_VECTOR_ELT(hap, 1, c_dMatrix_to_r_matrix(fromObject->hap->pBa));
    SET_VECTOR_ELT(hap, 2, c_dMatrix_to_r_matrix(fromObject->hap->pbA));
    SET_VECTOR_ELT(hap, 3, c_dMatrix_to_r_matrix(fromObject->hap->pba));
    // now input the hap ind the main list
    SET_VECTOR_ELT(ans,itemNumber++ , hap);
    SET_VECTOR_ELT(ans,itemNumber++ , c_dArray_to_r_array(fromObject->maf));
  }
  
  
  //now input the main names 
  
  itemNumber  = 0 ;
  if(verbose ==1)
    PROTECT(ans_name = allocVector(STRSXP, 26+extras));//added 5
  else
    PROTECT(ans_name = allocVector(STRSXP, 20+extras));//added 5
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kResult"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kLike"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("kr"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("a"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("uLike"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("LD"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("t"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("snp"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("position"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("double_recom"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("alim"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("poLike"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("back"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("chr"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("post"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("timesRun"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("timesConverged"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("convergenceInfo"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("usedSnps"));
  SET_STRING_ELT(ans_name,itemNumber++, mkChar("path"));  
  //extras
  if(!newVersion){
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("mea_global"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("pba_global"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("pBa_global"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("pbA_global"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("pBA_global"));    
    object.dinit_globals(pars->back);
  }


  if(verbose==1){
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("S"));      
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("choose"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("mea"));    
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("S1"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("hap"));
    SET_STRING_ELT(ans_name,itemNumber++, mkChar("maf"));
    //now input the hap names 
    PROTECT(hap_names = allocVector(STRSXP, 4));
    SET_STRING_ELT(hap_names, 0, mkChar("pBA"));
    SET_STRING_ELT(hap_names, 1, mkChar("pBa"));
    SET_STRING_ELT(hap_names, 2, mkChar("pbA"));
    SET_STRING_ELT(hap_names, 3, mkChar("pba"));
    setAttrib(hap, R_NamesSymbol, hap_names);	
  }
    

  //now install the names
  
  setAttrib(ans, R_NamesSymbol, ans_name);
  
  //now set up class name
  PROTECT(class_name = allocVector(STRSXP, 1));
  SET_STRING_ELT(class_name, 0, mkChar("pack.class"));
  classgets(ans, class_name);  
  
  setAttrib(ans, R_NamesSymbol, ans_name);
  if(verbose==0)
    UNPROTECT(3);
  else
    UNPROTECT(5);
  
  //clean up
 
  killFullResults(fromObject);
  killFunctionPars(pars);
  //printf("end of fun\n");
  return(ans) ;
}



