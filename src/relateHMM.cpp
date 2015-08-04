/*
  This is relateHMM
  last saved 22/2 2010  

*/

/*
  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of Relate 0.992
  
  HMMld is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Foobar is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Relate.  If not, see <http://www.gnu.org/licenses/>.
*/  



#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>//used for int to string conversion
#include <sys/stat.h> //used for file handles
#include <cfloat> //only used for getting system max double DBL_MAX
#include <vector> //used for handling joblists
#include <cstdlib>
#include <cstdio>
using namespace std;

#include "conf.h"//contains version number and SYS_MAX


//contains all the type definition
#ifndef _types_h
#define _types_h
#include "types.h"
#endif

//contains all the allocations
#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif


//this is the implementation of the ld calculation
#include "ld.h"
//this is simple variations of extracting from keeplists
#include "extractors.h"

#include "relateHMM.h"

//Optimization is written in ansi c 
extern "C" { 
#include "bfgs.h" 
}

// a std variabale 0=no info 3 alot of info
int  print_info;

/*
  since implementation is implemented as a 
  callback routine in std c
  the arguments must be static and defined in file
*/
dArray *relateHMM::Sk1;
dArray *relateHMM::Sk2;
dArray *relateHMM::Sk3;  
dArray *relateHMM::t;
dArray *relateHMM::alim;
double relateHMM::phi; //used only for calculating a,(so not in likelihood)
int relateHMM::double_recom;



//values used for branching different optimizations
int relateHMM::fixA;
int relateHMM::calcA;
int relateHMM::fixK2;
double relateHMM::fixA_val;
double relateHMM::fixK2_val;
double relateHMM::calcA_val;

//convergence info
int relateHMM::localTimesRun;
int relateHMM::localTimesConv;
dMatrix *relateHMM::convInfo;
pars *relateHMM::calculatedValues;


//these are used by fast all_pairs
iMatrix *relateHMM::genos_global;
dArray *relateHMM::maf1_global;
dArray *relateHMM::maf2_global;
dArray *relateHMM::pos_global;

hapStruct *relateHMM::hap_global;
bArray  *relateHMM::pre_calc_used_list;
#define MAX(x,y) ((x) > (y) ? (x):(y))
#define MIN(x,y) ((x) < (y) ? (x):(y))




/// Appends a matrix to a file,
void write_dMatrix_to_file(string str,const dMatrix *mat,int row,const iArray *keepList){
  /// @param str filename to use.
  /// @param mat The matrix to write.
  /// @param row The row of the matrix to write.
  /// @param keepList to use writing. If keeplist elem is zero then assume missing and write -1.
  ofstream myfile;
  myfile.open (str.c_str(),ios::app);
  int inPlace = 0;
  //  cout <<"size of keeplist:"<<res->keepList->x;
  for (int i=0;i < keepList->x ;i++){//loop through the keeplist
    if (keepList->array[i]==1){ //if snp is to be include write teh result from posterier
      myfile <<setprecision(POSTPRECISION)<< mat->matrix[row][inPlace]<<"\t";
      inPlace++;
    }
    else
      myfile << "-1" << "\t";
  }
  myfile << endl;
  myfile.close();
}

/// Will append an element to a file
void write_dArray_to_file(string str, const dArray *mat,int element){
  ofstream myfile;
  myfile.open (str.c_str(),ios::app);
  myfile <<setprecision(POSTPRECISION)<< mat->array[element] << endl;;
  myfile.close();
}

///apperantly no easy way of int to string conversion
string to_string1(int a){
  string str;
  stringstream ss;
  ss<<a;
  ss>>str;
  return str;
}


/// Returns the first filename that doesn't exist that is prefixed with a given string
string  update_filename(string filename){
  ///@param filename tre prefix filename.
  ///@param return If filename, filename1 and filename2 exists function will then return filename3
  struct stat buffer ;
  if (!stat( filename.c_str(), &buffer ) ){
    int i=0;
    for(;;){
      i++;
      string newname = filename + to_string1(i);
      if (stat( newname.c_str(), &buffer )){
	//print info that tells new filename
	printf("\t-> Filename: %s exists will instead use filename:%s\n",filename.c_str(),newname.c_str());
	return newname;
	break;
      }
    }
  }
    return filename;
}


/*
  myAbs and transpose
  is only used directly from or to the snp_pair_range

 */

void myAbs( dMatrix *matr){
  for (int i=0;i<matr->x;i++)
    for (int j=0;j<matr->y;j++)
      matr->matrix[i][j]=fabs(matr->matrix[i][j]);
}
iMatrix *transpose(const iMatrix *mat){
  iMatrix *retMat = allocIntMatrix(mat->y,mat->x);
  for (int i=0;i<retMat->x;i++)
    for (int j=0;j<retMat->y;j++)
      retMat->matrix[i][j]=mat->matrix[j][i];
  return retMat;
}


snpMatrix *snp_pair_range(iMatrix *v, int  depth){

  snpMatrix *retVal = new snpMatrix();
  int need_signed_r = 0;
  
  int rows = v->x;
  int width =v->y-1; //number of snps to process

  //alloc alle these with 0 as entry
  dMatrix *dprime = allocDoubleMatrix( depth,width,0);
  dMatrix *rmisc  = allocDoubleMatrix( depth,width,0);
  dMatrix *lod    = allocDoubleMatrix( depth,width,0);
  dMatrix *D      = allocDoubleMatrix( depth,width,0);
  dMatrix *pBA      = allocDoubleMatrix(  depth,width,0);
  dMatrix *pBa      = allocDoubleMatrix( depth,width,0);
  dMatrix *pbA      = allocDoubleMatrix(  depth,width,0);
  dMatrix *pba      = allocDoubleMatrix(  depth,width,0);
  iMatrix *transposed=transpose(v);
  //  fillup(rmisc);
  //  fillup(D);
  for(int idx_j = 0; idx_j < depth; idx_j++){
    for (int idx_i = 0  ; idx_i < v->y - 1 - idx_j; idx_i++) {

      geno_cptr res = get_geno_count(transposed->matrix[idx_i],transposed->matrix[idx_i+1+idx_j],rows);
      dprime->matrix[idx_j][idx_i] = res->dprime;
      D->matrix[idx_j][idx_i] = res->bigD/((res->total * 2.0)*(res->total * 2.0));
      pBA->matrix[idx_j][idx_i]  = res->u/(res->total * 2.0);
      pBa->matrix[idx_j][idx_i] = res->v/(res->total * 2.0);
      pbA->matrix[idx_j][idx_i] = res->w/(res->total * 2.0);
      pba->matrix[idx_j][idx_i] = res->x/(res->total * 2.0);
      
      if (need_signed_r) {
	if (res->rsq2 > 0) {
	  rmisc->matrix[idx_j][idx_i] =  res->sign_of_r * sqrt(res->rsq2);

	}else {
	  rmisc->matrix[idx_j][idx_i] = -2;
	}
      } else {

	rmisc->matrix[idx_j][idx_i] =  res->rsq2;

      }

      lod->matrix[idx_j][idx_i] = res->lod;
      free(res->expt);
      free(res);
  
    }
  }
  killMatrix(transposed);
  //collect result and return in a struct
  retVal->dprime=dprime;
  retVal->D=D;
  retVal->pBA=pBA;
  retVal->pBa=pBa;
  retVal->pbA=pbA;
  retVal->pba=pba;
  retVal->rmisc = rmisc;
  retVal->lod = lod;
  return retVal;
}


// return the 1-array, array
dArray *one_minus(const dArray *array){
  dArray *retVal = allocDoubleArray(array->x);
  for(int i=0;i<array->x;i++)
    retVal->array[i] =1.0 - array->array[i];
  return retVal;
}



/*
  sets the values in arg1 to randomvalues.
  ptr[0]:= from alim[0] to alim[1];
  sum(ptr[1],ptr[2],ptr[3])=1
  ptr[3] is a result holder
 */
void get_rand(double *ptr,const  dArray* alim,int fixK2,double fixK2_val){
  ptr[0] =alim->array[0]+alim->array[1]* ((double)rand()/(double)(RAND_MAX));
  if (fixK2)
    ptr[1] = fixK2_val;
  else
    ptr[1] = (double)rand()/(double)(RAND_MAX);
  double tmpSum = ptr[1];
  for (int i=2;i<4;i++){
    int outRand = rand();
    ptr[i] =(1-tmpSum)*( (double)outRand)/(double)(RAND_MAX);
    tmpSum +=ptr[i];
  }
  
  // last is result holder
  ptr[3] = - SYS_MAX;
}

/*
  returns a whichlist from a keeplist
*/
iArray *generateIndices(const iArray *keepList){
  int numTrue =0;
  for(int i=0;i<keepList->x;i++)
    if(keepList->array[i]==1)
      numTrue++;

  iArray *returnArray = allocIntArray(numTrue);
  int atPos = 0;
  for(int i=0;i<keepList->x;i++)
    if(keepList->array[i]==1){
      returnArray->array[atPos] = i;
      atPos++;
    }
  return returnArray;
	
}

//generates a truthtable
iMatrix *getPerm(){
  iMatrix *retVal = allocIntMatrix(81,4);

  
  int id=1;

  for(int i=0;i<3;i++){
    for(int j=0;j<27;j++)
      retVal->matrix[i*27+j][0]=id;
    id++;
  }
 
  id=1;
  for(int i=0;i<9;i++){
    for(int j=0;j<9;j++)
      retVal->matrix[i*9+j][1]=id;
    if(id!=3)
      id++;
    else
      id=1;
  } 

  id=1;
  for(int i=0;i<27;i++){
    for(int j=0;j<3;j++)
      retVal->matrix[i*3+j][2]=id;
    if(id!=3)
      id++;
    else
      id=1;
  } 

  id=1;
  for(int i=0;i<81;i++){
    retVal->matrix[i][3]=id;
    if(id!=3)
      id++;
    else
      id=1;
  }



  return retVal;
}

/*
  returns a keeplist that tells which snps to keep after maf checking
 */
iArray *getOKindices(iMatrix *matrix,dArray *mafArray,double min,int ind1, int ind2,iArray *keepList){
  //  int *keepList= new int[matrix->y];
  int numsToKeep=0;
  //printf("MINUMIMUN IN MAF STRIPPING IS:%f\n",min);
  for (int i=0;i<matrix->y;i++){
    //    if(matrix->matrix[ind1][i]==0 || matrix->matrix[ind2][i]==0 || mafArray->array[i]<min || mafArray->array[i]>(1-min) )
    if(matrix->matrix[ind1][i]!=0 && matrix->matrix[ind2][i]!=0 && mafArray->array[i]>min && mafArray->array[i]<(1-min) && mafArray->array[i]!=0 ){
      keepList->array[i]=1;
      numsToKeep++;
    }
    else
      keepList->array[i]=0;
  }
  
  iArray *returnArray = allocIntArray(numsToKeep);
  int inPlace=0;
  for (int i =0 ;i<matrix->y;i++)
    if(keepList->array[i]==1){
      returnArray->array[inPlace]=i;
      inPlace++;
    }
  //  delete [] keepList;
  return returnArray;
}

//returns the pairwise difference between elements
dArray *diff(dArray *var){
  dArray *retVal = allocDoubleArray(var->x+1);
  retVal->array[0]=0;
  retVal->array[retVal->x-1]=0;
  double myDifference =0;
  for (int i=1;i<var->x;i++){
    myDifference = var->array[i]-var->array[i-1];
    if (myDifference < 0)
      retVal->array[i]= 10000000.0 ;
    else
      retVal->array[i]= myDifference;
  }
  return retVal;
}

/*
  row_log returns the log of a matrix

 */
dArray *row_log(const dMatrix *var,int rowNum){
  dArray *retVal = allocDoubleArray(var->y);
  for (int i=0;i<var->y;i++)
    retVal->array[i] = log(var->matrix[rowNum][i]);
  return retVal;
}

dArray *row_log(const dArray *var){
  dArray *retVal = allocDoubleArray(var->x);
  for (int i=0;i<var->x;i++)
    retVal->array[i] = log(var->array[i]);
  return retVal;
}

//returns a keeplist of snps that won't have any NA
bArray *which_arent_NA(iMatrix *geno, int row1, int row2){
  bArray *retVal = allocBoolArray(geno->y);
  int numTrues= 0;
  for (int i=0;i<geno->y ; i++)
    if(geno->matrix[row1][i]!=0 && geno->matrix[row2][i]!=0 ){
      retVal->array[i] = 1;
      numTrues++;
    }
    else
      retVal->array[i] = 0;
  retVal->numTrue = numTrues;
  return retVal;
}



iMatrix *revCols(const iMatrix *matr){
  iMatrix *retVal = allocIntMatrix(matr->x,matr->y);
  
  for (int x=0;x<matr->x;x++){
    for (int i=0;i<matr->y;i++){
      retVal->matrix[x][i]=matr->matrix[x][matr->y-1-i];
    }
  }
  return retVal;
}

dMatrix *revCols(const dMatrix *matr){
  dMatrix *retVal = allocDoubleMatrix(matr->x,matr->y);
  
  for (int x=0;x<matr->x;x++){
    for (int i=0;i<matr->y;i++){
      retVal->matrix[x][i]=matr->matrix[x][matr->y-1-i];
    }
  }
  return retVal;
}

///This function will rbind a vector of zeroes, and will reverse the rest.
dMatrix *revCols_and_extend(dMatrix *matr){
  dMatrix *retVal = allocDoubleMatrix(matr->x,matr->y+1);
  
  //fill in zeroes as needed
  for (int i=0;i<matr->x;i++)
    retVal->matrix[i][0]=0;
  
  //now reverse order of snp
  for (int x=0;x<matr->x;x++){
    for (int i=0;i<matr->y;i++){
      retVal->matrix[x][i+1]=matr->matrix[x][matr->y-1-i];
    }
  }
  return retVal;
}

///Will return an vector array that tells which element in a column that has the heighest value
iArray *getHighestId(dMatrix *matr){
  iArray *retVal = allocIntArray(matr->y+1);
  retVal->array[0] = -1;
  if(matr->x==1){
    for(int i=1;i<retVal->x;i++)
      retVal->array[i]=0;
  }
  else
    for (int s=0;s<matr->y;s++){
      int highId =0;
      for (int d=1;d<matr->x;d++)
	if(matr->matrix[d][s] > matr->matrix[highId][s]) //changed in 0.91
	  highId=d;
      retVal->array[1+s]=highId;
    }
  return retVal;
}



dArray *dbl_recom(const dArray *t,const double a,const  double k,int lastFactor){
  dArray *retVal = allocDoubleArray(t->x);
  for (int i=0;i<t->x;i++)
    retVal->array[i]=(1.0-exp(-a*t->array[i]))*k+lastFactor*(exp(-a*t->array[i]));
  return retVal;
  
}

dArray *math_to(double factor1,const dArray *one,double factor2,const  dArray *two,double translate){
  dArray *retVal = allocDoubleArray(one->x);
  for (int i=0;i<retVal->x;i++)
    retVal->array[i]=one->array[i]*factor1+two->array[i]*factor2+translate;
  return retVal; 
}

dArray *copy(const dArray *in){
  dArray *retVal = allocDoubleArray(in->x);
  for (int i=0;i<retVal->x;i++)
    retVal->array[i]=in->array[i];
  return retVal;   
}

dArray *scale(dArray *in,double sc){
  dArray *retVal = allocDoubleArray(in->x);
  for (int i=0;i<retVal->x;i++)
    retVal->array[i]=sc*(in->array[i]);
  return retVal;   
}

dMatrix *remove_column(dMatrix *mat, int split){
  dMatrix *retVal = allocDoubleMatrix(mat->x,mat->y-1);
  for (int i=0;i<split;i++) //copy left values
    for(int j=0;j<mat->x;j++)
      retVal->matrix[j][i] = mat->matrix[j][i];
  for (int i=split;i<retVal->y;i++) //copy right values
    for(int j=0;j<mat->x;j++)
      retVal->matrix[j][i] = mat->matrix[j][i+1];
  killMatrix(mat);
  return retVal;
}

void fix_underflow(dArray *in){
  for(int i=0;i<in->x;i++)
    if(in->array[i]<0.000000000000001){
      in->array[i]=0;
    }
}

iArray *pruning(snpMatrix *ld,int ld_choose,int back,double prune_val){
  if(print_info>1){
    printf("\r");
    flush_print("pruning datastructures with prune_val:",prune_val);
  }

  dMatrix *mea;
  if(ld_choose)
    mea = revCols(ld->rmisc);
  else{
    mea = revCols(ld->D);
    myAbs(mea); 
  }
    

  
  /*
    the procedure is strange,
    if column sum is greater the pruneValue,
    then remove the entire SNP locus.
    then do a SHIFTUP operation on the next BACK-1 loci.
    SHIFT defined as removing the diagonal, and moving the rest up and appending zero
    EXAMPLE:
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]    1    6   11   16   21   26   31   36   41    46
[2,]    2    7   12   17   22   27   32   37   42    47
[3,]    3    8   13   18   23   28   33   38   43    48
[4,]    4    9   14   19   24   29   34   39   44    49
[5,]    5   10   15   20   25   30   35   40   45    50
   removing SNP 3 will turn the above to

     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
[1,]    1    6   16   21   26   31   36   41   46
[2,]    2    7   18   22   27   32   37   42   47
[3,]    3    8   19   24   28   33   38   43   48
[4,]    4    9   20   25   30   34   39   44   49
[5,]    5   10    0    0    0    0   40   45   50


  */
	      
  
  int snp = mea->y+1; //number of snps before process;
  iArray *retVal = allocIntArray(ld->rmisc->y+1); //snp length plus one
  retVal->array[0] = 1; // always keep first
  if(back>1) {
    for(int nsnp=1;nsnp<snp;nsnp++) {//iterate thourgh SNP's
      int doStripping=0;
      for(int j=0;j<back;j++)
	if( mea->matrix[j][nsnp-1] >prune_val){
	  doStripping=1;
	  break;//jump out of inner loop
	}
      if(doStripping) {
	//should update datastructure.
	retVal->array[nsnp] = 0; //col will be excluded from the keeplist
	if(back>2){
	  for (int p=0 ; p < back-1 ; p++){//number of colums to update
	    if(p+nsnp>=mea->y) //if no more colums exists, exit
	      break;
			    
	    for(int q=p+1;q<back-1;q++){
	      
	      mea->matrix[q][nsnp+p] = mea->matrix[q+1][p+nsnp];
	    }
	    
	    mea->matrix[back-1][nsnp+p] = 0;//input zero at bottom    
	  }
	}else //back must =2 so just remove 
	  
	  mea->matrix[1][nsnp] = 0;//shouldn't it be the next one?
      }
      else {
	
	//we shouldn't remove this snp so set 1 in keeplist
	retVal->array[nsnp] = 1;
      }
    }
  }
	

  else{
    //back is 1 so just check if less than prunevalue
    for (int i=1;i<retVal->x;i++)
      if(mea->matrix[0][i-1] <= prune_val)
	retVal->array[i] = 1;
  }

  if(print_info>1)
    flush_print("end of pruning\n");
  
  killMatrix(mea);
  return retVal;
}


/*
  returns the absolute value of the entries of a dMatrix
  Only used for setting mea=abs(mea)
*/


dMatrix *getMea(iMatrix *dat,int ld_choose,int back){
  dMatrix *mea;
  iMatrix *tmp = revCols(dat);
  snpMatrix *ld = snp_pair_range(tmp,back);
    
  if(ld_choose)
    mea = revCols(ld->rmisc);
  else{
    mea = revCols(ld->D);
    myAbs(mea); 
  }
  killMatrix(tmp);
  killSnpMatrix(ld);
  return mea;

}


bArray *pruning3(dMatrix *mea,int back,double prune_val){
  if(print_info>1){
    flush_print("\r");
    flush_print("Starting pruning of datastructures with prune_val:",prune_val);
  }
  int numTrue=0;
  
  /*
    the procedure is strange,
    if column sum is greater the pruneValue,
    then remove the entire SNP locus.
    then do a SHIFTUP operation on the next BACK-1 loci.
    SHIFT defined as removing the diagonal, and moving the rest up and appending zero
    EXAMPLE:
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]    1    6   11   16   21   26   31   36   41    46
[2,]    2    7   12   17   22   27   32   37   42    47
[3,]    3    8   13   18   23   28   33   38   43    48
[4,]    4    9   14   19   24   29   34   39   44    49
[5,]    5   10   15   20   25   30   35   40   45    50
   removing SNP 3 will turn the above to

     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
[1,]    1    6   16   21   26   31   36   41   46
[2,]    2    7   18   22   27   32   37   42   47
[3,]    3    8   19   24   28   33   38   43   48
[4,]    4    9   20   25   30   34   39   44   49
[5,]    5   10    0    0    0    0   40   45   50


  */
	      
  
  int snp = mea->y+1; //number of snps before process;
  bArray *retVal = allocBoolArray(mea->y+1); //snp length plus one
  retVal->array[0] = 1; // always keep first
  numTrue++;
  if(back>1) {
    for(int nsnp=1;nsnp<snp;nsnp++) {//iterate thourgh SNP's
      int doStripping=0;
      for(int j=0;j<back;j++)
	if( mea->matrix[j][nsnp-1] >prune_val){
	  doStripping=1;
	  break;//jump out of inner loop
	}
      if(doStripping) {
	//should update datastructure.
	retVal->array[nsnp] = 0; //col will be excluded from the keeplist
	if(back>2){
	  for (int p=0 ; p < back-1 ; p++){//number of colums to update
	    if(p+nsnp>=mea->y) //if no more colums exists, exit
	      break;
	    
	    for(int q=p+1;q<back-1;q++){
	      
	      mea->matrix[q][nsnp+p] = mea->matrix[q+1][p+nsnp];
	    }
	    
	    mea->matrix[back-1][nsnp+p] = 0;//input zero at bottom    
	  }
	}else //back must =2 so just remove 
	  
	  mea->matrix[1][nsnp] = 0;
      }
      else {
	//we shouldn't remove this snp so set 1 in keeplist
	retVal->array[nsnp] = 1;
	numTrue++;
      }
    }
  }
	

  else{
    //back is 1 so just check if less than prunevalue
    for (int i=1;i<retVal->x;i++)
      if(mea->matrix[0][i-1] <= prune_val){
	retVal->array[i] = 1;
	numTrue++;
      }
  }

  if(print_info>1)
    printf("\t end of pruning\n");
  retVal->numTrue = numTrue;
  killMatrix(mea);
  return retVal;
}


dMatrix *remove_cross_original(const dMatrix *mea,const bArray *keepList,int back){
  if(mea->y!=keepList->x-1){
    printf("keeplist in remove_cross should be of length plus 1");
    exit(0);
  }
  if(keepList->array[0]!=1){
    printf("First element in keeplist should be TRUE\n");
    exit(0);
  }

  dMatrix *retMat = allocDoubleMatrix(back,keepList->numTrue-1);
  dMatrix *tmp = allocDoubleMatrix(mea->x,mea->y);
  
  //copy mea into tmp
  for (int i=0;i<mea->x;i++)
    for(int j=0;j<mea->y;j++)
      tmp->matrix[i][j] = mea->matrix[i][j];
  
  if(print_info>2)
    printf("\t-> Start of remove_cross: dims of matrix:(%d,%d) length of keeplist:%d \n",mea->x,mea->y,keepList->x);
  
  /*
    the procedure is strange,
    Do a SHIFTUP operation on the next BACK-1 loci.
    SHIFT defined as removing the diagonal, and moving the rest up and appending zero
    EXAMPLE:
    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    [1,]    1    6   11   16   21   26   31   36   41    46
    [2,]    2    7   12   17   22   27   32   37   42    47
    [3,]    3    8   13   18   23   28   33   38   43    48
    [4,]    4    9   14   19   24   29   34   39   44    49
    [5,]    5   10   15   20   25   30   35   40   45    50
   removing SNP 3 will turn the above to

   [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
   [1,]    1    6   16   21   26   31   36   41   46
   [2,]    2    7   18   22   27   32   37   42   47
   [3,]    3    8   19   24   28   33   38   43   48
   [4,]    4    9   20   25   30   34   39   44   49
   [5,]    5   10    0    0    0    0   40   45   50
  */
	      
  if(mea->x>1) { // if back > 1
    int inPlace = 0;
    for(int nsnp=0;nsnp<mea->y;nsnp++) {//iterate thourgh SNP's ///COKE ROX
      
      if(keepList->array[nsnp+1]==1) {//we should keep this loci so copy it
	for(int d=0 ; d<back ; d++)
	  retMat->matrix[d][inPlace] = tmp->matrix[d][nsnp];
	//printf("inplace:%d\t nsnp:%d\n",inPlace,nsnp);
	inPlace++;
	
      }else{
	//printf("\t-> nsnp:%d\n",nsnp);
	//update structures
	
	if(mea->x >2){//back >2
	  for (int p=0 ; p <mea->x ; p++){//number of colums to update
	    if(p+1+nsnp>=mea->y) //if no more colums exists, exit
	      break;
	    
	    
	    for(int q=p;q<mea->x -1;q++){
	      //s	      printf("%d %d %d\n",q,p,nsnp);
	      tmp->matrix[q][1+nsnp+p] = tmp->matrix[q+1][1+p+nsnp];
	    }
	    tmp->matrix[mea->x-1][1+nsnp+p] = 0;//input zero at bottom    
	  }
	  // tmp->matrix[back-1][nsnp+back];
	}else
	  //back must =2 so just remove 
	  tmp->matrix[1][nsnp+1] = 0;
      }
    }
  }
  else{
    //back is 1, so we only need to copy values
    int inPlace = 0;
    for (int i=0 ; i<mea->y ; i++)
      if( keepList->array[i+1]) {
	retMat->matrix[0][inPlace] = tmp->matrix[0][i];
	inPlace++;
      }else
	retMat->matrix[0][inPlace] = 0;
  }

  //  if(print_info>1)
  // printf("\t end of remove_cross()\n");
  killMatrix(tmp);
  return retMat;
}
dMatrix *remove_cross(const dMatrix *mea,const bArray *keepList,int back){
  if(mea->y!=keepList->x-1){
    printf("keeplist in remove_cross should be of length plus 1");
    exit(0);
  }
  /*
    There are 2 cases that should be handled somewhat differently
    1. if first element in keeplist zero
    1. a if also the second element in keeplist zero
    2. otherwise
   */
  
  dMatrix *retMat = allocDoubleMatrix(back,keepList->numTrue-1);
  dMatrix *tmp = allocDoubleMatrix(mea->x,mea->y);
  int snpStart = 0;  //this should be the start loci
  //copy mea into tmp
  for (int i=0;i<mea->x;i++)
    for(int j=0;j<mea->y;j++)
      tmp->matrix[i][j] = mea->matrix[i][j];
  
  if(keepList->array[0]==0){
    /*
      first we will handle the the sick case.
      it shouldn't matter much for the speed since 
      this conditional will only happen once.
    */
    for (int p=1 ; p <mea->x ; p++){//number of colums to update
      if(p+1>=mea->y) //if no more colums exists, exit
	break;
      for(int q=p;q<mea->x -1;q++){
	tmp->matrix[q][p] = tmp->matrix[q+1][p];
      }
      tmp->matrix[mea->x-1][p] = 0;//input zero at bottom    
    }
    snpStart++;
    if(keepList->array[1]==0) { 
      for (int p=1 ; p <mea->x ; p++){//number of colums to update
	if(p+1>=mea->y) //if no more colums exists, exit
	  break;
	for(int q=p;q<mea->x -1;q++){
	  tmp->matrix[q][p+1] = tmp->matrix[q+1][p+1];
	}
	tmp->matrix[mea->x-1][p+1] = 0;//input zero at bottom    
      }
      snpStart++;
    } 
  }


  if(print_info>2)
    printf("\t-> Start of remove_cross: dims of matrix:(%d,%d) length of keeplist:%d \n",mea->x,mea->y,keepList->x);
  
  /*
    the procedure is strange,
    Do a SHIFTUP operation on the next BACK-1 loci.
    SHIFT defined as removing the diagonal, and moving the rest up and appending zero
    EXAMPLE:
    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    [1,]    1    6   11   16   21   26   31   36   41    46
    [2,]    2    7   12   17   22   27   32   37   42    47
    [3,]    3    8   13   18   23   28   33   38   43    48
    [4,]    4    9   14   19   24   29   34   39   44    49
    [5,]    5   10   15   20   25   30   35   40   45    50
   removing SNP 3 will turn the above to

   [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
   [1,]    1    6   16   21   26   31   36   41   46
   [2,]    2    7   18   22   27   32   37   42   47
   [3,]    3    8   19   24   28   33   38   43   48
   [4,]    4    9   20   25   30   34   39   44   49
   [5,]    5   10    0    0    0    0   40   45   50
  */

  if(mea->x>1) { // if back > 1
    int inPlace = 0;
    for(int nsnp=snpStart;nsnp<mea->y;nsnp++) {//iterate thourgh SNP's ///COKE ROX
      
      if(keepList->array[nsnp+1]==1) {//we should keep this loci so copy it
	for(int d=0 ; d<back ; d++)
	  retMat->matrix[d][inPlace] = tmp->matrix[d][nsnp];
	//printf("inplace:%d\t nsnp:%d\n",inPlace,nsnp);
	inPlace++;
	
      }else{
	//printf("\t-> nsnp:%d\n",nsnp);
	//update structures
	
	if(mea->x >2){//back >2
	  for (int p=0 ; p <mea->x ; p++){//number of colums to update
	    if(p+1+nsnp>=mea->y) //if no more colums exists, exit
	      break;
	    
	    
	    for(int q=p;q<mea->x -1;q++){
	      //s	      printf("%d %d %d\n",q,p,nsnp);
	      tmp->matrix[q][1+nsnp+p] = tmp->matrix[q+1][1+p+nsnp];
	    }
	    tmp->matrix[mea->x-1][1+nsnp+p] = 0;//input zero at bottom    
	  }
	  // tmp->matrix[back-1][nsnp+back];
	}else
	  //back must =2 so just remove 
	  tmp->matrix[1][nsnp+1] = 0;
      }
    }
  }
  else{
    //back is 1, so we only need to copy values
    int inPlace = 0;
    for (int i=0 ; i<mea->y ; i++)
      if( keepList->array[i+1]) {
	retMat->matrix[0][inPlace] = tmp->matrix[0][i];
	inPlace++;
      }else
	retMat->matrix[0][inPlace] = 0;
  }

  //  if(print_info>1)
  // printf("\t end of remove_cross()\n");
  killMatrix(tmp);
  return retMat;
}


double calculateA(double k0,double k1, double k2,double phi){
  double ma,mb,xa,xb,m,a,sq,pw;
  if(k2==0){
    ma = 1-log(k1)/log(2);
    mb = 0;
  }
  else{
    pw = pow(k1+2*k2,2);
    sq = sqrt(pw-4*k2);
    xa = (k1+2*k2+sq)/2.0;
    xb = k2/xa;
    ma = 1-log(xa)/log(2);
    mb = 1-log(xb)/log(2);
  }
  m = ma+mb;
  a = -m*log(1-phi);
  
  if(std::isnan(a)&& 0){//there is a bug in the bfgs algortihm, that makes the optim algo try outside of parameter space
    printf("m=%f , ma=%f ,xa=%f, sq=%f  mb=%f  , a=%f  k0=%f ,  k1=%f   , k2=%f \n",m,ma,xa,sq,mb,a,k0,k1,k2);
    printf("calc.a(k0=%f,k1=%f,k2=%f)\n",k0,k1,k2);
    printf("std::isnan in calculate.a\n");
  }
  /*
  else if(std::isinf(a)){
    printf("std::isinf in calculate.a -> k0=%f,k1=%f,k2=%f,phi=%f\n",k0,k1,k2,phi);
    if(k2==0)
      printf("k2=0,\t");
    printf("ma=%f,mb=%f,m=%f\n",ma,mb,m);
	

  }
  */
  return a;
}


// added in 0.987 20/2 2009 because of underflow error
//dMatrix* decode(decodePars *var){
dMatrix* decodek2(const dArray *ptr,const dArray *Sk1,const dArray *Sk2,const dArray *Sk3,const dArray *t,int double_recom){
  //  printf("will runk2 decode\n");
  int snp = Sk1->x;
  double a  = ptr->array[0];
  double k2 = ptr->array[1];
  double k1 = ptr->array[2];
  double k0 = 1-k2-k1;

  dArray *IBD00,*IBD01,*IBD02,*IBD10,*IBD11,*IBD12,*IBD20,*IBD21,*IBD22;
  if(double_recom){
    printf("assuming double recom\n");
    IBD00 = dbl_recom(t,a,k0,1);
    IBD01 = dbl_recom(t,a,k1,0);
    IBD02 = dbl_recom(t,a,k2,0);
    IBD10 = dbl_recom(t,a,k0,0);
    IBD11 = dbl_recom(t,a,k1,1);
    IBD12 = dbl_recom(t,a,k2,0);
    IBD20 = dbl_recom(t,a,k0,0);
    IBD21 = dbl_recom(t,a,k1,0);
    IBD22 = dbl_recom(t,a,k2,1);
  }
  else{
    IBD01 = dbl_recom(t,a,k1,0);//ok
    IBD02 = allocDoubleArray(t->x);

    if(k1!=1)
      for (int i=0;i<IBD02->x;i++)
	IBD02->array[i]= exp(-a*k1*t->array[i])*k2/(k1-1)+exp(-a*t->array[i])*k1+exp(-a*t->array[i])*k0*k1/(k1-1)+k2;//ok

    fix_underflow(IBD02);

    IBD00 = math_to(-1,IBD01,-1,IBD02,1);//ok

    if(k1==0){
      IBD10 = copy(IBD01);//ok
      IBD12 = copy(IBD01);//ok
    } else{
      IBD10 = scale(IBD01,k0/k1);//ok
      IBD12 = scale(IBD01,k2/k1);//ok
    }
    
    IBD11 = math_to(-1,IBD10,-1,IBD12,1);//ok
    IBD21 = copy(IBD01);
    IBD20 = allocDoubleArray(t->x);
    if(k1!=1)
      for (int i=0;i<t->x;i++)
	IBD20->array[i]= exp(-a*k1*t->array[i])*k0/(k1-1.0)+exp(-a*t->array[i])*k1+exp(-a*t->array[i])/(k1-1)*k2*k1-(1-k1)*k0/(k1-1);
    fix_underflow(IBD20);
    IBD22 = math_to(-1,IBD21,-1,IBD20,1);
  }

  //start Forward
  dArray *log_IBD0=allocDoubleArray(snp);
  dArray *log_IBD1=allocDoubleArray(snp);
  dArray *log_IBD2=allocDoubleArray(snp);

  log_IBD0->array[0] = log(k0)+ Sk3->array[0];
  log_IBD1->array[0] = log(k1)+ Sk2->array[0];
  log_IBD2->array[0] = log(k2)+ Sk1->array[0];

  double k=0;
  double logp0,logp1,logp2;
  for (int i=1;i<log_IBD0->x;i++){
    logp0 = exp(log_IBD0->array[i-1]-k);
    logp1 = exp(log_IBD1->array[i-1]-k);
    logp2 = exp(log_IBD2->array[i-1]-k);
    log_IBD0->array[i] = k+log(IBD00->array[i]*logp0+IBD10->array[i]*logp1+IBD20->array[i]*logp2)+Sk3->array[i];
    log_IBD1->array[i] = k+log(IBD01->array[i]*logp0+IBD11->array[i]*logp1+IBD21->array[i]*logp2)+Sk2->array[i];
    log_IBD2->array[i] = k+log(IBD02->array[i]*logp0+IBD12->array[i]*logp1+IBD22->array[i]*logp2)+Sk1->array[i];
    k = MAX(log_IBD0->array[i],MAX(log_IBD1->array[i],log_IBD2->array[i]));

  }
  logp0 = log_IBD0->array[snp-1];
  logp1 = log_IBD1->array[snp-1];
  logp2 = log_IBD2->array[snp-1];
  k = MAX(logp0, MAX(logp1, logp2));
  double l = k+log(exp(logp0-k) + exp(logp1-k)+ exp(logp2-k));
 
  dArray *bk_log_IBD0 = allocDoubleArray(snp);
  dArray *bk_log_IBD1 = allocDoubleArray(snp);
  dArray *bk_log_IBD2 = allocDoubleArray(snp);
  bk_log_IBD0->array[snp-1] = 0 + log(k0);
  bk_log_IBD1->array[snp-1] = 0 + log(k1);
  bk_log_IBD2->array[snp-1] = 0 + log(k2);
  k=0;
  for (int i=snp-2;i>=0;i--){
    logp0 = exp(Sk3->array[i+1] + bk_log_IBD0->array[i+1]-k);
    logp1 = exp(Sk2->array[i+1] + bk_log_IBD1->array[i+1]-k);
    logp2 = exp(Sk1->array[i+1] + bk_log_IBD2->array[i+1]-k);
    bk_log_IBD0->array[i] = k+log(IBD00->array[i+1]*logp0+IBD10->array[i+1]*logp1+IBD20->array[i+1]*logp2);
    bk_log_IBD1->array[i] = k+log(IBD01->array[i+1]*logp0+IBD11->array[i+1]*logp1+IBD21->array[i+1]*logp2);
    bk_log_IBD2->array[i] = k+log(IBD02->array[i+1]*logp0+IBD12->array[i+1]*logp1+IBD22->array[i+1]*logp2);
    k = MAX(bk_log_IBD2->array[i],MAX(bk_log_IBD1->array[i],bk_log_IBD0->array[i]));
  }
  logp0 =  bk_log_IBD0->array[1]+ Sk3->array[0];
  logp1 =  bk_log_IBD1->array[1]+ Sk2->array[0];
  logp2 =  bk_log_IBD2->array[1]+ Sk1->array[0];
  k = MAX(logp0, MAX(logp1, logp2));
  
  bk_log_IBD0->array[snp-1] = 0;
  bk_log_IBD1->array[snp-1] = 0;
  bk_log_IBD2->array[snp-1] = 0;
  k=0;
  double p0,p1,p2;
  for (int i=snp-2;i>=0;i--){

    p0 = exp(Sk3->array[i+1] + bk_log_IBD0->array[i+1]-k);
    p1 = exp(Sk2->array[i+1] + bk_log_IBD1->array[i+1]-k);
    p2 = exp(Sk1->array[i+1] + bk_log_IBD2->array[i+1]-k);
    bk_log_IBD0->array[i] = k+log(p2*IBD02->array[i+1]+p1*IBD01->array[i+1]+p0*IBD00->array[i+1]);
    bk_log_IBD1->array[i] = k+log(p2*IBD12->array[i+1]+p1*IBD11->array[i+1]+p0*IBD10->array[i+1]);
    bk_log_IBD2->array[i] = k+log(p2*IBD22->array[i+1]+p1*IBD21->array[i+1]+p0*IBD20->array[i+1]);
    k = MAX(bk_log_IBD2->array[i],MAX(bk_log_IBD1->array[i],bk_log_IBD0->array[i]));
  }

  logp0 =  bk_log_IBD0->array[0]+ Sk3->array[0]+log(k0);
  logp1 =  bk_log_IBD1->array[0]+ Sk2->array[0]+log(k1);
  logp2 =  bk_log_IBD2->array[0]+ Sk1->array[0]+log(k2);
  
  k = MAX(logp0, MAX(logp1, logp2));
  dMatrix *retMat = allocDoubleMatrix(3,snp);
  for (int i=0;i<snp;i++){
    retMat->matrix[0][i] = exp(log_IBD2->array[i]+bk_log_IBD2->array[i]-l);
    retMat->matrix[1][i] = exp(log_IBD1->array[i]+bk_log_IBD1->array[i]-l);
    retMat->matrix[2][i] = exp(log_IBD0->array[i]+bk_log_IBD0->array[i]-l);
  }
  
  //now clean up
  killArray(bk_log_IBD0);
  killArray(bk_log_IBD1);
  killArray(bk_log_IBD2);
  killArray(IBD00);
  killArray(IBD01);
  killArray(IBD02);
  killArray(IBD10);
  killArray(IBD11);
  killArray(IBD12);
  killArray(IBD20);
  killArray(IBD21);
  killArray(IBD22);
  killArray(log_IBD0);
  killArray(log_IBD1);
  killArray(log_IBD2);

  return retMat;
}




//added on 0.987 20/2 2009
dMatrix* decode(const dArray *ptr,const dArray *Sk1,const dArray *Sk2,const dArray *Sk3,const dArray *t,int double_recom){
  //  printf("will runk2=0 decode\n");
  int snp = Sk1->x;
  double a  = ptr->array[0];
  double k1 = ptr->array[2];
  double k0 = 1-0-k1;

  dArray *IBD00,*IBD01,*IBD10,*IBD11;
  if(double_recom){
    IBD00 = dbl_recom(t,a,k0,1);
    IBD01 = dbl_recom(t,a,k1,0);

    IBD10 = dbl_recom(t,a,k0,0);
    IBD11 = dbl_recom(t,a,k1,1);
  }
  else{
    IBD01 = dbl_recom(t,a,k1,0);//ok
    IBD00 = one_minus(IBD01);//math_to(-1,IBD01,-1,IBD02,1);//ok

    if(k1==0){
      IBD10 = copy(IBD01);//ok
    } else{
      IBD10 = scale(IBD01,k0/k1);//ok
    }
    
    IBD11 =one_minus(IBD10);// math_to(-1,IBD10,-1,IBD12,1);//ok
  }

  //start Forward
  dArray *log_IBD0=allocDoubleArray(snp);
  dArray *log_IBD1=allocDoubleArray(snp);
  
  log_IBD0->array[0] = log(k0)+ Sk3->array[0];
  log_IBD1->array[0] = log(k1)+ Sk2->array[0];
  
  double k=0;
  double logp0,logp1;//,logp2;
  for (int i=1;i<log_IBD0->x;i++){
    logp0 = exp(log_IBD0->array[i-1]-k);
    logp1 = exp(log_IBD1->array[i-1]-k);
    log_IBD0->array[i] = k+log(IBD00->array[i]*logp0+IBD10->array[i]*logp1)+Sk3->array[i];
    log_IBD1->array[i] = k+log(IBD01->array[i]*logp0+IBD11->array[i]*logp1)+Sk2->array[i];
    k = MAX(log_IBD0->array[i],log_IBD1->array[i]);//,log_IBD2->array[i]));

  }
  logp0 = log_IBD0->array[snp-1];
  logp1 = log_IBD1->array[snp-1];

  k = MAX(logp0,logp1);

  double l = k+log(exp(logp0-k) + exp(logp1-k));
 
  dArray *bk_log_IBD0 = allocDoubleArray(snp);
  dArray *bk_log_IBD1 = allocDoubleArray(snp);

  bk_log_IBD0->array[snp-1] = 0 + log(k0);
  bk_log_IBD1->array[snp-1] = 0 + log(k1);
  k=0;
  for (int i=snp-2;i>=0;i--){
    logp0 = exp(Sk3->array[i+1] + bk_log_IBD0->array[i+1]-k);
    logp1 = exp(Sk2->array[i+1] + bk_log_IBD1->array[i+1]-k);
    bk_log_IBD0->array[i] = k+log(IBD00->array[i+1]*logp0+IBD10->array[i+1]*logp1);
    bk_log_IBD1->array[i] = k+log(IBD01->array[i+1]*logp0+IBD11->array[i+1]*logp1);
    k = MAX(bk_log_IBD1->array[i],bk_log_IBD0->array[i]);
  }
  logp0 =  bk_log_IBD0->array[1]+ Sk3->array[0];
  logp1 =  bk_log_IBD1->array[1]+ Sk2->array[0];
  k = MAX(logp0, logp1);
  
  bk_log_IBD0->array[snp-1] = 0;
  bk_log_IBD1->array[snp-1] = 0;
  k=0;
  double p0,p1;
  for (int i=snp-2;i>=0;i--){
    p0 = exp(Sk3->array[i+1] + bk_log_IBD0->array[i+1]-k);
    p1 = exp(Sk2->array[i+1] + bk_log_IBD1->array[i+1]-k);
    bk_log_IBD0->array[i] = k+log(p1*IBD01->array[i+1]+p0*IBD00->array[i+1]);
    bk_log_IBD1->array[i] = k+log(p1*IBD11->array[i+1]+p0*IBD10->array[i+1]);
    k = MAX(bk_log_IBD1->array[i],bk_log_IBD0->array[i]);
  }

  logp0 =  bk_log_IBD0->array[0]+ Sk3->array[0]+log(k0);
  logp1 =  bk_log_IBD1->array[0]+ Sk2->array[0]+log(k1);
  k = MAX(logp0, logp1);
  dMatrix *retMat = allocDoubleMatrix(3,snp);
  for (int i=0;i<snp;i++){
    retMat->matrix[1][i] = exp(log_IBD1->array[i]+bk_log_IBD1->array[i]-l);
    retMat->matrix[2][i] = exp(log_IBD0->array[i]+bk_log_IBD0->array[i]-l);
  }
  
  //now clean up
  killArray(bk_log_IBD0);
  killArray(bk_log_IBD1);
  killArray(IBD00);
  killArray(IBD01);
  killArray(IBD10);
  killArray(IBD11);
  killArray(log_IBD0);
  killArray(log_IBD1);
  return retMat;
}



/*
  max3 and which_max are used only by viterbi function

 */
/// Extract's the max of 3 elements.
double max3(double a, double b, double c){
 if(a>b && a>c)
    return a;
  else if(b>a && b>c)
    return b;
  else 
    return c;
}
int which_max(double a,double b,double c){
  if(a>b && a>c)
    return 0;
  else if(b>a && b>c)
    return 1;
  return 2;
}


iArray* viterbi(const dArray *ptr,const dArray *Sk1,const dArray *Sk2,const dArray *Sk3,const dArray *t,int double_recom){
 
  int snp = Sk1->x;
  double a  = ptr->array[0];
  double k2 = ptr->array[1];
  double k1 = ptr->array[2];
  double k0 = 1-k2-k1;


  dArray *IBD00,*IBD01,*IBD02,*IBD10,*IBD11,*IBD12,*IBD20,*IBD21,*IBD22;
  if(double_recom){
    IBD00 = dbl_recom(t,a,k0,1);
    IBD01 = dbl_recom(t,a,k1,0);
    IBD02 = dbl_recom(t,a,k2,0);
    IBD10 = dbl_recom(t,a,k0,0);
    IBD11 = dbl_recom(t,a,k1,1);
    IBD12 = dbl_recom(t,a,k2,0);
    IBD20 = dbl_recom(t,a,k0,0);
    IBD21 = dbl_recom(t,a,k1,0);
    IBD22 = dbl_recom(t,a,k2,1);
  }
  else{
    IBD01 = dbl_recom(t,a,k1,0);//ok
    IBD02 = allocDoubleArray(t->x);

    if(k1!=1)
      for (int i=0;i<IBD02->x;i++)
	IBD02->array[i]= exp(-a*k1*t->array[i])*k2/(k1-1)+exp(-a*t->array[i])*k1+exp(-a*t->array[i])*k0*k1/(k1-1)+k2;//ok
    //
    //fix underflow
    fix_underflow(IBD02);

    IBD00 = math_to(-1,IBD01,-1,IBD02,1);//ok

    if(k1==0){
      IBD10 = copy(IBD01);//ok
      IBD12 = copy(IBD01);//ok
    } else{
      IBD10 = scale(IBD01,k0/k1);//ok
      IBD12 = scale(IBD01,k2/k1);//ok
    }
    
    IBD11 = math_to(-1,IBD10,-1,IBD12,1);//ok
    IBD21 = copy(IBD01);
    IBD20 = allocDoubleArray(t->x);
    if(k1!=1)
      for (int i=0;i<t->x;i++)
	IBD20->array[i]= exp(-a*k1*t->array[i])*k0/(k1-1.0)+exp(-a*t->array[i])*k1+exp(-a*t->array[i])/(k1-1)*k2*k1-(1-k1)*k0/(k1-1);
    fix_underflow(IBD20);
    IBD22 = math_to(-1,IBD21,-1,IBD20,1);
  }
  //same till here
  //start Forward
  dArray *log_IBD00 = row_log(IBD00);
  dArray *log_IBD10 = row_log(IBD10);
  dArray *log_IBD20 = row_log(IBD20);
  dArray *log_IBD01 = row_log(IBD01);
  dArray *log_IBD11 = row_log(IBD11);
  dArray *log_IBD21 = row_log(IBD21);
  dArray *log_IBD02 = row_log(IBD02);
  dArray *log_IBD12 = row_log(IBD12);
  dArray *log_IBD22 = row_log(IBD22);

  dArray *v0=allocDoubleArray(snp);
  dArray *v1=allocDoubleArray(snp);
  dArray *v2=allocDoubleArray(snp);

  double maxed = MAX(log(k0),MAX(log(k1),log(k2)));
  v0->array[0] = maxed + Sk3->array[0];
  v1->array[0] = maxed + Sk2->array[0];
  v2->array[0] = maxed + Sk1->array[0];

  iArray *ptr0 = allocIntArray(snp);
  iArray *ptr1 = allocIntArray(snp);
  iArray *ptr2 = allocIntArray(snp);
  
  ptr0->array[0] = -1;
  ptr1->array[0] = -1;
  ptr2->array[0] = -1;

  
  for (int i=1 ; i<snp ; i++){
    v0->array[i] = Sk3->array[i] + max3(log_IBD00->array[i]+v0->array[i-1],log_IBD10->array[i]+v1->array[i-1],log_IBD20->array[i]+v2->array[i-1] );
    v1->array[i] = Sk2->array[i] + max3(log_IBD01->array[i]+v0->array[i-1],log_IBD11->array[i]+v1->array[i-1],log_IBD21->array[i]+v2->array[i-1]);
    v2->array[i] = Sk1->array[i] + max3(log_IBD02->array[i]+v0->array[i-1],log_IBD12->array[i]+v1->array[i-1],log_IBD22->array[i]+v2->array[i-1]);
    ptr0->array[i] = which_max(v0->array[i-1]+log_IBD00->array[i],v1->array[i-1]+log_IBD10->array[i],v2->array[i-1]+log_IBD20->array[i]);
    ptr1->array[i] = which_max(v0->array[i-1]+log_IBD01->array[i],v1->array[i-1]+log_IBD11->array[i],v2->array[i-1]+log_IBD21->array[i]);
    ptr2->array[i] = which_max(v0->array[i-1]+log_IBD02->array[i],v1->array[i-1]+log_IBD12->array[i],v2->array[i-1]+log_IBD22->array[i]);
  }
  
  iArray *pi = allocIntArray(snp);
  int las = snp-1;

  pi->array[las] = which_max(v0->array[las],v1->array[las],v2->array[las]);
  for(int i=snp-1;i>0;i--){
    if(pi->array[i]==2)
      pi->array[i-1] = ptr2->array[i];
    else if(pi->array[i]==1)
      pi->array[i-1] = ptr1->array[i];
    else
      pi->array[i-1] = ptr0->array[i];

  }

  //now clean up
  killArray(IBD00);
  killArray(IBD01);
  killArray(IBD02);
  killArray(IBD10);
  killArray(IBD11);
  killArray(IBD12);
  killArray(IBD20);
  killArray(IBD21);
  killArray(IBD22);
  killArray(ptr0);
  killArray(ptr1);
  killArray(ptr2);
  killArray(log_IBD00);
  killArray(log_IBD01);
  killArray(log_IBD02);
  killArray(log_IBD10);
  killArray(log_IBD11);
  killArray(log_IBD12); 
  killArray(log_IBD20);
  killArray(log_IBD21);
  killArray(log_IBD22);
  killArray(v0);
  killArray(v1);
  killArray(v2);
  
  return pi;  
}





double like(const double *delta,const dArray *Sk1,const dArray *Sk2,const  dArray *Sk3,const dArray *t,int double_recom,const dArray *alim){

  double a= delta[0];
  double k2=delta[1];
  double k1=delta[2];
  

  double k0=1-(k1+k2);    
  dArray *IBD00,*IBD01,*IBD02,*IBD10,*IBD11,*IBD12,*IBD20,*IBD21,*IBD22;

  if(k2<0 || k1<0 || k0<0 || k2>1 || k2>1 || k2>1 || a>alim->array[1] )
    return SYS_MAX;
  
  //optim here 
  if(double_recom){
    IBD00 = dbl_recom(t,a,k0,1);
    IBD01 = dbl_recom(t,a,k1,0);
    IBD02 = dbl_recom(t,a,k2,0);
    IBD10 = dbl_recom(t,a,k0,0);
    IBD11 = dbl_recom(t,a,k1,1);
    IBD12 = dbl_recom(t,a,k2,0);
    IBD20 = dbl_recom(t,a,k0,0);
    IBD21 = dbl_recom(t,a,k1,0);
    IBD22 = dbl_recom(t,a,k2,1);
  }
  else{
    IBD01 = dbl_recom(t,a,k1,0);//ok
    IBD02 = allocDoubleArray(t->x);
	

    
    if(k1!=1){
      for (int i=0;i<IBD02->x;i++)
	IBD02->array[i]= exp(-a*k1*t->array[i])*k2/(k1-1)+exp(-a*t->array[i])*k1+exp(-a*t->array[i])*k0*k1/(k1-1)+k2;//ok
      //fix underflow
      fix_underflow(IBD02);
    }    
    else
      for (int i=0;i<IBD02->x;i++)
	IBD02->array[i]=0;
        
    IBD00 = math_to(-1,IBD01,-1,IBD02,1);//ok

    if(k1==0){
      IBD10 = copy(IBD01);//ok
      IBD12 = copy(IBD01);//ok
    } else{
      IBD10 = scale(IBD01,k0/k1);//ok
      IBD12 = scale(IBD01,k2/k1);//ok
    }
    
    IBD11 = math_to(-1,IBD10,-1,IBD12,1);//ok
    IBD21 = copy(IBD01);
    IBD20 = allocDoubleArray(t->x);
    if(k1!=1){
      for (int i=0;i<t->x;i++)
	IBD20->array[i]= exp(-a*k1*t->array[i])*k0/(k1-1.0)+exp(-a*t->array[i])*k1+exp(-a*t->array[i])/(k1-1)*k2*k1-(1-k1)*k0/(k1-1);
      fix_underflow(IBD20);
    }else
      for (int i=0;i<t->x;i++)
	IBD20->array[i]=0;
      
    
    IBD22 = math_to(-1,IBD21,-1,IBD20,1);
  }



  //  k0=0;
  //a=0.099275;
  //k1=0.979193;
  //k2=1-k1;

  double log_IBD0 = log(k0);
  double log_IBD1 = log(k1);
  double log_IBD2 = log(k2);
  double k =0;
  
  for (int i=1;i<t->x-1;i++){
    double logp0 = exp(log_IBD0+Sk3->array[i-1]-k);
    double logp1 = exp(log_IBD1+Sk2->array[i-1]-k);
    double logp2 = exp(log_IBD2+Sk1->array[i-1]-k);
    log_IBD0 = k+log(IBD00->array[i]*logp0+IBD10->array[i]*logp1+IBD20->array[i]*logp2);
    log_IBD1 = k+log(IBD01->array[i]*logp0+IBD11->array[i]*logp1+IBD21->array[i]*logp2);
    log_IBD2 = k+log(IBD02->array[i]*logp0+IBD12->array[i]*logp1+IBD22->array[i]*logp2);
    k = MAX(log(logp1)+k,MAX(log(logp2)+k,log(logp0)+k));

    if(std::isnan(k)||std::isnan(log_IBD0)||std::isnan(log_IBD1)||std::isnan(log_IBD2)){
      printf("error in likelihood, k %f,log_IBD-012,%f,%f,%f, a,k2,k1,k0)=:\t(%f,%f,%f,%f)\t logp: %f,%f,%f\n",k,log_IBD0,log_IBD1,log_IBD2,a,k2,k1,k0,logp0,logp1,logp2);
      printf("IBD00->array[i] %f,logp0 %f,IBD10->array[i] %f,logp1 %f,IBD20->array[i] %f,logp2 %f\n",IBD00->array[i],logp0,IBD10->array[i],logp1,IBD20->array[i],logp2);
      printf("IBD01->array[i] %f,logp0 %f,IBD11->array[i] %f,logp1 %f,IBD21->array[i] %f,logp2 %f\n",IBD01->array[i],logp0,IBD11->array[i],logp1,IBD21->array[i],logp2);
      printf("IBD02->array[i] %f,logp0 %f,IBD12->array[i] %f,logp1 %f,IBD22->array[i] %f,logp2 %f\n",IBD02->array[i],logp0,IBD12->array[i],logp1,IBD22->array[i],logp2);
      exit(0);
      break;
    }
  }

  double logp0 = Sk3->array[Sk3->x-1] + log_IBD0;
  double logp1 = Sk2->array[Sk2->x-1] + log_IBD1;
  double logp2 = Sk1->array[Sk1->x-1] + log_IBD2;

  k = MAX(logp0, MAX(logp1, logp2));
  
  double l =   k+log(exp(logp0-k) + exp(logp1-k)+ exp(logp2-k));

  
  //now clean up
  killArray(IBD00);
  killArray(IBD01);
  killArray(IBD02);
  killArray(IBD10);
  killArray(IBD11);
  killArray(IBD12);
  killArray(IBD20);
  killArray(IBD21);
  killArray(IBD22);
  
  return -l;
}

double* get_start(){
  double* ptr= new double[2];
  while(1){
    double k2= rand()/(double)RAND_MAX;
    double k1 =  ((1-k2)*rand())/(double)(RAND_MAX);
    if(!(pow((k1+2*k2),2.0)<4*k2)){
      ptr[0]=k2;
      ptr[1]=k1;
      break;
    }
  }
  
  return ptr;

}


//input is of lenth=2 (k2,k1)
double relateHMM::calcA_fun(double *var){

  double *parsToSend = new double[3];
  double k0,k1,k2,a,result;
  k2=var[0];
  k1= var[1];
  k0=1-k2-k1;
  
  if(pow((k1+2*k2),2.0)<4*k2){
    delete [] parsToSend;
    return SYS_MAX;
  }
  a =  calculateA(k0,k1,k2,phi); 
  if(std::isnan(a)||std::isinf(a)){
    delete [] parsToSend;
    return SYS_MAX;
  }
  parsToSend[0] = a;
  parsToSend[1] = var[0];
  parsToSend[2] = var[1];
  
  result = like(parsToSend,Sk1,Sk2,Sk3, t,double_recom,alim);
  //parsToSend[0] is a
  var[0]=parsToSend[1]; //k2
  var[1]=parsToSend[2]; //k1

  delete [] parsToSend;

  return result;
}



double *relateHMM::calcA_optim(double* var){
  double* fixed=get_start();
  double *lbd = new double[2];
  double *ubd = new double[2];
  lbd[0]=0;//these bounds are really not ok, but better than nothing
  ubd[0]=1;
  lbd[1]=0;
  ubd[1]=1;
  
  double *ptr = new double[2];
  double result;
  ptr[0]=fixed[0];
  ptr[1]=fixed[1];
  result = -findmax_bfgs(2,ptr,&calcA_fun,NULL,lbd,ubd,NULL,-1);

  var[0]=calculateA(1-ptr[1]-ptr[0],ptr[1],ptr[0],phi);
  var[1]=ptr[0];//k2
  var[2]=ptr[1]; //k1
  var[3]=result;
  delete [] ptr;
  delete [] lbd;
  delete [] ubd;
  delete[] fixed;
  return var;
}




double relateHMM::fixA_fun(double *var){
  double *parsToSend = new double[3];
  double k0,k1,k2,a,result;
  k2=var[0];
  k1= var[1];
  k0=1-k2-k1;
  if(pow((k1+2*k2),2.0)<4*k2){
    delete [] parsToSend;
    return SYS_MAX;
  }
  a =  var[0];
  if(std::isnan(a)||std::isinf(a)){
    delete [] parsToSend;
    return SYS_MAX;
  }
  parsToSend[0] = fixA_val;
  parsToSend[1] = var[0];
  parsToSend[2] = var[1];
  result = like(parsToSend,Sk1,Sk2,Sk3, t,double_recom,alim);
  var[0]=parsToSend[2];
  var[1]=parsToSend[1];
  delete [] parsToSend;
  return result;
}

double *relateHMM::fixA_optim(double* var){

  double *lbd = new double[2];
  double *ubd = new double[2];

  lbd[0] = 0;
  ubd[0] = 1;

  lbd[1]=0;
  ubd[1]=1;

  double *ptr = new double[2];
  double result;
  ptr[0] = var[1];//k1
  ptr[1] = var[2];//k2
  result = -findmax_bfgs(2,ptr,&fixA_fun,NULL,lbd,ubd,NULL,-1);
  
  var[0]= fixA_val; //a
  var[1]=ptr[0];
  var[2]=ptr[1];
  var[3]=result;
  delete [] ptr;
  delete [] lbd;
  delete [] ubd;
  
  return var;
}


double relateHMM::fixK2_fun(double *var){
  double *parsToSend = new double[3];
  double k0,k1,k2,a,result;
  k2=fixK2_val;
  a = var[0];
  k1= var[1];
  k0=1-k2-k1;
  if(pow((k1+2*k2),2.0)<4*k2){
    delete [] parsToSend;
    return SYS_MAX;
  }
  
  if(std::isnan(a)||std::isinf(a)){
    delete [] parsToSend;
    return SYS_MAX;
  }
  parsToSend[0] = var[0];
  parsToSend[1] = fixK2_val;
  parsToSend[2] = var[1];
  result = like(parsToSend,Sk1,Sk2,Sk3, t,double_recom,alim);
  var[0]=parsToSend[0];
  var[1]=parsToSend[2];
  delete [] parsToSend;
  return result;
}

double *relateHMM::fixK2_optim(double* var){

  double *lbd = new double[2];
  double *ubd = new double[2];

  lbd[0] = alim->array[0];
  ubd[0] = alim->array[1];

  lbd[1]=0;
  ubd[1]=1;

  double *ptr = new double[2];
  double result;
  ptr[0] = var[0];//a
  ptr[1] = var[2];//k2
  result = -findmax_bfgs(2,ptr,&fixK2_fun,NULL,lbd,ubd,NULL,-1);
  
  var[0]= ptr[0]; //a
  var[1]=fixK2_val;
  var[2]=ptr[1];
  var[3]=result;
  delete [] ptr;
  delete [] lbd;
  delete [] ubd;
  
  return var;
}



double relateHMM::fixK2_calcA_fun(double *var){
  double *parsToSend = new double[3];
  double k0,k1,k2,a,result;
  k2=fixK2_val;
  k1= var[0];
  k0=1-k2-k1;
  if(pow((k1+2*k2),2.0)<4*k2){
    delete [] parsToSend;
    return SYS_MAX;
  }
  a =  calculateA(k0,k1,k2,phi); 
  if(std::isnan(a)||std::isinf(a)){
    delete [] parsToSend;
    return SYS_MAX;
  }
  parsToSend[0] = a;
  parsToSend[1] = fixK2_val;
  parsToSend[2] = var[0];
  result = like(parsToSend,Sk1,Sk2,Sk3, t,double_recom,alim);
  var[0]=parsToSend[2];
  delete [] parsToSend;
  return result;
}

double *relateHMM::fixK2_calcA_optim(double* var){

  double *lbd = new double[1];
  double *ubd = new double[1];
  lbd[0]=-2*fixK2_val+0.5*sqrt(32*fixK2_val-16*fixK2_val);
  ubd[0]=1;

  double *ptr = new double[1];
  double result;
  ptr[0]=var[2];
  result = -findmax_bfgs(1,ptr,&fixK2_calcA_fun,NULL,lbd,ubd,NULL,-1);
  
  var[1]=fixK2_val; //k2
  var[2]=ptr[0]; //k1
  var[0]=calculateA(1-var[2]-var[1],var[2],var[1],phi);
  var[3]=result;
  delete [] ptr;
  delete [] lbd;
  delete [] ubd;
  
  return var;
}

double relateHMM::fixK2_fixA_fun(double *var){
  double *parsToSend = new double[3];
  double k0,k1,k2,result;
  k2=fixK2_val;
  k1= var[0];
  k0=1-k2-k1;
  if(pow((k1+2*k2),2.0)<4*k2){
    delete [] parsToSend;
    return SYS_MAX;
  }
  parsToSend[0] = fixA_val;
  
  parsToSend[1] = fixK2_val;
  parsToSend[2] = var[0];
  result = like(parsToSend,Sk1,Sk2,Sk3, t,double_recom,alim);
  var[0]=parsToSend[2];
  delete [] parsToSend;
  return result;
}

double *relateHMM::fixK2_fixA_optim(double* var){

  double *lbd = new double[1];
  double *ubd = new double[1];
  lbd[0]=0;
  ubd[0]=1;

  double *ptr = new double[1];
  double result;
  ptr[0]=var[2];
  
  result = -findmax_bfgs(1,ptr,&fixK2_fixA_fun,NULL,lbd,ubd,NULL,-1);
  
  var[0]=fixA_val;
  var[1]=fixK2_val; //k2
  var[2]=ptr[0]; //k1
  var[3]=result;
  delete [] ptr;
  delete [] lbd;
  delete [] ubd;
  
  return var;
}

double relateHMM::full_optim_fun(double *var){
  double result = like(var,Sk1,Sk2,Sk3, t,double_recom,alim);
  return result;
}


double *relateHMM::full_optim(double* var){
  double *lbd = new double[3];
  double *ubd = new double[3];
  double result = 
  lbd[0]=alim->array[0];
  lbd[1]=0;
  lbd[2]=0;
  ubd[0]=alim->array[1];
  ubd[1]=1;
  ubd[2]=1;
  result = -findmax_bfgs(3,var,&full_optim_fun,NULL,lbd,ubd,NULL,-1);
  var[3] = result;
  delete [] lbd;
  delete [] ubd;
  return var; 

}


dArray *relateHMM::run_optimization(int i,int j,double convTol,const dArray *alim){
  int seed = time(0); //used for randomness
  if(print_info)
    printf("\n\t->fixK2:%d\tcalcA:%d,fixA:%d\tfixA_val:%f\tseed:%d\n",fixK2,calcA,fixA,fixA_val,seed);
  dArray *retVal = allocDoubleArray(4);
  convInfo = allocDoubleMatrix(j,4);
  srand((unsigned) seed);

  double difference=0;
  double *tmpBest = new double [4];
  double *ptr = new double[4];
  localTimesRun=1;
  localTimesConv=1;
  get_rand(tmpBest,alim,fixK2,fixK2_val);
  //will now do optim branches

  if(fixK2 && calcA && !fixA){
    if(print_info)
      printf("\t->will do fixK2 && calcA optim firstRun\n");
    tmpBest = fixK2_calcA_optim(tmpBest);
  }
  else if(fixK2 && fixA && !calcA){
    if(print_info)
      printf("\t->will do fixk2 && fixA optim FirstrRun\n");
    tmpBest =  fixK2_fixA_optim(tmpBest);
  } 
  else if( fixK2 && !fixA && !calcA){
    if(print_info)
      printf("\t->will do fixk2 optim firstRun\n");
    tmpBest =  fixK2_optim(tmpBest);
  }
  else if( !fixK2 && fixA && !calcA){
    if(print_info)
      printf("\t->will do fixA optim firstRun\n");
    tmpBest =  fixA_optim(tmpBest);
  }
  else if(!fixK2 && !fixA && !calcA)  {
    if(print_info)
      printf("\t->will do full optim FirstRun \n");
    tmpBest = full_optim(tmpBest);
  }
  else if(!fixK2 && !fixA && calcA)  {
    if(print_info)
      printf("\t->will do calcA optim FirstRun \n");
    tmpBest = calcA_optim(tmpBest);
  }
  else{
    printf("\t->optimization not implemented\n");
    exit(0);
  }
  //added the alim check in 0.987 23/2 2009
  if(fabs((alim->array[1]-tmpBest[0])<PRECISION_COMPARISON))
    flush_print("The Estimated 'a' is defined on border\n");
  //copy point to list of convergened points
  for(int place=0;place<4;place++)
    convInfo->matrix[localTimesRun-1][place] =tmpBest[place];
  //  printf("timesrun=%d\ttimesConv=%d\n",timesRun,timesConv);
  while(localTimesConv<i &&localTimesRun<j){
    localTimesRun++;
    get_rand(ptr,alim,fixK2,fixK2_val);
    if(fixK2 && calcA && !fixA){
      if(print_info)
	printf("\t->will do fixK2 && calcA optim\n");
      ptr =  fixK2_calcA_optim(ptr);
    }
    else if(fixK2 && fixA && !calcA){
      if(print_info)
	printf("\t->will do fixK2 && fixA optim\n");
      ptr =  fixK2_fixA_optim(ptr);
      //print_array(ptr,3);
    }
    else if(fixK2 && !fixA && !calcA) {
      if(print_info)
	printf("\t->will do fixk2 optim\n");
      ptr  =  fixK2_optim(ptr);
    }
    else if(!fixK2 && fixA && !calcA) {
      if(print_info)
	printf("\t->will do fixA optim\n");
      ptr  =  fixA_optim(ptr);
    }
    else if(!fixK2 && !fixA && calcA) {
      if(print_info)
	printf("\t->will do calcA optim\n");
      ptr  =  calcA_optim(ptr);
    }
    else{
      if(print_info)
	printf("\t->will do full optim\n");
      ptr =  full_optim(ptr);
    }
    if(fabs((alim->array[1]-ptr[0])<PRECISION_COMPARISON))
      flush_print("The Estimated 'a' is defined on border\n");

    //copy point to list of convergened points
    for(int place=0;place<4;place++)
      convInfo->matrix[localTimesRun-1][place] =ptr[place];
    difference= fabs(ptr[3]-tmpBest[3]);
    if (difference < convTol){
      localTimesConv++;
      continue;
    }
    if(abs(ptr[3])<abs(tmpBest[3])){
      for (int i=0;i<4;i++)
	tmpBest[i]=ptr[i];
      localTimesConv=1;
    }
    
  }
  if(print_info)
    printf("\t->bestLike:%f \t a=%f,\t k =(%f , %f , %f)\t#%d/%d\n",tmpBest[3],tmpBest[0],tmpBest[1],tmpBest[2],1-tmpBest[1]-tmpBest[2],localTimesConv,localTimesRun);
  retVal->array[0] = tmpBest[0];
  retVal->array[1] = tmpBest[1];
  retVal->array[2] = tmpBest[2];
  retVal->array[3] = tmpBest[3];
  
  delete [] tmpBest;
  delete [] ptr;
  return retVal;
}


//used
dMatrix *LDemission( dArray *pBA,dArray *pbA,dArray *pBa,dArray *pba,iMatrix *dist,double epsilon){
  
  dMatrix *ep=allocDoubleMatrix(dist->x,dist->y);
  
  for(int i=0;i<dist->x;i++)
    for(int n=0;n<dist->y;n++)
      ep->matrix[i][n]=pow(epsilon,dist->matrix[i][n])*pow(1-epsilon,8-dist->matrix[i][n]);
  
  dMatrix *S=allocDoubleMatrix(3,pBa->x);
  
  for(int tal=0;tal<S->y;tal++) {

    //0
    S->matrix[0][tal] = (pow(pBA->array[tal],4.0))*ep->matrix[0][tal];
    S->matrix[1][tal] = (pow(pBA->array[tal],3.0))*ep->matrix[0][tal];
    S->matrix[2][tal] = (pow(pBA->array[tal],2.0))*ep->matrix[0][tal];
    //   continue;
    S->matrix[0][tal] += (4*pow(pBA->array[tal],3)*pbA->array[tal])*ep->matrix[1][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pbA->array[tal])*ep->matrix[1][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pbA->array[tal],2))*ep->matrix[2][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;
    //1
    S->matrix[0][tal] += (4*pow(pBA->array[tal],3)*pbA->array[tal])*ep->matrix[3][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pbA->array[tal])*ep->matrix[3][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pow(pbA->array[tal],2))*ep->matrix[4][tal];
    S->matrix[1][tal] += (pBA->array[tal]*pow(pbA->array[tal],2)+pow(pBA->array[tal],2)*pbA->array[tal])*ep->matrix[4][tal];
    S->matrix[2][tal] += (2*pBA->array[tal]*pbA->array[tal])*ep->matrix[4][tal];
    
    S->matrix[0][tal] += (4*pow(pbA->array[tal],3)*pBA->array[tal])*ep->matrix[5][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBA->array[tal])*ep->matrix[5][tal];
    S->matrix[2][tal] += 0;
    //2
    
    S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pbA->array[tal],2))*ep->matrix[6][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pbA->array[tal],3)*pBA->array[tal])*ep->matrix[7][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBA->array[tal])*ep->matrix[7][tal];
    S->matrix[2][tal] += 0;
       
    S->matrix[0][tal] += (pow(pbA->array[tal],4))*ep->matrix[8][tal];
    S->matrix[1][tal] += (pow(pbA->array[tal],3))*ep->matrix[8][tal];
    S->matrix[2][tal] += (pow(pbA->array[tal],2))*ep->matrix[8][tal];
    //3
    S->matrix[0][tal] += (4*pow(pBA->array[tal],3)*pBa->array[tal])*ep->matrix[9][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pBa->array[tal])*ep->matrix[9][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal]+\
      4*pow(pBA->array[tal],3)*pba->array[tal]+8*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal])*ep->matrix[10][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pba->array[tal]+2*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[10][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pBA->array[tal]*pBa->array[tal]*pow(pbA->array[tal],2)+4*pow(pBA->array[tal],2)*pba->array[tal]*pbA->array[tal])*ep->matrix[11][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    
    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal]+\
      4*pow(pBA->array[tal],3)*pba->array[tal]+8*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal])*ep->matrix[12][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pba->array[tal]+2*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[12][tal];
    S->matrix[2][tal] += 0;
    


    S->matrix[0][tal] += (8*pow(pBA->array[tal],2)*pbA->array[tal]*pba->array[tal]+8*pBA->array[tal]*pow(pbA->array[tal],2)*pBa->array[tal])*ep->matrix[13][tal];
    S->matrix[1][tal] += ((2*pBA->array[tal]*pbA->array[tal]*pba->array[tal]+2*pBA->array[tal]*pbA->array[tal]*pBa->array[tal]))*ep->matrix[13][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]+4*pow(pbA->array[tal],3)*\
      pBa->array[tal]+8*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal])*ep->matrix[14][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBa->array[tal]+2*pbA->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[14][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pBA->array[tal]*pBa->array[tal]*pow(pbA->array[tal],2)+4*pow(pBA->array[tal],2)*pba->array[tal]*pbA->array[tal])*ep->matrix[15][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]+4*pow(pbA->array[tal],3)*pBa->array[tal]\
			  +8*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal])*ep->matrix[16][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBa->array[tal]+2*pbA->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[16][tal];
    S->matrix[2][tal] += 0;
    

    S->matrix[0][tal] += (4*pow(pbA->array[tal],3)*pba->array[tal])*ep->matrix[17][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pba->array[tal])*ep->matrix[17][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[18][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pbA->array[tal]*pBA->array[tal])*ep->matrix[19][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+2*pow(pbA->array[tal],2)*pow(pBa->array[tal],2)) *ep->matrix[20][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pbA->array[tal]*pBA->array[tal] )*ep->matrix[21][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (8*pba->array[tal]*pBA->array[tal]*pBa->array[tal]*pbA->array[tal] )*ep->matrix[22][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;
      

    S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+4*pow(pba->array[tal],2)*pBA->array[tal]*pbA->array[tal] )*ep->matrix[23][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] +=(2*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+2*pow(pbA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[24][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] +=(4*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+4*pow(pba->array[tal],2)*pBA->array[tal]*pbA->array[tal])*ep->matrix[25][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (2*pow(pbA->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[26][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],3)*pBa->array[tal])*ep->matrix[27][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pBa->array[tal])*ep->matrix[27][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal]+4*pow(pBA->array[tal],3)*pba->array[tal]+8*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal])*ep->matrix[28][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pba->array[tal]+2*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[28][tal];
    S->matrix[2][tal] += 0;


    S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pbA->array[tal],2)+4*pow(pBA->array[tal],2)*pba->array[tal]*pbA->array[tal])*ep->matrix[29][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal]+4*pow(pBA->array[tal],3)*pba->array[tal]+8*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal])*ep->matrix[30][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pba->array[tal]+2*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[30][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (8*pow(pBA->array[tal],2)*pbA->array[tal]*pba->array[tal]+8*pBA->array[tal]*pow(pbA->array[tal],2)*pBa->array[tal])*ep->matrix[31][tal];
    S->matrix[1][tal] += (2*pBA->array[tal]*pbA->array[tal]*pba->array[tal]+2*pBA->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[31][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += ((4*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]+4*pow(pbA->array[tal],3)*pBa->array[tal]+8*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]))*ep->matrix[32][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBa->array[tal]+2*pbA->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[32][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pBA->array[tal]*pBa->array[tal]*pow(pbA->array[tal],2)+4*pow(pBA->array[tal],2)*pba->array[tal]*pbA->array[tal])*ep->matrix[33][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;


   S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]+4*pow(pbA->array[tal],3)*pBa->array[tal]+8*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal])*ep->matrix[34][tal];
   S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBa->array[tal]+2*pbA->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[34][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pbA->array[tal],3)*pba->array[tal])*ep->matrix[35][tal];
   S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pba->array[tal])*ep->matrix[35][tal];
   S->matrix[2][tal] += 0;


   S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[36][tal];
   S->matrix[1][tal] += (pow(pBA->array[tal],2)*pBa->array[tal]+pBA->array[tal]*pow(pBa->array[tal],2))*ep->matrix[36][tal];
   S->matrix[2][tal] += (2*pBA->array[tal]*pBa->array[tal])  *ep->matrix[36][tal];
   //
   S->matrix[0][tal] += (8*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+8*pBA->array[tal]*pbA->array[tal]*pow(pBa->array[tal],2))*ep->matrix[37][tal];
   S->matrix[1][tal] += (2*pBA->array[tal]*pba->array[tal]*pBa->array[tal]+2*pBA->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[37][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (8*pBA->array[tal]*pba->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[38][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

  

   S->matrix[0][tal] += (8*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+8*pBA->array[tal]*pbA->array[tal]*pow(pBa->array[tal],2))*ep->matrix[39][tal];
     S->matrix[1][tal] += (2*pBA->array[tal]*pba->array[tal]*pBa->array[tal]+2*pBA->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[39][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+8*pba->array[tal]*pbA->array[tal]*pBA->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pow(pbA->array[tal],2))*ep->matrix[40][tal];
   S->matrix[1][tal] += (pow(pBA->array[tal],2)*pba->array[tal]+pow(pba->array[tal],2)*pBA->array[tal]+pow(pBa->array[tal],2)*pbA->array[tal]+pow(pbA->array[tal],2)*pBa->array[tal])*ep->matrix[40][tal];
   S->matrix[2][tal] +=(2*pBA->array[tal]*pba->array[tal]+2*pBa->array[tal]*pbA->array[tal])*ep->matrix[40][tal];
   //

 
   S->matrix[0][tal] += (8*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+8*pbA->array[tal]*pBA->array[tal]*pow(pba->array[tal],2))*ep->matrix[41][tal];
   S->matrix[1][tal] += (2*pbA->array[tal]*pBa->array[tal]*pba->array[tal]+2*pbA->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[41][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (8*pBA->array[tal]*pba->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[42][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] +=    (8*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+8*pbA->array[tal]*pBA->array[tal]*pow(pba->array[tal],2))*ep->matrix[43][tal];
   S->matrix[1][tal] += (2*pbA->array[tal]*pBa->array[tal]*pba->array[tal]+2*pbA->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[43][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[44][tal];
   S->matrix[1][tal] += (pow(pbA->array[tal],2)*pba->array[tal]+pbA->array[tal]*pow(pba->array[tal],2))*ep->matrix[44][tal];
   S->matrix[2][tal] += (2*pbA->array[tal]*pba->array[tal])*ep->matrix[44][tal];
   
   
   S->matrix[0][tal] += (4*pow(pBa->array[tal],3)*pBA->array[tal])*ep->matrix[45][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pBA->array[tal])*ep->matrix[45][tal];
   S->matrix[2][tal] += 0;
   // is ok

   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal]+4*pow(pBa->array[tal],3)*pbA->array[tal]+8*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal])*ep->matrix[46][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pbA->array[tal]+2*pBa->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[46][tal];
   S->matrix[2][tal] += 0;
  
   S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pba->array[tal],2)+4*pow(pBa->array[tal],2)*pbA->array[tal]*pba->array[tal])*ep->matrix[47][tal];
     S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;
     
 
   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal]+4*pow(pBa->array[tal],3)*pbA->array[tal]+8*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal])*ep->matrix[48][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pbA->array[tal]+2*pBa->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[48][tal];
   S->matrix[2][tal] += 0;
   //is ok
   //   continue;       
   S->matrix[0][tal] += (8*pow(pBa->array[tal],2)*pba->array[tal]*pbA->array[tal]+8*pBa->array[tal]*pow(pba->array[tal],2)*pBA->array[tal])*ep->matrix[49][tal];
   S->matrix[1][tal] += (2*pBa->array[tal]*pba->array[tal]*pbA->array[tal]+2*pBa->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[49][tal];
   S->matrix[2][tal] += 0;
    
   S->matrix[0][tal] += (4*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal]+4*pow(pba->array[tal],3)*pBA->array[tal]+8*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal])*ep->matrix[50][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBA->array[tal]+2*pba->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[50][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pba->array[tal],2)+4*pow(pBa->array[tal],2)*pbA->array[tal]*pba->array[tal])*ep->matrix[51][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;


   S->matrix[0][tal] += (4*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal]+4*pow(pba->array[tal],3)*pBA->array[tal]+8*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal])*ep->matrix[52][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBA->array[tal]+2*pba->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[52][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pba->array[tal],3)*pbA->array[tal])*ep->matrix[53][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pbA->array[tal])*ep->matrix[53][tal];
   S->matrix[2][tal] += 0;
   
   S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[54][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pbA->array[tal]*pBA->array[tal])*ep->matrix[55][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+2*pow(pbA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[56][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pbA->array[tal]*pBA->array[tal])*ep->matrix[57][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (8*pba->array[tal]*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[58][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+4*pow(pba->array[tal],2)*pBA->array[tal]*pbA->array[tal])*ep->matrix[59][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] +=(2*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+2*pow(pbA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[60][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+4*pow(pba->array[tal],2)*pBA->array[tal]*pbA->array[tal])*ep->matrix[61][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (2*pow(pbA->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[62][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

 
   S->matrix[0][tal] += (4*pow(pBa->array[tal],3)*pBA->array[tal])*ep->matrix[63][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pBA->array[tal])*ep->matrix[63][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal]+4*pow(pBa->array[tal],3)*pbA->array[tal]+8*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal])*ep->matrix[64][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pbA->array[tal]+2*pBa->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[64][tal];

   S->matrix[2][tal] += 0;
   ////

   S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pba->array[tal],2)+4*pow(pBa->array[tal],2)*pbA->array[tal]*pba->array[tal])*ep->matrix[65][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal]+4*pow(pBa->array[tal],3)*pbA->array[tal]+8*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal])*ep->matrix[66][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pbA->array[tal]+2*pBa->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[66][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (8*pow(pBa->array[tal],2)*pba->array[tal]*pbA->array[tal]+8*pBa->array[tal]*pow(pba->array[tal],2)*pBA->array[tal])*ep->matrix[67][tal];
   S->matrix[1][tal] += (2*pBa->array[tal]*pba->array[tal]*pbA->array[tal]+2*pBa->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[67][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal]+4*pow(pba->array[tal],3)*pBA->array[tal]+8*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal])*ep->matrix[68][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBA->array[tal]+2*pba->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[68][tal];
   S->matrix[2][tal] += 0;
  ///////////////
   S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pba->array[tal],2)+4*pow(pBa->array[tal],2)*pbA->array[tal]*pba->array[tal])*ep->matrix[69][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal]+4*pow(pba->array[tal],3)*pBA->array[tal]+8*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal])*ep->matrix[70][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBA->array[tal]+2*pba->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[70][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pba->array[tal],3)*pbA->array[tal])*ep->matrix[71][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pbA->array[tal])*ep->matrix[71][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += pow(pBa->array[tal],4)*ep->matrix[72][tal];
   S->matrix[1][tal] += pow(pBa->array[tal],3)*ep->matrix[72][tal];
   S->matrix[2][tal] += pow(pBa->array[tal],2)*ep->matrix[72][tal];

   S->matrix[0][tal] += (4*pow(pBa->array[tal],3)*pba->array[tal])*ep->matrix[73][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pba->array[tal])*ep->matrix[73][tal];
   S->matrix[2][tal] += 0;
   
   S->matrix[0][tal] += (2*pow(pBa->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[74][tal];
   S->matrix[1][tal] +=0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBa->array[tal],3)*pba->array[tal])*ep->matrix[75][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pba->array[tal])*ep->matrix[75][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[76][tal];
   S->matrix[1][tal] += (pBa->array[tal]*pow(pba->array[tal],2)+pow(pBa->array[tal],2)*pba->array[tal])*ep->matrix[76][tal];
   S->matrix[2][tal] += (2*pBa->array[tal]*pba->array[tal])*ep->matrix[76][tal];

   S->matrix[0][tal] += (4*pow(pba->array[tal],3)*pBa->array[tal])*ep->matrix[77][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBa->array[tal])*ep->matrix[77][tal];
   S->matrix[2][tal] += 0; 

   S->matrix[0][tal] += (2*pow(pBa->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[78][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;


   S->matrix[0][tal] += (4*pow(pba->array[tal],3)*pBa->array[tal])*ep->matrix[79][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBa->array[tal])*ep->matrix[79][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += pow(pba->array[tal],4) * ep->matrix[80][tal];
   S->matrix[1][tal] += pow(pba->array[tal],3) *ep->matrix[80][tal];
   S->matrix[2][tal] += pow(pba->array[tal],2) *ep->matrix[80][tal];

  }
  killMatrix(ep);
  return S;
}
//used
dMatrix *LDemission_new( dArray *pBA,dArray *pbA,dArray *pBa,dArray *pba,dMatrix *ep){
  dMatrix *S=allocDoubleMatrix(3,pBa->x);
  
   for(int tal=0;tal<S->y;tal++) {

    //0
    S->matrix[0][tal] = (pow(pBA->array[tal],4.0))*ep->matrix[0][tal];
    S->matrix[1][tal] = (pow(pBA->array[tal],3.0))*ep->matrix[0][tal];
    S->matrix[2][tal] = (pow(pBA->array[tal],2.0))*ep->matrix[0][tal];
    //   continue;
    S->matrix[0][tal] += (4*pow(pBA->array[tal],3)*pbA->array[tal])*ep->matrix[1][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pbA->array[tal])*ep->matrix[1][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pbA->array[tal],2))*ep->matrix[2][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;
    //1
    S->matrix[0][tal] += (4*pow(pBA->array[tal],3)*pbA->array[tal])*ep->matrix[3][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pbA->array[tal])*ep->matrix[3][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pow(pbA->array[tal],2))*ep->matrix[4][tal];
    S->matrix[1][tal] += (pBA->array[tal]*pow(pbA->array[tal],2)+pow(pBA->array[tal],2)*pbA->array[tal])*ep->matrix[4][tal];
    S->matrix[2][tal] += (2*pBA->array[tal]*pbA->array[tal])*ep->matrix[4][tal];
    
    S->matrix[0][tal] += (4*pow(pbA->array[tal],3)*pBA->array[tal])*ep->matrix[5][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBA->array[tal])*ep->matrix[5][tal];
    S->matrix[2][tal] += 0;
    //2
    
    S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pbA->array[tal],2))*ep->matrix[6][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pbA->array[tal],3)*pBA->array[tal])*ep->matrix[7][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBA->array[tal])*ep->matrix[7][tal];
    S->matrix[2][tal] += 0;
       
    S->matrix[0][tal] += (pow(pbA->array[tal],4))*ep->matrix[8][tal];
    S->matrix[1][tal] += (pow(pbA->array[tal],3))*ep->matrix[8][tal];
    S->matrix[2][tal] += (pow(pbA->array[tal],2))*ep->matrix[8][tal];
    //3
    S->matrix[0][tal] += (4*pow(pBA->array[tal],3)*pBa->array[tal])*ep->matrix[9][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pBa->array[tal])*ep->matrix[9][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal]+\
      4*pow(pBA->array[tal],3)*pba->array[tal]+8*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal])*ep->matrix[10][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pba->array[tal]+2*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[10][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pBA->array[tal]*pBa->array[tal]*pow(pbA->array[tal],2)+4*pow(pBA->array[tal],2)*pba->array[tal]*pbA->array[tal])*ep->matrix[11][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    
    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal]+\
      4*pow(pBA->array[tal],3)*pba->array[tal]+8*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal])*ep->matrix[12][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pba->array[tal]+2*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[12][tal];
    S->matrix[2][tal] += 0;
    


    S->matrix[0][tal] += (8*pow(pBA->array[tal],2)*pbA->array[tal]*pba->array[tal]+8*pBA->array[tal]*pow(pbA->array[tal],2)*pBa->array[tal])*ep->matrix[13][tal];
    S->matrix[1][tal] += ((2*pBA->array[tal]*pbA->array[tal]*pba->array[tal]+2*pBA->array[tal]*pbA->array[tal]*pBa->array[tal]))*ep->matrix[13][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]+4*pow(pbA->array[tal],3)*\
      pBa->array[tal]+8*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal])*ep->matrix[14][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBa->array[tal]+2*pbA->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[14][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pBA->array[tal]*pBa->array[tal]*pow(pbA->array[tal],2)+4*pow(pBA->array[tal],2)*pba->array[tal]*pbA->array[tal])*ep->matrix[15][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]+4*pow(pbA->array[tal],3)*pBa->array[tal]\
			  +8*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal])*ep->matrix[16][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBa->array[tal]+2*pbA->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[16][tal];
    S->matrix[2][tal] += 0;
    

    S->matrix[0][tal] += (4*pow(pbA->array[tal],3)*pba->array[tal])*ep->matrix[17][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pba->array[tal])*ep->matrix[17][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[18][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pbA->array[tal]*pBA->array[tal])*ep->matrix[19][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+2*pow(pbA->array[tal],2)*pow(pBa->array[tal],2)) *ep->matrix[20][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pbA->array[tal]*pBA->array[tal] )*ep->matrix[21][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (8*pba->array[tal]*pBA->array[tal]*pBa->array[tal]*pbA->array[tal] )*ep->matrix[22][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;
      

    S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+4*pow(pba->array[tal],2)*pBA->array[tal]*pbA->array[tal] )*ep->matrix[23][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] +=(2*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+2*pow(pbA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[24][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] +=(4*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+4*pow(pba->array[tal],2)*pBA->array[tal]*pbA->array[tal])*ep->matrix[25][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (2*pow(pbA->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[26][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],3)*pBa->array[tal])*ep->matrix[27][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pBa->array[tal])*ep->matrix[27][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal]+4*pow(pBA->array[tal],3)*pba->array[tal]+8*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal])*ep->matrix[28][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pba->array[tal]+2*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[28][tal];
    S->matrix[2][tal] += 0;


    S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pbA->array[tal],2)+4*pow(pBA->array[tal],2)*pba->array[tal]*pbA->array[tal])*ep->matrix[29][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal]+4*pow(pBA->array[tal],3)*pba->array[tal]+8*pow(pBA->array[tal],2)*pbA->array[tal]*pBa->array[tal])*ep->matrix[30][tal];
    S->matrix[1][tal] += (2*pow(pBA->array[tal],2)*pba->array[tal]+2*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[30][tal];
    S->matrix[2][tal] += 0;
    
    S->matrix[0][tal] += (8*pow(pBA->array[tal],2)*pbA->array[tal]*pba->array[tal]+8*pBA->array[tal]*pow(pbA->array[tal],2)*pBa->array[tal])*ep->matrix[31][tal];
    S->matrix[1][tal] += (2*pBA->array[tal]*pbA->array[tal]*pba->array[tal]+2*pBA->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[31][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += ((4*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]+4*pow(pbA->array[tal],3)*pBa->array[tal]+8*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]))*ep->matrix[32][tal];
    S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBa->array[tal]+2*pbA->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[32][tal];
    S->matrix[2][tal] += 0;

    S->matrix[0][tal] += (4*pBA->array[tal]*pBa->array[tal]*pow(pbA->array[tal],2)+4*pow(pBA->array[tal],2)*pba->array[tal]*pbA->array[tal])*ep->matrix[33][tal];
    S->matrix[1][tal] += 0;
    S->matrix[2][tal] += 0;


   S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal]+4*pow(pbA->array[tal],3)*pBa->array[tal]+8*pow(pbA->array[tal],2)*pBA->array[tal]*pba->array[tal])*ep->matrix[34][tal];
   S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pBa->array[tal]+2*pbA->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[34][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pbA->array[tal],3)*pba->array[tal])*ep->matrix[35][tal];
   S->matrix[1][tal] += (2*pow(pbA->array[tal],2)*pba->array[tal])*ep->matrix[35][tal];
   S->matrix[2][tal] += 0;


   S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[36][tal];
   S->matrix[1][tal] += (pow(pBA->array[tal],2)*pBa->array[tal]+pBA->array[tal]*pow(pBa->array[tal],2))*ep->matrix[36][tal];
   S->matrix[2][tal] += (2*pBA->array[tal]*pBa->array[tal])  *ep->matrix[36][tal];
   //
   S->matrix[0][tal] += (8*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+8*pBA->array[tal]*pbA->array[tal]*pow(pBa->array[tal],2))*ep->matrix[37][tal];
   S->matrix[1][tal] += (2*pBA->array[tal]*pba->array[tal]*pBa->array[tal]+2*pBA->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[37][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (8*pBA->array[tal]*pba->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[38][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

  

   S->matrix[0][tal] += (8*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+8*pBA->array[tal]*pbA->array[tal]*pow(pBa->array[tal],2))*ep->matrix[39][tal];
     S->matrix[1][tal] += (2*pBA->array[tal]*pba->array[tal]*pBa->array[tal]+2*pBA->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[39][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+8*pba->array[tal]*pbA->array[tal]*pBA->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pow(pbA->array[tal],2))*ep->matrix[40][tal];
   S->matrix[1][tal] += (pow(pBA->array[tal],2)*pba->array[tal]+pow(pba->array[tal],2)*pBA->array[tal]+pow(pBa->array[tal],2)*pbA->array[tal]+pow(pbA->array[tal],2)*pBa->array[tal])*ep->matrix[40][tal];
   S->matrix[2][tal] +=(2*pBA->array[tal]*pba->array[tal]+2*pBa->array[tal]*pbA->array[tal])*ep->matrix[40][tal];
   //

 
   S->matrix[0][tal] += (8*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+8*pbA->array[tal]*pBA->array[tal]*pow(pba->array[tal],2))*ep->matrix[41][tal];
   S->matrix[1][tal] += (2*pbA->array[tal]*pBa->array[tal]*pba->array[tal]+2*pbA->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[41][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (8*pBA->array[tal]*pba->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[42][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] +=    (8*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+8*pbA->array[tal]*pBA->array[tal]*pow(pba->array[tal],2))*ep->matrix[43][tal];
   S->matrix[1][tal] += (2*pbA->array[tal]*pBa->array[tal]*pba->array[tal]+2*pbA->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[43][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[44][tal];
   S->matrix[1][tal] += (pow(pbA->array[tal],2)*pba->array[tal]+pbA->array[tal]*pow(pba->array[tal],2))*ep->matrix[44][tal];
   S->matrix[2][tal] += (2*pbA->array[tal]*pba->array[tal])*ep->matrix[44][tal];
   
   
   S->matrix[0][tal] += (4*pow(pBa->array[tal],3)*pBA->array[tal])*ep->matrix[45][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pBA->array[tal])*ep->matrix[45][tal];
   S->matrix[2][tal] += 0;
   // is ok

   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal]+4*pow(pBa->array[tal],3)*pbA->array[tal]+8*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal])*ep->matrix[46][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pbA->array[tal]+2*pBa->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[46][tal];
   S->matrix[2][tal] += 0;
  
   S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pba->array[tal],2)+4*pow(pBa->array[tal],2)*pbA->array[tal]*pba->array[tal])*ep->matrix[47][tal];
     S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;
     
 
   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal]+4*pow(pBa->array[tal],3)*pbA->array[tal]+8*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal])*ep->matrix[48][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pbA->array[tal]+2*pBa->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[48][tal];
   S->matrix[2][tal] += 0;
   //is ok
   //   continue;       
   S->matrix[0][tal] += (8*pow(pBa->array[tal],2)*pba->array[tal]*pbA->array[tal]+8*pBa->array[tal]*pow(pba->array[tal],2)*pBA->array[tal])*ep->matrix[49][tal];
   S->matrix[1][tal] += (2*pBa->array[tal]*pba->array[tal]*pbA->array[tal]+2*pBa->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[49][tal];
   S->matrix[2][tal] += 0;
    
   S->matrix[0][tal] += (4*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal]+4*pow(pba->array[tal],3)*pBA->array[tal]+8*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal])*ep->matrix[50][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBA->array[tal]+2*pba->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[50][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pba->array[tal],2)+4*pow(pBa->array[tal],2)*pbA->array[tal]*pba->array[tal])*ep->matrix[51][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;


   S->matrix[0][tal] += (4*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal]+4*pow(pba->array[tal],3)*pBA->array[tal]+8*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal])*ep->matrix[52][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBA->array[tal]+2*pba->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[52][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pba->array[tal],3)*pbA->array[tal])*ep->matrix[53][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pbA->array[tal])*ep->matrix[53][tal];
   S->matrix[2][tal] += 0;
   
   S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[54][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pbA->array[tal]*pBA->array[tal])*ep->matrix[55][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (2*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+2*pow(pbA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[56][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBA->array[tal],2)*pba->array[tal]*pBa->array[tal]+4*pow(pBa->array[tal],2)*pbA->array[tal]*pBA->array[tal])*ep->matrix[57][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (8*pba->array[tal]*pBA->array[tal]*pBa->array[tal]*pbA->array[tal])*ep->matrix[58][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+4*pow(pba->array[tal],2)*pBA->array[tal]*pbA->array[tal])*ep->matrix[59][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] +=(2*pow(pBA->array[tal],2)*pow(pba->array[tal],2)+2*pow(pbA->array[tal],2)*pow(pBa->array[tal],2))*ep->matrix[60][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pbA->array[tal],2)*pBa->array[tal]*pba->array[tal]+4*pow(pba->array[tal],2)*pBA->array[tal]*pbA->array[tal])*ep->matrix[61][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (2*pow(pbA->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[62][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

 
   S->matrix[0][tal] += (4*pow(pBa->array[tal],3)*pBA->array[tal])*ep->matrix[63][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pBA->array[tal])*ep->matrix[63][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal]+4*pow(pBa->array[tal],3)*pbA->array[tal]+8*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal])*ep->matrix[64][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pbA->array[tal]+2*pBa->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[64][tal];

   S->matrix[2][tal] += 0;
   ////

   S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pba->array[tal],2)+4*pow(pBa->array[tal],2)*pbA->array[tal]*pba->array[tal])*ep->matrix[65][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal]+4*pow(pBa->array[tal],3)*pbA->array[tal]+8*pow(pBa->array[tal],2)*pba->array[tal]*pBA->array[tal])*ep->matrix[66][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pbA->array[tal]+2*pBa->array[tal]*pBA->array[tal]*pba->array[tal])*ep->matrix[66][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (8*pow(pBa->array[tal],2)*pba->array[tal]*pbA->array[tal]+8*pBa->array[tal]*pow(pba->array[tal],2)*pBA->array[tal])*ep->matrix[67][tal];
   S->matrix[1][tal] += (2*pBa->array[tal]*pba->array[tal]*pbA->array[tal]+2*pBa->array[tal]*pba->array[tal]*pBA->array[tal])*ep->matrix[67][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal]+4*pow(pba->array[tal],3)*pBA->array[tal]+8*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal])*ep->matrix[68][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBA->array[tal]+2*pba->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[68][tal];
   S->matrix[2][tal] += 0;
  ///////////////
   S->matrix[0][tal] += (4*pBa->array[tal]*pBA->array[tal]*pow(pba->array[tal],2)+4*pow(pBa->array[tal],2)*pbA->array[tal]*pba->array[tal])*ep->matrix[69][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal]+4*pow(pba->array[tal],3)*pBA->array[tal]+8*pow(pba->array[tal],2)*pBa->array[tal]*pbA->array[tal])*ep->matrix[70][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBA->array[tal]+2*pba->array[tal]*pbA->array[tal]*pBa->array[tal])*ep->matrix[70][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pba->array[tal],3)*pbA->array[tal])*ep->matrix[71][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pbA->array[tal])*ep->matrix[71][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += pow(pBa->array[tal],4)*ep->matrix[72][tal];
   S->matrix[1][tal] += pow(pBa->array[tal],3)*ep->matrix[72][tal];
   S->matrix[2][tal] += pow(pBa->array[tal],2)*ep->matrix[72][tal];

   S->matrix[0][tal] += (4*pow(pBa->array[tal],3)*pba->array[tal])*ep->matrix[73][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pba->array[tal])*ep->matrix[73][tal];
   S->matrix[2][tal] += 0;
   
   S->matrix[0][tal] += (2*pow(pBa->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[74][tal];
   S->matrix[1][tal] +=0;
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBa->array[tal],3)*pba->array[tal])*ep->matrix[75][tal];
   S->matrix[1][tal] += (2*pow(pBa->array[tal],2)*pba->array[tal])*ep->matrix[75][tal];
   S->matrix[2][tal] += 0;

   S->matrix[0][tal] += (4*pow(pBa->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[76][tal];
   S->matrix[1][tal] += (pBa->array[tal]*pow(pba->array[tal],2)+pow(pBa->array[tal],2)*pba->array[tal])*ep->matrix[76][tal];
   S->matrix[2][tal] += (2*pBa->array[tal]*pba->array[tal])*ep->matrix[76][tal];

   S->matrix[0][tal] += (4*pow(pba->array[tal],3)*pBa->array[tal])*ep->matrix[77][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBa->array[tal])*ep->matrix[77][tal];
   S->matrix[2][tal] += 0; 

   S->matrix[0][tal] += (2*pow(pBa->array[tal],2)*pow(pba->array[tal],2))*ep->matrix[78][tal];
   S->matrix[1][tal] += 0;
   S->matrix[2][tal] += 0;


   S->matrix[0][tal] += (4*pow(pba->array[tal],3)*pBa->array[tal])*ep->matrix[79][tal];
   S->matrix[1][tal] += (2*pow(pba->array[tal],2)*pBa->array[tal])*ep->matrix[79][tal];
   S->matrix[2][tal] += 0;

   double tmp1 = pow(pba->array[tal],4);
   double tmp2 =  ep->matrix[80][tal];
   S->matrix[0][tal] += tmp1  * tmp2;
   S->matrix[1][tal] += pow(pba->array[tal],3) *ep->matrix[80][tal];
   S->matrix[2][tal] += pow(pba->array[tal],2) *ep->matrix[80][tal];

  }
   
  return S;
}

void fillup(dMatrix *mat){
  for (int i=0; i<mat->x;i++)
    for(int j=0; j<mat->y;j++)
      mat->matrix[i][j] = 0;
}








dMatrix *emission2(dArray *maf,dArray *maf2,double epsilon,dMatrix *dist){
  //make ep matrix
  dMatrix *ep=allocDoubleMatrix(9,maf->x);
  for(int i=0;i<maf->x;i++){
    for (int j=0;j<9;j++){
      //  cout <<"dist->matrix[j][i]:\t"<<dist->matrix[j][i]<<endl;
      //cout <<"4.0-dist->matrix[j][i]:\t"<<4.0-dist->matrix[j][i]<<endl;
      ep->matrix[j][i]=pow(epsilon,0.0+dist->matrix[j][i])*pow(1.0-epsilon,4.0-dist->matrix[j][i]);
    }
  }

  dMatrix *retVal=allocDoubleMatrix(3,maf->x);

  for(int tal=0;tal<maf->x;tal++){
    retVal->matrix[2][tal] = pow(maf->array[tal],2.0)*ep->matrix[0][tal];
    retVal->matrix[1][tal] = pow(maf->array[tal],3)*ep->matrix[0][tal];
    retVal->matrix[0][tal] = pow(maf->array[tal],4)*ep->matrix[0][tal];
    
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 2*pow(maf->array[tal],2)*maf2->array[tal]*ep->matrix[1][tal];
    retVal->matrix[0][tal] += 4*pow(maf->array[tal],3)*maf2->array[tal]*ep->matrix[1][tal];
  
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 0;
    retVal->matrix[0][tal] += 2*pow(maf->array[tal],2)*pow(maf2->array[tal],2)*ep->matrix[2][tal];
    
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 2*pow(maf->array[tal],2)*maf2->array[tal]*ep->matrix[3][tal];
    retVal->matrix[0][tal] += 4*pow(maf->array[tal],3)*maf2->array[tal]*ep->matrix[3][tal];
  
    retVal->matrix[2][tal] += 2*maf2->array[tal]*maf->array[tal]*ep->matrix[4][tal];
    retVal->matrix[1][tal] += maf2->array[tal]*maf->array[tal]*(maf2->array[tal]+maf->array[tal])*ep->matrix[4][tal];
    retVal->matrix[0][tal] += 4*pow(maf2->array[tal],2)*pow(maf->array[tal],2)*ep->matrix[4][tal];
  
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 2*pow(maf2->array[tal],2)*maf->array[tal]*ep->matrix[5][tal];
    retVal->matrix[0][tal] += 4*pow(maf2->array[tal],3)*maf->array[tal]*ep->matrix[5][tal];
  
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 0;
    retVal->matrix[0][tal] += 2*pow(maf2->array[tal],2)*pow(maf->array[tal],2)*ep->matrix[6][tal];
    
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 2*pow(maf2->array[tal],2)*maf->array[tal]*ep->matrix[7][tal];
    retVal->matrix[0][tal] += 4*pow(maf2->array[tal],3)*maf->array[tal]*ep->matrix[7][tal];
  
    retVal->matrix[2][tal] += pow(maf2->array[tal],2.0)*ep->matrix[8][tal];
    retVal->matrix[1][tal] += pow(maf2->array[tal],3)*ep->matrix[8][tal];
    retVal->matrix[0][tal] += pow(maf2->array[tal],4)*ep->matrix[8][tal];
    
  }
  killMatrix(ep);
  return retVal;
}



dMatrix *emission2_new(dArray *maf,dArray *maf2, dMatrix *dist){
  //make ep matrix
  dMatrix *ep = dist;

  dMatrix *retVal=allocDoubleMatrix(3,maf->x);

  for(int tal=0;tal<maf->x;tal++){
    retVal->matrix[2][tal] = pow(maf->array[tal],2.0)*ep->matrix[0][tal];
    retVal->matrix[1][tal] = pow(maf->array[tal],3)*ep->matrix[0][tal];
    retVal->matrix[0][tal] = pow(maf->array[tal],4)*ep->matrix[0][tal];
    
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 2*pow(maf->array[tal],2)*maf2->array[tal]*ep->matrix[1][tal];
    retVal->matrix[0][tal] += 4*pow(maf->array[tal],3)*maf2->array[tal]*ep->matrix[1][tal];
  
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 0;
    retVal->matrix[0][tal] += 2*pow(maf->array[tal],2)*pow(maf2->array[tal],2)*ep->matrix[2][tal];
    
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 2*pow(maf->array[tal],2)*maf2->array[tal]*ep->matrix[3][tal];
    retVal->matrix[0][tal] += 4*pow(maf->array[tal],3)*maf2->array[tal]*ep->matrix[3][tal];
  
    retVal->matrix[2][tal] += 2*maf2->array[tal]*maf->array[tal]*ep->matrix[4][tal];
    retVal->matrix[1][tal] += maf2->array[tal]*maf->array[tal]*(maf2->array[tal]+maf->array[tal])*ep->matrix[4][tal];
    retVal->matrix[0][tal] += 4*pow(maf2->array[tal],2)*pow(maf->array[tal],2)*ep->matrix[4][tal];
  
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 2*pow(maf2->array[tal],2)*maf->array[tal]*ep->matrix[5][tal];
    retVal->matrix[0][tal] += 4*pow(maf2->array[tal],3)*maf->array[tal]*ep->matrix[5][tal];
  
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 0;
    retVal->matrix[0][tal] += 2*pow(maf2->array[tal],2)*pow(maf->array[tal],2)*ep->matrix[6][tal];
    
    retVal->matrix[2][tal] += 0;
    retVal->matrix[1][tal] += 2*pow(maf2->array[tal],2)*maf->array[tal]*ep->matrix[7][tal];
    retVal->matrix[0][tal] += 4*pow(maf2->array[tal],3)*maf->array[tal]*ep->matrix[7][tal];
  
    retVal->matrix[2][tal] += pow(maf2->array[tal],2.0)*ep->matrix[8][tal];
    retVal->matrix[1][tal] += pow(maf2->array[tal],3)*ep->matrix[8][tal];
    retVal->matrix[0][tal] += pow(maf2->array[tal],4)*ep->matrix[8][tal];
    
  }

  return retVal;
}




//geno is in this context the genotypes for an individual
dMatrix *error(iArray *geno, double epsilon, int snp){
  dMatrix *ep = allocDoubleMatrix(3,snp);
  for (int s=0;s<snp;s++){
    //pos 1
    if(geno->array[s]==1)
      ep->matrix[0][s] = pow((1-epsilon),2.0);
    if(geno->array[s]==2)
      ep->matrix[0][s] = (1-epsilon)*epsilon;
    if(geno->array[s]==3)
      ep->matrix[0][s] = pow(epsilon,2.0);
    //pos2
    if(geno->array[s]==1)
      ep->matrix[1][s] = 2*(1-epsilon)*epsilon;
    if(geno->array[s]==2)
      ep->matrix[1][s] = pow((1-epsilon),2.0) + pow(epsilon,2.0);
    if(geno->array[s]==3)
      ep->matrix[1][s] = 2*(1-epsilon)*epsilon;
    //pos 3
    if(geno->array[s]==1)
      ep->matrix[2][s] = pow((epsilon),2.0);
    if(geno->array[s]==2)
      ep->matrix[2][s] = (1-epsilon)*epsilon;
    if(geno->array[s]==3)
      ep->matrix[2][s] = pow((1-epsilon),2.0);
  }
	
  return ep;
}

/*
  return and matrix of dim 9 x snp
 */

dMatrix *errorNoLD(iArray *ind1,iArray *ind2,int snp,double epsilon){
  dMatrix *e1 = error(ind1,epsilon,snp); 
  dMatrix *e2 = error(ind2,epsilon,snp);
  
  dMatrix *ep = allocDoubleMatrix(9,snp);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++){
      for (int s=0;s<snp;s++){
	ep->matrix[i*3+j][s] = e1->matrix[i][s] * e2->matrix[j][s];
      }
    }
  killMatrix(e1);
  killMatrix(e2);
  return ep;
}

iArray *shiftLeft(int i, iArray *tmp){
  iArray *res = allocIntArray(tmp->x-i);
  for (int p=i;p<tmp->x;p++)
    res->array[p-i] = tmp->array[p]; 
  return res;
}

dMatrix *errorLD(iArray *ind1,iArray *ind2,iArray *ind1t,iArray *ind2t,int snp,double epsilon){
  iArray *leftE1 = shiftLeft(1,ind1);
  iArray *leftE2 = shiftLeft(1,ind2);
  dMatrix *e1 = error(leftE1 ,epsilon,snp); 
  dMatrix *e2 = error(leftE2,epsilon,snp);
  dMatrix *e1t = error(ind1t,epsilon,snp); 
  dMatrix *e2t = error(ind2t,epsilon,snp);
  
  dMatrix *ep = allocDoubleMatrix(81,snp);
  
  
  int num = 0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	for (int l=0;l<3;l++){
	  for (int s=0;s<snp;s++){
	    ep->matrix[num][s] =  e1->matrix[i][s] * e2->matrix[j][s]*  e1t->matrix[k][s] * e2t->matrix[l][s];
	    
	  }
	  num++;
	}
      }
    }
  }
  killArray(leftE1);
  killArray(leftE2);
  killMatrix(e1);
  killMatrix(e2);
  killMatrix(e1t);
  killMatrix(e2t);
  return ep;
}




dMatrix *LD(iArray* ind1,iArray* ind2,dArray *maf, dArray *maf2,double epsilon){
  dMatrix *distMatrix = allocDoubleMatrix(9,maf->x);
  int perm[9][2] = {{1,1},{1,2},{1,3}, {2,1},{2,2},{2,3},   {3,1}, {3,2},{3,3} } ;


  for(int i=0;i<maf->x;i++){
    for(int j=0;j<9;j++){
      distMatrix->matrix[j][i] = fabs(perm[j][0] - ind1->array[i]) + fabs(perm[j][1]-ind2->array[i]);
    }
  }
  dMatrix *returnMatrix=emission2(maf,maf2,epsilon,distMatrix);


  killMatrix(distMatrix);
  
  return returnMatrix;
}

dMatrix *LD2_new(iArray *choose,dArray *maf,iArray *ind1,iArray *ind2,dMatrix *pba,dMatrix *pBa,dMatrix* pbA, dMatrix *pBA,dMatrix *S1,double epsilon){

  //extract function skips first elem
  dArray *pBAt = extract(pBA,choose);//OK
  dArray *pbAt = extract(pbA,choose);//OK
  dArray *pBat = extract(pBa,choose);//OK
  dArray *pbat = extract(pba,choose);//OK
  iArray *ind1t = extract(ind1,choose);//ok
  iArray *ind2t = extract(ind2,choose);//ok
  dArray *maft = math_to(1,pBAt,1,pBat,0);
  dArray *maft2=one_minus(maft);
  dMatrix *ep2 = errorLD(ind1,ind2,ind1t,ind2t,ind1->x-1,epsilon);
  dMatrix *S2 = LDemission_new(pBAt,pbAt,pBat, pbat,ep2); //o
  dMatrix *ept = errorNoLD(ind1t,ind2t,ind1t->x,epsilon);
  dMatrix *S1t =emission2_new(maft,maft2,ept);//////
  dMatrix *result = allocDoubleMatrix(3,maf->x);
  result->matrix[0][0] = S1->matrix[0][0];
  result->matrix[1][0] = S1->matrix[1][0];
  result->matrix[2][0] = S1->matrix[2][0];
  for (int i=0;i<S2->y;i++)
    for (int j=0;j<3;j++){
      result->matrix[j][1+i] = S2->matrix[j][i]/ S1t->matrix[j][i];
    }
  killMatrix(ep2);
  killMatrix(ept);
  killMatrix(S1t);
  killMatrix(S2);
  killArray(pBAt);
  killArray(pBat);
  killArray(pbAt);
  killArray(pbat);
  killArray(ind1t);
  killArray(ind2t);
  killArray(maft);
  killArray(maft2);
  return result;
}


iArray *getRow(iMatrix *matrix,int n_row){
  iArray *retVal = allocIntArray(matrix->y);
  for(int i=0;i<matrix->y;i++)
    retVal->array[i]=matrix->matrix[n_row][i];
  return retVal;
}

dArray *calcMaf(iMatrix *matrix){
  dArray* returnArray= allocDoubleArray(matrix->y);

  //iterate through all rows
  
  for (int i=0;i<matrix->y;i++){
    int tmpSum=0;
    int numElems=0;
    int tmp;
    for (int j=0;j<matrix->x;j++){
      tmp = matrix->matrix[j][i];
      if (tmp==0)
	continue;
      else{
	tmpSum+=(tmp-1);
	numElems++;
      }
    }
    returnArray->array[i]=(((double)tmpSum)/numElems)/2;
  }
  return returnArray;
  
}

//merge the second into the first
void mergeLists(iArray* keepList,iArray *pruneList){
  
  int atPos = 0;
  //printf("keeplist len:%d\tprunelist len:%d\n",keepList->x,pruneList->x);

  for (int i=0;i<keepList->x;i++)
    if(keepList->array[i]==1){
      if(pruneList->array[atPos]==1)
	keepList->array[i] = 1;
      else{
	keepList->array[i] = 0;
      }
      atPos++;
    }
}

/*

 */

bArray *merge_bArrays(bArray* first,bArray *second) {
  bArray *retVal = allocBoolArray(first->x);
  retVal->numTrue = second->numTrue;
 
  int atPos=0;
  for (int i=0 ; i<first->x ; i++)
    if(first->array[i]==1){
      if(second->array[atPos]==1)
	retVal->array[i] = 1;
      else
	retVal->array[i] = 0;
      atPos++;
    }
  return retVal;
}




pars *getPars(iMatrix *data, dArray *position,int pair1,int pair2,double min \
	      ,double epsilon,int back,int ld_choose,int ld_adj,int doPrune,int double_recom,double prune_val) {
  
  //first check if ind1 SNP's, is non NA, on atleast on loci with ind2;
  //bugfix for 0.71 version
  for(int i=0;i<data->y;i++){
    if(data->matrix[pair1][i]!=0 && data->matrix[pair2][i]!=0){
      break;
    }
    if(data->y-1==i){
      //if we haven't breaked out of the loop efter the last loci
      //then we should exit
      printf("\t-> Problem in genotype dataset. pair=(%d,%d) doesn't share 1 non-missing SNP. will exit\n",pair1,pair2);
      exit(0);
    }
  }


  dMatrix *S=NULL;
  int ind,snp;
  //these are used by both prune and non prune calculations
  dArray *mafarray = calcMaf(data);// ok //remember R version uses missing=NA

  iArray *keepList = allocIntArray(mafarray->x);
  iArray *okList=getOKindices(data,mafarray,min,pair1,pair2,keepList);//using maf
  if(print_info>1)
    printf("\t-> NumSnps after stripping, mafmin AND NA's:%d\n",okList->x);

  iMatrix *strippedData=extractOK(okList,data);
  
  //these are specific to either prune or non prune
  dArray *maf2=NULL;
  iArray *choose=NULL;
  iArray *ind1=NULL;
  iArray *ind2=NULL;
  
  dArray *maf=NULL;
  dArray *pos=NULL;
  
  //placeholder user for substitution
  dArray  *dArrayTmp=NULL;
  iArray  *iArrayTmp=NULL;
  iMatrix *iMatrixTmp=NULL;

  //return value structures;
  pars *returnVal = new pars();
  returnVal->choose=NULL;
  hapStruct *hap = new hapStruct(); 
  hap->pBA=NULL;  hap->pBa=NULL;  hap->pbA=NULL;  hap->pba=NULL;

  //do maf strip
  
  maf2=extractOK(okList,mafarray);
  maf=one_minus(maf2);
  pos=extractOK(okList,position);
  
  ind1 = getRow(strippedData,pair1);//ok
  ind2 = getRow(strippedData,pair2);//ok
    


  if(doPrune){
    if(print_info>1)
      flush_print("Will begin pruning procedure\n");
    iMatrix *tmp = revCols(strippedData);
    if(print_info>1)
      flush_print("Calculating ld patterns might take some time...");
    snpMatrix *ld = snp_pair_range(tmp,back);
    
    iArrayTmp = pruning(ld,ld_choose,back,prune_val);
    
    killArray(okList);
    okList = generateIndices(iArrayTmp);
    if(print_info>1)
      printf("\t-> Number of snps after pruning:%d\n",okList->x);
    
    dArrayTmp = extractOK(okList,mafarray);

    killArray(mafarray);
    mafarray = dArrayTmp;
    
    killArray(maf2);
    killArray(maf);
    
    //copy array
    maf2=allocDoubleArray(mafarray);
    //maf2=mafarray;
    maf = one_minus(mafarray);

    iMatrixTmp = extractOK(okList,strippedData);
    killMatrix(strippedData);
    strippedData = iMatrixTmp;
    killArray(ind1);
    killArray(ind2);
    ind1=getRow(strippedData,pair1);
    ind2=getRow(strippedData,pair2);
    
    
    
    dArrayTmp = extractOK(okList,pos);
    

    mergeLists(keepList, iArrayTmp);
    killArray(pos);
    pos = dArrayTmp;
    killSnpMatrix(ld);
    killMatrix(tmp);
    killArray(iArrayTmp);
  }

  ind = strippedData->x;
  snp = strippedData->y;

  dMatrix *ep1 = errorNoLD(ind1,ind2,strippedData->y,epsilon);//ok to the 10^-7 digit

  dMatrix *S1 = emission2_new(maf,maf2,ep1);

  if (ld_adj&&back>0) {
    
    iMatrix *tmp = revCols(strippedData);
    snpMatrix *ld = snp_pair_range(tmp,back);
    dMatrix *pba=revCols_and_extend(ld->pba);
    dMatrix *pBa=revCols_and_extend(ld->pbA);
    dMatrix *pbA=revCols_and_extend(ld->pBa);
    dMatrix *pBA=revCols_and_extend(ld->pBA);
    dMatrix *mea;
    if(ld_choose)
      mea = revCols(ld->rmisc);
    else{
      mea = revCols(ld->D);
      myAbs(mea);
    }
    
    
    //remove after use
    killMatrix(tmp);
    killSnpMatrix(ld);
    
    choose = getHighestId(mea) ;
    
    //new ld calculation
    S=LD2_new(choose,maf,ind1,ind2,pba,pBa,pbA,pBA,S1,epsilon);
    hap->mea=mea;
    hap->pba=pba;
    hap->pBa=pBa;
    hap->pbA=pbA;
    hap->pBA=pBA;
    returnVal->choose=choose;
   
   
  }
  else
    S=S1;
  //
  returnVal->hap=hap;
  //initialize function pars
  // PIK PIK PIK these are global for now
  returnVal->t  = diff(pos);//ok
  returnVal->S = S;
  returnVal->S1 = S1;
  returnVal->maf = maf;
  returnVal->maf2 = maf2;
  
  returnVal->ind1 = ind1;
  returnVal->ind2 = ind2;
  returnVal->pos = pos;
  returnVal->keepList = keepList;
  //now clean up
  killMatrix(ep1);
  //killArray(maf2);
  killArray(mafarray);
  killArray(okList);
  killMatrix(strippedData);
  return returnVal;
}
 

iMatrix *extract_cols_of_iMatrix(iMatrix *inData,bArray *sel){
  iMatrix *outData = allocIntMatrix(inData->x,sel->numTrue);
  int pos = 0;
  for(int i=0;i<inData->y;i++)
    if(sel->array[i]){
      for (int j=0;j<inData->x;j++)
	outData->matrix[j][pos] = inData->matrix[j][i];
      pos++;
    }
  return outData;
}

bArray *maf_stripper(dArray *mafs,double min){
  bArray *retVal = allocBoolArray(mafs->x);
  int numTrues = 0;
  for(int i=0 ; i<mafs->x ; i++)
    if(mafs->array[i]!=0 && mafs->array[i] > min && mafs->array[i] < (1-min)){
      retVal->array[i]=1;
      numTrues++;
    }
    else
      retVal->array[i] = 0;
  retVal->numTrue = numTrues;
  return retVal;
}


void print_bArray(char *st,bArray *ba){
  printf("\n %s : %d out of %d\n",st,ba->numTrue,ba->x);
  for (int i=0;i<ba->x;i++)
    cout << ba->array[i] << " ";
  cout << endl;
       
}

void print_iMatrix(char *st,iMatrix *dm,char *file){
  if(file==NULL){
    printf("\n %s : (%d,%d)\n",st,dm->x,dm->y);
    for (int i=0;i<dm->x;i++){
      for (int j=0;j<dm->y;j++)
	cout << dm->matrix[i][j] << " ";
      cout << endl;
    }
  } else{
    ofstream myfile;
    myfile.open(file);
    printf("\n %s : (%d,%d)\n",st,dm->x,dm->y);
    for (int i=0;i<dm->x;i++){
      for (int j=0;j<dm->y;j++)
	myfile << dm->matrix[i][j] << " ";
      myfile << endl;
    }
    myfile.close();
  }
}

bArray *extend_bArray_with_one(const bArray * in){
  bArray *out = allocBoolArray(in->x+1);
  out->numTrue=in->numTrue+1;
  out->array[0] = 1;
  for(int i=0 ; i<in->x ; i++)
    out->array[i+1] = in->array[i];
  return out;

}


void relateHMM::pre_calc(const functionPars *pars) {

  if(print_info>1)
    flush_print("Will begin pre-calculation\n");


  /*
    First we calculate the maf for all loci
   */
  dArray *tmpDarray = calcMaf(pars->data);
  bArray *mafList = maf_stripper(tmpDarray,pars->min);
  killArray(tmpDarray);
  iMatrix *tmpiMatrix= extract_cols_of_iMatrix(pars->data,mafList);
    
  if(print_info>1)
    printf("\t-> NumSnps after mafstripping:%d\n",mafList->numTrue);
    
  /*
    Should we also strip some loci because of ld
   */
  if(pars->doPrune){
    if(print_info>1){
      flush_print("Will begin pruning\n");
      flush_print("Calculating ld patterns, might take some time...");
    }
    dMatrix *meas = getMea(tmpiMatrix,pars->LD,pars->back);
    bArray *pruneList = pruning3(meas,pars->back,pars->prune_val);
    if(print_info>1)
      printf("\t-> Number of snps after pruning:%d\n",pruneList->numTrue);
    genos_global = extract_cols_of_iMatrix(tmpiMatrix,pruneList);
    pre_calc_used_list = merge_bArrays(mafList,pruneList);
    killArray(mafList);
    killArray(pruneList);
    killMatrix(tmpiMatrix);
    //killMatrix(meas); 
  }else{
    genos_global = tmpiMatrix;
    pre_calc_used_list = mafList;
  }
  
  /*
    We don't want to calculate things more than once.
    So we'l keep a global maf1, maf2 and haplos/measure
  */
  maf1_global = calcMaf(genos_global);
  maf2_global = one_minus(maf1_global);
  pos_global = extractOK(pre_calc_used_list,pars->position);
  if(pars->ld_adj && pars->back>0){
    iMatrix *reversed = revCols(genos_global);
    snpMatrix *ld = snp_pair_range(reversed,pars->back2);
    hap_global = new hapStruct();
  //lets also input mea directly
    if(pars->LD)
      hap_global->mea  = revCols(ld->rmisc);
    else{
      hap_global->mea = revCols(ld->D); //prepare
      myAbs(hap_global->mea); 
    }
    hap_global->pba = revCols_and_extend(ld->pba);
    hap_global->pBa = revCols_and_extend(ld->pbA); //this is correct
    hap_global->pbA =  revCols_and_extend(ld->pBa); //this is correct
    hap_global->pBA = revCols_and_extend(ld->pBA);
    
    killMatrix(reversed);
    killSnpMatrix(ld);
    
  }
  if(print_info>1)
    puts("\t-> pre-calculation is over");
}

iArray *bArray_to_iArray(bArray *in){
  iArray *out = allocIntArray(in->x);
  for(int i=0;i<in->x;i++)
    out->array[i] = in->array[i];
  return out;
}

void relateHMM::init_globals(const functionPars *pars){

  hap_global = NULL;  
  print_info = pars->print_results;
  if(print_info){
    
    printFunctionPars(pars);

  }
  alim=pars->alim;
  
  
  
  fixA=pars->fixA;
  fixA_val = pars->fixA_val;
  fixK2=pars->fixK;
  fixK2_val=pars->fixK_val;
  double_recom = pars->double_recom;
  calcA = pars->calcA;
  phi = pars->phi;
  pre_calc(pars);
}

void relateHMM::dinit_globals(int back){
  killMatrix(genos_global);  
  if(back>0)
    killHapStruct(hap_global);

  killArray(maf1_global);
  killArray(maf2_global);
  killArray(pre_calc_used_list);
  killArray(pos_global);
 
}



fullResults *relateHMM::sgl(const functionPars *pars) {
  
  //now everything is prepared
  double delta[3] = {0,0,0};

  
  //check first check if individuals share atleast one non-NA loci
  for(int i=0;i<genos_global->y;i++) {
    if(genos_global->matrix[pars->pair->array[0]][i]!=0 && genos_global->matrix[pars->pair->array[1]][i]!=0){
      break;
    } if(genos_global->y-1==i){
      printf("\t-> Problem in genotype dataset. pair=(%d,%d) doesn't share 1 non-missing SNP. will exit\n",pars->pair->array[0],pars->pair->array[1]);
      exit(0);
    }
  }
  //now check which loci are NA
  bArray *NA_keepList = which_arent_NA(genos_global,pars->pair->array[0],pars->pair->array[1]);

  if(print_info >1)
    printf("\t-> Number of snps after NA stripped: %d\n",NA_keepList->numTrue);
  
  //now check if we should remove more snp's due to ld
  //if(pars->adj&&pars->back>0){

  
  iMatrix *tmpiMatrix = extract_cols_of_iMatrix(genos_global,NA_keepList);
  iArray *ind1 = getRow(tmpiMatrix,pars->pair->array[0]);//ok
  iArray *ind2 = getRow(tmpiMatrix,pars->pair->array[1]);//ok
  bArray *COKE_ROX = merge_bArrays(pre_calc_used_list,NA_keepList);
  dArray *pos_local = extractOK(NA_keepList,pos_global); 
  t = diff(pos_local);
  dMatrix *ep1 = errorNoLD(ind1,ind2,tmpiMatrix->y,pars->epsilon);
  dArray *maf2_local=extractOK(NA_keepList,maf1_global);
  dArray *maf1_local  = one_minus(maf2_local); 
  dMatrix *S1 = emission2_new(maf1_local,maf2_local,ep1);
  killMatrix(tmpiMatrix);

  iArray *choose = NULL;
  dMatrix* S = NULL;
  dMatrix *mea_local = NULL; dMatrix *pba_local = NULL;
  dMatrix *pbA_local = NULL; dMatrix *pBa_local = NULL;
  dMatrix *pBA_local = NULL;
  if(pars->back >0 && pars->ld_adj){
    if(print_info>1)
      printf("\t-> Will now remove these loci and update ld table\n");
    bArray *extended_NA_keepList = extend_bArray_with_one(NA_keepList);
    mea_local = remove_cross(hap_global->mea,NA_keepList,pars->back);
    pba_local = remove_cross(hap_global->pba,extended_NA_keepList,pars->back);
    pbA_local = remove_cross(hap_global->pbA,extended_NA_keepList,pars->back);
    pBa_local = remove_cross(hap_global->pBa,extended_NA_keepList,pars->back);
    pBA_local = remove_cross(hap_global->pBA,extended_NA_keepList,pars->back);
    
    //  dMatrix *S1=LD(ind1,ind2,maf,maf2,epsilon);// ok
    //ok to the 10^-8 digit
    killArray(extended_NA_keepList);
    choose =  getHighestId(mea_local) ;
  
    S = LD2_new(choose,maf1_local,ind1,ind2,pba_local,pBa_local,pbA_local,pBA_local,S1,pars->epsilon);
   
  }else
    S = S1;
  
  Sk1 = row_log(S,2);
  Sk2 = row_log(S,1);
  Sk3 = row_log(S,0);  


   //now calculate the u.like and po.like
 

  delta[0]=alim->array[0];
  double uLike = like(delta,Sk1,Sk2,Sk3,t,pars->double_recom,alim);//,0,0,0,0,0,NULL);
  delta[2]=1;
  double poLike = like(delta,Sk1,Sk2,Sk3,t,pars->double_recom,alim);//,0,0,0,0,0,NULL);
  if(print_info)
    printf("\t->uLike:%f\t poLike:%f",uLike,poLike);
  if(pars->doPar){
    delta[0]=pars->par->array[0];
    delta[1]=pars->par->array[1];
    delta[2]=pars->par->array[2];
  }  



  /*
    :now do optimization:
    first  par: is number of times to converge to same value before exiting
    second par: is the max number of times to run
    second par: are the bounds of the first parameter 2 optimize
  */

  dArray *bestLike;
  //dArray *likeVal;//used if want the point likelihood.
  //return NULL;
  if(!pars->doPar){
    bestLike = run_optimization(pars->timesToConverge,pars->timesToRun,pars->convTol,alim);
  }
  else{
    bestLike=allocDoubleArray(4);
    for(int i=0;i<3;i++)
      bestLike->array[i]=delta[i];
    bestLike->array[3] = like(delta,Sk1,Sk2,Sk3, t,double_recom,alim);//,0,0,0,0,0,NULL);
  }

  //now run decode and find the viterbi path
  dMatrix *decodeResult;
  if(bestLike->array[1]!=0)//if k2!=0
    decodeResult = decodek2(bestLike,Sk1,Sk2,Sk3,t,pars->double_recom);
  else
    decodeResult = decode(bestLike,Sk1,Sk2,Sk3,t,pars->double_recom);

  iArray *path = viterbi(bestLike,Sk1,Sk2,Sk3,t,pars->double_recom);
  
  //now prepare return values, and input them to R
  dArray *kResults=allocDoubleArray(3);
  kResults->array[0]=bestLike->array[1];
  kResults->array[1]=bestLike->array[2];
  kResults->array[2]=1-kResults->array[0]-kResults->array[1];
  
  fullResults *res = new fullResults();
  hapStruct *hap = new hapStruct(); 
  hap->pba=pba_local;
  hap->pBa=pBa_local;
  hap->pbA=pbA_local;
  hap->pBA=pBA_local;
  hap->mea=mea_local;
  res->timesRun = localTimesRun;
  res->timesConverged = localTimesConv;
  res->kLike = bestLike->array[3];
  res->S = S;
  res->a = bestLike->array[0];
  res->uLike = uLike;
  res->LD = pars->LD;
  res->t = t;
  res->snp = S->y;
  res->position = pos_local;
  res->poLike = poLike;
  res->post = decodeResult;
  res->kResult = kResults;
  res->kr = kResults->array[0]/4.0+kResults->array[1]/2.0;
  res->double_recom = pars->double_recom;
  res->alim=alim;
  res->choose=choose;
  res->back = pars->back;
  res->hap = hap;
  res->S1 = S1;
  res->maf = maf1_local;
  
  res->path= path;
  //res->keepList = bArray_to_iArray(COKE_ROX);
  res->keepList = bArray_to_iArray(NA_keepList);
  res->chr = extractOK(COKE_ROX,pars->chr);
  
  if(pars->back>0){
    res->mea_global = hap_global->mea;
    res->pba_global = hap_global->pba;
    res->pBa_global = hap_global->pBa;
    res->pbA_global = hap_global->pbA;
    res->pBA_global = hap_global->pBA;
  }
  //if we have haven't done any optimization, then we don't have an optimization list
  if(pars->doPar)
    res->convInfo=NULL;
  else
    res->convInfo = convInfo;
  if(print_info>1)
    convInfo->print("convergenceinfo is:",NULL);
  // clean up
  //  killMatrix(decodeResult);
  killArray(bestLike);
  //  delete pars;
  killArray(Sk1);
  killArray(Sk2);
  killArray(Sk3);
  killMatrix(ep1);
  killArray(NA_keepList);
  killArray(ind1);
  killArray(ind2);
  killArray(COKE_ROX);
  killArray(maf2_local);
  
  return res;

}



fullResults *relateHMM::single_pair(const functionPars *pars) {


  print_info = pars->print_results;
  if(print_info){
    printFunctionPars(pars);
  }
  alim=pars->alim;
  
  double delta[3] = {0,0,0};
  
  fixA=pars->fixA;
  fixA_val = pars->fixA_val;
  fixK2=pars->fixK;
  fixK2_val=pars->fixK_val;
  double_recom = pars->double_recom;
  calcA = pars->calcA;
  phi = pars->phi;
  

  /*
    will now do all the prelim calculations
    check relateHMM.h for full list of included values
  */

  //
  
  calculatedValues = getPars(pars->data,pars->position,pars->pair->array[0],pars->pair->array[1], \
	    pars->min,pars->epsilon,pars->back,pars->LD,pars->ld_adj,pars->doPrune,pars->double_recom,pars->prune_val);


  Sk1 = row_log(calculatedValues->S,2);
  Sk2 = row_log(calculatedValues->S,1);
  Sk3 = row_log(calculatedValues->S,0);
  t = calculatedValues->t;
  
   //now calculate the u.like and po.like
  delta[0]=alim->array[0];
  double uLike = like(delta,Sk1,Sk2,Sk3,t,pars->double_recom,alim);//,0,0,0,0,0,NULL);
  delta[2]=1;
  double poLike = like(delta,Sk1,Sk2,Sk3,t,pars->double_recom,alim);//,0,0,0,0,0,NULL);
  if(print_info)
    printf("\t->uLike:%f\t poLike:%f",uLike,poLike);
  if(pars->doPar){
    delta[0]=pars->par->array[0];
    delta[1]=pars->par->array[1];
    delta[2]=pars->par->array[2];
  }  

  /*
    :now do optimization:
    first  par: is number of times to converge to same value before exiting
    second par: is the max number of times to run
    second par: are the bounds of the first parameter 2 optimize
  */

  dArray *bestLike;
  //dArray *likeVal;//used if want the point likelihood.
  //return NULL;
  if(!pars->doPar){
    bestLike = run_optimization(pars->timesToConverge,pars->timesToRun,pars->convTol,alim);
  }
  else{
    bestLike=allocDoubleArray(4);
    for(int i=0;i<3;i++)
      bestLike->array[i]=delta[i];
    bestLike->array[3] = like(delta,Sk1,Sk2,Sk3, t,double_recom,alim);//,0,0,0,0,0,NULL);
  }

  //now run decode and find the viterbi path
  dMatrix *decodeResult;
  if(bestLike->array[1]!=0)//if k2!=0
    decodeResult = decodek2(bestLike,Sk1,Sk2,Sk3,t,pars->double_recom);
  else
    decodeResult = decode(bestLike,Sk1,Sk2,Sk3,t,pars->double_recom);

  iArray *path = viterbi(bestLike,Sk1,Sk2,Sk3,t,pars->double_recom);
  
  //now prepare return values, and input them to R
  dArray *kResults=allocDoubleArray(3);
  kResults->array[0]=bestLike->array[1];
  kResults->array[1]=bestLike->array[2];
  kResults->array[2]=1-kResults->array[0]-kResults->array[1];
  
  fullResults *res = new fullResults();
  res->timesRun = localTimesRun;
  res->timesConverged = localTimesConv;
  res->kLike = bestLike->array[3];
  res->S = calculatedValues->S;
  res->a = bestLike->array[0];
  res->uLike = uLike;
  res->LD = pars->LD;
  res->t = t;
  res->snp = calculatedValues->maf->x;
  res->position = calculatedValues->pos;
  res->poLike = poLike;
  res->post = decodeResult;
  res->kResult = kResults;
  res->kr = kResults->array[0]/4.0+kResults->array[1]/2.0;
  res->double_recom = pars->double_recom;
  res->alim=alim;
  res->choose=calculatedValues->choose;
  res->back = pars->back;
  res->hap = calculatedValues->hap;
  res->S1 = calculatedValues->S1;
  res->maf = calculatedValues->maf;

  

  iArray *tskTmp = generateIndices(calculatedValues->keepList);
  res->chr = extractOK(tskTmp,pars->chr);
  killArray(tskTmp);
  
  res->path= path;
  res->keepList = calculatedValues->keepList;
  //if we have haven't done any optimization, then we don't have an optimization list
  if(pars->doPar)
    res->convInfo=NULL;
  else
    res->convInfo = convInfo;
  if(print_info>1&&!pars->doPar)
    convInfo->print("convergenceinfo is:",NULL);


  // clean up
  //  killMatrix(decodeResult);
  killArray(bestLike);
  //  delete pars;
  killArray(Sk1);
  killArray(Sk2);
  killArray(Sk3);



  
  //  delete decPars;

  killArray(calculatedValues->ind1);
  killArray(calculatedValues->ind2);
  killArray(calculatedValues->maf2);
  delete calculatedValues;
  return res;
}

int relateHMM::fast_all_pairs(functionPars *pars){
   printf("\t-> Will dump these data, that are used for analysis after the maf/prune removal\n");

   printf("\t\tgeno:\t\"%s\"\n\t\tchromo:\t\"%s\"\n\t\tpos:\t\"%s\"\n\t\tkeepList:\"%s\"\n",\
	  pars->stripped_geno.c_str(),pars->stripped_chr.c_str(),pars->stripped_pos.c_str(),pars->stripped_keeplist.c_str());
  
   genos_global->print(NULL,pars->stripped_geno.c_str());
   pos_global->print(NULL,pars->stripped_pos.c_str());
  iArray *tmp = extractOK(pre_calc_used_list,pars->chr);
  tmp->print(NULL,pars->stripped_chr.c_str());
  killArray(tmp);
  pre_calc_used_list->print(NULL,pars->stripped_keeplist.c_str());


  int numIndividuals = pars->data->x;
  fullResults *res = NULL;
  

  pars->postFilename = update_filename(pars->postFilename);
  pars->kFilename = update_filename(pars->kFilename);
  int numberOfPairsToRun = numIndividuals*(numIndividuals-1)/2;
  printf("\t-> number of individuals:%d\n",numIndividuals);
  printf("\t-> number of pairs to run:%d\n",numberOfPairsToRun);
  fflush(stdout);
  int numTimes=0;
  for (int p1=0;p1<numIndividuals-1;p1++)
    for(int p2=p1+1 ; p2<numIndividuals; p2++ ){
      printf("\r\t-> ind1:%d\tind2:%d   %d/%d ",p1,p2,numTimes,numberOfPairsToRun); //space after final int, is important otherwise screen output wont look nice
      fflush(stdout);
      pars->pair->array[0] = p1;
      pars->pair->array[1] = p2;
      res = sgl(pars);
      if(pars->fixK==1 && pars->fixK_val==0){
	//printf("The special case...");
	
	write_dMatrix_to_file(pars->postFilename,res->post,2,res->keepList);//this is the third column,will append
	write_dArray_to_file(pars->kFilename,res->kResult,1);//this is the second elment,will append
      }else{
	//printf("The general case...");
	write_dMatrix_to_file(pars->postFilename,res->post,2,res->keepList);//this is the third column,will append
	write_dMatrix_to_file(pars->postFilename,res->post,1,res->keepList);//this is the sec column,will append
	write_dArray_to_file(pars->kFilename,res->kResult,1);//this is the second elment,will append
	write_dArray_to_file(pars->kFilename,res->kResult,0);//this is the first elment,will append
      }
      killFullResults(res);
      numTimes++;
    }
  
  
  printf("\n");
  return 0;
}

int relateHMM::all_pairs(functionPars *pars) {
  // print_info = 0; // COKE ROX pars->print_results;
  int numIndividuals = pars->data->x;
  //  numIndividuals = 4;
  fullResults *res = NULL;
  printf("\t-> number of individuals:%d\n",numIndividuals);
  pars->postFilename = update_filename(pars->postFilename);
  pars->kFilename = update_filename(pars->kFilename);

  for (int p1=0;p1<numIndividuals-1;p1++)
    for(int p2=p1+1 ; p2<numIndividuals; p2++ ){
      printf("\r\t->ind1:%d\tind2:%d ",p1,p2); //space after final int, is important otherwise
      fflush(stdout);
      pars->pair->array[0] = p1;
      pars->pair->array[1] = p2;
      res = single_pair(pars);
      if(pars->fixK==1 && pars->fixK_val==0){
	//printf("The special case...");
	write_dMatrix_to_file(pars->postFilename,res->post,2,res->keepList);//this is the third column,will append
	write_dArray_to_file(pars->kFilename,res->kResult,1);//this is the second elment,will append
      }else{
	//printf("The general case...");
	write_dMatrix_to_file(pars->postFilename,res->post,2,res->keepList);//this is the third column,will append
	write_dMatrix_to_file(pars->postFilename,res->post,1,res->keepList);//this is the sec column,will append
	write_dArray_to_file(pars->kFilename,res->kResult,1);//this is the second elment,will append
	write_dArray_to_file(pars->kFilename,res->kResult,0);//this is the first elment,will append
      }
      

      killFullResults(res);
      
    }

  printf("\n");
  return 0;
}


int relateHMM::joblist_fast_all_pairs(functionPars *pars,vector<iArray*> pa){
   printf("\t-> Will dump these data, that are used for analysis after the maf/prune removal\n");
  printf("\t\tgeno:\t\"%s\"\n\t\tchromo:\t\"%s\"\n\t\tpos:\t\"%s\"\n\t\tkeepList:\"%s\"\n",GENOTYPE_FILE,CHROMOSOMES_FILE,POSITIONS_FILE,KEEPLIST_FILE);
  
  genos_global->print(NULL,GENOTYPE_FILE);
  pos_global->print(NULL,POSITIONS_FILE);
  iArray *tmp = extractOK(pre_calc_used_list,pars->chr);
  tmp->print(NULL,CHROMOSOMES_FILE);
  killArray(tmp);
  pre_calc_used_list->print(NULL,KEEPLIST_FILE);




  flush_print("Will run a specified joblist\n");
  //int numIndividuals = pars->data->x;
  fullResults *res = NULL;
  flush_print("number of pairs to run: ",(int)pa.size());
  printf("\n");
  pars->postFilename = update_filename(pars->postFilename);
  pars->kFilename = update_filename(pars->kFilename);
  for(unsigned int i=0;i<pa.size();i++){
    killArray(pars->pair);
    pars->pair = pa[i];
    printf("\r\t-> ind1:%d\tind2:%d  %d/%d",pars->pair->array[0],pars->pair->array[1],i+1,(int)pa.size()); //space after final int, is important otherwise
    fflush(stdout);
    res = sgl(pars);
    if(pars->fixK==1 && pars->fixK_val==0){
      //printf("The special case...");
      write_dMatrix_to_file(pars->postFilename,res->post,2,res->keepList);
      //this is the third column,will append
      write_dArray_to_file(pars->kFilename,res->kResult,1);
      //this is the second elment,will append
    }else{
      //printf("The general case...");
      write_dMatrix_to_file(pars->postFilename,res->post,2,res->keepList);
      //this is the third column,will append
      write_dMatrix_to_file(pars->postFilename,res->post,1,res->keepList);
      //this is the sec column,will append
      write_dArray_to_file(pars->kFilename,res->kResult,1);
      //this is the second elment,will append
      write_dArray_to_file(pars->kFilename,res->kResult,0);
      //this is the first elment,will append
    }
      killFullResults(res);
  }
  
  
  printf("\n");
  return 0;
}


