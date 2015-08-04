/*
  Copyright (C) 2009 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.987
  
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

#include <iostream>
#include <fstream>
#include "relateHMM.h"
#include <string>
#include <iomanip>//used for setting number of decimals in posterior
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <vector>
#include <sys/stat.h>
#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif

#include "asort.h"

#include "conf.h"
#include "extractors.h"
#include "filereader_and_conversions.h"

using namespace std;



void wrt_dMatrix_to_file(string str, dMatrix *mat,int row,fullResults *res){
  ofstream myfile;
  myfile.open (str.c_str(),ios::app);
  int inPlace = 0;
  //  cout <<"size of keeplist:"<<res->keepList->x;
  for (int i=0;i < res->keepList->x ;i++){//loop through the keeplist
    if (res->keepList->array[i]==1){ //if snp is to be include write teh result from posterier
      myfile <<setprecision(POSTPRECISION)<< mat->matrix[row][inPlace]<<"\t";
      inPlace++;
    }
    else
      myfile << "-1" << "\t";
  }
  myfile << endl;
  myfile.close();
}
//will append an element to a file
void wrt_dArray_to_file(string str, dArray *mat,int element){
  ofstream myfile;
  myfile.open (str.c_str(),ios::app);
  myfile << mat->array[element] << endl;;
  myfile.close();
}



void post_processing_of_option_pars(functionPars *pars,int myChoose){
  if(pars->calcA && pars->fixK && pars->fixK_val==0 && pars->timesToRun!=1 ){
    printf("\t-> Calc.A and fixK=0 is a one-dim optimization\n");
    printf("\t-> So will only need to run optimization once.\n");
    pars->timesToRun = 1;
  }
  /*//0.96
  if(pars->back2 < pars->back){
    printf("\t-> back2 should be greater than back\n");
    printf("\t-> back2 will now be set twice the value of back\n");
    pars->back2 = pars->back*2;
    }*/
  pars->back2= pars->back*2;
  if(pars->doPrune && pars->back==0){
    printf("\t-> Must supply a back > 0, when you want to prune\n");
    exit(0);

  }
  /* //0.96    
  if(pars->print_results)
    if(myChoose)
      printf("\t-> Will use the old/slow method                                        \n");
    else
      printf("\t-> Will use the new/fast method                                        \n");
  */
}


void setDefaultParsRun(functionPars *pars){

  iArray *pair = allocIntArray(2);
  dArray *alim = allocDoubleArray(2);
  dArray *par = allocDoubleArray(3);
  pars->data=NULL;
  pars->pair= pair;
  pars->chr=NULL;
  pars->position=NULL;
  pars->par=par;
  pars->alim=alim;

  pars->alim->array[0]=0.001;
  pars->alim->array[1]=0.15;
  pars->pair->array[0]=0;
  pars->pair->array[1]=1;
  pars->doPar=0;
  pars->min=0;
  pars->double_recom=0;
  pars->LD=1; //ld=1=rsq, ld=0=D
  pars->ld_adj=1; 
  pars->epsilon=0.01;
  pars->back=5;
  pars->doPrune = 0;
  pars->fixA=0;
  pars->fixK=1;
//  pars->abs=0;
  pars->calcA=1;
  pars->fixA_val=0.0;
  pars->fixK_val=0;
  pars->prune_val = 0.2;
  pars->phi = 0.013;
  pars->convTol = 0.1;
  pars->timesToConverge=2;
  pars->timesToRun=3;
  pars->print_results=1;
  pars->back2= 10;// COKE_ROX

  //set default stripped filenames, these are defined in conf.h
  pars->stripped_geno = GENOTYPE_FILE;
  pars->stripped_pos = POSITIONS_FILE;
  pars->stripped_chr = CHROMOSOMES_FILE;
  pars->stripped_keeplist = KEEPLIST_FILE;
}

void info(){
  printf("\t--------------\n");
  printf("%s\n",FLATFILES);
  printf("%s\n",PLINKFILES);
  printf("%s",REQ);
  printf("%s\n",PREINFO);
  printf("%s\n",INPUTFILES);
  printf("%s\n",ALLPAIRS);
  printf("%s\n",OUTPUTFILES);
  printf("%s\n",CONTACT);

}


int main(int argc, char *argv[]){

  for(int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  string geno= "";
  string pos = "";
  string chr = "";
  

  string option = "options.txt";
  string postFilename = "postoutput.post";
  string kFilename = "koutput.k";


  int myChoose =1;
  int argPos=1;

  //added in 0.83
  int genoTaken=0;
  int posTaken=0;
  int chrTaken=0;
  int individualsToTestList=0;
  int doPlink=0;//added in 0.97;
  string plink_fam;//added in 0.97;
  string plink_bim;//added in 0.97;
  string plink_bed;//added in 0.97;
  bArray *plinkKeep = NULL; //added in 0.97;
  vector<iArray*> re;//added in 0.93
  if(argc==1){
    info();
    return 0;
  }
  functionPars *pars =  new functionPars();
  setDefaultParsRun(pars);
  while(argPos <argc){
    if (strcmp(argv[argPos],"-o")==0){
      option  = argv[argPos+1]; 
    }
    else if(strcmp(argv[argPos],"-g")==0){
      geno  = argv[argPos+1]; 
      genoTaken=1;
    }
    else if(strcmp(argv[argPos],"--help")==0){
      info();
      return 0;
    }
    else if(strcmp(argv[argPos],"-p")==0){
      pos  = argv[argPos+1];
      posTaken=1;
    }
    else if(strcmp(argv[argPos],"-c")==0){
      chr  = argv[argPos+1];
      chrTaken=1;
    }    
    else if(strcmp(argv[argPos],"-post")==0){
       postFilename = argv[argPos+1]; 
    }    
    else if(strcmp(argv[argPos],"-k")==0){
       kFilename = argv[argPos+1]; 
    }
    else if(strcmp(argv[argPos],"-v")==0){
      pars->print_results = atoi(argv[argPos+1]);
     
    }
    else if(strcmp(argv[argPos],"-plink-bed")==0){
      plink_bed = (argv[argPos+1]);
      doPlink++;
    }
    else if(strcmp(argv[argPos],"-plink-bim")==0){
      plink_bim = argv[argPos+1];
      doPlink++;
    }else if(strcmp(argv[argPos],"-out")==0){
      std::string prefix = std::string(argv[argPos+1]);
      std::cout<<"prefix:" << prefix<<std::endl;
      postFilename = prefix+".post";
      kFilename = prefix+".k";
      pars->stripped_geno = prefix + "_stripped.geno";
      pars->stripped_chr = prefix + "_stripped.chr";
      pars->stripped_pos = prefix + "_stripped.pos";
      pars->stripped_keeplist = prefix + "_stripped.keep";

    }
    else if(strcmp(argv[argPos],"-file")==0){
      std::string p_str =string( argv[argPos+1]);
      if(p_str.length()>4){
	std::string ext = p_str.substr(p_str.length()-4,p_str.length());
	if (!ext.compare(".chr")||!ext.compare(".geno")||!ext.compare(".pos")){
	  std::string front = p_str.substr(0,p_str.length()-4);
	  pos = (front+".pos");
	  chr = (front+".chr");
	  geno = (front+".geno");
	}else{
	  pos = (p_str+".pos");
	  chr = (p_str+".chr");
	  geno = (p_str+".geno");

	}}else{
	  pos = (p_str+".pos");
	  chr = (p_str+".chr");
	  geno = (p_str+".geno");

     	}
      posTaken=1;
      genoTaken=1;
      chrTaken=1;
    }
    else if(strcmp(argv[argPos],"-plink")==0){
      std::string p_str =string( argv[argPos+1]);
      if(p_str.length()>4){
	std::string ext = p_str.substr(p_str.length()-4,p_str.length());
	if (!ext.compare(".bed")||!ext.compare(".bim")||!ext.compare(".fam")){
	  std::string front = p_str.substr(0,p_str.length()-4);
	  plink_bim = (front+".bim");
	  plink_fam = (front+".fam");
	  plink_bed = (front+".bed");
	}else{
	  plink_bim = (p_str+".bim");
	  plink_fam = (p_str+".fam");
	  plink_bed = (p_str+".bed");	
	}}else{
	plink_bim = (p_str+".bim");
	plink_fam = (p_str+".fam");
	plink_bed = (p_str+".bed");	
     	}

	doPlink+=3;
    }
    else if(strcmp(argv[argPos],"-plink-fam")==0){
      plink_fam = argv[argPos+1];
      doPlink++;
    }

    else if(strcmp(argv[argPos],"-choose")==0){
      myChoose = atoi(argv[argPos+1]); 
     
    } else if(strcmp(argv[argPos],"-d")==0){
      re = getTestIndividuals(argv[argPos+1]);
      individualsToTestList = 1;
    }
    else{
      printf("\nArgument unknown will exit: %s \n",argv[argPos]);
      info();
      return 0;
    }
    argPos+=2;
  }
  if(pars->print_results){
    printf("\t-> relateHMM version: %.3f\t",VERSION);
    time_t t; //added time output in 0.985
    time(&t);
    std::cout<<ctime(&t);
  }
  //   Changed in 0.83
  
  if(!doPlink){
    if(pars->print_results){
      printf("\t-> We will use text file as input...\n");
      printf("\t\tgeno:\t\"%s\"\n\t\tchromo:\t\"%s\"\n\t\tpos:\t\"%s\"\n",geno.c_str(),chr.c_str(),pos.c_str());

    }
    //We are given solely -geno -chromo -pos args
    if(!genoTaken){
      printf("\t-> You must supply atleast a genotype file -g , will exit\n");
      exit(0);
    }
    
    if(!posTaken){
      printf("\t-> You must supply atleast a positionfile file -p , will exit\n");
      exit(0); // remember
    }
    pars->position = getPos(pos.c_str());
    
    if(!chrTaken){
      if(pars->print_results!=0)
	printf("\t-> no chromosomefile given, will assume all loci is on same chromosome\n");
      
      pars->chr = allocIntArray(pars->position->x);
      for(int i=0;i<pars->chr->x;i++)
	pars->chr->array[i] = -1;
      
    }
    else
      pars->chr = getChr(chr.c_str());

    if(pars->print_results!=0)
      flush_print("Reading genotype data, depending on size this might take som time...\n");
    pars->data = getData(geno.c_str());
  }else{
    if(genoTaken||chrTaken||posTaken){
      printf("\t-> You must use either plink binary file or -g,-p file as input, will exit\n");
      return 0;
    }
    //We are using plink binary files
    printf("\t-> Will assume these are the plink files:\n\t\tbed: %s\n\t\tbim: %s\n\t\tfam: %s\n",plink_bed.c_str(),plink_bim.c_str(),plink_fam.c_str());
    if(doPlink!=3){
      printf("If you want to run relateHMM with plink binary files\n");
      printf("you need to supply all .bim .bed .fam f\n");
      return 0;
    }
    int numInds = numberOfLines(plink_fam.c_str())-1;//the number of individuals is just the number of lines
    plinkKeep = doBimFile(pars,plink_bim.c_str()," \t");

    if(pars->print_results){
      printf("\t-> Plink file contains %d autosomale SNP's\n",plinkKeep->numTrue);
      flush_print("Will now read entire plink binary file into memory...");
    }
    iMatrix *tmp = bed_to_iMatrix(plink_bed.c_str(),numInds,plinkKeep->x);
    if(pars->print_results)
      flush_print("\r\t-> Will now read entire plink binary file into memory... done\n");
    if(tmp->x==plinkKeep->numTrue)
      
      pars->data = tmp;
    else{
      pars->data = extractOK(plinkKeep,tmp);
      killMatrix(tmp);
    }
    killArray(plinkKeep);
    mysort(pars,pars->print_results); //0.92
  }
  //    printf("\r                                                                               ");
    if(pars->print_results!=0)
      printf("\r");
  
  if(pars->data->y != pars->chr->x || pars->position->x != pars->data->y){
    printf("Dimension of data input doesn't have compatible dimensions, program will exit\n");
    printf("Dimension of genodata:=(%d,%d), positions:=%d, chromosomes:=%d\n",pars->data->x,pars->data->y,pars->position->x,pars->chr->x);
    return 0;
  }

  pars->kFilename = kFilename;
  pars->postFilename = postFilename;
  getOptions(option.c_str(),pars);
  post_processing_of_option_pars(pars,myChoose);
  //0.96
  if(pars->doAllPairs||individualsToTestList)//always run fast version
    myChoose=0;

  if(chrTaken&&!doPlink)
    mysort(pars,pars->print_results); //0.92

  relateHMM *myObject = new relateHMM();
  if(pars->print_results)
    printf("\t-> Dimension of genodata:=(%d,%d), positions:=%d, chromosomes:=%d \n",pars->data->x,pars->data->y,pars->position->x,pars->chr->x);
  //doErrorCheck
  int i1 = pars->pair->array[0];
  int i2 = pars->pair->array[1];
  int ulim=pars->data->x-1;
  if(i1<0 || i1>ulim || i2<0 || i2>ulim){
    printf("You want to run analysis on pair=(%d,%d), genotype dataset only contains (%d...%d), will exit\n",i1,i2,0,pars->data->x);
    exit(0);
  }
  if(!individualsToTestList){
    if(!pars->doAllPairs){
      //    pars->print_results = 2;
      fullResults *res;
      if(myChoose)
	res = myObject->single_pair(pars);
      else{
	myObject->init_globals(pars);
	res = myObject->sgl(pars);
	myObject->dinit_globals(pars->back);
      }   
      if(pars->kFilename.compare("")) { 
	if(fexists(pars->kFilename))
	  printf("\t-> Filename: %s exists will append new results to file\n",pars->kFilename.c_str());
	
	if(pars->fixK==1 && pars->fixK_val==0){
	  //printf("The special k case...");
	  wrt_dArray_to_file(pars->kFilename,res->kResult,1);
	  //this is the second elment,will append
	}
	else{
	  //printf("The general k  case...");
	  wrt_dArray_to_file(pars->kFilename,res->kResult,1);
	  //this is the second elment,will append
	  wrt_dArray_to_file(pars->kFilename,res->kResult,0);
	  //this is the first elment,will append
	}
      }
      if(pars->postFilename.compare("")){
	if(fexists(pars->postFilename))
	  printf("\t-> Filename: %s exists will append new results to file\n",pars->postFilename.c_str());
	if(pars->fixK==1 && pars->fixK_val==0){
	  //printf("The special post case...");
	  wrt_dMatrix_to_file(pars->postFilename,res->post,2,res);
	  //this is the third column,will append
	}
	else{
	  //printf("The general post case...");
	  wrt_dMatrix_to_file(pars->postFilename,res->post,2,res);
	  //this is the third column,will append
	  wrt_dMatrix_to_file(pars->postFilename,res->post,1,res);
	  //this is the sec column,will append
	}
      }
      killFullResults(res);
    }
    else{
      // pars->print_results=0;
      //int returnValue = myObject.all_pairs(pars);
      if(myChoose)
	myObject->all_pairs(pars);
      else{
	myObject->init_globals(pars);
	myObject->fast_all_pairs(pars);
	myObject->dinit_globals(pars->back);
      }
      
    }
  }else{
    myObject->init_globals(pars);
    myObject->joblist_fast_all_pairs(pars,re);
    myObject->dinit_globals(pars->back);
  }
  delete myObject;
  if(pars->print_results){
    time_t t; //added time output in 0.985
    time(&t);
    std::cout<<"\t-> "<<ctime(&t);
  }
  killFunctionPars(pars); 
} 
