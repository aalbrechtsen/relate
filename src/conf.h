/*
  generic constants used at compile time

 */
#define VERSION 0.992
#define SYS_MAX 10000000.0 //apparantly theres a problem with DBL_MAX
#define POSTPRECISION 4 //set the number of decimals printet in posterior file
#define PLINK_POS_SCALING 1000000 //used for scaling the position from a plink format

//used for comparing floating points differences smaller will be assumed to be the same
#define PRECISION_COMPARISON 0.00001 
#define MAX_ELEMS_PER_LINE 5000000


//These are the default filenames used when dumping the stripped results
#define CHROMOSOMES_FILE "stripped.chr"
#define GENOTYPE_FILE "stripped.geno"
#define POSITIONS_FILE "stripped.pos"
#define KEEPLIST_FILE "keep.index"


#define PREINFO "\t--------------------------------\n\t\tThis is just a small help box, for Relate version 0.992! \n"
#define INPUTFILES "\tInputfiles:\n\t\tInputfiles can be tab-seperated whitespace or any of {,.:;}. \n\tFor one-dimensional inputdata \
newlines are also accepted as seperators.\n\tThe genotypes files should have columns as SNP's and rows as individuals.\n\tThe -d file \
should contain the number of individuals, with zero or one\n\tindicating if the individual should be included.\n"
#define ALLPAIRS "\tAllpairs:\n\t\tTo run allpairs you must change options file. If the -d argument\n\tis given at runtime the program\
 will run all combinations of pairs. \n"
#define OUTPUTFILES "\tOutputfiles:\n\tIf some files already exists these will NOT be overwritten,\n\tbut new ones will be created by appending an integer to the filesnames.\n"
#define CONTACT "\tContact:\n\t\tThis program is copyrighted and are under the gpl license.\n\
\tThis mean that there are no warranty, and the authors take no\n\tresponsibility. But the authors \
will of cause take any kind honour\n\tor credit.\n\tBugs or improvement can be reported to the following emails\
\n\t\tmain author: albrecht@binf.ku.dk\n\t\tprogrammer : thorfinn@binf.ku.dk\
\n\t\twebsite: http://staff.pubhealth.ku.dk/~ande/software\n"

#define REQ "\t(You need to supply atleast a genotype file, and a positionfile,\n\tor all three plink files.)\n\
\tSorting will be performed on all input data, so the ordering of SNP's\n \
\tshould not matter as long as all input files share the same order.\n\
\t(For more elaborate info consult the manual or check the website)\n"

#define FLATFILES "\tSYNTAX: flatfiles\n\t./relateHMM -o OPTION -g GENOTYPES \
 -p POSITIONS\n\t-c CHROMO -d TESTINDIVIDUALS -v VERBOSE\n\t-post POSTOUTPUT -k KOUTPUT\n"

#define PLINKFILES "\tSYNTAX: plink binary files\n\t./relateHMM -o OPTION -plink-bed test.bed \
 -plink-fam test.fam \n\t-plink-bim test.bim -d TESTINDIVIDUALS -v VERBOSE\n\t-post POSTOUTPUT -k KOUTPUT\n"
