/*Input: A set of WMs (PhyloGibbs/TRANSFAC format). Optionally a prior weight for each WM */
/*       A set of multiple alignments. Maybe use input where each sequence group has own phylogenetic tree.*/
/*       A phylogenetic tree                       */
/*       Background models fr each of the species. Think of it as lengh 1 WM*/
/*       Column probabilities according to Nacho's model (unknown WM). Prior for unknown WMs of different length*/
/*       A parameter file that says how/what calculation to do */
/*       also says which species is reference***/
/*                                                            */
/*Calculations:                                                */
/*             Use given priors and WMs to calculate posteriors of WM occurrences **/
/*             Fit the priors to maximize likelihood of the data                  **/
/*             Use EM to optimize the WMs (local maximum likelihood of data)      **/
/*             Calculate Z scores for particular segments of the input (enhancer finding) **/
/*             These tasks can be combined                                                **/

/*Modules needed:                                                            */
/*               Read multiple alignments                                    */
/*               Read phylogenetic tree                                      */
/*               Read background models for each species                     */
/*               Set background models for internal nodes                    */
/*               A module for calculating the score of a alignment column given the tree and WM (or bg) **/
/*               A module for calculating Z value of a segment               */
/*               Calculating posteriors for each position each WM            */
/*               Optimizing priors                                           */
/*               EM on WMs (expected counts from posteriors)                 */


/***Data structures needed****/
/**Input multiple alignments with respect to reference species.**/
/***   -set of reference species sequences***/
/***   -for each position in reference columns of bases in other species***/
/***    -each column of bases has species in same order***/
/***    -when we  ***/
/***   -for each 'tile length' group of species that are ungapped wrt reference at that tile-length ***/

/***Calculating sum of parses:             ***************/
/***  For getting Z(n) we need to sum over all WMs we can put*****/
/***  For each WM w of length l_w we need to score the whole segment of the alignment that corresponds to bases
      n-l_w+1 through n in the reference sequence           ****/
/***  First we need to know which sequences are gaplessly aligned with the reference ***/
/**** The others are background per default *******/
/**** Then we score each of the sequences individually for the WM -> all positive go into foreground**/
/**** Then we go through all columns in the alignment segment  ****/
/**** For each column some species will be scored as bg and some as fg. We calculate score for each column***/
/***  The column scores are calculated recursively: input fg/bg pattern + WM column + bases + global tree ***/
/**** NOTE: we only need to score the columns where the reference species has a base ****/
/****       all other columns the bases in the other species will always be scored bg (independent of configuration)****/
/***        and so they will cancel out of the score                                                                ****/
/*** Gaps are scored as 'missing data', i.e. unknown bases... must just be factor 1/4?? ******************************/

/***Fill in matrix with every position/every sequence/every matrix, probability of putting side for WM there****/
/***After FW and BW Zs are calculated for a sequence we can use this to fill in posteriors *********************/
/***The set of posteriors would immediately allows us to update priors ****************************************/
/***Updating WM entries is more difficult.. i.e. see Saurabh's paper ****************************************/
/***However if we ignore contribution from other species then it might be more easy **************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>
#include "evomodel.h"
#include <stdlib.h>

void dummy();

/***** structure definitions *****/
const char digits[11]="0123456789";
const int treelength = 8192;

typedef struct {
  char name[1024];
  int len;
  /***put in a single array of length 4*len ****/
  double **mat;
  double **lmat;
  double **matrev;
  double **lmatrev;
  double *colsum;
  double prior;
  double priorinit;
} wm;


typedef struct blubber {
  char name[1024]; /***is empty for internal nodes****/
  int index; /**index in array of species **/
  int level; /***level in the tree***/
  struct blubber *parent;
  double q;
} treenode;

typedef struct baba {
  double q;/***distance to parent***/
  int parent;/***parent index in full treestructure****/
  int index;/******index in full tree structure****/
} reducedtreenode;


typedef struct {
  int numspec;
  int len;
  char **specnames;
  int *specindex;
  int type; /** 0 = A region, 1 = B left gene, 2 = B right gene, 3 = C region ***/
  /***make single array of length numspec * len ****/
  char **columns; /*** bases A,C,G,T,- mapped to 0,1,2,3,4   [pos][species] ****/
} alignment;


/*********** global variables **************/
int numnodes;
int numspec;
int numwms;
int numgroups;
int maxlen;
int otherwmlen;
double bg[4];
const double selec_cutoff = 1.1; /****To get rid of sites that are all N*****/
reducedtreenode **reducedtree;



/**************** FUNCTION DECLARATIONS *************/
int read_param_file(char *filename,char *refspecies,char *ttreestring,int *em_prior,int *em_wm,int *do_segment,int *seglen,int *bgorder,double *otherwmprior,double *bgprior,int *restrictparses,char *otherwmfile,char *sitefile,char *priorfile,char *loglikfile,char *column_statfile,char *cons_statfile,char *mutcountfile,char *segfile,double *minposterior,int *printsiteals,char *otherwm_proffile,int *proflen, int *do_ep, double *priordiff, double *wmdiff, int *winlen, int *steplen, int *do_evo, int *markov_order,char *epfilename, int *useproximities,char *mode, double *minposteriorWM, int *UFEprint, char *mybgfile);
void mychomp(char *s);
int parsetree(char *ttreestring,treenode **nodelist,char *refspecies, int prox);
int read_bg_filenames(char *filename,char **bgfiles,int numnodes,treenode **nodelist);
int countwms(char *filename);
int readwms(char *filename,wm *wms,int *num_with_prior,double *wm_prior_sum);
int myisnumber(char s);
void printwm(wm *thiswm);
double getinf(double *n);
char consensus(double na, double nc, double ng, double nt);
char* twodigitstr(int m);
int readsequences(alignment *alignments,char *filename,int numgroups);
int parsenames(alignment *alignments,treenode **nodelist);
void setleafprob(double *prob,char letter);
int getseqsizes(char *filename,int *numgroups,int *groupsizes,int *seqlengths);
double scorewm(alignment *al,int pos,int len, double **mat,int spec);
double scorewmforward(alignment *al,int pos,int len, double **mat,int spec);
double scorebg(char *column,treenode **nodelist);
void forwardbackward(alignment *al,double **prob,double **probrev,treenode **nodelist,wm *wms,int restrictparses,double *otherwm,double *Z,int numwms);
void printsites(alignment *al,double **prob,double **probrev,wm *wms,FILE *sitefile,double minposterior,int printsiteals, int printUFE);
void getpostsum(alignment *al,double **prob,double **probrev,double *postsums,wm *wms, int numwms);
void print_priors(double *postsums,FILE *priorfile,double totpostsum,wm *wms);
double update_priors(double *postsums,double *totpostsum,wm *wms, int numwms);
double scorerefwm(alignment *al,int pos,int len,double **mat);
double scorerefwmforward(alignment *al,int pos,int len,double **mat);
int read_otherwmfile(char *otherwmfile,double *otherwm,treenode **nodelist);
int getindex(char *column,int *selec);
int gettotindex(char *column);
void printals(alignment *alignments);
void indextocol(int i, char *tmpcol);
int indextototindex(int i);
void colstats(alignment *alignments,treenode **nodelist,double *otherwm,char *column_statfile);
int get_region_type(char *name);
void cons_inf(alignment *alignments,treenode **nodelist,char *cons_statfile, char *mutcountfile);
int basenum(char let);
void print_tree(treenode **nodelist);
double alrefscore(alignment *al,int pos,int len,double **mat,char *tmpseq);
void printsegments(alignment *al,double *Z,FILE *segfp,int seglen);
char revbase(char base);
void otherwm_profs(alignment *alignments,treenode **nodelist,double *otherwm,char *otherwm_proffile,int proflen);
void get_selection_vector(treenode **nodelist,int *selec);
reducedtreenode **init_reducedtrees(treenode **nodelist);
double ratiocolumn(char *column,double *w,reducedtreenode *thistree);
reducedtreenode *tree_from_selec(int *selec);

// for the wm refinement
void update_wms(alignment *al,double **prob,double **probrev,wm *wms,FILE *sitefile,double minposterior,int printsiteals, double ***postcounts);
double diffwms(double ***postcounts,wm *wms, int numwms);

// for the enhancer finding
int extractaln(alignment *input, alignment *output, int start, int len) ;
int n2ng(alignment *input, int start, int step);

// for the nth order markov model
int markovorder = 0;
double *bgfwd, *bgbwd;

void getbgmodelfromseqs(alignment *al);
void read_cond_bg_probs(double *bgmok_tmp, int mo, char *filename);

inline int getprevseqs(alignment *al, int spec, int pos, char *seqret);
inline int getprevseqs_fwd(alignment *al, int spec, int pos, char *seqret);

inline int seq2index(char *seq);
inline void index2seq(int index, char* tmp);
inline char int2base(int number);
inline int base2int(char base);

double scorerefwm_nth(alignment *al,int pos,int len,double **mat);
double scorerefwmforward_nth(alignment *al,int pos,int len,double **mat);
double scorewm_nth(alignment *al,int pos,int len, double **mat,int spec);
double scorewmforward_nth(alignment *al,int pos,int len, double **mat,int spec);
double scorebg_nth(char *column,treenode **nodelist);
double ratiocolumn_nth(alignment *al,int pos,char *column,double *w,reducedtreenode *thistree, int refpos);
double ratiocolumnforward_nth(alignment *al,int pos, char *column,double *w,reducedtreenode *thistree, int refpos);
void forwardbackward_nth(alignment *al,double **prob,double **probrev,treenode **nodelist,wm *wms,int restrictparses,double *otherwm,double *Z,int numwms,double *Y);

int gaplesswithref(alignment *al,int pos,int len,int spec);
int gaplesswithref_fwd(alignment *al,int pos,int len,int spec);

inline int doesMOfit(alignment *al, int spec, int pos, int wmlen);
inline int doesMOfit_fwd(alignment *al, int spec, int pos, int wmlen);

inline void flip(char *in, char* out, int len);

void node2kids(treenode **nodelist);
void node2kidsreduced(reducedtreenode *thistree, int **kids);

inline double myround(double nbr, int digits);

void check_sequence(alignment *al);


double EXPratiocolumn_nth(alignment *al,int pos, char *column,double *w,reducedtreenode *thistree, double *fg, double *bg);
double EXPratiocolumnforward_nth(alignment *al,int pos, char *column,double *w,reducedtreenode *thistree, double *fg, double *bg);
void EXPforwardbackward(alignment *al,double **prob,double **probrev,treenode **nodelist,wm *wms,int restrictparses,double *otherwm,double *Z,int numwms, double *Y);

/********** MAIN STARTS*************/

/***TODO: ****/
/**make option to only test parses with sites fg>bg ***/
/***distance dependent prior for each WM***/
/***  for this we need one prob dist going downstream from a gene and one going upstream ***/
/***  Then at each point we can put a site down for the left gene or for the right gene with separate probs*****/
/***  Somehow we will also need some space smoothing. *****/
/**WM independently for each branch in the tree***/
/**   For this we need a temporary WM to hold the WMs for each branch**/
/**   This is specified by the WMs at the leafs and taking averages going up..maybe weighted by q for each branch***/
/***updating the WM counts***/


/* nth order background model 

- Functions involved:

  1)  scorerefwm
  2)  scorerefwmforward
  3)  scorewm(alignment *al,int pos,int len, double **mat,int spec);
  4)  scorewmforward(alignment *al,int pos,int len, double **mat,int spec);
  5)  scorebg(char *column,treenode **nodelist);
  6)  forwardbackward

*/



int main(int argc, char *argv[]) {

  
  FILE *wmfile,*priorfp,*sitefp,*segfp,*loglikfp;
  treenode **nodelist;
  char refspecies[1024];
  char epfilename[1024],otherwmfile[1024],sitefile[1024],priorfile[1024],loglikfile[1024],column_statfile[1024],cons_statfile[1024],mutcountfile[1024],segfile[1024],otherwm_proffile[1024];
  char treestring[8192];
  int useproximities = 0;
  char **bgfiles;
  int spec,em_wm, em_prior, seglen, do_segment,i,j,k,seqnum,pos,bgorder,num_with_prior,*selec,*selecrev,prior_done,wm_done,first;
  alignment *alignments;
  wm *wms, *wmsbg;
  double *stored_priors;
  double wm_prior_sum = 0;
  double av_wm_prior,totpostsum,priordif,minposterior,minposteriorWM,wmdif;
  double otherwmprior,bgprior,**prob,**probrev,*postsums,***postcounts, ***orgpostcounts, *otherwm,*Z,*Y,*Zbg,*Ybg;
  int *seqlengths,*groupsizes, restrictparses, printsiteals,proflen;

  int do_ep = 0; // 1: enhancer prediction, 0: wmref or bsp
  double diffprior = 0.000075;
  double diffwm = 0.0075;
  int winlen = 400;
  int steplen = 50;
  int do_evo = 0; // creating UFE model: yes or no

  char mode[1024];
  int UFEprint = 0;
  char mybgfile[1024];

  /**set some default values***/
  strcpy(epfilename,"");
  strcpy(otherwmfile,"");
  strcpy(sitefile,"");
  strcpy(priorfile,"");
  strcpy(loglikfile,"");
  strcpy(column_statfile,"");
  strcpy(cons_statfile,"");
  strcpy(mutcountfile,"");
  strcpy(otherwm_proffile,"");
  strcpy(mode,"bsp");
  strcpy(mybgfile,"");

  postcounts = NULL;
  orgpostcounts = NULL;
  em_wm = 0;
  em_prior = 0;
  seglen = 0;
  do_segment = 0;
  otherwmprior = 0;
  bgprior = 0;
  proflen = 1;
  /**set default background model***/
  bg[0] = 0.25;
  bg[1] = 0.25;
  bg[2] = 0.25;
  bg[3] = 0.25;
  
  if(argc != 4) {
    fprintf(stderr,"  Usage: motevo input_sequences param_file wm_file\n");
    return (0);
  }
  
  /** read parameter file ****/
  if(read_param_file(argv[2],refspecies,treestring,&em_prior,&em_wm,&do_segment,&seglen,&bgorder,&otherwmprior,&bgprior,&restrictparses,otherwmfile,sitefile,priorfile,loglikfile,column_statfile,cons_statfile,mutcountfile,segfile,&minposterior,&printsiteals,otherwm_proffile,&proflen,&do_ep,&diffprior,&diffwm,&winlen,&steplen,&do_evo,&markovorder,epfilename,&useproximities,mode,&minposteriorWM,&UFEprint,mybgfile))
    {
      fprintf(stderr,"Error reading parameter file\n");
      return 1;
    }

  // do some other things with the parameters
  do_ep=0;
  em_wm=0;

  if( strcmp(mode,"TFBS")== 0 )     {  do_ep=0; em_wm=0; }
  else if( strcmp(mode,"ENH")==0 ) {  do_ep=1; em_wm=0; }
  else if( strcmp(mode,"WMREF")== 0 ) {  do_ep=0; em_wm=1; }
  else {
    fprintf(stderr,"Parameter 'Mode' in the parameter file must be set to TFBS, ENH, or WMREF!\nWill do TFBS by default.");
    strcpy(mode,"TFBS");
    do_ep=0; em_wm=0;
  }

  printf("Running in %s mode\n",mode);

  /***now we parse the tree****/
  /***First count nodes and species****/
  pos = 0;
  numnodes = 1;
  while(treestring[pos] != ';' && treestring[pos] != '\0')
    {
      if(treestring[pos] == ':')
	{
	  ++numnodes;
	}
      ++pos;
    }
  nodelist = (treenode **) calloc(numnodes,sizeof(treenode *));
  bgfiles = (char **) calloc(numnodes,sizeof(char *));
  for(i=0;i<numnodes;++i)
    {
      bgfiles[i] = (char *) calloc(1024,sizeof(char));
      bgfiles[i][0] = '\0';
    }
  if(parsetree(treestring,nodelist,refspecies,useproximities))
    {
      fprintf(stderr,"Error in parsing treesting\n");
      return 1;
    }
  /***print out tree for checking***/
  print_tree(nodelist);

  /******data structure for reduced trees when gaps occur********/
  reducedtree = (reducedtreenode **) init_reducedtrees(nodelist);

  /***now read the bg file names ****/
  if(read_bg_filenames(argv[2],bgfiles,numnodes,nodelist))
    {
      fprintf(stderr,"Error in reading background files\n");
      return 1;
    }
  
  /***read the motifs from the motif file****/
  wmfile = (FILE *) fopen(argv[3],"r");

  /**first count the number of WMs***/
  numwms = countwms(argv[3]);
  wms = (wm *) calloc(numwms+2,sizeof(wm));
  wmsbg = (wm *) calloc(2,sizeof(wm)); // 0: bg, 1: otherwm
  // to store the priors
  stored_priors = (double *) calloc(numwms+2,sizeof(double));

  num_with_prior = 0;
  if(numwms > 0)
    {
      if(readwms(argv[3],wms,&num_with_prior,&wm_prior_sum))
	{
	  fprintf(stderr,"Error reading WMs\n");
	  return 1;
	}
    }

  /***print the WMs for checking****/
  for(i=0;i<numwms;++i)
    {
      printwm(wms+i);
    }

  /****read values for nacho model****/
  if(strcmp(otherwmfile,"") != 0)
    {
      if(otherwmprior == 0)
	{
	  fprintf(stderr,"Warning: otherwmfile specified but prior is zero\n");
	}
      /**hold probability ratios for all 5^numspec possible letter combinations***/
      otherwm = (double *) calloc((int) pow(5.0,(double) numspec),sizeof(double));
      read_otherwmfile(otherwmfile,otherwm,nodelist);
    }
  /**check that file was specified when positive prior**/
  if(otherwmprior > 0 && strcmp(otherwmfile,"") == 0)
    {
      fprintf(stderr,"Error: otherwm.prior is %lf but no otherwmfile was specified\n",otherwmprior);
      return 1;
    }
 
  /***wm priors. Add a pseudocount which is the average of the existing ones****/
  if(num_with_prior > 0)
    {
      av_wm_prior = wm_prior_sum/num_with_prior;
    }
  else
    {
      av_wm_prior = 1.0;
    }
  wm_prior_sum = 0;
  for(i=0;i<numwms;++i)
    {
      if(wms[i].prior == 0)
	wms[i].prior = av_wm_prior;
      wm_prior_sum += wms[i].prior;
    }
  wm_prior_sum += otherwmprior;
  for(i=0;i<numwms;++i)
    {
      wms[i].prior *= ((1.0-bgprior)/wm_prior_sum);
      wms[i].priorinit =  wms[i].prior;
      printf("final prior %s is %lf %lf\n",wms[i].name,wms[i].prior,wms[i].priorinit);
    }
  otherwmprior *= (1-bgprior)/wm_prior_sum;

  /**set otherwm and bg as normal wms***/
  wms[numwms].prior = bgprior;
  wms[numwms].priorinit = bgprior;
  wms[numwms+1].prior = otherwmprior;
  wms[numwms+1].priorinit = otherwmprior;
  wms[numwms+1].len = otherwmlen;
  wms[numwms].len = 1;
  strcpy(wms[numwms+1].name,"UFEwm");
  strcpy(wms[numwms].name,"background");

  printf("UFE wm prior %lf %lf\n",otherwmprior,wms[numwms+1].priorinit);
  printf("bg prior %lf %lf\n",bgprior,wms[numwms].priorinit);

  /***otherwm matrix***/
  wms[numwms+1].mat = (double **) calloc(otherwmlen,sizeof(double *));
  for(pos=0;pos<otherwmlen;++pos)
    {
      wms[numwms+1].mat[pos] = (double *) calloc(4,sizeof(double));
      for(i = 0;i<4;++i)
	{
	  wms[numwms+1].mat[pos][i] = 1.0;
	}
    }


  /**print some info on parameters read from parameter file***/
  printf("ref species is %s\n",refspecies);
  printf("treestring is %s\n",treestring);
  printf("em_prior is %d\n",em_prior);
  printf("em_wm is %d\n",em_wm);
  printf("do_segment is %d\n",do_segment);
  printf("seglen is %d\n",seglen);
  /*printf("bgorder %d\n",bgorder);
  for(i=0;i<numspec;++i)
    {
      printf("bgfile %d name %s file %s\n",i,nodelist[i]->name,bgfiles[i]);
    }
  */
  for(i=0;i<4;++i)
    printf("bg[%d] = %lf\n",i,bg[i]);
 
  /***now read in the input file***/
  numgroups = 0;
  seqnum = 0;
  /**assuming at most 50000 groups and 250000 sequences***/
  groupsizes = (int *) calloc(150000,sizeof(int));
  seqlengths = (int *) calloc(750000,sizeof(int));
  if(getseqsizes(argv[1],&numgroups,groupsizes,seqlengths))
    {
      fprintf(stderr,"Error reading input file\n");
      return 1;
    }

  /***now get the memory for the alignment structure**/
  alignments = (alignment *) calloc(numgroups,sizeof(alignment));
  for(i=0;i<numgroups;++i)
    {
      alignments[i].numspec = groupsizes[i];
      alignments[i].len = seqlengths[seqnum];
      ++seqnum;
      /**check all aligmments have same length*/
      for(j=1;j<groupsizes[i];++j)
	{
	  if(seqlengths[seqnum] != alignments[i].len)
	    {
	      fprintf(stderr,"Error: sequence number %d in group %d has length %d which is different from length %d of first sequence in group\n",(j+1),i+1,seqlengths[seqnum],alignments[i].len);
	      return 1;
	    }
	  ++seqnum;
	}
      /***memory for the species/sequence names and indices****/
      alignments[i].specnames = (char **) calloc(alignments[i].numspec,sizeof(char *));
      alignments[i].specindex = (int *) calloc(alignments[i].numspec,sizeof(int));
      for(j=0;j<alignments[i].numspec;++j)
	{
	  alignments[i].specnames[j] = (char *) calloc(1024,sizeof(char));
	}
      /**memory for the columns***/
      alignments[i].columns = (char **) calloc(alignments[i].len,sizeof(char *));
      for(j=0;j<alignments[i].len;++j)
	{
	  alignments[i].columns[j] = (char *) calloc(alignments[i].numspec,sizeof(char));
	}
    }
  /***no longer need seqlengths and groupsizes***/
  free(seqlengths);
  free(groupsizes);

  printf("allocated memory\n");

  /***now read in the sequences****/
  if(readsequences(alignments,argv[1],numgroups))
    {
      fprintf(stderr,"Error in reading sequences from file\n");
      return 1;
    }
  printf("read sequences\n");

  /***parse the sequences and get everything in the format we want it in***/
  if(parsenames(alignments,nodelist))
    {
      fprintf(stderr,"Error in parsing sequences\n");
      return 1;
    }
  
  printf("parsed sequences\n");

  /***print out all alignments as test****/
  /*printals(alignments);*/

  /* init nth order background model */
  printf("mybgfile is set to %s\n",mybgfile);

  if( strcmp(mybgfile,"")==0 ) {
    printf("Calculating background frequencies from input sequence\n");
    getbgmodelfromseqs(alignments);
  }
  else {
    printf("Reading background frequencies from input file\n");
    read_cond_bg_probs(bgfwd, markovorder, mybgfile);
  }


  /***statistics on observed frequency different columns and expected frequencies under bg and fg****/
  
  if((strcmp(otherwmfile,"") != 0) && (strcmp(column_statfile,"") != 0))
    colstats(alignments,nodelist,otherwm,column_statfile);
  else if((strcmp(otherwmfile,"") == 0) && (strcmp(column_statfile,"") != 0))
    {
      fprintf(stderr,"Error: you asked for a column statistics file but didn't give an UFE-model file \n");
      return 1;
    }
  /***conservation statistics per region and mutation counts, if filenames specified****/
  if((strcmp(cons_statfile,"") != 0) || (strcmp(mutcountfile,"") != 0))
    cons_inf(alignments,nodelist,cons_statfile,mutcountfile);
  
  /***otherwm model profiles if filename is specified****/
  if((strcmp(otherwmfile,"") != 0) && (strcmp(otherwm_proffile,"") != 0))
    otherwm_profs(alignments,nodelist,otherwm,otherwm_proffile,proflen);
  else if((strcmp(otherwmfile,"") == 0) && (strcmp(otherwm_proffile,"") != 0))
    {
      fprintf(stderr,"Error: you asked for an 'conservation' profile statistics file but didn't give an UFE-model file\n");
      return 1;
    }


  /***memory for forward and backward partition sums**/
  maxlen = 0;
  for(i=0;i<numgroups;++i)
    {
      if(alignments[i].len > maxlen)
	maxlen = alignments[i].len;
    }
  /***memory for matrix of posterior probabilities site to occur***/
  prob = (double **) calloc(maxlen,sizeof(double *));
  /***memory on for sites on opposite strand*****/
  probrev = (double **) calloc(maxlen,sizeof(double *));
  Z = (double *) calloc(maxlen,sizeof(double));
  Y = (double *) calloc(maxlen,sizeof(double));
  Zbg = (double *) calloc(maxlen,sizeof(double));
  Ybg = (double *) calloc(maxlen,sizeof(double));
  for(pos=0;pos<maxlen;++pos)
    {
      /***for each position number of wms + bg + nacho's model***/
      prob[pos] = (double *) calloc(numwms+2,sizeof(double));
      probrev[pos] = (double *) calloc(numwms+2,sizeof(double));
    }
  /***memory for 'selection tree'***/
  selec = (int *) calloc(numnodes,sizeof(int));
  selecrev = (int *) calloc(numnodes,sizeof(int));

  /***if em_prior = 1 we need space for the sum of posteriors for all factors***/
  postsums = (double *) calloc(numwms+2,sizeof(double));
  if(em_prior)
    prior_done = 0;
  else
    prior_done = 1;
  
  /***if em_wm = 1 we need also space for counts of bases at each position of each wm***/
  if(em_wm)  {
    postcounts = (double ***) calloc(numwms,sizeof(double **));
    orgpostcounts = (double ***) calloc(numwms,sizeof(double **));
    for(j=0;j<numwms;++j)
      {
	postcounts[j] = (double **) calloc(wms[j].len,sizeof(double *));
	orgpostcounts[j] = (double **) calloc(wms[j].len,sizeof(double *));
	for(pos=0;pos<wms[j].len;++pos) {
	  postcounts[j][pos] = (double *) calloc(4,sizeof(double));
	  orgpostcounts[j][pos] = (double *) calloc(4,sizeof(double));
	}
      }
    wm_done = 0;
  }
  else
    {
      wm_done = 1;
    }

  
  int win = winlen;
  int step = steplen;

  /*** enhancer predictions ***/
  /* Remarks:
    - If the score log [P(D|wms,o,bg)/P(D|o,bg) ] is < 0, then decrease diffprior, because this should never happen!
    - To do:
       - make sure it also works on single sequence
  */
  if(do_ep == 1) {

    int max_len = 0;
    for(i=0;i<numgroups;++i) 
      if(alignments[i].len > max_len)
	max_len = alignments[i].len;
    
    printf("Enhancer predictions:\n");
    printf("  window length is %i\n",win);
    printf("  step length is %i\n",step);
    printf("  em prior difference is %f\n",diffprior);

    alignment *alntmp = NULL;
    alntmp = (alignment *) malloc(sizeof(alignment));
    alntmp->len = max_len;
    
    alntmp->numspec = numspec;
    alntmp->specnames = (char **) calloc(alntmp->numspec,sizeof(char *));
    alntmp->specindex = (int *) calloc(alntmp->numspec,sizeof(int));
    
    for(j=0;j<alntmp->numspec;++j)
      alntmp->specnames[j] = (char *) calloc(1024,sizeof(char));
  
    alntmp->columns = (char **) calloc(alntmp->len,sizeof(char *));
    for(j=0;j<alntmp->len;++j)
      alntmp->columns[j] = (char *) calloc(alntmp->numspec,sizeof(char));

    char *epfile = epfilename;
    FILE *epfp = fopen(epfile, "w");

    //header for file
    fprintf(epfp,"identifier\twin_start\twinstop\twinscore");
    for(j=0;j<(2+numwms);++j) 
      fprintf(epfp,"\t%s",wms[j].name);
    fprintf(epfp,"\n");

      
    for(i=0;i<numgroups;++i) {

      int pos;
      for(pos=0;pos<alignments[i].len;pos+=step) {

	if( extractaln(alignments+i,alntmp,pos,win) == -1) { break; } 
	

	fprintf(epfp,"%s\t%i\t%i\t",alignments[i].specnames[0],pos,pos+win-1);


	/*** P(D|otherwm,bg,wms) ***/      
	// reset priors
	int j;
	for(j=0;j<(2+numwms);++j) {
	  wms[j].prior = wms[j].priorinit;      
	  //printf("%i %i %f\n",j,numwms,wms[j].prior);
	}
	
	int counter=0;
	if(em_prior == 0) 
	  prior_done = 1;
	else
	  prior_done = 0;

	while(!prior_done) {
	  
	  // prior updating
	  for(j=0;j<(numwms+2);++j)
	    postsums[j] = 0.0;

	 
	  if(markovorder==0)
	    forwardbackward(alntmp,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms);
	  else
	    forwardbackward_nth(alntmp,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms,Y);


	  /***sums of posteriors ***/	  
	  getpostsum(alntmp,prob,probrev,postsums,wms,numwms);
	  priordif = update_priors(postsums,&totpostsum,wms,numwms);
	  //printf("\tP(D|o,wm,bg): priordiff = %f\n",priordif);
	  counter++;

	  if(priordif < diffprior)
	    prior_done = 1;

	}
	// save priors
	for(j=0;j<(2+numwms);++j) {
	  stored_priors[j] = wms[j].prior;
	}
	

	if(markovorder==0)
	  forwardbackward(alntmp,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms);
	else 
	  forwardbackward_nth(alntmp,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms,Y);
	//forwardbackward(alignments+i,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms);

	for(j=0;j<(2+numwms);++j) {
	  wms[j].prior = wms[j].priorinit;      
	  //	  printf("\t prior after: P(D|bg,o,wms)\t%i %i %f\n",j,numwms,wms[j].prior);
	}

	/*** P(D|otherwm,bg) ***/
	// create wms for only otherwm and bg
	wmsbg[0] = wms[numwms]; 
	wmsbg[0].prior = wmsbg[0].priorinit;
	
	if(wms[numwms+1].priorinit > 0) {
	  wmsbg[1] = wms[numwms+1]; 
	  wmsbg[1].prior = 1-wmsbg[0].prior;
	}
	
	// prior updating
	for(j=0;j<(numwms+2);++j)
	  postsums[j] = 0.0;
	
	for(j=0;j<(numwms+2);++j) {
	  //	  printf("%i %f (%f)\n",j,wms[j].prior,wms[j].priorinit);
	}

	//printf("prior for otherwm P(D|bg,o): %f\n\n",wms[numwms+1].priorinit);

	counter = 0;
	if(em_prior == 0) 
	  prior_done = 1;
	else
	  prior_done = 0;

	while(!prior_done && wms[numwms+1].priorinit > 0) {
	  
	  if(markovorder==0)
	    forwardbackward(alntmp,prob,probrev,nodelist,wmsbg,restrictparses,otherwm,Zbg,0);
	  else
	    forwardbackward_nth(alntmp,prob,probrev,nodelist,wmsbg,restrictparses,otherwm,Zbg,0,Ybg);
	  //forwardbackward(alignments+i,prob,probrev,nodelist,wmsbg,restrictparses,otherwm,Zbg,0);

	  /***sums of posteriors ***/
	  getpostsum(alntmp,prob,probrev,postsums,wmsbg,0);
	  priordif = update_priors(postsums,&totpostsum,wmsbg,0);

	  //printf("\tP(D|o,bg): priordiff[%i] = %f\n",counter,priordif);
	  counter++;
	  if(priordif < diffprior)
	    prior_done = 1;
	}

	if(wms[numwms+1].priorinit > 0) {
	  if(markovorder==0)
	    forwardbackward(alntmp,prob,probrev,nodelist,wmsbg,restrictparses,otherwm,Zbg,0);
	  else 
	    forwardbackward_nth(alntmp,prob,probrev,nodelist,wmsbg,restrictparses,otherwm,Zbg,0,Ybg);
	  //forwardbackward(alignments+i,prob,probrev,nodelist,wmsbg,restrictparses,otherwm,Zbg,0);
	}

 	double fgZ = 0;
	double bgZ = 0;
	for(j=0;j<alntmp->len;++j) {
	  fgZ += log(Z[j]);
	  //printf("%i: Z[%i]=%f  fgZ=%f\n",j,j,Z[j],fgZ);
	  if(wms[numwms+1].priorinit > 0) 
	    bgZ += log(Zbg[j]);
	}

	char *id = alignments[i].specnames[0];
		
	if(wms[numwms+1].priorinit > 0) {
	  //fprintf(epfp,"%f\t%f\t%f\n",fgZ,bgZ,fgZ-bgZ);
	  fprintf(epfp,"%f",fgZ-bgZ);

	  for(j=0;j<(2+numwms);++j) {
	    fprintf(epfp,"\t%f",stored_priors[j]);
	  }
	  fprintf(epfp,"\n");
	  
	}
	else {
	  fprintf(epfp,"%f",fgZ);
	  
	  for(j=0;j<(2+numwms);++j) {
	    fprintf(epfp,"\t%f",stored_priors[j]);
	  }
	  fprintf(epfp,"\n");

	}

      }
    }

    fclose(epfp);
    printf("Enhancer predictions done ...\n");
    return 1;
  }
  
  else { /*** wm refinement or binding site predictions ***/

    // To do:
    // - wm refinement: constraint min posterior
    // - wmdiff and priordiff - specify in file
    
    /***the updating of prior and/or WM****/
    while(!wm_done || !prior_done)
      {
	first = 0;
	/***reset to zero the counts of the number of sites for each factor***/
	for(j=0;j<(numwms+2);++j)
	  postsums[j] = 0.0;
	
	/***clear counts of number of bases in each column of each WM****/
	if(em_wm) {
	  for(j=0;j<numwms;++j) {
	    for(pos=0;pos<wms[j].len;++pos) {
	      for(k=0;k<4;++k) {
		postcounts[j][pos][k] = 0.0;
		orgpostcounts[j][pos][k] = 0.0;
	      }
	    }
	  }
	}
	
	/***forward backward algorithm ****/
	for(i=0;i<numgroups;++i) {

	  
	  /**calculate partition sums with all WMs and put posteriors into prob and probrev matrices***/
	  //forwardbackward(alignments+i,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms);
	  if(markovorder==0)
	    forwardbackward(alignments+i,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms);
	  else 
	    forwardbackward_nth(alignments+i,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms,Y);
	  
	  /***sums of posteriors ***/
	  if(em_prior)
	    getpostsum(alignments+i,prob,probrev,postsums,wms,numwms);
	  
	  /***sums of posteriors for wms ***/
	  if(em_wm) {
	    update_wms(alignments+i,prob,probrev,wms,sitefp,minposteriorWM,printsiteals,postcounts);
	    
	    const double PC = 0.5;

	    for(j=0;j<numwms;++j) {
	      for(pos=0;pos<wms[j].len;++pos) {
		for(k=0;k<4;++k)
		  orgpostcounts[j][pos][k] = postcounts[j][pos][k];
	      }
	    }
	  }

	} 
	/*
	for(j=0;j<numwms;++j) {
	  for(pos=0;pos<wms[j].len;++pos) {
	    for(k=0;k<4;++k)
	      printf("numwms=%i, wm=%i, pos=%i, value[%i]=%f\n",numwms,j,pos,k,postcounts[j][pos][k]);
	  }
	}
	*/
	
	/***update the wms given the current posteriors and calculate difference***/
	if(em_wm) {
	  
	  double wm_diff = diffwms(postcounts,wms,numwms);
	  
	  if(wm_diff < diffwm) {
	    wm_done = 1;
	    
	    
	    /* print wms and priors */
	    /* write weight matrix to file  */
	    char pwmfile[1024];
	    strcpy(pwmfile,"wms.updated");
	    FILE *pwmfp = fopen(pwmfile, "w");
	    for (j=0; j<numwms; ++j) {
	      
	      
	      fprintf(pwmfp, "//\n");
	      fprintf(pwmfp, "NA %s\n", wms[j].name);
	      fprintf(pwmfp, "P0\tA\tC\tG\tT\tcons\tinf\n");
	      for (pos=0; pos<wms[j].len; ++pos) {
		if ((pos>=0)&&(pos<9))
		  fprintf(pwmfp, "0%d\t", pos+1);
		if ((pos>=9))
		  fprintf(pwmfp, "%d\t", pos+1);

		double cs[4];
		for (k=0; k<4; ++k)
		  cs[k] = orgpostcounts[j][pos][k];

		// consensus and information
		char cons = consensus(cs[0],cs[1],cs[2],cs[3]);
		double inform = getinf(cs);
		
		for (k=0; k<4; ++k)
		  fprintf(pwmfp, "%.3f\t", orgpostcounts[j][pos][k]);
		fprintf(pwmfp,"%c\t%.3f\n",cons,inform);
		
	      }
	      fprintf(pwmfp, "//\n");
	    }
	    fclose(pwmfp);
	    
	    
	  }
	  else { 
	    
	    // update wms
	    int i;
	    for (j=0; j<(numwms); ++j)
	      for (pos=0; pos<wms[j].len; ++pos)
		for (i=0; i<4; ++i) {
		  wms[j].mat[pos][i]=postcounts[j][pos][i];
		  wms[j].lmat[pos][i] = log(wms[j].mat[pos][i]);
		}
	    
	    for (i=0; i<(numwms); ++i)
	      for (j=0; j<wms[i].len; ++j)
		for (k=0; k<4; ++k) {
		  wms[i].matrev[j][k] = wms[i].mat[wms[i].len-1-j][3-k];
		  wms[i].lmatrev[j][k]  = wms[i].lmat[wms[i].len-1-j][3-k];
		}
	  }
	}
	
	/***update the priors given the current posteriors and calculate difference***/
	if(em_prior)
	  {
	    priordif = update_priors(postsums,&totpostsum,wms,numwms);
	    printf("total prior dif %lf\n",priordif);
	    if(priordif < diffprior)
	      {
		prior_done = 1;
	      }
	  }
      }
    /***finallly print out the results***/
    if(strcmp(sitefile,"") != 0)
      {
	sitefp = fopen(sitefile,"w");
      }
    if(strcmp(loglikfile,"") != 0)
      {
	loglikfp = fopen(loglikfile,"w");
      }
    if(strcmp(segfile,"") != 0)
      {
	segfp = fopen(segfile,"w");
      }
    for(j=0;j<(numwms+2);++j)
      postsums[j] = 0.0;
    
    for(i=0;i<numgroups;++i)
      {

	check_sequence(alignments+i);

	if(markovorder==0)
	  forwardbackward(alignments+i,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms);
	else 
	  forwardbackward_nth(alignments+i,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms,Y);
	//forwardbackward(alignments+i,prob,probrev,nodelist,wms,restrictparses,otherwm,Z,numwms);


	if(strcmp(loglikfile,"") != 0)
	  {
	    /**calculate log-likelihood of this sequence***/
	    double tF = 0.0;
	    for(int b=0;b<alignments[i].len;b++) {
	      tF += log(Z[b]);
	    }
	    /***now print line in the file***/
	    fprintf(loglikfp,"%s %lf\n",alignments[i].specnames[0],tF);
	  }
	 
	if(strcmp(sitefile,"") != 0)
	  {
	    printsites(alignments+i,prob,probrev,wms,sitefp,minposterior,printsiteals,UFEprint);
	  }
	if(strcmp(segfile,"") != 0)
	  {
	    printsegments(alignments+i,Z,segfp,seglen);
	  }
	getpostsum(alignments+i,prob,probrev,postsums,wms,numwms);
	
      }
    if(strcmp(sitefile,"") != 0)
      {
	fclose(sitefp);
      }
    if(strcmp(segfile,"") != 0)
      {
	fclose(segfp);
      }
    if(strcmp(loglikfile,"") != 0)
      {
	fclose(loglikfp);
      }
    
    /**print out priors if requested****/
    if(strcmp(priorfile,"") != 0)
      {
	priordif = update_priors(postsums,&totpostsum,wms,numwms);
	priorfp = fopen(priorfile,"w");
	print_priors(postsums,priorfp,totpostsum,wms);
	fclose(priorfp);
      }
    
  }
  return (0);
}

/***************** FUNCTIONS **********************************/


void print_tree(treenode **nodelist)
{
  int i;
  for(i=0;i<numnodes;++i)
    {
      printf("node %d name %s qval %g parent %s\n",i,nodelist[i]->name,nodelist[i]->q,(nodelist[i]->parent)->name);
    }
  return;
}
  


int basenum(char let)
{
  if(let == 'A' || let == 'a')
    return 0;
  else if(let == 'C' || let == 'c')
    return 1;
  else if(let == 'G' || let == 'g')
    return 2;
  else if(let == 'T' || let == 't')
    return 3;
  /***return a random number***/
  else 
    return ((int) (4.0*rand()/(RAND_MAX+1.0)));
}

void cons_inf(alignment *alignments,treenode **nodelist,char *cons_statfile, char *mutcountfile)
{
  FILE *consfile;
  char reflet, let;
  int *conscount,i,spec,reflen,pos;
  int mutcount[4][4];
  int totmut = 0;
  conscount = (int *) calloc(numspec,sizeof(int));
  consfile = NULL;
  if(strcmp(cons_statfile,"") != 0)
    {
      consfile = (FILE *) fopen(cons_statfile,"w");
    }

  for(let=0;let<4;++let)
    {
      for(reflet=0;reflet<4;++reflet)
	{
	  mutcount[reflet][let] = 0;
	}
    }

  for(i=0;i<numgroups;++i)
    {
      for(spec=0;spec<numspec;++spec)
	conscount[spec] = 0;
      reflen = 0;
      for(pos=0;pos<alignments[i].len;++pos)
	{
	  reflet = alignments[i].columns[pos][0];
	  if(reflet != '-')
	    {
	      ++reflen;
	      for(spec=1;spec<numspec;++spec)
		{
		  let = alignments[i].columns[pos][spec];
		  if(reflet == let)
		    ++conscount[spec];
		  else if(let != '-')
		    {
		      ++totmut;
		      ++mutcount[basenum(reflet)][basenum(let)];
		    }
		}
	    }
	}
      if(consfile != NULL)
	{
	  fprintf(consfile,"%d %d ",i,reflen);
	  for(spec=1;spec<numspec;++spec)
	    fprintf(consfile,"%g ",((double) conscount[spec])/((double) reflen));
	  fprintf(consfile,"\n");
	}
    }
  if(consfile != NULL)
    fclose(consfile);


  consfile = NULL;
  if(strcmp(mutcountfile,"") != 0)
    {
      consfile = fopen("mutcounts","w");
      for(reflet=0;reflet<4;++reflet)
	{
	  for(let=0;let<4;++let)
	    {
	      fprintf(consfile,"%d %d %g\n",reflet,let,((double) mutcount[reflet][let])/((double) totmut));
	    }
	}
      fclose(consfile);
    }
  return;
}


void otherwm_profs(alignment *alignments,treenode **nodelist,double *otherwm,char *otherwm_proffile,int proflen)
{
  FILE *proffile;
  int i,index,totindex,refpos,pos,spec,count,reflen, shift;
  double *scores,scoresum,avscore;
  char reflet;
  static int *selec = NULL;
  if(selec == NULL)
    {
      selec = (int *) calloc(numspec,sizeof(int));
    }
  for(spec=0;spec<numspec;++spec)
    {
      selec[spec] = 1;
    }

  proffile = NULL;
  if(strcmp(otherwm_proffile,"") != 0)
    {
      proffile = (FILE *) fopen(otherwm_proffile,"w");
    }
  
  
  for(i=0;i<numgroups;++i)
    {
      refpos = 0;
      scores = (double *) calloc(alignments[i].len,sizeof(double));
      for(pos=0;pos<alignments[i].len;++pos)
	{
	  reflet = alignments[i].columns[pos][0];
	  if(reflet != '-')
	    {
	      index = getindex(alignments[i].columns[pos],selec);
	      scores[refpos] = otherwm[index];
	      ++refpos;
	     }
	}
      reflen = refpos;
     
      scoresum = 0;
      /**get first proflen positions***/
      for(pos=0;pos<proflen;++pos)
	scoresum += scores[pos];
      for(pos=proflen;pos<reflen;++pos)
	{
	  avscore = (scoresum/((double) proflen));
	  fprintf(proffile,"%d %d %lf %d %s\n",pos-proflen,pos-1,avscore,i,alignments[i].specnames[0]);
	  scoresum += scores[pos];
	  scoresum -= scores[pos-proflen];
	}
      free(scores);
    }
  fclose(proffile);
  return;

}




void colstats(alignment *alignments,treenode **nodelist,double *otherwm,char *column_statfile)
{
  FILE *colfile;
  char tmpcol[1024];
  int i, spec,pos,index,totindex,totcols;
  int *colcountA,*colcountB,*colcountC;
  int *totcountA,*totcountB,*totcountC;
  int numcols = (int) pow(2.0,(double) numspec);
  totcountA = (int *) calloc(numcols,sizeof(int));
  totcountB = (int *) calloc(numcols,sizeof(int));
  totcountC = (int *) calloc(numcols,sizeof(int));
  numcols = (int) pow(5.0,(double) numspec);
  colcountA = (int *) calloc(numcols,sizeof(int));
  colcountB = (int *) calloc(numcols,sizeof(int));
  colcountC = (int *) calloc(numcols,sizeof(int));
  double score;
  /****selection on species set to 1 by default***/
  static int *selec = NULL;
  if(selec == NULL)
    {
      selec = (int *) calloc(numspec,sizeof(int));
    }
  for(spec=0;spec<numspec;++spec)
    {
      selec[spec] = 1;
    }


  for(i=0;i<numgroups;++i)
    {
      for(pos=0;pos<alignments[i].len;++pos)
	{
	  index = getindex(alignments[i].columns[pos],selec);
	  totindex = gettotindex(alignments[i].columns[pos]);
	  /*printf("column %s index %d totindex %d\n",alignments[i].columns[pos],index,totindex);*/
	  if(alignments[i].type == 0)
	    {
	      ++colcountA[index];
	      ++totcountA[totindex];
	    }
	  else if(alignments[i].type == 3)
	    {
	      ++colcountC[index];
	      ++totcountC[totindex];
	    }
	  else
	    {
	      ++colcountB[index];
	      ++totcountB[totindex];
	    }
	}
    }
  colfile = (FILE *) fopen(column_statfile,"w");
  for(i=0;i<numcols;++i)
    {
      indextocol(i,tmpcol);
      totindex = indextototindex(i);
      score = scorebg(tmpcol,nodelist);
      if(totcountA[totindex] == 0)
	totcountA[totindex] = 1;
      if(totcountB[totindex] == 0)
	totcountB[totindex] = 1;
      if(totcountC[totindex] == 0)
	totcountC[totindex] = 1;
      fprintf(colfile,"%s countA %g countB %g countC %g bg %g fg %g sumA %d sumB %d sumC %d\n",tmpcol,((double) colcountA[i])/((double) totcountA[totindex]),((double) colcountB[i])/((double) totcountB[totindex]),((double) colcountC[i])/((double) totcountC[totindex]),score,score*otherwm[i],totcountA[totindex],totcountB[totindex],totcountC[totindex]);
    }

  free(totcountA);
  free(totcountB);
  free(totcountC);
  free(colcountA);
  free(colcountB);
  free(colcountC);
  return;
}

void indextocol(int i, char *tmpcol)
{
  int spec,thisnum;
  for(spec=0;spec<numspec;++spec)
    {
      thisnum = (i % 5);
      if(thisnum == 0)
	tmpcol[spec] = 'A';
      else if(thisnum == 1)
	tmpcol[spec] = 'C';
      else if(thisnum == 2)
	tmpcol[spec] = 'G';
      else if(thisnum == 3)
	tmpcol[spec] = 'T';
      else
	tmpcol[spec] = '-';
      i -= thisnum;
      i /= 5;
    }
  tmpcol[numspec] = '\0';
  return;
}

int indextototindex(int i)
{
  int spec,thisnum;
  int totindex = 0;
  int offset = 1;
  for(spec=0;spec<numspec;++spec)
    {
      thisnum = (i % 5);
      if(thisnum == 4)
	totindex += offset;
      i -= thisnum;
      i /= 5;
      offset *= 2;
    }
  return totindex;
}

void printals(alignment *alignments)
{
  int i, spec,pos;
  for(i=0;i<numgroups;++i)
    {
      printf("group %d\n",i);
      for(spec=0;spec<numspec;++spec)
        {
          printf("%s\n",alignments[i].specnames[spec]);
          for(pos=0;pos<alignments[i].len;++pos)
            {
              printf("%c",alignments[i].columns[pos][spec]);
            }
          printf("\n");
        }
    }
}

int read_otherwmfile(char *otherwmfile,double *otherwm,treenode **nodelist)
{
  int i,found,spec,index,offset;
  FILE *infile;
  char s[1024],tmpname[1024],orname[1024];
  int *nodemap;
  double dum1,dum2,ratio,score;
   
  
  nodemap = (int *) calloc(numspec,sizeof(int));
  mychomp(otherwmfile);
  infile = (FILE *) fopen(otherwmfile,"r");
  if(infile == NULL)
    {
      fprintf(stderr,"Otherwmfile file %s does not exist or cannot be opened.\nExiting.\n",otherwmfile);
      return 1;
    }
  /**reading the file****/
  spec = 0;
  while (fgets(s,1024,infile)) /**get line from the file***/
    {
      mychomp(s);
      if(s[0] == '>') /***species name*****/
	{
	  i = 1;
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(tmpname,s+i);
	  /*printf("read species name %s\n",tmpname);*/
	  found = 0;
	  i=0;
	  while(found == 0 && i < numspec)
	    {
	      if(strcmp(tmpname,nodelist[i]->name) == 0)
		{
		  found = 1;
		  nodemap[spec] = i;
		  ++spec;
		}
	      ++i;
	    }
	  if(found == 0)
	    {
	      fprintf(stderr,"Error: cannot find species in tree with name %s from otherwmfile %s\n",tmpname,otherwmfile);
	      return 1;
	    }
	  /*printf("name %s maps to node %d with name %s\n",tmpname,nodemap[spec-1],nodelist[nodemap[spec-1]]->name);*/
	}
      /**all species need to have been read by now***/
      else
	{
	  if(spec != numspec)
	    {
	      fprintf(stderr,"Error: there were %d species in file %s but %d in tree\n",spec,otherwmfile,numspec);
	      return 1;
	    }
	  if(sscanf(s,"%s %lf %lf %lf",tmpname,&dum1,&dum2,&ratio) == 0)
	    {
	      fprintf(stderr,"Error: Couldn't scan letters and three doubles from %s\n",s);
	      return 1;
	    }
	  if(strlen(tmpname) != numspec)
	    {
	      fprintf(stderr,"Error: number of letters in %s doesn't match number of species %d\n",tmpname,numspec);
	      return 1;
	    }
	  index = 0;
	  /*we also need to make an ordered set of bases for calc. the bg model***/
	  orname[numspec] = '\0';
	  for(i=0;i<numspec;++i)
	    {
	      offset = (int) pow(5.0,(double) nodemap[i]);
	      orname[nodemap[i]] = tmpname[i];
	      if(tmpname[i] == 'A')
		index += int(0.0);
  	      else if(tmpname[i] == 'C')
		index += int(1.0*offset);
	      else if(tmpname[i] == 'G')
		index += int(2.0*offset);
	      else if(tmpname[i] == 'T')
		index += int(3.0*offset);
	      else if(tmpname[i] == '-')
		index += int(4.0*offset);
	      else
		{
		  fprintf(stderr,"Error: cannot parse character %d position %d of line %s\n",tmpname[i],i,s);
		  return 1;
		}
	    }
	  /****only do columns where there is a letter in the reference****/
	  if(index % 5 != 4){
	    score = scorebg(orname,nodelist);
	    if((dum1-score)*(dum1-score)/((dum1+score)*(dum1+score)) > 0.01)
	      {
		fprintf(stderr,"DEVIATION %s %g mybg %g\n",tmpname,dum1,score);
	      }
	    otherwm[index] = ratio;
	  }
	  else
	    {
	      otherwm[index] = 0;
	    }
	  
	}
    }
  fclose(infile);
  return 0;
}

void print_priors(double *postsums,FILE *priorfile,double totpostsum,wm *wms)
{
  int j;
  double lensum = 0.0;
  for(j=0;j<(2+numwms);++j)
    {
      lensum += postsums[j]*((double) wms[j].len);
    }

  fprintf(priorfile,"WM_name final_prior nr_of_sites density\n");

  for(j=0;j<(2+numwms);++j)
    {
      fprintf(priorfile,"%s %g %g %g\n",wms[j].name,postsums[j],(totpostsum*postsums[j]),postsums[j]*((double) wms[j].len)/lensum);
    }
  return;
}

double update_priors(double *postsums,double *totpostsum,wm *wms,int numwms)
{
  int j;
  *totpostsum = 0.0;
  double priordif = 0.0;
  for(j=0;j<(2+numwms);++j)
    {
      if(postsums[j] >= 0.01)
	*totpostsum += postsums[j];
      else
	postsums[j] = 0.0;
    }
  /***normalize****/
  for(j=0;j<(2+numwms);++j)
    {
      /**if less than 0.01 sites in whole data.. set to zero****/
      postsums[j] /= *totpostsum;
      if(postsums[j] > 0)
	priordif += (postsums[j]-wms[j].prior)*(postsums[j]-wms[j].prior)/((postsums[j]+wms[j].prior)*(postsums[j]+wms[j].prior));
      wms[j].prior = postsums[j];
    }
  priordif /= ((double) (numwms+2));

  return (priordif);
}

void getpostsum(alignment *al,double **prob,double **probrev,double *postsums,wm *wms,int numwms)
{
  int pos,j,refpos;

  refpos = 0;
  for(pos=0;pos<al->len;++pos)
    {
      if(al->columns[pos][0] != '-')
	{
	  ++refpos;
	  for(j=0;j<numwms;++j)
	    {
	      /**check wm can fit and has nonzero prior****/
	      if(wms[j].len<=refpos && wms[j].prior > 0)
		{
		  postsums[j] += prob[pos][j]+probrev[pos][j];
		}
	    }
	  /***count posterior for the background model***/
	  postsums[numwms] += prob[pos][numwms];
	  /***otherwm model****/
	  postsums[numwms+1] += prob[pos][numwms+1];
	}
    }
  return;
}

void printsegments(alignment *al,double *Z,FILE *segfp,int seglen)
{
  int pos;
  double sum, mean;;
  /***if too short to fit a segment just get the whole thing****/
  if(al->len < seglen)
    {
      sum += 0;
      for(pos=0;pos<al->len;++pos)
	{
	  sum += log(Z[pos]);
	}
      mean = sum/((double) al->len);
      fprintf(segfp,"0 %d %g %g %s\n",al->len,sum,mean,al->specnames[0]);
    }
  else
    {
      /***get the first segment****/
      for(pos=0;pos<seglen;++pos)
	{
	  sum += log(Z[pos]);
        }
      mean = sum/((double) seglen);
      fprintf(segfp,"0 %d %g %g %s\n",seglen,sum,mean,al->specnames[0]);
      pos = seglen;
      while(pos< al->len)
	{
	  sum += log(Z[pos]);
	  sum -= log(Z[pos-seglen]);
	  fprintf(segfp,"%d %d %g %g %s\n",pos,seglen,sum,mean,al->specnames[0]);
	  ++pos;
	}
    }
  return;
}

void printsites(alignment *al,double **prob,double **probrev,wm *wms,FILE *sitefile,double minposterior,int printsiteals, int printUFE)
{
  int pos, refpos,j,basecount,thispos,spec,startpos;
  static double *selec;
  double totsum = 0.0;
  double score,alscore;
  char tmpseq[100];
  selec = (double *) calloc(numspec,sizeof(double));

  refpos = 0;
  for(pos=0;pos<al->len;++pos)
    {
      if(al->columns[pos][0] != '-')
	{
	  ++refpos;
	  /**wm models***/
	  for(j=0;j<(numwms+2);++j)
	    {
	      if( printUFE==0 && j==(numwms+1) ) { continue; }
	      /**check wm can fit and has nonzero prior****/
	      if(wms[j].len<=refpos && wms[j].prior > 0 && j != numwms)
		{
		  totsum += prob[pos][j]*((double) wms[j].len);
		  totsum += probrev[pos][j]*((double) wms[j].len);
		  if(prob[pos][j] >= minposterior)
		    {
		      fprintf(sitefile,"%d-%d + %lf %s %s\n",refpos-wms[j].len+1,refpos,prob[pos][j],wms[j].name,al->specnames[0]);
		      /****print the bases at the positions in site if wanted****/
		      if(printsiteals != 0)
			{
			  startpos = pos;
			  basecount = 1;
			  while(basecount < wms[j].len)
			    {
			      --startpos;
			      if(al->columns[startpos][0] != '-')
				{
				  ++basecount;
				}
			    }
			  selec[0] = scorerefwm(al,pos,wms[j].len,wms[j].mat);
			  for(spec=1;spec<numspec;++spec)
			    {
			      selec[spec] = scorewm(al,pos,wms[j].len,wms[j].mat,spec);
			    }
			  /***now print out the alignment****/
			  for(spec=0;spec<numspec;++spec)
			    {
			      if(selec[spec] > selec_cutoff || spec == 0)
				{
				  for(thispos=startpos;thispos<=pos;++thispos)
				    {
				      if(al->columns[thispos][0] != '-')
					fprintf(sitefile,"%c",al->columns[thispos][spec]);
				    }
				  fprintf(sitefile," %g %s\n",log(selec[spec]),al->specnames[spec]);
				}
			    }
			}
		    }
		  if(probrev[pos][j] >= minposterior)
		    {
		      fprintf(sitefile,"%d-%d - %lf %s %s\n",refpos-wms[j].len+1,refpos,probrev[pos][j],wms[j].name,al->specnames[0]);
		      /****print the bases at the positions in site if wanted****/
		      if(printsiteals != 0)
			{
			  startpos = pos;
			  basecount = 1;
			  while(basecount < wms[j].len)
			    {
			      --startpos;
			      if(al->columns[startpos][0] != '-')
				{
				  ++basecount;
				}
			    }
			  selec[0] =  scorerefwm(al,pos,wms[j].len,wms[j].matrev);
			  for(spec=1;spec<numspec;++spec)
			    {
			      selec[spec] = scorewm(al,pos,wms[j].len,wms[j].matrev,spec);
			    }
			  /***now print out the alignment****/
			  for(spec=0;spec<numspec;++spec)
			    {
			      if(selec[spec] > selec_cutoff || spec == 0)
				{
				  for(thispos=pos;thispos>=startpos;--thispos)
				    {
				      if(al->columns[thispos][0] != '-')
					fprintf(sitefile,"%c",revbase(al->columns[thispos][spec]));
				    }
				  fprintf(sitefile," %g %s\n",log(selec[spec]),al->specnames[spec]);
				}
			    }
			}
		    }
		}
	    }
	  totsum += prob[pos][numwms];
	}
    }
  /*printf("positions %d totsum %lf\n",refpos,totsum);*/
  free(selec);
  return;
}


/**get selection vector (1=node under selection,0=node is under bg) for all nodes in the tree**/
void get_selection_vector(treenode **nodelist,int *selec)
{
  int spec,thispat;
  
  /**lookup table from selection pattern of the species to whole tree using max-parsimony***/
  /***only fill in this table once*****/
  static int *lookup_sel_mapping = NULL;
  if(lookup_sel_mapping == NULL)
    {
      /***Note we assume that the number of species is relatively small compared to the number of columns in the alignments * numwms****/
      /***Otherwise this is a relative costly to calculate****/
      
      int pat,min,tot_pat;
      int *tmp_selec, *gainlosswm,*gainlossbg;
      gainlosswm = (int *) calloc(numnodes,sizeof(int));
      gainlossbg = (int *) calloc(numnodes,sizeof(int));
      tmp_selec = (int *) calloc(numnodes,sizeof(int));

      /*** tot_pat = 2^(numspec-1) possible patterns, because reference is always 1*********/
      tot_pat = (1 << (numspec-1));
      lookup_sel_mapping = (int *) calloc(tot_pat,sizeof(int));

      for(pat=0;pat<tot_pat;++pat)
	{
	  /***initialize output integer to zero******/
	  lookup_sel_mapping[pat] = 0;
	  /***Reference is under selection per definition***/
	  tmp_selec[0] = 1;
	  gainlosswm[0] = 0;
	  gainlossbg[0] = 1;
	  /**read selection of other species from the pattern and set gains and losses***/
	  for(spec=1;spec<numspec;++spec)
	    {
	      /***get the bit shifted spec-1 to the right****/
	      tmp_selec[spec] = ((pat>>(spec-1)) & 1);
	      if(tmp_selec[spec])
		{
		  gainlosswm[spec] = 0;
		  gainlossbg[spec] = 1;
		}
	      else
		{
		  gainlosswm[spec] = 1;
		  gainlossbg[spec] = 0;
		}
	    }

	  /**Calculate gain and losses for the rest of the tree***/
	  for(spec=0;spec<numnodes;++spec)
	    {
	      /***only for non-root nodes****/
	      if(nodelist[spec]->parent != 0)
		{
		  if(gainlosswm[spec] < (gainlossbg[spec]+1))
		    {
		      min = gainlosswm[spec];
		    }
		  else
		    {
		      min = gainlossbg[spec]+1;
		    }
		  gainlosswm[nodelist[spec]->parent->index] += min;
		  
		  if(gainlossbg[spec] < (gainlosswm[spec]+1))
		    {
		      min = gainlossbg[spec];
		    }
		  else
		    {
		      min = gainlosswm[spec]+1;
		    }
		  gainlossbg[nodelist[spec]->parent->index] += min;
		}
	    }
	  /***Now go from the root down****/
	  spec = numnodes-1;
	  
	  /***If WM set selection output pattern, which has only the internal nodes***/
	  if(gainlosswm[spec] <= gainlossbg[spec])
	    lookup_sel_mapping[pat] += (1 << (spec-numspec));
	  
	  /***now do the rest of the internal nodes****/
	  for(spec=numnodes-2;spec>=numspec;--spec)
	    {
	      /***parent has the site****/
	      if(selec[nodelist[spec]->parent->index])
		{
		  if(gainlosswm[spec] <= (gainlossbg[spec]+1))
		    lookup_sel_mapping[pat] += (1 << (spec-numspec));
		}
	      /***parent does not have the site***/
	      else
		{
		  if(gainlossbg[spec] > (gainlosswm[spec]+1))
		   lookup_sel_mapping[pat] += (1 << (spec-numspec));
		}
	    }
	}
      /***free the arrays we used to fill in the table****/
      free(tmp_selec);
      free(gainlosswm);
      free(gainlossbg);
    }

  /**Set selection pattern for the other leafs of the tree***/
  thispat = 0;
  for(spec=1;spec<numspec;++spec)
    {
      if(selec[spec])
	thispat += (1 << (spec-1));
    }
  /***output selection vector for the internal nodes***/
  thispat = lookup_sel_mapping[thispat];
  for(spec=numspec;spec<numnodes;++spec)
    selec[spec] = ((thispat >> (spec-numspec)) & 1);

  return;
}

void forwardbackward(alignment *al,double **prob,double **probrev,treenode **nodelist,wm *wms,int restrictparses,double *otherwm,double *Z,int numwms)
{
  int pos,j,refpos,relpos,abspos,spec,index,allselec;
  double fac, facrev,score,postsum;
  static double *F = NULL;
  if(F == NULL)
    F = (double *) calloc(maxlen,sizeof(double));
  static double *R = NULL;
  if(R == NULL)
    R = (double *) calloc(maxlen+1,sizeof(double));
 
  double bgprior = wms[numwms].prior;
  double otherwmprior = wms[numwms+1].prior;

  static int *selec = NULL;
  if(selec == NULL)
    selec = (int *) calloc(numnodes,sizeof(int));
  reducedtreenode *thistree;

  /***run over all positions***/
  /**at each position run over all wms***/
  /***check which species are background vs. foreground by scoring each species for the wm**/
  /***calculate score for wm by taking this 'selection' array****/
  refpos = 0;
  for(pos=0;pos<al->len;++pos)
    {
      /***forward***/
      if(al->columns[pos][0] != '-')
	{
	  ++refpos;
	  /***score and record column under background only****/
	  /***we store 1/probability****/
	  F[pos] = bgprior;
       
	  /**wm models***/
	  for(j=0;j<numwms;++j)
	    {
	      /**check wm can fit and has nonzero prior****/
	      if(wms[j].len<=refpos && wms[j].prior > 0)
		{
		  /***if restrictparses is set we only put a site when the reference scores better than background***/
		  /***otherwise we consider all parses for the reference****/
		  if(restrictparses)
		    {

		      score = scorerefwm(al,pos,wms[j].len,wms[j].mat);

		      if(score > selec_cutoff)
			selec[0] = 1;
		      else
			selec[0] = 0;
		    }
		  else
		    {
		      selec[0] = 1;
		    }
		  
		  /***only score when reference has selection***/
		  if(selec[0]> 0)
		    {
		      /**get selection for the other leafs****/
		      for(spec=1;spec<numspec;++spec)
			{
			  if(scorewm(al,pos,wms[j].len,wms[j].mat,spec) > selec_cutoff)
			    selec[spec] = 1;
			  else
			    selec[spec] = 0;
			}
		      /**Get structure of the reduced tree that runs only on these species***/
		      thistree = tree_from_selec(selec);

		      /***get full score going backwards***/
		      relpos = wms[j].len-1;
		      
		      fac = 0.5*wms[j].prior * ratiocolumn(al->columns[pos],wms[j].mat[relpos],thistree);
		      --relpos;
		      abspos = pos-1;
		      while(relpos >= 0)
			{
			  /***only positions with base in reference****/
			  if(al->columns[abspos][0] != '-')
			    {
			      /**check that it is not negative***/
			      if(abspos < 0)
				{
				  fprintf(stderr,"HELP NEGATIVE abspos\n");
				}
			      fac *= F[abspos]*ratiocolumn(al->columns[abspos],wms[j].mat[relpos],thistree);
			      --relpos;
			    }
			  --abspos;
			}
		      F[pos] += fac;
		      prob[pos][j] = fac;
		    }
		  else
		    prob[pos][j] = 0;
		    
		  /***selection for the reverse-strand*******/
		  if(restrictparses)
		    {
		      score = scorerefwm(al,pos,wms[j].len,wms[j].matrev);
		      if(score > selec_cutoff)
			selec[0] = 1;
		      else
			selec[0] = 0;
		    }
		  else
		      selec[0] = 1;

		  if(selec[0] > 0)
		    {
		      /**get selection for the other leafs****/
		      for(spec=1;spec<numspec;++spec)
			{
			  if(scorewm(al,pos,wms[j].len,wms[j].matrev,spec) > selec_cutoff)
			    selec[spec] = 1;
			  else
			    selec[spec] = 0;
			}
		      /**Get structure of the reduced tree that runs only on these species***/
		      thistree = tree_from_selec(selec);

		      /***get full score going backwards***/
		      relpos = wms[j].len-1;
		      fac = 0.5*wms[j].prior * ratiocolumn(al->columns[pos],wms[j].matrev[relpos],thistree);
		      --relpos;
		      abspos = pos-1;
		      while(relpos >= 0)
			{
			  /***only positions with base in reference****/
			  if(al->columns[abspos][0] != '-')
			    {
			      /**check that it is not negative***/
			      if(abspos < 0)
				{
				  fprintf(stderr,"HELP NEGATIVE abspos\n");
				}
			      fac *= F[abspos]*ratiocolumn(al->columns[abspos],wms[j].matrev[relpos],thistree);
			      --relpos;
			    }
			  --abspos;
			}
		      F[pos] += fac;
		      probrev[pos][j] = fac;
		    }
		  else
		    {
		      probrev[pos][j] = 0;
		    }
		}
	    }
	  /****Other WM****/
	  prob[pos][numwms+1] = 0;
	  probrev[pos][numwms+1] = 0; 
	  if(otherwmlen <= refpos && otherwmprior > 0)
	    {
	      selec[0] = 1;
	      /**check which species are ungapped aligned wrt reference****/
	      for(spec=1;spec<numspec;++spec)
		{
		  score = scorewm(al,pos,wms[numwms+1].len,wms[numwms+1].mat,spec);
		  if(score > selec_cutoff)
		    selec[spec] = 1;
		  else
		    selec[spec] = 0;
		}
	      /***get full score going backwards***/
	      relpos = otherwmlen-1;
	      /**get number that this column corresponds to****/
	      index = getindex(al->columns[pos],selec);
	      fac = otherwmprior * otherwm[index];
	      --relpos;
	      abspos = pos-1;
	      while(relpos >= 0)
		{
		  /***only positions with base in reference****/
		  if(al->columns[abspos][0] != '-')
		    {
		      index = getindex(al->columns[abspos],selec);
		      fac *= otherwm[index]*F[abspos];
		      --relpos;
		    }
		  --abspos;
		}
	      F[pos] += fac;
	      prob[pos][numwms+1] = fac;
	    }
	  /***store 1/F[pos] ****/
	  Z[pos] = F[pos];
	  F[pos] = 1.0/F[pos];
	}
      /**column has gap in reference score ratio = 1***/
      else
	{
	  F[pos] = 1.0;
	  Z[pos] = 1.0;
	}
    }

  /***backwards calculation R*****/
  refpos = 0;
  for(pos=(al->len)-1;pos>=0;--pos)
    {
      if(al->columns[pos][0] != '-')
	{
	  ++refpos;
	  /**first bg model**/
	  R[pos] = bgprior;
	  /**wm models***/
	  for(j=0;j<numwms;++j)
	    {
	      /**check wm can fit and has nonzero prior****/
	      if(wms[j].len<=refpos && wms[j].prior > 0)
		{
		  if(restrictparses)
		    {
		      score = scorerefwmforward(al,pos,wms[j].len,wms[j].mat);
		      if(score > selec_cutoff)
			selec[0] = 1;
		      else 
			selec[0] = 0;
		    }
		  else
		    selec[0] = 1;

		  if(selec[0]>0)
		    {
		      /**get selection for the other leafs****/
		      for(spec=1;spec<numspec;++spec)
			{
			  if(scorewmforward(al,pos,wms[j].len,wms[j].mat,spec) > selec_cutoff)
			    selec[spec] = 1;
			  else
			    selec[spec] = 0;
			}
		     
		      /**Get structure of the reduced tree that runs only on these species***/
		      thistree = tree_from_selec(selec);

		      /***get full score going forwards***/
		      relpos = 0;
		      fac = 0.5 * wms[j].prior * ratiocolumn(al->columns[pos],wms[j].mat[relpos],thistree);
		      ++relpos;
		      abspos = pos+1;
		      while(relpos < wms[j].len)
			{
			  /***only positions with base in reference****/
			  if(al->columns[abspos][0] != '-')
			    {
			      fac *= ratiocolumn(al->columns[abspos],wms[j].mat[relpos],thistree)*R[abspos];
			      ++relpos;
			    }
			  ++abspos;
			}
		      R[pos] += fac;
		    }
		  /***now reverse strand motif******/
		  if(restrictparses)
		    {
		      score = scorerefwmforward(al,pos,wms[j].len,wms[j].matrev);
		      if(score > selec_cutoff)
			selec[0] = 1;
		      else 
			selec[0] = 0;
		    }
		  else
		      selec[0] = 1;
		  
		  if(selec[0] > 0)
		    {
		      /**get selection for the other leafs****/
		      for(spec=1;spec<numspec;++spec)
			{
			  if(scorewmforward(al,pos,wms[j].len,wms[j].matrev,spec) > selec_cutoff)
			    selec[spec] = 1;
			  else
			    selec[spec] = 0;
			}
		      /**Get structure of the reduced tree that runs only on these species***/
		      thistree = tree_from_selec(selec);


		      /***get full score going forwards***/
		      relpos = 0;
		      facrev = 0.5*wms[j].prior * ratiocolumn(al->columns[pos],wms[j].matrev[relpos],thistree);
		      /***now all the other columns****/
		      ++relpos;
		      abspos = pos+1;
		      while(relpos < wms[j].len)
			{
			  /***only positions with base in reference****/
			  if(al->columns[abspos][0] != '-')
			    {
			      facrev *= ratiocolumn(al->columns[abspos],wms[j].matrev[relpos],thistree)*R[abspos];
			      ++relpos;
			    }
			  ++abspos;
			}
		      R[pos] += facrev;
		    }
		}
	    }
	  /****Other WM****/
	  if(otherwmlen<=refpos && otherwmprior > 0)
	    {
	      selec[0] = 1;
	      /**check which species are ungapped aligned wrt reference****/
	      for(spec=1;spec<numspec;++spec)
		{
		  score = scorewmforward(al,pos,wms[numwms+1].len,wms[numwms+1].mat,spec);
		  if(score > selec_cutoff)
		    selec[spec] = 1;
		  else
		    selec[spec] = 0;
		}

	      /***get full score going forwards***/
	      relpos = 0;
	      /**get number that this column corresponds to****/
	      index = getindex(al->columns[pos],selec);
	      fac = otherwmprior * otherwm[index];
	      ++relpos;
	      abspos = pos+1;
	      while(relpos < otherwmlen)
		{
		  /***only positions with base in reference****/
		  if(al->columns[abspos][0] != '-')
		    {
		      index = getindex(al->columns[abspos],selec);
		      fac *= otherwm[index]*R[abspos];
		      ++relpos;
		    }
		  ++abspos;
		}
	      R[pos] += fac;
	    }
	  R[pos] = 1.0/R[pos];
	}
      /**column has gap in reference score ratio = 1***/
      else
	{
	  R[pos] = 1.0;
	}
    }
  /**calculate ratio of R/F from pos to end for each pos (needed for posteriors)***/
  R[al->len] = 1.0;
  for(pos=(al->len)-1;pos>=0;--pos)
    {
      R[pos] = F[pos]*R[pos+1]/R[pos];
    }
  /***finally get all the posteriors***/
  refpos = 0;
  for(pos=0;pos<al->len;++pos)
    {
      if(al->columns[pos][0] != '-')
	{

	  ++refpos;
	  /**wm models***/
	  for(j=0;j<numwms;++j)
	    {
	      /**check wm can fit and has nonzero prior****/
	      if(wms[j].len<=refpos && wms[j].prior > 0)
		{
		  score = R[pos+1]*prob[pos][j]*F[pos];
		  prob[pos][j] = score;
		  score = R[pos+1]*probrev[pos][j]*F[pos];
		  probrev[pos][j] = score;
		}
	    }
	  /***count posterior for the background model***/
	  score = R[pos+1]*bgprior*F[pos];
	  prob[pos][numwms] = score;
	  /***other wm****/
	  score = R[pos+1]*prob[pos][numwms+1]*F[pos];
	  prob[pos][numwms+1] = score;
	}
    }
  return;
}

int gettotindex(char *column)
{
  int spec,index;
  char letter;
  int offset = 1;
  index = 0;
  for(spec=0;spec<numspec;++spec)
    {
      letter = column[spec];
      if(letter == '-')
	index += offset;
      offset *= 2;
    }
  return index;
}

int getindex(char *column,int *selec)
{
  int spec,index;
  char letter;
  int offset = 1;
  index = 0;
  for(spec=0;spec<numspec;++spec)
    {
      /***put gap if not under selection***/
      if(selec[spec] == 0)
	{
	  index += 4*offset;
	}
      /***else put the letter****/
      else
	{
	  letter = column[spec];
	  if(letter == 'A' || letter == 'a')
	    index += 0;
	  else if(letter == 'C' || letter == 'c')
	    index += offset;
	  else if(letter == 'G' || letter == 'g')
	    index += 2*offset;
	  else if(letter == 'T' || letter == 't')
	    index += 3*offset;
	  else if(letter == '-')
	    index += 4*offset;
	  else
	    index += 4*offset;/**this is default, i.e. treated as a gap***/
	}
      offset *= 5;
    }
  return index;
}


double scorerefwmforward(alignment *al,int pos,int len,double **mat)
{
  int relpos;
  char refletter;
  double score,fac,bfac;
  static double w[4];
  /***run backwards over reference counting each base***/

  relpos = 0;
  score = 1.0;
  while(relpos < len)
    {
      /***letter in reference species****/
      if(pos >= al->len)
	{
	  return 0;
	}
      refletter = al->columns[pos][0];
      if(refletter != '-')
	{
	  setleafprob(w,refletter);
	  fac = 0;
	  bfac = 0;
	  fac += mat[relpos][0]*w[0];
	  bfac += w[0]*bg[0];
	  fac += mat[relpos][1]*w[1];
	  bfac += w[1]*bg[1];
	  fac += mat[relpos][2]*w[2];
	  bfac += w[2]*bg[2];
	  fac += mat[relpos][3]*w[3];
	  bfac += w[3]*bg[3];
	  score *= (fac/bfac);
	  
	  ++relpos;
	}
      ++pos;
    }

  return score;
}


double scorewmforward(alignment *al,int pos,int len, double **mat,int spec) {

  int relpos;
  char letter,refletter;
  double score,fac,bfac;
  static double w[4];
  /***run backwards over reference counting each base***/

  relpos = 0;
  score = 1.0;
  while(relpos < len)
    {
      /***letter in reference species****/
      if(pos >= al->len)
	{
	  return 0;
	}
      refletter = al->columns[pos][0];
      letter = al->columns[pos][spec];
      if(refletter != '-')
	{
	  /**letter in reference aligned to gap -> return zero***/
	  if(letter == '-')
	    {
	      return 0;
	    }
	  else
	    {
	      setleafprob(w,letter);
	      fac = 0;
	      bfac = 0;
	      fac += mat[relpos][0]*w[0];
	      bfac += w[0]*bg[0];
	      fac += mat[relpos][1]*w[1];
	      bfac += w[1]*bg[1];
	      fac += mat[relpos][2]*w[2];
	      bfac += w[2]*bg[2];
	      fac += mat[relpos][3]*w[3];
	      bfac += w[3]*bg[3];
	      score *= (fac/bfac);
	    }
	  ++relpos;
	}
      else
	{
	  /**gap in reference aligned to letter -> return zero ***/
	  if(letter != '-')
	    {
	      return 0;
	    }
	} 
      ++pos;
    }
  return score;
}


double alrefscore(alignment *al,int pos,int len,double **mat,char *tmpseq)
{
  int relpos;
  char refletter;
  double score;
  /***run backwards over reference counting each base***/
  
  relpos = len-1;
  score = 0.0;
  tmpseq[len] = '\0';
  while(relpos >= 0)
    {
      /***letter in reference species****/
      if(pos < 0)
	{
	  return 0;
	}
      refletter = al->columns[pos][0];
      if(refletter != '-')
	{
	  tmpseq[relpos] = refletter;
	  if(refletter == 'A' || refletter == 'a')
	    score += log(mat[relpos][0]/bg[0]);
	  else if(refletter == 'C' || refletter == 'c')
	    score += log(mat[relpos][1]/bg[1]);
	  else if(refletter == 'G' || refletter == 'g')
	    score += log(mat[relpos][2]/bg[2]);
	  else if(refletter == 'T' || refletter == 't')
	    score += log(mat[relpos][3]/bg[3]);
	  else
	    fprintf(stderr,"cannot get score for letter %c\n",refletter);
	  --relpos;
	}
      --pos;
    }
  return score;
}


double scorerefwm(alignment *al,int pos,int len,double **mat)
{
  int relpos;
  char refletter;
  double score,fac,bfac;
  static double w[4];
  /***run backwards over reference counting each base***/

  relpos = len-1;
  score = 1.0;
  while(relpos >= 0)
    {
      /***letter in reference species****/
      if(pos < 0)
	{
	  return 0;
	}
      refletter = al->columns[pos][0];
      if(refletter != '-')
	{
	  setleafprob(w,refletter);
	  fac = 0;
	  bfac = 0;
	  fac += mat[relpos][0]*w[0];
	  bfac += w[0]*bg[0];
	  fac += mat[relpos][1]*w[1];
	  bfac += w[1]*bg[1];
	  fac += mat[relpos][2]*w[2];
	  bfac += w[2]*bg[2];
	  fac += mat[relpos][3]*w[3];
	  bfac += w[3]*bg[3];
	  score *= (fac/bfac);
	  --relpos;
	}
      --pos;
    }
  return score;
}
  

double scorewm(alignment *al,int pos,int len, double **mat,int spec)
{
  int relpos;
  char letter,refletter;
  double score,fac,bfac;
  static double w[4];
  /***run backwards over reference counting each base***/

  relpos = len-1;
  score = 1.0;
  while(relpos >= 0)
    {
      /***If pos runs off the alignment on the left then we cannot fit window****/
      if(pos < 0)
	{
	  return 0;
	}
      /***letter in reference species****/
      refletter = al->columns[pos][0];
      letter = al->columns[pos][spec];
      if(refletter != '-')
	{
	  /**letter in reference aligned to gap -> return zero***/
	  if(letter == '-')
	    {
	      return 0;
	    }
	  else
	    {
	      setleafprob(w,letter);
	      fac = 0;
	      bfac = 0;
	      fac += mat[relpos][0]*w[0];
	      bfac += w[0]*bg[0];
	      fac += mat[relpos][1]*w[1];
	      bfac += w[1]*bg[1];
	      fac += mat[relpos][2]*w[2];
	      bfac += w[2]*bg[2];
	      fac += mat[relpos][3]*w[3];
	      bfac += w[3]*bg[3];
	      score *= (fac/bfac);
	    }
	  --relpos;
	}
      else
	{
	  /**gap in reference aligned to letter -> return zero ***/
	  if(letter != '-')
	    {
	      return 0;
	    }
	} 
      --pos;
    }
  return score;
}



void setleafprob(double *prob,char letter)
{
  switch (letter)
    {
    case '-': /**gap***/
      prob[0] = -1; /***this signals that this leaf is not to count***/
      prob[1] = -1;
      prob[2] = -1;
      prob[3] = -1;
      break;
    case 'a':
      prob[0] = 1.0;
      prob[1] = 0.0;
      prob[2] = 0.0;
      prob[3] = 0.0;
      break;
    case 'A':
      prob[0] = 1.0;
      prob[1] = 0.0;
      prob[2] = 0.0;
      prob[3] = 0.0;
      break;
    case 'c':
      prob[0] = 0.0;
      prob[1] = 1.0;
      prob[2] = 0.0;
      prob[3] = 0.0;
      break;
    case 'C':
      prob[0] = 0.0;
      prob[1] = 1.0;
      prob[2] = 0.0;
      prob[3] = 0.0;
      break;
    case 'g':
      prob[0] = 0.0;
      prob[1] = 0.0;
      prob[2] = 1.0;
      prob[3] = 0.0;
      break;
    case 'G':
      prob[0] = 0.0;
      prob[1] = 0.0;
      prob[2] = 1.0;
      prob[3] = 0.0;
      break;
    case 't':
      prob[0] = 0.0;
      prob[1] = 0.0;
      prob[2] = 0.0;
      prob[3] = 1.0;
      break;
    case 'T':
      prob[0] = 0.0;
      prob[1] = 0.0;
      prob[2] = 0.0;
      prob[3] = 1.0;
      break;
    case 'r':
      prob[0] = 0.5;
      prob[1] = 0.0;
      prob[2] = 0.5;
      prob[3] = 0.0;
      break;
    case 'R':
      prob[0] = 0.5;
      prob[1] = 0.0;
      prob[2] = 0.5;
      prob[3] = 0.0;
      break;
    case 'y':
      prob[0] = 0.0;
      prob[1] = 0.5;
      prob[2] = 0.0;
      prob[3] = 0.5;
      break;
    case 'Y':
      prob[0] = 0.0;
      prob[1] = 0.5;
      prob[2] = 0.0;
      prob[3] = 0.5;
      break;
    case 'm':
      prob[0] = 0.5;
      prob[1] = 0.5;
      prob[2] = 0.0;
      prob[3] = 0.0;
      break;
    case 'M':
      prob[0] = 0.5;
      prob[1] = 0.5;
      prob[2] = 0.0;
      prob[3] = 0.0;
      break;
    case 'k':
      prob[0] = 0.0;
      prob[1] = 0.0;
      prob[2] = 0.5;
      prob[3] = 0.5;
      break;
    case 'K':
      prob[0] = 0.0;
      prob[1] = 0.0;
      prob[2] = 0.5;
      prob[3] = 0.5;
      break;
    case 's':
      prob[0] = 0.0;
      prob[1] = 0.5;
      prob[2] = 0.5;
      prob[3] = 0.0;
      break;
   case 'S':
      prob[0] = 0.0;
      prob[1] = 0.5;
      prob[2] = 0.5;
      prob[3] = 0.0;
      break;
    case 'w':
      prob[0] = 0.5;
      prob[1] = 0.0;
      prob[2] = 0.0;
      prob[3] = 0.5;
      break;
    case 'W':
      prob[0] = 0.5;
      prob[1] = 0.0;
      prob[2] = 0.0;
      prob[3] = 0.5;
      break;
    case 'h':
      prob[0] = 0.33333333;
      prob[1] = 0.33333333;
      prob[2] = 0.0;
      prob[3] = 0.33333333;
      break;
     case 'H':
      prob[0] = 0.33333333;
      prob[1] = 0.33333333;
      prob[2] = 0.0;
      prob[3] = 0.33333333;
      break;
    case 'b':
      prob[0] = 0.0;
      prob[1] = 0.33333333;
      prob[2] = 0.33333333;
      prob[3] = 0.33333333;
      break;
    case 'B':
      prob[0] = 0.0;
      prob[1] = 0.33333333;
      prob[2] = 0.33333333;
      prob[3] = 0.33333333;
      break;
    case 'v':
      prob[0] = 0.33333333;
      prob[1] = 0.33333333;
      prob[2] = 0.33333333;
      prob[3] = 0.0;
      break;
    case 'V':
      prob[0] = 0.33333333;
      prob[1] = 0.33333333;
      prob[2] = 0.33333333;
      prob[3] = 0.0;
      break;
    case 'd':
      prob[0] = 0.33333333;
      prob[1] = 0.0;
      prob[2] = 0.33333333;
      prob[3] = 0.33333333;
      break;
    case 'D':
      prob[0] = 0.33333333;
      prob[1] = 0.0;
      prob[2] = 0.33333333;
      prob[3] = 0.33333333;
      break;
    case 'n':
      prob[0] = 0.25;
      prob[1] = 0.25;
      prob[2] = 0.25;
      prob[3] = 0.25;
      break;
    case 'N':
      prob[0] = 0.25;
      prob[1] = 0.25;
      prob[2] = 0.25;
      prob[3] = 0.25;
      break;
    default:
      prob[0] = 0.25;
      prob[1] = 0.25;
      prob[2] = 0.25;
      prob[3] = 0.25;
    }
  return;
}

/*******create the data-structure for all the reduced columns******/
reducedtreenode **init_reducedtrees(treenode **nodelist)
{
  double qval;
  int gapstruc,spec,nodecount,node,parindex;
  /***number of possible structures = 2^(numspec-1) because reference is always present****/
  int numstruc = (1 << (numspec-1));
  reducedtreenode **reducedtree;
  reducedtree = (reducedtreenode **) calloc(numstruc,sizeof(reducedtreenode *));
  /**count input letters to each node****/
  int *letters = (int *) calloc(numnodes,sizeof(int));
  
  /****run over all possible gap structures****/
  for(gapstruc=0;gapstruc<numstruc;++gapstruc)
    {
      /****initialize input letters to zero*****/
      for(spec=0;spec<numnodes;++spec)
	letters[spec] = 0;

      /***reference always has a letter****/
      letters[0] = 1;
      /***input from reference to its parent***/
      letters[nodelist[0]->parent->index] += 1;
      nodecount = 1;/***keep track of the number of nodes in this tree**/
      for(spec=1;spec<numspec;++spec)
	{
	  /***read whether this species has a letter in current gap structure***/
	  letters[spec] = ((gapstruc >> (spec-1)) & 1);
	  if(letters[spec])
	    {
	      letters[nodelist[spec]->parent->index] += 1;
	      ++nodecount;/***one node for every species with a letter***/
	    }
	}
      /***go over internal nodes of the reference tree**/
      for(spec=numspec;spec<numnodes;++spec)
	{
	  /***add input to parent if current is with inputs***/
	  if(letters[spec] > 0 && nodelist[spec]->parent != 0)
	    letters[nodelist[spec]->parent->index] += 1;
	}
      /****count internal nodes that are used*****/
      for(spec=numspec;spec<numnodes;++spec)
	{
	  if(letters[spec] > 1)
	    ++nodecount;/***another node when at least 2 input*****/
	}
      /***allocate memory for this tree. One extra to indicate when we 'ran off' the tree (= for error checking)*****/
      reducedtree[gapstruc] = (reducedtreenode *) calloc(nodecount+1,sizeof(reducedtreenode));
      /******Now fill in the information for the nodes of this tree****/
      node = 0;
      for(spec=0;spec<numspec;++spec)
	{
	  if(letters[spec] > 0)
	    {
	      /****index of this node******/
	      reducedtree[gapstruc][node].index = spec;
	      /***go find parent index and qvalue******/
	      parindex = nodelist[spec]->parent->index;
	      qval = nodelist[spec]->q;
	      while(letters[parindex] < 2 && nodelist[parindex]->parent != 0)
		{
		  /**multiply by q for the branch going from parent***/
		  qval *= nodelist[parindex]->q;
		  /***step up to the next parent****/
		  parindex = nodelist[parindex]->parent->index;
		}
	      /***parent found *****/
	      if(letters[parindex] >= 2)
		{
		  reducedtree[gapstruc][node].parent = parindex;
		  reducedtree[gapstruc][node].q = qval;
		}
	      /***otherwise this node becomes the root***/
	      else
		{
		  reducedtree[gapstruc][node].parent = -1;
		  reducedtree[gapstruc][node].q = -1;
		}
	      ++node;
	    }
	}
      /*******now go and do this for the internal nodes*********/
       for(spec=numspec;spec<numnodes;++spec)
	 {
	   /***check this node is used******/
	   if(letters[spec] >= 2)
	     {
	      /****index of this node******/
	      reducedtree[gapstruc][node].index = spec;
	      /***If this was the original root***/
	      if(nodelist[spec]->parent == 0)
		{
		  reducedtree[gapstruc][node].parent = -1;
		  reducedtree[gapstruc][node].q = -1;
		}
	      /***Else this was an internal node******/
	      else
		{
		  parindex = nodelist[spec]->parent->index;
		  qval = nodelist[spec]->q;
		  while(letters[parindex] < 2 && nodelist[parindex]->parent != 0)
		    {
		      /**multiply by q for the branch going from parent***/
		      qval *= nodelist[parindex]->q;
		      /***step up to the next parent****/
		      parindex = nodelist[parindex]->parent->index;
		    }
		  /***parent found *****/
		  if(letters[parindex] >= 2)
		    {
		      reducedtree[gapstruc][node].parent = parindex;
		      reducedtree[gapstruc][node].q = qval;
		    }
		  else
		    {
		      reducedtree[gapstruc][node].parent = -1;
		      reducedtree[gapstruc][node].q = -1;
		    }
		}
	      ++node;
	     }
	 }
       /***indicate we have run off the tree****/
       reducedtree[gapstruc][node].index = -1;
    }  
  return reducedtree;
}


double scorebg(char *column,treenode **nodelist)
{
  int i,base,parbase,k,parent;
  char letter;
  double qval,thissum;
  static double *mat = NULL;
  if(mat == NULL)
    mat = (double *) calloc(4*numnodes,sizeof(double));
  /***set all probabilities to 1********/
  for(i=0;i<(4*numnodes);++i)
    {
      mat[i] = 1.0;
    }
  /*****set leaf probs and calc gapstructure index*******/
    int gapstruc = 0;
  letter = column[0];
  base = 0;
  setleafprob(&(mat[0]),letter);
  for(i=1;i<numspec;++i)
    {
      letter = column[i];
      base += 4;
      if(letter != '-')
	{
	  setleafprob(&(mat[base]),letter);
	  gapstruc += (1 << (i-1));
	}
    }
  
  
  /**get the reduced tree for this gap pattern*****/
  reducedtreenode *thistree = reducedtree[gapstruc];  
    
  i = 0;
  /**for all nodes below root******/
  while(thistree[i].index != -1)
    {
      /***set the probability vector for this node****/
      base = 4*(thistree[i].index);
      thissum = 0;
      for(k=0;k<4;++k)
	{
	  thissum += bg[k] * mat[base+k];
	}
      /***not the root*****/
      if(thistree[i].parent >=0)
	{
	  parbase = 4*(thistree[i].parent);
	  qval = thistree[i].q;
	  for(k=0;k<4;++k)
	    {
	      /**multiply contribution from this branch***/
	      mat[parbase+k] *= (thissum*(1.0-qval)+qval*mat[base+k]);
	    }
	  /***next node in this tree***/
	  ++i;
	}
      /***the root of the tree****/
      else
	{
	  return thissum;
	}
    }
  fprintf(stderr,"ERROR, in scorebg function we ran off the tree without encountering a root\n");
  return 1;
}

reducedtreenode *tree_from_selec(int *selec)
{
  int i;
  int gapstruc = 0;
  for(i=1;i<numspec;++i)
    {
      if(selec[i])
	gapstruc += (1 << (i-1));
    }
  return reducedtree[gapstruc];
}


/***Calculate ratio fg/bg on reduced tree of the species which have the site****/
double ratiocolumn(char *column,double *w,reducedtreenode *thistree)
{
  int i,cur,base,parbase,k;
  char letter;
  double qval,sumfg,sumbg;
  static double *matfg = NULL;
  if(matfg == NULL)
    matfg = (double *) calloc(4*numnodes,sizeof(double));
  static double *matbg = NULL;
  if(matbg == NULL)
    matbg = (double *) calloc(4*numnodes,sizeof(double));

  /***set all probabilities to 1********/
  for(i=0;i<4*numnodes;++i)
    {
      matfg[i] = 1.0;
      matbg[i] = 1.0;
    }

  i = 0;
  /**for all nodes******/
  while(thistree[i].index != -1)
    {
      cur = thistree[i].index;
      base = 4*cur;
      /***if it is a leaf, put its probability***/
      if(cur<numspec)
	{
	  letter = column[cur];
	  setleafprob(&(matfg[base]),letter);
	  setleafprob(&(matbg[base]),letter);
	}
      sumfg = 0;
      sumbg = 0;
      for(k=0;k<4;++k)
	{
	  sumfg += w[k]*matfg[base+k];
	  sumbg += bg[k]*matbg[base+k];
	}
      /***not the root*****/
      if(thistree[i].parent >=0)
	{
	  qval = thistree[i].q;
	  parbase = 4*(thistree[i].parent);
	  for(k=0;k<4;++k)
	    {
	      /**multiply contribution from this branch***/
	      matfg[parbase+k] *= (sumfg*(1.0-qval)+qval*matfg[base+k]);
	      matbg[parbase+k] *= (sumbg*(1.0-qval)+qval*matbg[base+k]);
	    }
	  /**next node in the tree***/
	  ++i;
	}
      /***the root of the tree****/
      else
	{
	  return (sumfg/sumbg);
	}
    }
  fprintf(stderr,"ERROR, in ratiocolumn function we ran off the tree without encountering a root\n");
  return 1;
}


/***Read the names of the sequences in the alignments and map them to the tree*/
int parsenames(alignment *alignments,treenode **nodelist)
{
  int g,pos,spec,refspec,found,refnum,foundrefspec,reflen,refindex;
  char *thisname,*refname,*occur, **tmpcolumns, **tmpnames;

  /***stop compiler warning**/
  refnum = -1;
  /*printf("IN function PARSENAMES number of groups %d\n",numgroups);*/
  /***go over all groups****/
  for(g=0;g<numgroups;++g)
    {
      /*printf("group %d\n",g);*/
      foundrefspec = 0;
      /***find the indices for the species***/
      for(spec=0;spec<alignments[g].numspec;++spec)
	{
	  /*printf("name %s\n",alignments[g].specnames[spec]);*/
	  thisname = alignments[g].specnames[spec];
	  found = 0;
	  for(refspec=0;refspec<numspec;++refspec)
	    {
	      /***only check leafs***/
	      refname = nodelist[refspec]->name;
	      /**check if refname occurs in thisname****/
	      occur = strstr(thisname,refname);
	      /**first occurrence of substring***/
	      if(occur != NULL && found == 0)
		{
		  found = 1;
		  refnum = refspec;
		}
	      /**found species name more than once****/
	      else if(occur != NULL && found == 1)
		{
		  fprintf(stderr,"Error: the name %s contains the name of more than one species identifier from the tree\nincluding %s and %s\n",thisname,refname,nodelist[refnum]->name);
		  return 1;
		}
	    }
	  /***found no species name****/
	  if(found == 0)
	    {
	      fprintf(stderr,"Error: found no species identifier in sequence name %s group %d spec %d\n",thisname,g,spec);
	      return 1;
	    }
	  /**set reference number***/
	  alignments[g].specindex[spec] = refnum;
	  /*printf("seq %s has specnum %d with name %s\n",alignments[g].specnames[spec],refnum,nodelist[refnum]->name);*/
	  /***check if this is the reference species***/
	  if(refnum == 0)
	    {
	      foundrefspec = 1;
	      reflen = 0;
	      refindex = spec;
	    }
	  /* printf("species %d name %s mapped to num %d\n",spec,alignments[g].specnames[spec],alignments[g].specindex[spec]);*/
	}
    
      /**check we found the reference species***/
      if(foundrefspec == 0)
	{
	  fprintf(stderr,"Error: did not find the reference species in group %d\n",g);
	  return 1;
	}
     
      /***now refill the columns taking into account the right order****/
      tmpcolumns = (char **) calloc(alignments[g].len,sizeof(char *));
      for(pos=0;pos<alignments[g].len;++pos)
	{
	  tmpcolumns[pos] = (char *) calloc(numspec,sizeof(char));
	  /***initialize all to zero***/
	  for(spec=0;spec<numspec;++spec)
	    {
	      tmpcolumns[pos][spec] = 0;
	    }
	}
      /***Go fill in all columns with species in right order*****/
      for(pos=0;pos<alignments[g].len;++pos)
	{
	  /*printf("pos %d: ",pos);*/
	  /****go over all species****/
	  for(spec=0;spec<alignments[g].numspec;++spec)
	    {
	      tmpcolumns[pos][alignments[g].specindex[spec]] = alignments[g].columns[pos][spec];
	    }
	  /***Fill in dashes for absent species***/
	  for(spec=0;spec<numspec;++spec)
	    {
	      if(tmpcolumns[pos][spec] == 0)
		{
		  tmpcolumns[pos][spec] = '-';
		}
	      /*printf("spec %d: %c ",spec,tmpcolumns[pos][spec]);*/
	    }
	  /*printf("\n");*/
	}
      /***Now we clear all the memory associates with columns****/
      for(pos=0;pos<alignments[g].len;++pos)
	{
	  free(alignments[g].columns[pos]);
	}
      free(alignments[g].columns);
      /***now let the columns point to tmpcolumns***/
      alignments[g].columns = tmpcolumns;

      /***udate the names***/
      tmpnames = (char **) calloc(numspec,sizeof(char *));
      for(spec=0;spec<numspec;++spec)
	{
	  tmpnames[spec] = (char *) calloc(1024,sizeof(char));
	  strcpy(tmpnames[spec],"");
	}
      for(spec=0;spec<alignments[g].numspec;++spec)
	{
	  strcpy(tmpnames[alignments[g].specindex[spec]],alignments[g].specnames[spec]);
	}
      for(spec=0;spec<alignments[g].numspec;++spec)
	{
	  free(alignments[g].specnames[spec]);
	}
      free(alignments[g].specnames);
      alignments[g].specnames = tmpnames;
      alignments[g].type = get_region_type(alignments[g].specnames[0]);
      free(alignments[g].specindex);
    }
  
  return 0;
}

int get_region_type(char *name)
{
  int occur;
  if(strstr(name,"|A|") != NULL)
    return 0;
  else if(strstr(name,"|C|") != NULL)
    return 3;
  /****for now we are assuming everything else is B of type 1****/
  else
    return 1;
}



int readsequences(alignment *alignments,char *filename,int numgroups)
{
  FILE *inputfile;
  int c;
  int curgroup = 0;
  int specnum = 0;
  int pos = 0;
  char line[1024], *retval;
  inputfile = fopen(filename,"r");
  if(inputfile == NULL)
    {
      fprintf(stderr,"input file %s doesn't exist or cannot be opened\n",filename);
      return 1;
    }

  for(curgroup=0;curgroup<numgroups;++curgroup)
    {
      retval = fgets(line,1024,inputfile);
      /***read until one gets a line where first two are > > ****/
      while( retval != NULL && !(line[0] == '>' && line[1] == '>'))
	retval = fgets(line,1024,inputfile);

      if(retval == NULL)
	{
	  fprintf(stderr,"Error: group %d doesn't start with >>\n",curgroup);
	  return 1;
	}
      mychomp(line);
      strcpy(alignments[curgroup].specnames[0],line+2);
      /*printf("name first sequence:%s\n",alignments[curgroup].specnames[0]);*/
      /***read the letters****/
      pos = 0;
      while(pos<alignments[curgroup].len)
	{
	  c = fgetc(inputfile);
	  if(!isspace(c))
	    {
	      if( c=='N' ) 
		c = '-';
	      
	      alignments[curgroup].columns[pos][0] = c;
	      ++pos;
	    }
	}
      /**read to end of line****/
      c = fgetc(inputfile);
      while(c != '\n')
	{
	  c = fgetc(inputfile);
	}
      /***read in the names and sequences of the other species****/
      for(specnum=1;specnum<alignments[curgroup].numspec;++specnum)
	{
	  retval = fgets(line,1024,inputfile);
	  /**keep reading until one finds line that starts with > ****/
	  /***read until one gets a line where first two are > > ****/
	  while( retval != NULL && !(line[0] == '>'))
	    retval = fgets(line,1024,inputfile);
	  if(retval == NULL)
	    {
	      fprintf(stderr,"Error: group %d specnum %d doesn't start with >\n",curgroup,specnum);
	      return 1;
	    }
	  mychomp(line);
	  strcpy(alignments[curgroup].specnames[specnum],line+1);
	  /*printf("name sequence %d is:%s\n",specnum,alignments[curgroup].specnames[specnum]);*/
	  pos = 0;
	  while(pos<alignments[curgroup].len)
	    {
	      c = fgetc(inputfile);
	      if(!isspace(c))
		{
		  if( c == 'N' ) { c = '-'; }

		  alignments[curgroup].columns[pos][specnum] = c;
		  ++pos;
		}
	    }
	  /**read to end of line****/
	  c = fgetc(inputfile);
	  while(c != '\n')
	    {
	      c = fgetc(inputfile);
	    }
	}
    }
  fclose(inputfile);

  return 0;
}


int getseqsizes(char *filename,int *numgroups,int *groupsizes,int *seqlengths)
{
  FILE *inputfile;
  int c;
  int seqnum = 0;
  inputfile = fopen(filename,"r");
  if(inputfile == NULL)
    {
      fprintf(stderr,"input file %s doesn't exist or cannot be opened\n",filename);
      return 1;
    }
  
  c = fgetc(inputfile);
  while(c != EOF)
    {
      /*printf("read characeter %c\n",c);*/
      if(c == '>')
	{
	  ++seqnum;
	  seqlengths[seqnum-1] = 0;
	  c = fgetc(inputfile);
	  /***new group ***/
	  if(c == '>')
	    {
	      /**new group**/
	      ++(*numgroups);
	      groupsizes[*numgroups-1] = 1;
	    }
	  /**new sequence in group***/
	  else
	    {
	      ++groupsizes[*numgroups-1];
	    }
	  /**now read to the end of the nameline****/
	    while(c != '\n')
	      {
		c = fgetc(inputfile);
	      }
	}
      else if(!isspace(c))
	{
	  ++seqlengths[seqnum-1];
	}
      c = fgetc(inputfile);
    }
  fclose(inputfile);
  return 0;
}


char consensus(double na, double nc, double ng, double nt) 
{
    double ntot;
    char c;
    
    ntot=na+nc+ng+nt;
    /* changed 0.5 to 0.51 because of pathology that with
       equally-divided bases roundoff error determines which is picked */
    if (na > 0.51*ntot)
        c='A';
    else if (nc > 0.51*ntot)
        c='C';
    else if (ng > 0.51*ntot)
        c='G';
    else if (nt > 0.51*ntot)
        c='T';
    else if ((na+nc) > 3.0*ntot/4.0)
        c='M';
    else if ((na+ng) > 3.0*ntot/4.0)
        c='R';
    else if ((na+nt) > 3.0*ntot/4.0)
        c='W';
    else if ((nc+ng) > 3.0*ntot/4.0)
        c='S';
    else if ((nc+nt) > 3.0*ntot/4.0)
        c='Y';
    else if ((ng+nt) > 3.0*ntot/4.0)
        c='K';
    else if (na+nc+ng > 7.0*ntot/8.0)
        c='V';
    else if (na+nc+nt > 7.0*ntot/8.0)
        c='H';
    else if (na+ng+nt > 7.0*ntot/8.0)
        c='D';
    else if (nc+ng+nt > 7.0*ntot/8.0)
        c='B';
    else
        c='N';
    return c;
}


double getinf(double *n)
{
  /***information score = log(prob) - log(background prob) = ****/
  /*** gamma(2.0) - gamma(sum+2.0) + sum_i gamma(ni+0.5)-gamma(0.5) ***/
  /*** -sum *log(1.0/4.0)***/

  int i;
  double inf;
  double sum = 0;
  for(i=0;i<4;++i)
    {
      sum += n[i];
    }
  inf = sum * log(4.0);
  inf += lgamma(2.0);
  inf -= lgamma(sum+2.0);
  for(i=0;i<4;++i)
    {
      if(n[i]>0)
	{
	  inf += (lgamma(n[i]+0.5)-lgamma(0.5));
	}
    }
  inf /= (sum*log(2.0));
  return inf;
}


char* twodigitstr(int m)
{
    char *c;

    c = (char *) malloc(4);
    if ((m<0)||(m>99)) {
        c[0]='E';
        c[1]='E';
    }
    else {
        c[0]=digits[m/10];
        c[1]=digits[m%10];
    }
    c[2]='\0';
    return c;
}


void printwm(wm *thiswm)
{
  char cons;
  int pos,letter;
  double inf;
  double n[4];
  printf("//\n");
  printf("NA %s\n",thiswm->name);
  printf("PW %lf\n",thiswm->prior);
  printf("%2s    %6c     %6c     %6c     %6c     %6s      %6s\n","PO",'A','C','G','T',"cons","inf");
  for(pos=0;pos<thiswm->len;++pos)
    {
      for(letter=0;letter<4;++letter)
	{
	  inf = -0.5+(thiswm->colsum[pos]+2.0)*(thiswm->mat[pos][letter]);
	  if(inf<0)
	    {
	      inf = 0;
	    }
	  n[letter] = inf;
	}
      inf = getinf(n);
      cons = consensus(n[0],n[1],n[2],n[3]);
      printf("%2s    %6.2f     %6.2f     %6.2f     %6.2f     %6c      %6.2f\n",twodigitstr(pos), n[0],n[1],n[2],n[3],cons,inf);
    }
  printf("//\n");

  return;
}


 int readwms(char *filename,wm *wms,int *num_with_prior,double *wm_prior_sum)
{
  FILE *wmfile;
  char s[1024];
  double tmpmat[100][4],sum; /***temporary matrix for storing the WM entries. Max length 100 ****/
  int inmotif,lastpos,thislen,thispos,j,k,i;
  int wmnum = 0;
  int matpos = 0;
  int pos,letnum;
  double na,nc,ng,nt,weight;
  wmfile = fopen(filename,"r");
  if(wmfile == NULL)
    {
      fprintf(stderr,"WM file %s doesn't exist or cannot be opened\n",filename);
      return 1;
    }
  inmotif = 0;
  thislen = 0;
  lastpos = 0;

  while (fgets(s,1024,wmfile)) /**get line from the file***/
    {
      /***line with start of motif indicator***/
      if(inmotif == 0 && s[0] == '/' && s[1] == '/')
	{
	  inmotif = 1;
	  lastpos = 0;
	  inmotif = 1;
	  thislen = 0;
	  matpos = 0;
	  wms[wmnum].prior = 0.0;
	}
      /**name line, already in motif***/
      else if((s[0] == 'N' && s[1] == 'A' && inmotif == 1) || (s[0] == 'I' && s[1] == 'D' && inmotif == 1))
	{
	  /**set name line***/
	  mychomp(s);
	  pos = 2;
	  while(isspace(s[pos]) && s[pos] != '\0'){
	    ++pos;
	  }
	  if(s[pos] == '\0'){
	    --pos;
	  }
	  if(pos < 2){
	    pos = 2;
	  }
	  strcpy(wms[wmnum].name,s+pos);
	}
      /**name line starting a new motif****/
      else if((s[0] == 'N' && s[1] == 'A' && inmotif == 0) || (s[0] == 'I' && s[1] == 'D' && inmotif == 0))
	{
	  lastpos = 0;
	  inmotif = 1;
	  thislen = 0;
	  matpos = 0;
	  wms[wmnum].prior = 0.0;
	  /**set name line***/
	  mychomp(s);
	  pos = 2;
	  while(isspace(s[pos]) && s[pos] != '\0'){
	    ++pos;
	  }
	  if(s[pos] == '\0'){
	    --pos;
	  }
	  if(pos < 2){
	    pos = 2;
	  }
	  strcpy(wms[wmnum].name,s+pos);
	}
      /***line with counts***/
      else if(myisnumber(s[0]) &&( myisnumber(s[1]) || isblank(s[1])))
	{
	  if(inmotif == 0)
	    {
	      fprintf(stderr, "error: encountered a line with motif counts without having read ID (identifier) or NA (name) for motif %d\n",wmnum);
	      return 1;
	    }
	  sscanf(s,"%d %lf %lf %lf %lf",&thispos,&na,&nc,&ng,&nt);
	  if(thispos-lastpos != 1 && !(thispos == 0 && lastpos == 0))
	    {
	      fprintf(stderr,"Warning: position numbering went from %d to %d in input motif %s number %d\n",lastpos,thispos,wms[wmnum].name,wmnum);
	    }
	  lastpos = thispos;
	  tmpmat[matpos][0] = na;
	  tmpmat[matpos][1] = nc;
	  tmpmat[matpos][2] = ng;
	  tmpmat[matpos][3] = nt;
	  ++matpos;
	  ++thislen;
	}
      /**end of a WM***/
      else if( s[0] == '/' && s[1] == '/' && inmotif == 1)
	{
	  wms[wmnum].len = thislen;
	  wms[wmnum].mat = (double **) calloc(thislen,sizeof(double *));
	  for(matpos=0;matpos<thislen;++matpos)
	    {
	      wms[wmnum].mat[matpos] = (double *) calloc(4,sizeof(double));
	      for(letnum = 0;letnum<4;++letnum)
		{
		  wms[wmnum].mat[matpos][letnum] = tmpmat[matpos][letnum];
		}
	    }
	  inmotif = 0;
	  ++wmnum;
	}
      /***line prior weight****/
      else if(s[0] == 'P' && s[1] == 'W' && inmotif == 1)
        {
          sscanf(s,"PW %lf\n",&weight);
          wms[wmnum].prior = weight;
        }

    }
  fclose(wmfile);

  /*motif ended without a // **/
  if(inmotif)
    {
      wms[wmnum].len = thislen;
      wms[wmnum].mat = (double **) calloc(thislen,sizeof(double *));
	  for(matpos=0;matpos<thislen;++matpos)
	    {
	      wms[wmnum].mat[matpos] = (double *) calloc(4,sizeof(double));
	      for(letnum = 0;letnum<4;++letnum)
		{
		  wms[wmnum].mat[matpos][letnum] = tmpmat[matpos][letnum];
		}
	    }
	  inmotif = 0;
	  ++wmnum;
    }
  /*printf("there were %d wms\n",wmnum);*/
  /***add pseudocounts, normalize, and create the matrix of log values****/
  for(i=0;i<wmnum;++i)
    {
      /*printf("WM %s\n",wms[i].name);*/
      *wm_prior_sum += wms[i].prior;
      if(wms[i].prior > 0)
	{
	  ++(*num_with_prior);
	}
      /***memory for log probability***/
      wms[i].lmat = (double **) calloc(wms[i].len,sizeof(double *));
      /***memory for the reverse complement matrix****/
      wms[i].matrev = (double **) calloc(wms[i].len,sizeof(double *));
      wms[i].lmatrev = (double **) calloc(wms[i].len,sizeof(double *));
      /***total counts in each column****/
      wms[i].colsum = (double *) calloc(wms[i].len,sizeof(double));
      /***go add pseudo-counts to all entries and normalize****/
      for(j=0;j<wms[i].len;++j)
	{
	  /*printf("pos %d ",j);*/
	  sum = 0;
	  for(k=0;k<4;++k)
	    {
	      sum += (wms[i].mat[j][k]+0.5);
	    }
	  wms[i].lmat[j] = (double *) calloc(4,sizeof(double));
	  wms[i].matrev[j] = (double *) calloc(4,sizeof(double));
	  wms[i].lmatrev[j] = (double *) calloc(4,sizeof(double));
	  wms[i].colsum[j] = sum-2.0;
	  for(k=0;k<4;++k)
	    {
	      wms[i].mat[j][k] = (wms[i].mat[j][k]+0.5)/sum;
	      /*printf("%g ",wms[i].mat[j][k]);*/
	      wms[i].lmat[j][k] = log(wms[i].mat[j][k]);
	    }
	  /*printf("\n");*/
	}
      /***set the reverse WMs****/
      for(j=0;j<wms[i].len;++j)
	{
	  /*printf("reverse pos %d ",j);*/
	  for(k=0;k<4;++k)
	    {
	      wms[i].matrev[j][k] = wms[i].mat[wms[i].len-1-j][3-k];
	      /*printf("%g ",wms[i].matrev[j][k]);*/
	      wms[i].lmatrev[j][k] = wms[i].lmat[wms[i].len-1-j][3-k];
	    }
	  /*printf("\n");*/
	}

      /*printwm(&(wms[i]));*/
    }

  return (0);
}

int myisnumber(char s)
{
  if(s == '0' || 
     s == '1' || 
     s == '2' || 
     s == '3' ||
     s == '4' ||
     s == '5' ||
     s == '6' || 
     s == '7' ||
     s == '8' ||
     s == '9'   )
    {
      return 1;
    }
  else{
    return 0;
  }
}



int countwms(char *filename)
{
  FILE *wmfile;
  char s[1024];
  int numwms = 0;
  int inmotif;
  wmfile = fopen(filename,"r");
  if(wmfile == NULL)
    {
      fprintf(stderr,"WM file %s doesn't exist or cannot be opened\n",filename);
      return 1;
    }
  inmotif = 0;
   while (fgets(s,1024,wmfile)) /**get line from the file***/
    {
      mychomp(s);
      if(inmotif == 0 && s[0] == '/' && s[1] == '/')
	{
	  ++numwms;
	  inmotif = 1;
	}
       else if((s[0] == 'N' && s[1] == 'A' && inmotif == 0) || (s[0] == 'I' && s[1] == 'D' && inmotif == 0))
	 {
	   ++numwms;
	   inmotif = 1;
	 }
      else if( s[0] == '/' && s[1] == '/' && inmotif == 1)
	{
	  inmotif = 0;
	}
    }
   fclose(wmfile);
   return (numwms);
}


int read_param_file(char *filename,char *refspecies,char *ttreestring,int *em_prior,int *em_wm,int *do_segment,int *seglen,int *bgorder,double *otherwmprior,double *bgprior,int *restrictparses,char *otherwmfile,char *sitefile,char *priorfile,char *loglikfile,char *column_statfile,char *cons_statfile,char *mutcountfile,char *segfile,double *minposterior,int *printsiteals,char *otherwm_proffile, int *proflen, int *do_ep, double *priordiff, double *wmdiff, int *winlen, int *steplen, int *do_evo,int *markov_order,char *epfilename, int *useproximities,char *mode, double *minposteriorWM, int *UFEprint, char *mybgfile)
{

  FILE *paramfile;
  char s[treelength],keyword[256],letter;
  int i;
  int gotrefspecies = 0;
  double val;

  paramfile = (FILE *) fopen(filename,"r");
  if(paramfile == NULL)
    {
      fprintf(stderr,"Parameter file %s does not exist or cannot be opened.\nExiting.\n",filename);
      return 1;
    }
  /**reading the file****/
  while (fgets(s,treelength,paramfile)) /**get line from the file***/
    {
      mychomp(s);
      /**read up to first space***/
      i =0;
      while(!isspace(s[i]) && s[i] != '\0')
	{
	  keyword[i] = s[i];
	  ++i;
	}
      keyword[i] = '\0';
      
      /***now read parameter depending on keyword****/
      if(strcmp(keyword,"refspecies") == 0)
	{
	  /**find next nonwhite space***/
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  /**copy name into refspecies***/
	  strcpy(refspecies,s+i);
	  gotrefspecies = 1;
	}
      else if(strcmp(keyword,"EMprior") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  if(s[i] == '1')
	    {
	      *em_prior = 1;
	    }
	  else if(s[i] == '0')
	    {
	      *em_prior = 0;
	    }
	  else
	    {
	      fprintf(stderr,"Error: EMprior should be 1 (do EM) or 0 (no EM)\n");
	      return 1;
	    }
	}
      else if(strcmp(keyword,"xxxEMwm") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  if(s[i] == '1')
	    {
	      *em_wm = 1;
	    }
	  else if(s[i] == '0')
	    {
	      *em_wm = 0;
	    }
	  else
	    {
	      fprintf(stderr,"Error: EMwm should be 1 (do EM) or 0 (no EM)\n");
	      return 1;
	    }
	}
      else if(strcmp(keyword,"UseProximities") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  if(s[i] == '1')
	    {
	      *useproximities = 1;
	    }
	  else if(s[i] == '0')
	    {
	      *useproximities = 0;
	    }
	  else
	    {
	      fprintf(stderr,"Error: UseProximities be 1 (use proximities) or 0 (use time)\n");
	      return 1;
	    }
	}
      else if(strcmp(keyword,"xxxDoUFE") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  if(s[i] == '1')
	    {
	      *do_evo = 1;
	    }
	  else if(s[i] == '0')
	    {
	      *do_evo = 0;
	    }
	  else
	    {
	      fprintf(stderr,"Error: DoUFE should be 1 (create file for UFE model) or 0 (don't create file for UFE model)\n");
	      return 1;
	    }
	}

      else if(strcmp(keyword,"xxxDoEP") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  if(s[i] == '1')
	    {
	      *do_ep = 1;
	    }
	  else if(s[i] == '0')
	    {
	      *do_ep = 0;
	    }
	  else
	    {
	      fprintf(stderr,"Error: DoEP should be 1 (do EP) or 0 (no EP)\n");
	      return 1;
	    }
	}
      else if(strcmp(keyword,"TREE") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(ttreestring,s+i);
	}
      else if(strcmp(keyword,"Mode") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(mode,s+i);
	}

      else if(strcmp(keyword,"DOsegment") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  if(s[i] == '1')
	    {
	      *do_segment = 1;
	    }
	  else if(s[i] == '0')
	    {
	      *do_segment = 0;
	    }
	  else
	    {
	      fprintf(stderr,"Error: DOsegment should be 1 (get segment scores) or 0 (no segment scores)\n");
	      return 1;
	    }
	 
	}
      else if(strcmp(keyword,"UFEprint") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  if(s[i] == '1')
	    {
	      *UFEprint = 1;
	    }
	  else if(s[i] == '0')
	    {
	      *UFEprint = 0;
	    }
	  else
	    {
	      fprintf(stderr,"Error: UFEprint should be 1 (print UFE sites) or 0 (don't print it)\n");
	      return 1;
	    }
	 
	}
      else if(strcmp(keyword,"segmentlength") == 0)
	{
	  *seglen = atoi(s+i);
	}
      else if(strcmp(keyword,"winlen") == 0)
	{
	  *winlen = atoi(s+i);
	}
      else if(strcmp(keyword,"markovorderBG") == 0)
	{
	  *markov_order = atoi(s+i);
	}
      else if(strcmp(keyword,"steplen") == 0)
	{
	  *steplen = atoi(s+i);
	}
      else if(strcmp(keyword,"bgorder") == 0)
	{
	  *bgorder = atoi(s+i);
	}
      else if(strcmp(keyword,"UFEwmprior") == 0)
	{
	  *otherwmprior = atof(s+i);	  
	}
      else if(strcmp(keyword,"bgprior") == 0)
	{
	  *bgprior = atof(s+i);
	}
      else if(strcmp(keyword,"priordiff") == 0)
	{
	  *priordiff = atof(s+i);
	}
      else if(strcmp(keyword,"wmdiff") == 0)
	{
	  *wmdiff = atof(s+i);
	}
      else if(strcmp(keyword,"restrictparses") == 0)
	{
	  *restrictparses = atoi(s+i);
	} 
      else if(strcmp(keyword,"UFEwmfile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(otherwmfile,s+i);
	}
      else if(strcmp(keyword,"sitefile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(sitefile,s+i);
	}
      else if(strcmp(keyword,"priorfile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(priorfile,s+i);
	}
      else if(strcmp(keyword,"loglikfile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(loglikfile,s+i);
	}
      else if(strcmp(keyword,"CRMfile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(epfilename,s+i);
	}
      else if(strcmp(keyword,"mybgfile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(mybgfile,s+i);
	}
      else if(strcmp(keyword,"segmentfile") == 0)
        {
          while(isspace(s[i]) && s[i] != '\0')
            {
              ++i;
            }
          strcpy(segfile,s+i);
        }
      else if(strcmp(keyword,"column_statfile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(column_statfile,s+i);
	}
      else if(strcmp(keyword,"cons_statfile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(cons_statfile,s+i);
	}
      else if(strcmp(keyword,"mutcountfile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(mutcountfile,s+i);
	}
      else if(strcmp(keyword,"UFEwmproffile") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  strcpy(otherwm_proffile,s+i);
	}
      else if(strcmp(keyword,"UFEwmlen") == 0)
	{
	  otherwmlen = atoi(s+i);
	}
      else if(strcmp(keyword,"minposterior") == 0){
	*minposterior = atof(s+i);
      }
      else if(strcmp(keyword,"minposteriorWM") == 0){
	*minposteriorWM = atof(s+i);
      }
      else if(strcmp(keyword,"printsiteals") == 0){
	*printsiteals = atoi(s+i);
      }
      else if(strcmp(keyword,"bg") == 0)
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	  sscanf(s+i,"%c %lf\n",&letter,&val);
	  printf("letter %c val is %lf\n",letter,val);
	  if(letter == 'A' || letter == 'a')
	    bg[0] = val;
	  else if(letter == 'C' || letter == 'c')
	    bg[1] = val;
	  else if(letter == 'G' || letter == 'g')
	    bg[2] = val;
	  else if(letter == 'T' || letter == 't')
	    bg[3] = val;
	  else
	    {
	      fprintf(stderr,"Error: cannot recognize letter %c\n",letter);
	      return 1;
	    }
	}
    }
  fclose(paramfile);
  return (0);
}

int read_bg_filenames(char *filename,char **bgfiles,int numnodes,treenode **nodelist)
{
  FILE *paramfile;
  char s[1024],keyword[256],tmpstr[1024];
  int i,pos,index,found;
  paramfile = (FILE *) fopen(filename,"r");
  if(paramfile == NULL)
    {
      fprintf(stderr,"Parameter file %s does not exist or cannot be opened.\nExiting.\n",filename);
      return 1;
    }

  while (fgets(s,1024,paramfile)) /**get line from the file***/
    {
      mychomp(s);
      /**read up to first space***/
      i =0;
      while(!isspace(s[i]) && s[i] != '\0')
	{
	  keyword[i] = s[i];
	  ++i;
	}
      keyword[i] = '\0';
      if(strcmp(keyword,"bgfile") == 0) /**read a background file***/
	{
	  while(isspace(s[i]) && s[i] != '\0')
	    {
	      ++i;
	    }
	   pos = 0;
	   while(!isspace(s[i]) && s[i] != '\0'){
	     tmpstr[pos] = s[i];
	     ++pos;
	     ++i;
	   }
	   tmpstr[pos] = '\0';
	   /**now find the species that goes with this name**/
	   found = 0;
	   for(index=0;index<numnodes;++index)
	     {
	       if(strcmp(nodelist[index]->name,tmpstr) == 0)
		 {
		   found = 1;
		   break;
		 }
	     }
	   if(found == 0)
	     {
	       return (1);
	     }
	   else
	     {
	       while(isspace(s[i]) && s[i] != '\0')
		 {
		   ++i;
		 }
	       strcpy(bgfiles[index],s+i);
	     }
	}
    }
  fclose(paramfile);
  return (0);
}


/*********parse the input tree******************/
int parsetree(char *ttreestring,treenode **nodelist,char *refspecies, int prox)
{
  treenode *thisnode, *descendant;
  treenode **rawlist;
  char numstring[treelength],tmpstr[treelength];
  int numpos,absindex,nodeindex;
  int level = 0;
  int maxlevel =0;
  int pos = 1;
  int found = 0;

  if(ttreestring[0] != '(')
    {
      fprintf(stderr,"Error: treestring should start with (\n");
      return 1;
    }
  thisnode = (treenode *) malloc(sizeof(treenode));
  thisnode->parent = NULL;
  thisnode->level = level;
  strcpy(thisnode->name,"");

  rawlist = (treenode **) calloc(numnodes,sizeof(treenode *));
  nodeindex = 0;
  rawlist[nodeindex] = thisnode;
  ++nodeindex;

  pos = 1;
  while(ttreestring[pos] != ';' && ttreestring[pos] != '\0')
    {
      if(ttreestring[pos] == '(') /***new subtree****/
	{
	  ++level;
	  descendant = (treenode *) malloc(sizeof(treenode));
	  descendant->parent = thisnode;
	  strcpy(descendant->name,"");
	  descendant->level = level;
	  thisnode = descendant;
	  rawlist[nodeindex] = thisnode;
	  ++nodeindex;
	  ++pos;
	}
      else if(ttreestring[pos] == ')') /***end subtree ***/
	{
	  thisnode = thisnode->parent;
	  --level;
	  ++pos;
	}
      else if(ttreestring[pos] == ':') /** proximity parent current node ***/
	{
	  numpos = 0;
	  ++pos;
	  while(isdigit(ttreestring[pos]) || ttreestring[pos] == '.' || isspace(ttreestring[pos]))
	    {
	      numstring[numpos] = ttreestring[pos];
	      ++numpos;
	      ++pos;
	    }
	  numstring[numpos] = '\0';
	  double qtmp = atof(numstring);
	  if( prox == 0) // convert time to proximities
	    qtmp = exp(-qtmp);

	  thisnode->q = qtmp;
	}
      else if(ttreestring[pos] == ',') /** comma separates nodes at same level, return to parent ***/
	{
	  thisnode = thisnode->parent; /***move to the parent***/
	  --level;
	  ++pos;
	}
      else if(!isspace(ttreestring[pos]))/**any other nonspace character is considered name of string ***/
	{
	  /**we are getting a leaf****/
	  descendant = (treenode *) malloc(sizeof(treenode));
	  ++level;
	  ++numspec;
	  descendant->parent = thisnode;
	  descendant->level = level;
	  if(level >= maxlevel)
	    maxlevel = level;
	  thisnode=descendant;
	  rawlist[nodeindex] = thisnode;
	  ++nodeindex;
	  numpos = 0;
	  while(ttreestring[pos] != '(' && ttreestring[pos] != ')' && ttreestring[pos] != ':' && ttreestring[pos] != ',')
	    {
	      thisnode->name[numpos] = ttreestring[pos];
	      ++pos;
	      ++numpos;
	    }

	  /***empty name is not allowed****/
	  if(numpos == 0)
	    {
	      fprintf(stderr,"Eror: leaf in tree has empty name\n");
	      return 1;
	    }
	  thisnode->name[numpos] = '\0';
	}
      else
	{
	  ++pos;
	}
    }

  /**First get the reference********/
  for(nodeindex=0;nodeindex<numnodes;++nodeindex)
    {
      thisnode = rawlist[nodeindex];
      if(strcmp(thisnode->name,refspecies) == 0)
	{
	  nodelist[0] = thisnode;
	  thisnode->index = 0;
	  found = 1;
	  break;
	}
    }
  if(found == 0)
    {
      fprintf(stderr,"no nodename matches the reference species name %s\n",refspecies);
      return 1;
    }

  /***now get all the leafs of the tree****/
  absindex = 1;
  for(nodeindex=0;nodeindex<numnodes;++nodeindex)
    {
      thisnode = rawlist[nodeindex];
      /**leaf that is not reference****/
      if(strcmp(thisnode->name,"") != 0 && strcmp(thisnode->name,refspecies) != 0)
	{
	  nodelist[absindex] = thisnode;
	  thisnode->index = absindex;
	  ++absindex;
	}
    }
  
  /***Go through internal nodes starting at deepest level***/
  for(level=maxlevel;level>=0;--level)
    {
      for(nodeindex=0;nodeindex<numnodes;++nodeindex)
	{
	  thisnode = rawlist[nodeindex];
	  if(thisnode->level == level && strcmp(thisnode->name,"") == 0)
	    {
	      nodelist[absindex] = thisnode;
	      thisnode->index = absindex;
	      ++absindex;
	    }
	}
    }
  /***Finally fill in names for the internal nodes***/
  for(absindex=0;absindex<numnodes;++absindex)
    {
      thisnode = nodelist[absindex];
      if(thisnode->parent != NULL)
	{
	  strcpy(tmpstr,thisnode->name);
	  strcat(tmpstr,"_");
	  strcat(thisnode->parent->name,tmpstr);
	}
    }
  free(rawlist);
  return 0;
}





void mychomp(char *s) 
{
  if (s[strlen(s)-1]=='\n')
    s[strlen(s)-1]='\0';
  if (s[strlen(s)-1]=='\r')
    s[strlen(s)-1]='\0';
  
  return;
}

char revbase(char base)
{
  if(base == 'A')
    return 'T';
  else if(base == 'C')
    return 'G';
  else if(base == 'G')
    return 'C';
  else if(base == 'T')
    return 'A';
  else if(base == 'N')
    return 'N';
  else if(base == 'W')
    return 'W';
  else if(base == 'S')
    return 'S';
  else if(base == 'Y')
    return 'R';
  else if(base == 'R')
    return 'Y';
  else if(base == 'a')
    return 't';
  else if(base == 'c')
    return 'g';
  else if(base == 'g')
    return 'c';
  else if(base == 't')
    return 'a';
  else if(base == 'n')
    return 'n';
  else if(base == 'w')
    return 'w';
  else if(base == 's')
    return 's';
  else if(base == 'y')
    return 'r';
  else if(base == 'r')
    return 'y';
  else
    return 'N';
}



/* for the wm refinement */
void update_wms(alignment *al,double **prob,double **probrev,wm *wms,FILE *sitefile,double minposterior,int printsiteals, double ***postcounts) {
  
  int pos, refpos,j,basecount,thispos,spec,startpos;
  static double *selec;
  double totsum = 0.0;
  double score,alscore;
  char tmpseq[100];
  selec = (double *) calloc(numspec,sizeof(double));
  
  refpos = 0;
  for(pos=0;pos<al->len;++pos) {
    
    if(al->columns[pos][0] != '-') {
      
      ++refpos;
      /**wm models***/
      for(j=0;j<(numwms);++j) {
	
	/**check wm can fit and has nonzero prior****/
	if(wms[j].len<=refpos && wms[j].prior > 0) {
	  if(prob[pos][j] >= minposterior) {
	    
	    /****print the bases at the positions in site if wanted****/
	    
	    startpos = pos;
	    basecount = 1;
	    while(basecount < wms[j].len) {
	      
	      --startpos;
	      if(al->columns[startpos][0] != '-')
		++basecount;
	    }
	    selec[0] = scorerefwm(al,pos,wms[j].len,wms[j].mat);
	    spec = 0;
	    
	    /***now print out the alignment****/
	    if(selec[spec] > selec_cutoff) {
	      
	      int counter = -1;
	      for(thispos=startpos;thispos<=pos;++thispos) {
		
		if(al->columns[thispos][0] != '-') {
		  counter++;
		  if (al->columns[thispos][spec]=='A' || al->columns[thispos][spec]=='a') { postcounts[j][counter][0] += prob[pos][j]; }
		  if (al->columns[thispos][spec]=='C' || al->columns[thispos][spec]=='c') { postcounts[j][counter][1] += prob[pos][j]; }
		  if (al->columns[thispos][spec]=='G' || al->columns[thispos][spec]=='g') { postcounts[j][counter][2] += prob[pos][j]; }
		  if (al->columns[thispos][spec]=='T' || al->columns[thispos][spec]=='t') { postcounts[j][counter][3] += prob[pos][j]; }
		}
	      }
	    }
	  }
	  
	  if(probrev[pos][j] >= minposterior) {
	    
	    /****Print the bases at the positions in site if wanted****/
	    startpos = pos;
	    basecount = 1;
	    while(basecount < wms[j].len) {
	      --startpos;
	      if(al->columns[startpos][0] != '-')
		++basecount;
	    }
	    selec[0] =  scorerefwm(al,pos,wms[j].len,wms[j].matrev);
	    spec = 0;
	    
	    /***now print out the alignment****/
	    for(spec=0;spec<numspec;++spec) {
	      if(selec[spec] > selec_cutoff) {
		
		int counter = -1;
		for(thispos=pos;thispos>=startpos;--thispos) {
		  if(al->columns[thispos][0] != '-') {
		    counter++;
		    if (al->columns[thispos][spec]=='A' || al->columns[thispos][spec]=='a') { postcounts[j][counter][3] += probrev[pos][j]; }
		    if (al->columns[thispos][spec]=='C' || al->columns[thispos][spec]=='c') { postcounts[j][counter][2] += probrev[pos][j]; }
		    if (al->columns[thispos][spec]=='G' || al->columns[thispos][spec]=='g') { postcounts[j][counter][1] += probrev[pos][j]; }
		    if (al->columns[thispos][spec]=='T' || al->columns[thispos][spec]=='t') { postcounts[j][counter][0] += probrev[pos][j]; }

		    
		  }
		}
	      }
	    }
	  }
	}
      }
    }    
  } 
  /*printf("positions %d totsum %lf\n",refpos,totsum);*/
  free(selec);
  
  return;

}
  
  
double diffwms(double ***postcounts,wm *wms, int numwms) {
  
  // calc difference (old wm vs refined wm)

  int j,pos,i;
  const double PC = 0.5;
  /***** normalize postcounts *****/
  for (j=0; j<(numwms); ++j)
    for (pos=0; pos<wms[j].len;++pos) {
      double tot;
      tot=4*PC + postcounts[j][pos][0]+postcounts[j][pos][1]+postcounts[j][pos][2]+postcounts[j][pos][3];
      
      postcounts[j][pos][0] += PC;
      postcounts[j][pos][1] += PC;
      postcounts[j][pos][2] += PC;
      postcounts[j][pos][3] += PC;

      postcounts[j][pos][0]/=tot;
      postcounts[j][pos][1]/=tot;
      postcounts[j][pos][2]/=tot;
      postcounts[j][pos][3]/=tot;
    }
  
  double diff=0;
  
  for (j=0; j<(numwms); ++j)
    for (pos=0; pos<wms[j].len; ++pos)
      for(i=0; i<4; ++i)
	if (wms[j].mat[pos][i] > 0 || postcounts[j][pos][i]  > 0)
	  diff +=(wms[j].mat[pos][i]-postcounts[j][pos][i])*(wms[j].mat[pos][i]-postcounts[j][pos][i])/
	    (wms[j].mat[pos][i]+postcounts[j][pos][i])*(wms[j].mat[pos][i]+postcounts[j][pos][i]);
  
  printf("wmdiff %f\n", diff);
  
  return diff;


}



// enhancer predictions


int extractaln(alignment *input, alignment *output, int start, int len) {

  int i = 0;

  output->numspec = input->numspec;

  // get real start
  int pos = 0;
  int posrel = 0;

  while(posrel < start && pos < input->len) {
    if( (input->columns[pos][0] != '-') && (input->columns[pos][0] != 'N') ) {
      posrel++;
    }

    pos++;
  }

  start = pos;

  // check if you can fit a window of length len
  pos = start;
  posrel = 0;
  int reallen = 0;

  while(posrel < len && pos < input->len) {
    if( (input->columns[pos][0] != '-') && (input->columns[pos][0] != 'N') ) {
      posrel++;
    }

    pos++;
    reallen++;
  }

  if(posrel < (len-5)) { return -1; }
  output->len = reallen;
  
  // copy names
  for(i=0;i<numspec;i++) 
    output->specnames[i] = input->specnames[i];
  
  // copy sequences
  pos = start;
  posrel = 0;
  int pos0 = 0;
  while(posrel < len && pos < input->len) {

    for(i=0;i<numspec;i++) {
      output->columns[pos0][i] = input->columns[pos][i];
    }

    if( (input->columns[pos][0] != '-') && (input->columns[pos][0] != 'N') ) {
      posrel++;
    }
    pos++;
    pos0++;
  }

  return 0;
}


/***** HERE ARE ALL THE EXTENDED FUNCTIONS TO SUPPORT THE NTH ORDER MARKOV ORDER BACKGROUND MODEL *****/
/***** The original functions are still there. The name of the new functions is the old name followed by _nth. ******/

// extended background model

inline int base2int(char base) {

  if(base == 'A' || base == 'a')
    return 0; 
  if(base == 'C' || base == 'c')
    return 1;
  if(base == 'G' || base == 'g')
    return 2;
  if(base == 'T' || base == 't')
    return 3;

  return(-1);
}


inline char int2base(int number) {

  if(number==0)
    return 'A';
  if(number==1)
    return 'C';
  if(number==2)
    return 'G';
  if(number==3)
    return 'T';


  return(-1);
}


void getbgmodelfromseqs(alignment *al) {

  printf("Calculating background model (markov order = %i) ... \n",markovorder);
  
  const int pc = 1;

  // count combinations
  int i;
  int *counts = (int *) calloc( int(pow(4,markovorder+1)), sizeof(int) );
  for(i=0;i<(pow(4,markovorder+1));++i)
    counts[i]=0;

  // extract ref sequence and delete remove all gaps
  for(i=0;i<numgroups;++i) {

    int pos;
    int bases = 0;
    for(pos=0;pos<al[i].len;++pos) {
      if( al[i].columns[pos][0] != 'N' && al[i].columns[pos][0] != '-'  && al[i].columns[pos][0] != ' ') 
	bases++;
    }
    //printf("nbr of bases: %i\n",bases);

    char *refseq = (char *) calloc(bases, sizeof(char));
    bases=0;
    for(pos=0;pos<al[i].len;++pos) {
      if( al[i].columns[pos][0] != 'N' && al[i].columns[pos][0] != '-' && al[i].columns[pos][0] != ' ') 
	refseq[bases++] = al[i].columns[pos][0];
    }
    //printf("%s\n",refseq);
    
    bases--;
    for(pos=markovorder;pos<bases;++pos) {
      int index = 0;

      int pos2;
      for(pos2=(pos-markovorder);pos2<=pos;++pos2) 
	index += base2int(refseq[pos2]) * ( (int) pow(4,(pos2-pos+markovorder)) );

      if(index<0) {
	printf("Something is very wrong with the extended background model!\n");
	exit(-1);
      }

      counts[index]++;
    }

    free(refseq);
  }

  // allocate memory for background models
  bgfwd = (double *) calloc( int(pow(4,markovorder+1)), sizeof(double) );

  if(markovorder==0) {
    int k;
    for(k=0;k<4;++k)
      bgfwd[k] = bg[k];
    return;
  }

  // calculate conditional background frequencies: P(A|C) = #(A|C) / sum_i #(i|C)
  for(i=0;i<(pow(4,markovorder+1));++i) {
    int curtot = 0;
    int prev = (i % (int)pow(4, markovorder) );
    int j;
    for (j=0; j<4; ++j)
      curtot += counts[(int)pow(4,markovorder)*j+prev];

    bgfwd[i]= ( (double)(counts[i]+pc) / (double)(curtot+4*pc) );

    char *tmp = (char *) calloc(markovorder,sizeof(char));
    index2seq(i,tmp);
    if(markovorder < 4)
      printf("%s: %i = %lf\n",tmp,i,bgfwd[i]);
    //printf("%s: %i = %lf, counts=%i, countstot=%i,pc=%i\n",tmp,i,bgfwd[i],counts[i],curtot,pc);
    free(tmp);
  }

  printf("done\n");
  free(counts);
}


int shift_right(alignment *al, int pos, int len) {

  int abs = pos;
  int by = 0;
  int shift = 0;


  while( by<len ) {

    if(abs>=al->len) {
      printf("shift order gapless doesMOfit: something went wrong: %i\n",pos);
      return -1;
    }

    if( al->columns[abs][0] == '-' )
      ;
    else {
      by++; 
    }

    shift++;
    abs++;
  }
  
  shift--;
  //printf("\tby = %i\tabs = %i\n",shift,abs);
  return shift;
}


// "original" functions

void forwardbackward_nth(alignment *al,double **prob,double **probrev,treenode **nodelist,wm *wms,int restrictparses,double *otherwm,double *Z,int numwms, double *Y) {

  const int verbose = 0;
  double vbg=-99;
  double vfg=-99;
  double vF=-99;
  double vR=-99;
  
  int pos,k,j,refpos,relpos,abspos,spec,index,allselec;
  double fac, facrev,score,postsum;
  static double *F = NULL;
  if(F == NULL)
    F = (double *) calloc(maxlen,sizeof(double));
  static double *R = NULL;
  if(R == NULL)
    R = (double *) calloc(maxlen+1,sizeof(double));
 
  double bgprior = wms[numwms].prior;
  double otherwmprior = wms[numwms+1].prior;

  static int *selec = NULL;
  if(selec == NULL)
    selec = (int *) calloc(numnodes,sizeof(int));
  reducedtreenode *thistree;


  int **nth = (int **) calloc(al->len,sizeof(int *));
  for(j=0;j<al->len;++j) 
    nth[j] = (int *) calloc(numwms,sizeof(int ));



  /*
  int h;
  for(h=0;h<maxlen;++h) {
    for(j=0;j<numwms;++j) {
      prob[h][j] = 0.0;
      probrev[h][j] = 0.0;
    }

    Z[h] = 0.0;
    Y[h] = 0.0;
    F[h] = 1.0;
    R[h] = 1.0;
  }
  */


  /*
  # selection pattern the same for fwd and bwd?
    - problem at position 245 (fwd) and 224 (bwd)
    - selection not the same!
       - compare scorewm versus scorewmforward!
  # check one position:
      - is sumfg or sumbg the same for the leafs?
      - still the same after summing over the nodes?
  */

  double fs, fsrev;

  /***run over all positions***/
  /**at each position run over all wms***/
  /***check which species are background vs. foreground by scoring each species for the wm**/
  /***calculate score for wm by taking this 'selection' array****/
  refpos = 0;
  for(pos=0;pos<al->len;++pos)
    {
      /***forward***/
      if(al->columns[pos][0] != '-')
	{
	  ++refpos;
	  /***score and record column under background only****/
	  /***we store 1/probyability****/
	  F[pos] = bgprior;
	  

	  /**********************************************************************
	     - check if nth background model can be used
	     - store this information for each position and wm
     	     -> most consistent: use the same wm model for all wms at a given position!
	  **************************************************************/

	  int use_nth = 1;
	  
	  for(j=0;j<numwms;++j) {

	    nth[pos][j] = 1;

	    if(wms[j].len<=refpos) {

	      for(int s=0;s<numspec;++s) {

		if( gaplesswithref(al,pos,wms[j].len,s)==1 ) { // species gaplessly aligned -> will be included in our selection

		  if( doesMOfit(al,s,pos,wms[j].len)!=1 ) { // if for a gaplessly species the MO model does not fit, only use 0th order
		    use_nth = 0;
		    nth[pos][j] = 0;
		    }
		}
	      }
	    }
	    else { // wm does not fit
	      nth[pos][j] = 0;
	    }
	  }
	
	  
	  for(j=0;j<numwms;++j)
	    use_nth *= nth[pos][j];

	  for(j=0;j<numwms;++j)
	    nth[pos][j] = use_nth;
	  

	  
	  
	  //printf("MO_order_fwd = %i\trefpos = %i\tpos = %i\tnth[0] = %i\tnth[1] = %i\n",use_nth,refpos,pos,nth[pos][0],nth[pos][1]);
	  

	  /* end checks for nth order bg model */


	  /**wm models***/
	  for(j=0;j<numwms;++j)
	    {
	      //xxx = j;
	      /**check wm can fit and has nonzero prior****/
	      if(wms[j].len<=refpos && wms[j].prior > 0)
		{
		  /***if restrictparses is set we only put a site when the reference scores better than background***/
		  /***otherwise we consider all parses for the reference****/
		  if(restrictparses) {

		    if(use_nth)
		      score = scorerefwm_nth(al,pos,wms[j].len,wms[j].mat); 
		    else
		      score = scorerefwm(al,pos,wms[j].len,wms[j].mat); 
		    
		    if(score > selec_cutoff)
		      selec[0] = 1;
		    else
		      selec[0] = 0;
		  }
		  else {
		    selec[0] = 1;
		  }



		  /***only score when reference has selection***/
		  if(selec[0]> 0)
		    {
		      /**get selection for the other leafs****/
		      for(spec=1;spec<numspec;++spec)
			{
			  if(use_nth) {
			    //printf("scorewm_nth: spec = %i\tpos = %i\trefpos = %i\tvalue = %g\n",spec,pos,refpos,scorewm_nth(al,pos,wms[j].len,wms[j].mat,spec));
			    if(scorewm_nth(al,pos,wms[j].len,wms[j].mat,spec) > selec_cutoff)
			    //if(2.0> selec_cutoff)
			      selec[spec] = 1;
			    else
			      selec[spec] = 0;
			  }
			  else {
			    //			    printf("scorewm_0th: spec = %i\tpos = %i\trefpos = %i\tvalue = %g\n",spec,pos,refpos,scorewm(al,pos,wms[j].len,wms[j].mat,spec) );
			    if(scorewm(al,pos,wms[j].len,wms[j].mat,spec) > selec_cutoff)
			      selec[spec] = 1;
			    else
			      selec[spec] = 0;
			  }
			}
		      /**Get structure of the reduced tree that runs only on these species***/
		      thistree = tree_from_selec(selec);

		      /*
		      printf("selection_fwd:\trefpos = %i\t",refpos);
		      for(int x=0;x<numspec;++x) {
			printf("s[%i] = %i\t",x,selec[x]);
		      }
		      printf("\n");
		      */

		      /***get full score going backwards***/		      		     
		      relpos = wms[j].len-1;
		      if(use_nth) {
			fac = 0.5*wms[j].prior * ratiocolumn_nth(al,pos,al->columns[pos],wms[j].mat[relpos],thistree,refpos);
			//fs = ratiocolumn_nth(al,pos,al->columns[pos],wms[j].mat[relpos],thistree,refpos);
			//printf("ratiocolumn [nth]\tfwd\tmat\tpos = %i\trefpos = %i\tabspos = %i\tvalue = %g\t%c\twm = %i\n",pos,refpos,pos,ratiocolumn_nth(al,pos,al->columns[pos],wms[j].mat[relpos],thistree,refpos),al->columns[pos][0],xxx);
		      }
		      else {
			fac = 0.5*wms[j].prior * ratiocolumn(al->columns[pos],wms[j].mat[relpos],thistree);
			//fs = ratiocolumn(al->columns[pos],wms[j].mat[relpos],thistree);
			//printf("ratiocolumn [0th]\tfwd\tmat\tpos = %i\trefpos = %i\tabspos = %i\tvalue = %g\t%c\twm = %i\n",pos,refpos,pos,ratiocolumn(al->columns[pos],wms[j].mat[relpos],thistree),al->columns[pos][0],xxx);
		      }

		      --relpos;
		      abspos = pos-1;
		      while(relpos >= 0)
			{
			  /***only positions with base in reference****/
			  if(al->columns[abspos][0] != '-')
			    {
			      /**check that it is not negative***/
			      if(abspos < 0)
				{
				  fprintf(stderr,"HELP NEGATIVE abspos\n");
				}

			      if(use_nth) {
				fac *= F[abspos]*ratiocolumn_nth(al,abspos,al->columns[abspos],wms[j].mat[relpos],thistree,refpos);
				//fs *= ratiocolumn_nth(al,abspos,al->columns[abspos],wms[j].mat[relpos],thistree,refpos);
				//double v = ratiocolumn_nth(al,abspos,al->columns[abspos],wms[j].mat[relpos],thistree,refpos);
				//printf("ratiocolumn [nth]\tfwd\tmat\tpos = %i\trefpos = %i\tabspos = %i\tvalue = %g\t%c\twm = %i\n",pos,refpos,abspos,v,al->columns[abspos][0],xxx);
		      
			      }
			      else {
				fac *= F[abspos]*ratiocolumn(al->columns[abspos],wms[j].mat[relpos],thistree);
				//fs *= ratiocolumn(al->columns[abspos],wms[j].mat[relpos],thistree);
				//double v = ratiocolumn(al->columns[abspos],wms[j].mat[relpos],thistree);
				//printf("ratiocolumn [0th]\tfwd\tmat\tpos = %i\trefpos = %i\tabspos = %i\tvalue = %g\t%c\twm = %i\n",pos,refpos,abspos,v,al->columns[abspos][0],xxx);
			      }

			      --relpos;
			    }
			  --abspos;
			}
		      F[pos] += fac;

		      //printf("combinations\tfwd\tmat\tpos = %i\trefpos = %i\tvalue = %g\twm = %i\n",pos,refpos,fs,xxx);

		      prob[pos][j] = fac;
		    }
		  else
		    prob[pos][j] = 0;

		  /***selection for the reverse-strand*******/
		  if(restrictparses)
		    {
		      if(use_nth)
			score = scorerefwm_nth(al,pos,wms[j].len,wms[j].matrev);
		      else
			score = scorerefwm(al,pos,wms[j].len,wms[j].matrev);
		      
		      if(score > selec_cutoff)
			selec[0] = 1;
		      else
			selec[0] = 0;
		    }
		  else
		      selec[0] = 1;

		  if(selec[0] > 0)
		    {
		      /**get selection for the other leafs****/
		      for(spec=1;spec<numspec;++spec)
			{
			  if(use_nth) {
			    //printf("scorewm_nth: spec = %i\trefpos = %i\tvalue = %g\n",spec,refpos,scorewm_nth(al,pos,wms[j].len,wms[j].matrev,spec));
			    if(scorewm_nth(al,pos,wms[j].len,wms[j].matrev,spec) > selec_cutoff)
			      selec[spec] = 1;
			    else
			      selec[spec] = 0;
			  }
			  else {
			    //printf("scorewm_0th: spec = %i\trefpos = %i\tvalue = %g\n",spec,refpos,scorewm(al,pos,wms[j].len,wms[j].matrev,spec));
			    if(scorewm(al,pos,wms[j].len,wms[j].matrev,spec) > selec_cutoff)
			      selec[spec] = 1;
			    else
			      selec[spec] = 0;
			  }

			  
			}
		      /**Get structure of the reduced tree that runs only on these species***/
		      thistree = tree_from_selec(selec);


		      /***get full score going backwards***/
		      relpos = wms[j].len-1;
		      if(use_nth) {
			fac = 0.5*wms[j].prior * ratiocolumn_nth(al,pos,al->columns[pos],wms[j].matrev[relpos],thistree,refpos);
			if( verbose==1)
			  EXPratiocolumn_nth(al,pos,al->columns[pos],wms[j].matrev[relpos],thistree,&vfg,&vbg);
		      }
		      else
			fac = 0.5*wms[j].prior * ratiocolumn(al->columns[pos],wms[j].matrev[relpos],thistree);

		      if( verbose==1) {
			
			printf("forward, minus, one, pos=%i, wm=%i, use_nth=%i\n",pos,j,use_nth);
			printf("\t wm[%i].prior=%f\n",j,wms[j].prior);
			printf("\t fg=%f, bg=%f, fac=%f\n",vfg,vbg,fac);
		      }

		      --relpos;
		      abspos = pos-1;
		      

		      while(relpos >= 0)
			{
			  /***only positions with base in reference****/
			  if(al->columns[abspos][0] != '-')
			    {
			      /**check that it is not negative***/
			      if(abspos < 0)
				{
				  fprintf(stderr,"HELP NEGATIVE abspos\n");
				}
			      if(use_nth) {
				fac *= F[abspos]*ratiocolumn_nth(al,abspos,al->columns[abspos],wms[j].matrev[relpos],thistree,refpos);
				if( verbose==1)
				  EXPratiocolumn_nth(al,abspos,al->columns[abspos],wms[j].matrev[relpos],thistree,&vfg,&vbg);
			      }
			      else
				fac *= F[abspos]*ratiocolumn(al->columns[abspos],wms[j].matrev[relpos],thistree);

			      if( verbose==1) {
				
				printf("forward, minus, two, pos=%i, abspos=%i, wm=%i\n",pos,abspos,j);
				printf("\t fg=%f, bg=%f, fac=%f, F=%f\n",vfg,vbg,fac,F[abspos]);
			      }
			      
			      --relpos;
			    }
			  --abspos;
			}
		      F[pos] += fac;
		      
		      probrev[pos][j] = fac;
		    }
		  else
		    {
		      probrev[pos][j] = 0;
		    }
		}
	    }
	  /****Other WM****/
	  prob[pos][numwms+1] = 0;
	  probrev[pos][numwms+1] = 0; 
	  if(otherwmlen <= refpos && otherwmprior > 0)
	    {
	      selec[0] = 1;
	      /**check which species are ungapped aligned wrt reference****/
	      for(spec=1;spec<numspec;++spec)
		{
		  score = scorewm(al,pos,wms[numwms+1].len,wms[numwms+1].mat,spec);
		  if(score > selec_cutoff)
		    selec[spec] = 1;
		  else
		    selec[spec] = 0;
		}
	      /***get full score going backwards***/
	      relpos = otherwmlen-1;
	      /**get number that this column corresponds to****/
	      index = getindex(al->columns[pos],selec);
	      fac = otherwmprior * otherwm[index];
	      --relpos;
	      abspos = pos-1;
	      while(relpos >= 0)
		{
		  /***only positions with base in reference****/
		  if(al->columns[abspos][0] != '-')
		    {
		      index = getindex(al->columns[abspos],selec);
		      fac *= otherwm[index]*F[abspos];
		      --relpos;
		    }
		  --abspos;
		}
	      F[pos] += fac;
	      prob[pos][numwms+1] = fac;
	      
	    }
	  /***store 1/F[pos] ****/
	  Z[pos] = F[pos];
	  F[pos] = 1.0/F[pos];

	}
      /**column has gap in reference score ratio = 1***/
      else
	{
	  F[pos] = 1.0;
	  Z[pos] = 1.0;
	}
    }

  int goingdown = refpos;
  if(verbose==1) {
    double sum, partsum;
    for(int x=0;x<al->len;++x) {
      
      sum += log(F[x]);
      if(x>25 && x<870) 
	partsum += log(F[x]);
    }
    
    printf("forward, partsum, sum=%f [%f]\n",sum,partsum);
  }

  //printf("backwards calculation R\n");
  

  /***backwards calculation R*****/
  refpos = 0;
  for(pos=(al->len)-1;pos>=0;--pos)
    {
      if(al->columns[pos][0] != '-')
	{
	  ++refpos;
	  goingdown--;
	  /**first bg model**/
	  R[pos] = bgprior;

	  // check for each wm if you can use a nth order bg model
	  int use_nth = 1;
	  
	  int whichnth[numwms];

	  for(j=0;j<numwms;++j) {
	    
	    if(wms[j].len<=refpos) {
	      
	      // get for each position, position correspoding to the (left_2_right) direction
	      int shift = shift_right(al,pos,wms[j].len);

	      if( shift == -1 ) { // can't shift to the right, alignment too short
		use_nth = 0;
		break;
	      }

	      if( (pos+shift)<al->len ) {
		
		whichnth[j]= nth[pos+shift][j];
	      }
	    }
	  else
	    whichnth[j]= 0;	    
	  }


	  
	  //printf("MO_order_bwd = %i\trefpos = %i\tpos = %i\tnth[%i] = %i\tnth[%i] = %i\n",use_nth,goingdown,pos,0,whichnth[0],1,whichnth[1]);
	  
	  
	  /* end checks for nth order bg model */

	  
	  /**wm models***/
	  for(j=0;j<numwms;++j)
	    {
	      //xxx = j;
	      /**check wm can fit and has nonzero prior****/
	      if(wms[j].len<=refpos && wms[j].prior > 0)
		{
		  if(restrictparses)
		    {
		      if(whichnth[j]) {
			score = scorerefwmforward_nth(al,pos,wms[j].len,wms[j].mat);
		      }
		      else
			score = scorerefwmforward(al,pos,wms[j].len,wms[j].mat);

		      if(score > selec_cutoff)
			selec[0] = 1;
		      else 
			selec[0] = 0;
		    }
		  else
		    selec[0] = 1;

		  //printf("can wm fit: bwd: refpos = %i\tgoingdown = %i\tselec[0] = %i\n",refpos,goingdown,selec[0]);


		  if(selec[0]>0)
		    {
		      /**get selection for the other leafs****/
		      for(spec=1;spec<numspec;++spec)
			{
			  //printf("scorewmforward_nth: spec = %i\tpos = %i\trefpos = %i\tvalue = %g\n",spec,pos,goingdown,scorewmforward_nth(al,pos,wms[j].len,wms[j].mat,spec));
			  if(whichnth[j]) {

			    if(scorewmforward_nth(al,pos,wms[j].len,wms[j].mat,spec) > selec_cutoff)
			      selec[spec] = 1;
			    else
			      selec[spec] = 0;
			  }
			  else {
			    //printf("scorewmforward_0th: spec = %i\tpos = %i\trefpos = %i\tvalue = %g\n",spec,pos,goingdown,scorewmforward(al,pos,wms[j].len,wms[j].mat,spec));
			    if(scorewmforward(al,pos,wms[j].len,wms[j].mat,spec) > selec_cutoff)
			      selec[spec] = 1;
			    else
			      selec[spec] = 0;
			  }
			}
		     
		      /**Get structure of the reduced tree that runs only on these species***/
		      thistree = tree_from_selec(selec);
		    
		      /*
		      printf("selection_bwd:\trefpos = %i\t",goingdown);
                      for(int x=0;x<numspec;++x) {
                        printf("s[%i] = %i\t",x,selec[x]);
                      }
                      printf("\n");
		      */



		      /***get full score going forwards***/


		      relpos = 0;

		      if(whichnth[j]) {
			fac = 0.5 * wms[j].prior * ratiocolumnforward_nth(al,pos,al->columns[pos],wms[j].mat[relpos],thistree,goingdown);
			//fs = ratiocolumn_nth(al,pos,al->columns[pos],wms[j].mat[relpos],thistree,goingdown);
			//printf("ratiocolumn [nth]\tbwd\tmat\tpos = %i\trefpos = %i\tabspos = %i\tvalue = %g\t%c\twm = %i\n",pos,goingdown,pos,ratiocolumn_nth(al,pos,al->columns[pos],wms[j].mat[relpos],thistree,goingdown),al->columns[pos][0],xxx);
		      }
		      else {
			fac = 0.5 * wms[j].prior * ratiocolumn(al->columns[pos],wms[j].mat[relpos],thistree);
			//fs = ratiocolumn(al->columns[pos],wms[j].mat[relpos],thistree);
			//printf("ratiocolumn [0th]\tbwd\tmat\tpos = %i\trefpos = %i\tabspos = %i\tvalue = %g\t%c\twm = %i\n",pos,goingdown,pos,ratiocolumn(al->columns[pos],wms[j].mat[relpos],thistree),al->columns[pos][0],xxx);
		      }



		      ++relpos;
		      abspos = pos+1;
		      while(relpos < wms[j].len)
			{
			  /***only positions with base in reference****/
			  if(al->columns[abspos][0] != '-')
			    {
			      if(whichnth[j]) {
				fac *= ratiocolumnforward_nth(al,abspos,al->columns[abspos],wms[j].mat[relpos],thistree,goingdown)*R[abspos];
				//fs *= ratiocolumnforward_nth(al,abspos,al->columns[abspos],wms[j].mat[relpos],thistree,goingdown);
				//double v = ratiocolumnforward_nth(al,abspos,al->columns[abspos],wms[j].mat[relpos],thistree,goingdown);
				//printf("ratiocolumn [nth]\tbwd\tmat\tpos = %i\trefpos = %i\tabspos = %i\tvalue = %g\t%c\twm = %i\n",pos,goingdown,abspos,v,al->columns[abspos][0],xxx);
			      }
			      else {
				fac *= ratiocolumn(al->columns[abspos],wms[j].mat[relpos],thistree)*R[abspos];
				//fs *= ratiocolumn(al->columns[abspos],wms[j].mat[relpos],thistree);
				//double v = ratiocolumn(al->columns[abspos],wms[j].mat[relpos],thistree);
				//printf("ratiocolumn [0th]\tbwd\tmat\tpos = %i\trefpos = %i\tabspos = %i\tvalue = %g\t%c\twm = %i\n",pos,goingdown,abspos,v,al->columns[abspos][0],xxx);
				
			      }

			      ++relpos;
			    }
			  ++abspos;
			}
		      R[pos] += fac;

		      //printf("combinations\tbwd\tmat\tpos = %i\trefpos = %i\tvalue = %g\twm = %i\n",pos,goingdown,fs,xxx);
		      //printf("pos=%i, R=%f\n",pos,fac);
		    }
		  /***now reverse strand motif******/
		  if(restrictparses)
		    {
		      if(whichnth[j])
			score = scorerefwmforward_nth(al,pos,wms[j].len,wms[j].matrev);
		      else
			score = scorerefwmforward(al,pos,wms[j].len,wms[j].matrev);

		      if(score > selec_cutoff)
			selec[0] = 1;
		      else 
			selec[0] = 0;
		    }
		  else
		      selec[0] = 1;
		  
		  if(selec[0] > 0)
		    {
		      /**get selection for the other leafs****/
		      for(spec=1;spec<numspec;++spec)
			{
			  if(whichnth[j]) {
			    if(scorewmforward_nth(al,pos,wms[j].len,wms[j].matrev,spec) > selec_cutoff)
			      selec[spec] = 1;
			    else
			      selec[spec] = 0;
			  }
			  else {
			    if(scorewmforward(al,pos,wms[j].len,wms[j].matrev,spec) > selec_cutoff)
			      selec[spec] = 1;
			    else
			      selec[spec] = 0;
			  }
			}
		      /**Get structure of the reduced tree that runs only on these species***/
		      thistree = tree_from_selec(selec);


		      /***get full score going forwards***/
		      relpos = 0;
		      if(whichnth[j]) {
			facrev = 0.5*wms[j].prior * ratiocolumnforward_nth(al,pos,al->columns[pos],wms[j].matrev[relpos],thistree,goingdown);
			if( verbose==1)
			  EXPratiocolumnforward_nth(al,pos,al->columns[pos],wms[j].matrev[relpos],thistree,&vfg,&vbg);
		      }
		      else
			facrev = 0.5*wms[j].prior * ratiocolumn(al->columns[pos],wms[j].matrev[relpos],thistree);

		      if( verbose==1) {
			
			printf("backward, minus, one, pos=%i, wm=%i, use_nth=%i\n",pos,j,use_nth);
			printf("\t wm[%i].prior=%f\n",j,wms[j].prior);
			printf("\t fg=%f, bg=%f, fac=%f\n",vfg,vbg,facrev);
		      }

		      /***now all the other columns****/
		      ++relpos;
		      abspos = pos+1;
		      while(relpos < wms[j].len)
			{
			  /***only positions with base in reference****/
			  if(al->columns[abspos][0] != '-')
			    {
			      if(whichnth[j]) {
				facrev *= ratiocolumnforward_nth(al,abspos,al->columns[abspos],wms[j].matrev[relpos],thistree,goingdown)*R[abspos];
				if( verbose==1)
				  EXPratiocolumnforward_nth(al,abspos,al->columns[abspos],wms[j].matrev[relpos],thistree,&vfg,&vbg);
			      }
			      else 
				facrev *= ratiocolumn(al->columns[abspos],wms[j].matrev[relpos],thistree)*R[abspos];

			      if( verbose==1) {
			       printf("backward, minus, two, pos=%i, abspos=%i, wm=%i, use_nth=%i\n",pos,abspos,j,use_nth);
			       printf("\t fg=%f, bg=%f, fac=%f, R=%f\n",vfg,vbg,fac,R[abspos]);
			      }

			      ++relpos;
			    }
			  ++abspos;
			}
		      R[pos] += facrev;
		      
		    }
		}
	    }

	  /****Other WM****/
	  if(otherwmlen<=refpos && otherwmprior > 0)
	    {
	      selec[0] = 1;
	      /**check which species are ungapped aligned wrt reference****/
	      for(spec=1;spec<numspec;++spec)
		{
		  score = scorewmforward(al,pos,wms[numwms+1].len,wms[numwms+1].mat,spec);
		  if(score > selec_cutoff)
		    selec[spec] = 1;
		  else
		    selec[spec] = 0;
		}

	      /***get full score going forwards***/
	      relpos = 0;
	      /**get number that this column corresponds to****/
	      index = getindex(al->columns[pos],selec);
	      fac = otherwmprior * otherwm[index];
	      ++relpos;
	      abspos = pos+1;
	      while(relpos < otherwmlen)
		{
		  /***only positions with base in reference****/
		  if(al->columns[abspos][0] != '-')
		    {
		      index = getindex(al->columns[abspos],selec);
		      fac *= otherwm[index]*R[abspos];
		      ++relpos;
		    }
		  ++abspos;
		}
	      R[pos] += fac;
	    }
	  Y[pos] = R[pos];
	  R[pos] = 1.0/R[pos];
	}
      /**column has gap in reference score ratio = 1***/
      else
	{
	  Y[pos] = 1.0;
	  R[pos] = 1.0;
	}
    }


  
  if(verbose==1) {
    double sum, partsum;
    for(int x=0;x<al->len;++x) {
      
      sum += log(R[x]);
      if(x>25 && x<870) 
	partsum += log(R[x]);
    }
    
    printf("backward, partsum, sum=%f [%f]\n",sum,partsum);
  }


  int t;
  for(t=0;t<al->len;++t) 
    free(nth[t]);
  free(nth);
  
  
  /**calculate ratio of R/F from pos to end for each pos (needed for posteriors)***/
  R[al->len] = 1.0;
  for(pos=(al->len)-1;pos>=0;--pos)
    {
      R[pos] = F[pos]*R[pos+1]/R[pos];
    }
  /***finally get all the posteriors***/
  refpos = 0;
  for(pos=0;pos<al->len;++pos)
    {
      if(al->columns[pos][0] != '-')
	{

	  ++refpos;
	  /**wm models***/
	  for(j=0;j<numwms;++j)
	    {
	      /**check wm can fit and has nonzero prior****/
	      if(wms[j].len<=refpos && wms[j].prior > 0)
		{
		  score = R[pos+1]*prob[pos][j]*F[pos];
		  prob[pos][j] = score;
		  score = R[pos+1]*probrev[pos][j]*F[pos];
		  probrev[pos][j] = score;
		}
	    }
	  /***count posterior for the background model***/
	  score = R[pos+1]*bgprior*F[pos];
	  prob[pos][numwms] = score;
	  /***other wm****/
	  score = R[pos+1]*prob[pos][numwms+1]*F[pos];
	  prob[pos][numwms+1] = score;
	}
    }

  return;
}



int seq2index(char *seq) { 

  // actgg|a  ... the _last_ position is the reference - always!

  int j;
  int index=0;
    
  for(j=markovorder;j>=0;j--) {
    index += (int) pow(4,j) * base2int(seq[j]);
  }

  return index;
}


int seq2indexcond(char *seq) { 

  // actgg|a  ... the _last_ position is the reference - always!
  // *****  (only gives you the index of the conditional part)

  int j;
  int index=0;

  if((markovorder-1)<0) {

    printf("There's a problem in [seq2indexcond]\n");
    exit(-1);
  }
    
  for(j=(markovorder-1);j>=0;j--) {
    index += (int) pow(4,j) * base2int(seq[j]);
  }

  return index;
}


void index2seq(int index, char *tmp) {

  int j;
  for(j=(markovorder);j>=0;j--) {

    int rest = index % (int) pow(4,j);
    int coeff = (int) (index-rest) / (int) pow(4,j);
    tmp[j] = int2base(coeff);
    index = rest;
  }

}


double scorerefwm_nth(alignment *al,int pos,int len,double **mat)
{
  int relpos;
  char refletter;
  double score,fac,bfac;
  static double w[4];

  static char *seqret = NULL;
  if(seqret == NULL)
    seqret = (char *) calloc(markovorder+1,sizeof(char));

  /***run backwards over reference counting each base***/

  relpos = len-1;
  score = 1.0;
  while(relpos >= 0)
    {
      /***letter in reference species****/
      if(pos < 0)
	{
	  return 0;
	}
      refletter = al->columns[pos][0];
      if(refletter != '-')
	{
	  if( getprevseqs(al,0,pos,seqret)<0 ) {
	    //printf("Something is wrong. Couldn't fit bg model. [scorerefwm_nth] {pos=%i, sec=%i, rel=%i}\n",pos,0,relpos);
	    //exit(-1);
	  }

	  int index = seq2indexcond(seqret);

	  setleafprob(w,refletter);
	  fac = 0;
	  bfac = 0;
	  fac += mat[relpos][0]*w[0];
	  bfac += w[0]*bgfwd[index];
	  fac += mat[relpos][1]*w[1];
	  bfac += w[1]*bgfwd[index + (int)pow(4,markovorder)];
	  fac += mat[relpos][2]*w[2];
	  bfac += w[2]*bgfwd[index + 2*(int)pow(4,markovorder)];
	  fac += mat[relpos][3]*w[3];
	  bfac += w[3]*bgfwd[index + 3*(int)pow(4,markovorder)];
	  score *= (fac/bfac);
	  --relpos;
	}
      --pos;
    }
  //free(seqret);
  return score;
}

double scorewm_nth(alignment *al,int pos,int len, double **mat,int spec) {

  int relpos;
  char letter,refletter;
  double score,fac,bfac;
  static double w[4];

  static char *seqret = NULL;
  if( seqret==NULL )
    seqret = (char *) calloc(markovorder+1,sizeof(char));

  /***run backwards over reference counting each base***/

  int START = pos;
  
  relpos = len-1;
  score = 1.0;
  while(relpos >= 0)
    {
      /***If pos runs off the alignment on the left then we cannot fit window****/
      if(pos < 0)
	{
	  //free(seqret);
	  return 0;
	}
      /***letter in reference species****/
      refletter = al->columns[pos][0];
      letter = al->columns[pos][spec];

      

      if(refletter != '-')
	{
	  /**letter in reference aligned to gap -> return zero***/
	  if(letter == '-')
	    {
	      //free(seqret);
	      return 0;
	    }
	  else
	    {
	      if( getprevseqs(al,spec,pos,seqret)<0 ) {
		//printf("Something is wrong. Couldn't fit bg model. [scorewm_nth] {pos=%i, sec=%i, rel=%i, name=%s}\n",pos,spec,relpos,al->specnames[spec]);
		//free(seqret);
		return 0;
	      }
	      
	      int index = seq2indexcond(seqret);


	      setleafprob(w,letter);
	      fac = 0;
	      bfac = 0;
	      fac += mat[relpos][0]*w[0];
	      bfac += w[0]*bgfwd[index];
	      fac += mat[relpos][1]*w[1];
	      bfac += w[1]*bgfwd[index + (int)pow(4,markovorder)];
	      fac += mat[relpos][2]*w[2];
	      bfac += w[2]*bgfwd[index + 2*(int)pow(4,markovorder)];
	      fac += mat[relpos][3]*w[3];
	      bfac += w[3]*bgfwd[index + 3*(int)pow(4,markovorder)];
	      score *= (fac/bfac);
	    }
	  --relpos;
	}
      else
	{
	  /**gap in reference aligned to letter -> return zero ***/
	  if(letter != '-')
	    {
	      //free(seqret);
	      return 0;
	    }
	} 
      --pos;
    }
  //free(seqret);
  return score;
}

double scorewmforward_nth(alignment *al,int pos,int len, double **mat,int spec) {

  int relpos;
  char letter,refletter;
  double score,fac,bfac;
  static double w[4];

  static char *seqret = NULL;
  if( seqret==NULL )
    seqret = (char *) calloc(markovorder+1,sizeof(char));

  /***run backwards over reference counting each base***/
  relpos = 0;
  score = 1.0;
  while(relpos < len)
    {
      /***letter in reference species****/
      if(pos >= al->len)
	{
	  //free(seqret);
	  return 0;
	}
      refletter = al->columns[pos][0];
      letter = al->columns[pos][spec];

      if(refletter != '-')
	{
	  /**letter in reference aligned to gap -> return zero***/
	  if(letter == '-')
	    {
	      //free(seqret);
	      return 0;
	    }
	  else
	    {
	      if( getprevseqs(al,spec,pos,seqret)<0 ) {

		int pp = 1;
		while ( getprevseqs(al,spec,pos+pp,seqret)<0 ) {
		  pp++;
		}
		//printf("Something is wrong. Couldn't fit bg model. [scorerefwmforward_nth], {pos=%i, relpos=%i, len=%i, seq=%s}\n",pos,relpos,len,seqret);
		//return 0;
	      }

	      int index = seq2indexcond(seqret);
	      
	      setleafprob(w,letter);
	      fac = 0;
	      bfac = 0;
	      fac += mat[relpos][0]*w[0];
	      bfac += w[0]*bgfwd[index];
	      fac += mat[relpos][1]*w[1];
	      bfac += w[1]*bgfwd[index + (int)pow(4,markovorder)];
	      fac += mat[relpos][2]*w[2];
	      bfac += w[2]*bgfwd[index + 2*(int)pow(4,markovorder)];
	      fac += mat[relpos][3]*w[3];
	      bfac += w[3]*bgfwd[index + 3*(int)pow(4,markovorder)];
	      score *= (fac/bfac);
	    }
	  ++relpos;
	}
      else
	{
	  /**gap in reference aligned to letter -> return zero ***/
	  if(letter != '-')
	    {
	      //free(seqret);
	      return 0;
	    }
	} 
      ++pos;
    }
  //free(seqret);
  return score;
}





double scorerefwmforward_nth(alignment *al,int pos,int len,double **mat) {

  int relpos;
  char refletter;
  double score,fac,bfac;
  static double w[4];

  static char *seqret = NULL;
  if( seqret==NULL)
    seqret = (char *) calloc(markovorder+1,sizeof(char));

  /***run backwards over reference counting each base***/

  relpos = 0;
  score = 1.0;
  while(relpos < len)
    {
      /***letter in reference species****/
      if(pos >= -1+al->len)
	{
	  //free(seqret);
	  return 0;
	}
      refletter = al->columns[pos][0];
      if(refletter != '-')
	{
	  if( getprevseqs(al,0,pos,seqret)<0 ) {

	    int pp = 1;
	    while ( getprevseqs(al,0,pos+pp,seqret)<0 ) {
	      pp++;
	    }
	    //printf("Something is wrong. Couldn't fit bg model. [scorerefwmforward_nth], {pos=%i, relpos=%i, len=%i, seq=%s}\n",pos,relpos,len,seqret);
	    //exit(-1);
	  }

	  //flip(seqret,seqretflip,markovorder+1);
	  //int index = seq2indexcond(seqretflip);
	  int index = seq2indexcond(seqret);

	  setleafprob(w,refletter);
	  fac = 0;
	  bfac = 0;
	  fac += mat[relpos][0]*w[0];
	  bfac += w[0]*bgfwd[index];
	  fac += mat[relpos][1]*w[1];
	  bfac += w[1]*bgfwd[index + (int)pow(4,markovorder)];
	  fac += mat[relpos][2]*w[2];
	  bfac += w[2]*bgfwd[index + 2*(int)pow(4,markovorder)];
	  fac += mat[relpos][3]*w[3];
	  bfac += w[3]*bgfwd[index + 3*(int)pow(4,markovorder)];
	  score *= (fac/bfac);
	  
	  ++relpos;
	}
      ++pos;
    }

  //free(seqret);
  return score;
}


inline int getprevseqs(alignment *al, int spec, int pos, char *seqret) {

  int len = 0;
  int relpos = pos;
  
  while(len<=markovorder && relpos>=0) {
    
    if(al->columns[relpos][spec] != '-') {

      seqret[markovorder-len] = al->columns[relpos][spec];
      ++len;
    }
    --relpos;
  }


  if(len != (1+markovorder)) {
    //printf("relpos=%i, pos=%i, exit len=%i, name=%s\n",relpos,pos,len,al->specnames[spec]);
    return(-1);
  }
  return len;
}

inline void flip(char *in, char* out, int len) {

  int i;
  for(i=0;i<len;++i)
    out[i] = in[len-1-i];

}


inline int getprevseqs_fwd(alignment *al, int spec, int pos, char *seqret) {

  int len = 0;
  int relpos = pos;
  
  while(len<=markovorder && relpos<al->len) {
    
    if(al->columns[relpos][spec] != '-') {
      seqret[markovorder-len] = al->columns[relpos][spec];
      ++len;
    }
    ++relpos;
  }

  if(len != (1+markovorder)) {
    return(-1);
  }
  return len;
}



inline int doesMOfit(alignment *al, int spec, int pos, int wmlen) {

  int len = 0;
  int relpos = pos;

  relpos++;
  while(relpos>=0 && len < wmlen) {

    relpos--;
    if(al->columns[relpos][spec] != '-')
      len++;
  }

  if(relpos<=0)
    return 0;

 
  len = 0;
  while(len<=markovorder && relpos>=0) {
    
    if(al->columns[relpos][spec] != '-')
      ++len;
    
    --relpos;
  }
  
  if(len != (1+markovorder)) {
    return 0;
  }
  return 1;
}


inline int doesMOfit_fwd(alignment *al, int spec, int pos, int wmlen) {

  int len = 0;
  int relpos = pos;

  while(relpos<al->len && len < wmlen) {
    
    if(al->columns[relpos][spec] != '-')
      len++;
    relpos++;
    
  }
  
  if(relpos<0)
    return 0;
 
  len = 0;
  while(len<=markovorder && relpos<al->len) {
    

    if(al->columns[relpos][spec] != '-')
      ++len;
    
    ++relpos;
  }
  
  if(len != (1+markovorder)) {
    return 0;
  }
  return 1;
}



int gaplesswithref(alignment *al,int pos,int len,int spec) {

  int relpos;
  char letter,refletter;

  //  printf("in_gaplesswithref_fwd\tpos = %i\tspec = %i\n",pos,spec);

  /***run backwards over reference counting each base***/
  relpos = len-1;
  while(relpos >= 0) {
      /***If pos runs off the alignment on the left then we cannot fit window****/
      if(pos < 0)
	return 0;

      /***letter in reference species****/
      refletter = al->columns[pos][0];
      letter = al->columns[pos][spec];
      if(refletter != '-') {
	/**letter in reference aligned to gap -> return zero***/
	if(letter == '-')
	  {
	    return 0;
	  }
	
	--relpos;
      }
      else {
	/**gap in reference aligned to letter -> return zero ***/
	if(letter != '-')
	  
	  return 0;
	
      } 
      --pos;
  }
  return 1;
}

int gaplesswithref_fwd(alignment *al,int pos,int len,int spec) {

  int relpos;
  char letter,refletter;

  //  printf("in_gaplesswithref_bwd\tpos = %i\tspec = %i\n",pos,spec);

  /***run backwards over reference counting each base***/
  relpos = len-1;
  while(relpos >= 0) {
      /***If pos runs off the alignment on the left then we cannot fit window****/
      if(pos < 0)
	return 0;

      /***letter in reference species****/
      refletter = al->columns[pos][0];
      letter = al->columns[pos][spec];
      if(refletter != '-') {
	/**letter in reference aligned to gap -> return zero***/
	if(letter == '-')
	  {
	    return 0;
	  }
	
	--relpos;
      }
      else {
	/**gap in reference aligned to letter -> return zero ***/
	if(letter != '-')
	  
	  return 0;
	
      } 
      --pos;
  }
  return 1;
}



double ratiocolumn_nth(alignment *al,int pos, char *column,double *w,reducedtreenode *thistree, int refpos) {


  int j,i,cur,base,parbase,k;
  char letter;
  double qval,sumfg,sumbg;
  static double *matfg = NULL;
  if(matfg == NULL)
    matfg = (double *) calloc(4*numnodes,sizeof(double));
  static double *matbg = NULL;
  if(matbg == NULL)
    matbg = (double *) calloc(4*numnodes,sizeof(double));

  static char *seqret = NULL;
  if( seqret==NULL)
    seqret = (char *) calloc(markovorder+1,sizeof(char));


  int **kids; // to get for each node its leafs (the bases of the species)
  kids = (int **) calloc(numnodes,sizeof(int *)); // we know that the root has all species as leafs/children
  for(j=0;j<numnodes;++j)
    kids[j] = (int *)calloc(numspec,sizeof(int ));

  for(j=0;j<numnodes;++j)
    for(i=0;i<numspec;++i)
      kids[j][i]=0;

  node2kidsreduced(thistree,kids);

  int POSf = 155; // 

  //  printf("reduced_fwd: 0: %i, 1: %i\trefpos = %i\tabspos = %i\n",kids[0][0],kids[1][0],refpos,pos);
  /*
  i=0;
  printf("TREE_fwd:\trefpos = %i\tpos = %i\t",refpos,pos);
  while(thistree[i].index != -1) {

    printf("%i\t",thistree[i].index);
    ++i;
  }
  printf("\n");
  */

  /***set all probabilities to 1********/
  for(i=0;i<4*numnodes;++i)
    {
      matfg[i] = 1.0;
      matbg[i] = 1.0;
    }

  i = 0;
  /**for all nodes******/
  while(thistree[i].index != -1)
    {
      cur = thistree[i].index;
      base = 4*cur;
      /***if it is a leaf, put its probability***/
      if(cur<numspec) {
	letter = column[cur];
	setleafprob(&(matfg[base]),letter);
	setleafprob(&(matbg[base]),letter);

	
      }
      sumfg = 0;
      sumbg = 0;
      for(k=0;k<4;++k)
	{
	  sumfg += w[k]*matfg[base+k];

	  // sumbg += bg[k]*matbg[base+k]; // original
	  if(cur<numspec) { 

	    int spec = cur;
	    if( getprevseqs(al,spec,pos,seqret)<0 ) {
	      //printf("Something is wrong. Couldn't fit bg model. [ratiocolumn_nth,leafs] {pos=%i, spec=%i, name=%s}\n",pos,spec,al->specnames[spec]);
	      //exit(-1);
	    }
	      
	    int index = seq2indexcond(seqret);
	    sumbg += matbg[base+k] * bgfwd[index + k*(int)pow(4,markovorder)];
	    //if( refpos == POSf  && matbg[base+k]>0) 
	      {
		//printf("inratio_fwd_bgfwd_leafs\t%s\trefpos = %i\tabspos = %i\t fg = %g\tbg = %g\n",seqret,refpos,pos,sumfg,sumbg);
	    }

	  }
	  else { // for nodes: do average of cond bases of children
	    cur = thistree[i].index;
	    //printf("notkids_bwd\trefpos = %i\tabspos\n",refpos,pos);

	    double bgtmp = 0;
	    int counter = 0;
	    
	    for(int x=0;x<numspec;++x)
	      if(kids[cur][x]==1) {
		if( getprevseqs(al,x,pos,seqret)<0 ) {
		  //printf("Something is wrong. Couldn't fit bg model. [ratiocolumn_nth,inner] {pos=%i, spec=%i, name=%s}\n",pos,x,al->specnames[x]);
		  //exit(-1);
		}
		
		int index = seq2indexcond(seqret);

		

		bgtmp += bgfwd[index + k*(int)pow(4,markovorder)];
		counter++;
	      }
	    bgtmp /= counter;
	    
	    sumbg += bgtmp*matbg[base+k]; // original


	    //if( refpos == POSf  && matbg[base+k]>0) 
	      {
		//printf("more than reference_fwd\trefpos = %i\n",refpos);
		//printf("inratio_fwd_bgfwd_node\t%s\trefpos = %i\tabspos = %i\t fg = %g\tbg = %g\n","xx",refpos,pos,sumfg,sumbg);
	    }

	  }
	}
      
      /***not the root*****/
      if(thistree[i].parent >=0)
	{
	  qval = thistree[i].q;
	  parbase = 4*(thistree[i].parent);
	  for(k=0;k<4;++k)
	    {
	      /**multiply contribution from this branch***/
	      matfg[parbase+k] *= (sumfg*(1.0-qval)+qval*matfg[base+k]);
	      matbg[parbase+k] *= (sumbg*(1.0-qval)+qval*matbg[base+k]);
	    }
	  /**next node in the tree***/
	  ++i;
	}
      /***the root of the tree****/
      else
	{
	  for(j=0;j<numnodes;++j)
	    free(kids[j]);
	  free(kids);
	  
	  
	  
	  //free(seqret);

	  //if( refpos == POSf ) 
	  {
	    
	    //printf("inratio_fwd_END\t%c\trefpos = %i\tabspos = %i\t fg = %g\tbg = %g\n",column[0],refpos,pos,sumfg,sumbg);
	  }

	  return (sumfg/sumbg);
	  

	}
    }
  fprintf(stderr,"ERROR, in ratiocolumn function we ran off the tree without encountering a root\n");
  return 1;
}



double ratiocolumnforward_nth(alignment *al,int pos, char *column,double *w,reducedtreenode *thistree, int refpos) {

  int j,i,cur,base,parbase,k;
  char letter;
  double qval,sumfg,sumbg;
  static double *matfg = NULL;
  if(matfg == NULL)
    matfg = (double *) calloc(4*numnodes,sizeof(double));
  static double *matbg = NULL;
  if(matbg == NULL)
    matbg = (double *) calloc(4*numnodes,sizeof(double));

  static char *seqret = NULL;
  if( seqret==NULL)
    seqret = (char *) calloc(markovorder+1,sizeof(char));


  int **kids; // to get for each node its leafs (the bases of the species)
  kids = (int **) calloc(numnodes,sizeof(int *)); // we know that the root has all species as leafs/children
  for(j=0;j<numnodes;++j)
    kids[j] = (int *)calloc(numspec,sizeof(int ));

  for(j=0;j<numnodes;++j)
    for(i=0;i<numspec;++i)
      kids[j][i]=0;

  node2kidsreduced(thistree,kids);

  int POSb = 134;

  //printf("reduced_bwd: 0: %i, 1: %i\trefpos = %i\tabspos = %i\n",kids[0][0],kids[1][0],refpos,pos);
  /*
  i=0;
  printf("TREE_bwd:\trefpos = %i\tpos = %i\t",refpos,pos);
  while(thistree[i].index != -1) {

    printf("%i\t",thistree[i].index);
    ++i;
  }
  printf("\n");
  */

  /***set all probabilities to 1********/
  for(i=0;i<4*numnodes;++i)
    {
      matfg[i] = 1.0;
      matbg[i] = 1.0;
    }

  i = 0;
  /**for all nodes******/
  while(thistree[i].index != -1)
    {
      cur = thistree[i].index;
      base = 4*cur;
      /***if it is a leaf, put its probability***/
      if(cur<numspec)
	{
	  letter = column[cur];
	  setleafprob(&(matfg[base]),letter);
	  setleafprob(&(matbg[base]),letter);
	}
      sumfg = 0;
      sumbg = 0;
      for(k=0;k<4;++k)
	{
	  sumfg += w[k]*matfg[base+k];

	  if(cur<numspec) { 

	    int spec = cur;
	    if( getprevseqs(al,spec,pos,seqret)<0 ) {
	      //printf("Something is wrong. Couldn't fit bg model. [ratiocolumnforward_nth,leafs] {pos=%i, spec=%i, name=%s}\n",pos,spec,al->specnames[spec]);

	      int pp = 1;
	      while ( getprevseqs(al,spec,pos+pp,seqret)<0 ) {
		pp++;
	      }
	    }
	
	    int index = seq2indexcond(seqret);
	    sumbg += matbg[base+k] * bgfwd[index + k*(int)pow(4,markovorder)];

	    //if( refpos == POSb  && matbg[base+k]>0 ) 
	      {
		//printf("inratio_bwd_bgfwd_leafs\t%s\trefpos = %i\tabspos = %i\t fg = %g\tbg = %g\n",seqret,refpos,pos,sumfg,sumbg);
	    }

	  }
	  else { // for nodes: do average of cond bases of children
	    cur = thistree[i].index;
	    
	    //printf("notkids_bwd\trefpos = %i\tabspos = %i\t",refpos,pos);
	    /*
	    int o=0;
	    printf("TREE_bwd:\trefpos = %i\tpos = %i\t",refpos,pos);
	    while(thistree[o].index != -1) {
	      
	      printf("%i\t",thistree[o].index);
	      ++o;
	    }
	    printf("\n");
	    */

	    double bgtmp = 0;
	    int counter = 0;
	    
	    for(int x=0;x<numspec;++x)
	      if(kids[cur][x]==1) {
		if( getprevseqs(al,x,pos,seqret)<0 ) {
		  //printf("Something is wrong. Couldn't fit bg model. [ratiocolumnforward_nth,inner] {pos=%i, spec=%i, name=%s}\n",pos,x,al->specnames[x]);

		  int pp = 1;
		  while ( getprevseqs(al,x,pos+pp,seqret)<0 ) {
		    pp++;
		  }
		}

		
		int index = seq2indexcond(seqret);

		

		bgtmp += bgfwd[index + k*(int)pow(4,markovorder)];
		counter++;
	      }
	    bgtmp /= counter;

	    sumbg += bgtmp*matbg[base+k]; // original

	    //if( refpos == POSb  && matbg[base+k]>0 ) 
	      {

		//printf("inratio_bwd_bgfwd_node\t%s\trefpos = %i\tabspos = %i\t fg = %g\tbg = %g\n","xx",refpos,pos,sumfg,sumbg);
	    }

	  }
	}
      
      /***not the root*****/
      if(thistree[i].parent >=0)
	{
	  qval = thistree[i].q;
	  parbase = 4*(thistree[i].parent);
	  for(k=0;k<4;++k)
	    {
	      /**multiply contribution from this branch***/
	      matfg[parbase+k] *= (sumfg*(1.0-qval)+qval*matfg[base+k]);
	      matbg[parbase+k] *= (sumbg*(1.0-qval)+qval*matbg[base+k]);
	    }
	  /**next node in the tree***/
	  ++i;
	}
      /***the root of the tree****/
      else
	{
	  for(j=0;j<numnodes;++j)
	    free(kids[j]);
	  free(kids);

	  //free(seqret);

	  //	  printf("\t fg=%f\n",sumfg);
	  //printf("\t bg=%f\n",sumbg);

	  //if( refpos == POSb ) 
	  {
	    //printf("more than reference_bwd\trefpos = %i\n",refpos);
	    //printf("inratio_bwd_END\t%c\trefpos = %i\tabspos = %i\t fg = %g\tbg = %g\n",column[0],refpos,pos,sumfg,sumbg);
	  }
	  return (sumfg/sumbg);
	}
    }
  fprintf(stderr,"ERROR, in ratiocolumn function we ran off the tree without encountering a root\n");
  return 1;
}



void node2kidsreduced(reducedtreenode *thistree, int **kids) {

  int j, i;

  i=0;
  while(thistree[i].index != -1) {

    int cur = thistree[i].index;

    if(cur<numspec) {
	kids[cur][cur]=1;
    }
    
    if(thistree[i].parent >=0) {

      for(int k=0;k<numspec;++k)
	if(kids[cur][k]==1)
	kids[thistree[i].parent][k]=1;
    }

    ++i;
  }

  /*  
  printf("node=\t0\t1\t2\t3\t4\t5\t6\n");
  for(j=0;j<numnodes;++j) {
    printf("node=%i\t",j);
    for(i=0;i<numspec;++i)
      printf("%i\t",kids[j][i]);
    printf("\n");
  }
   printf("\n\n");
  */
  
}


inline double myround(double nbr, int digits) {

  return (pow(0.1,digits) * round( nbr * pow(10,digits) ) );
}











int n2ng(alignment *input, int start,int step) {


  int pos = start;
  int posrel = start;
  
  while(posrel < (step+start) && pos < input->len) {
    if( (input->columns[pos][0] != '-') && (input->columns[pos][0] != 'N') ) {
      posrel++;
    }
    pos++;
  }

  printf("%i %i\n",pos,posrel);

  return pos;

}




void read_cond_bg_probs(double *bgmok_tmp, int mo, char *filename) {

  printf("Reading background %s\n",filename);


  ifstream bgfile;
  bgfile.open(filename, ios::out);
  if (bgfile.fail() == true) {

    fprintf(stderr,"Could not read background file. Exiting.\n");
    exit(-1);
  }


  bgfwd = (double *) calloc( int(pow(4,markovorder+1)), sizeof(double) );
  
  
  
  while (!bgfile.eof() ) {

    string line;
    getline(bgfile, line);
    
    int pos = line.find('\t');
    if(pos < 0) { continue; }

    string id, value;

    id=line.substr(0, pos);
    value=line.substr(pos+1, line.length()-pos-1);

    char *seq =  new char[id.size()+1];
    strcpy(seq,id.c_str());
    int index = seq2index(seq);

    cout << "id = " << id << "\t" << "value = " << value << "\t" << "index = " << index << endl;


    double condprob = atof(value.c_str());

    //cout << "id: " << id << "  value: " << value  "  index: " << endl; //index << "  condprob " << condprob << endl;                                                                             
    //cout << "bgmok_tmp" << " " << mo << " " << index << " " << condprob << endl;                                                                                                                 
        //    int index = atoi(id.c_str());
    //    double condprob = atof(value.c_str());
    bgfwd[index]=condprob;
  }
  //exit(-1);                                                                                                                                                                                            
  

  bgfile.close();
}




void check_sequence(alignment *al) {

  
  for(int spec=1;spec<numspec;++spec) {

    int nbr_notgaps = 0;

    for(int pos=0;pos<al->len;++pos) {

      if( al->columns[pos][spec] != '-' ) { nbr_notgaps++; }
    }

    if( nbr_notgaps <= markovorder ) {

      for(int pos=0;pos<al->len;++pos) 
	al->columns[pos][spec] = '-';
    }
  }
}
