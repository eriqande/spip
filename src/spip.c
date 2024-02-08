/*
 *  spip_m2.c
 *  PRwP
 *
 *  Created by Eric Anderson on 6/6/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

/*
 
 
 
 SPIP SOURCE CODE WRITTEN BY ERIC C ANDERSON
 
 eric.anderson@noaa.gov
 
 Public domain source code developed by Federal employee
 
 
 
 */

#define DEFINED 1
#define UN_EXTERN

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "MathStatRand.h"
#include "ECA_MemAlloc.h"
#include "ranlib.h"
#include "MCTypesEtc.h"
#include "ECA_Opt2.h"


#define REALLOC_BLOCK 5

/* SOME ENUM'S */
enum dsn_type {
	NEG_BINOM,
	POISSON,
	BINARY
};

enum sex_enum {
	MALE,
	FEMALE
};

enum samtype { /* an enum to say when an adult was sampled */
	BEF_KILL, /* before the episode of killing */
	AFT_KILL, /* after the episode of killing */
	WHILE_REPRO  /* while it was in the act of reproducing */
};


/* SOME TYPEDEFS */

typedef struct {  /* a structure to hold all the population-specific parameters */
	int MaxAge;  /* the oldest that an individual can be and still reproduce */
	
	double *femsurv;  /* survival rates between age classes for females */
	
	double *malesurv;  /* survival rates between age classes for males */
	
	double *Fasrf;  /* Female age-specific relative fecundities */
	
	double *Masrp;  /* Male age-specific relative reproductive prowess */
	
	double *FPR;  /* Female probability of reproducing at all at a given age */
	double *MPR;  /* Male probability of reproducing at all at a given age */
	/* the two above are useful for salmonid-like life histories */
	
	/* these two take the place of the SemelPar parameter that I had before
	 because I think that it might be nice to have sex-specific semelparity 
	 and also because the probability of dying after reproduction might be different
	 for different ages.  Notice that this death is separate from the death that 
	 occurs from the failure to survive */
	double *FPDPR;  /* Female prob of death post-reproduction */
	double *MPDPR; /* Male prob of death post-reproduction */
	
	double Fcv;    /* Female coefficient of variation in offspring number */
	
	int ReproInhib;  /* age of a living offspring that still inhibits reproduction.  Only for Binary */
	
	double MCVrs;  /* male variance in reproductive success */
	
	double sticky_m_alpha;
	
	double MateSpec;  /* the mate specificity parameter */
	
	/*	int YearlyVarSwitch;  whether the male reproductive rates are "sticky" across years
	 or not */
	
	double SexRatio;  /* prob that a newborn is a male */
	
	enum dsn_type OffsDsn;
	
	int CohortSizesRandom;  /* this is a flag that says whether the number of offspring in each cohort is fixed at a certain
	 value, or is allowed to be the random sum of the number of offspring each individual makes.  If
	 fixed gCohortSizesRandom==0, and if not it is 1.  If fixed, we have fewer options for the offspring number
	 distribution --- namely it can be poisson of negative binomial (because the conditional branching
	 process model can be dealt with in those cases).  In other words, gCohortSizesRandom==0 means fixed cohort
	 size and ==1 means random cohort size. */
	
	int DiscardBefore;  /* year before which any pedigree information is not traced back through */
	int TossAll;        /* if this is 1, then all individuals born between DiscardBefore and DiscardBefore+
	 MaxAge-1 inclusive will have NULL parents and their PEDIGREE will say 0 for both
	 parents.  If it is not 1, then those same individuals will only have null parents
	 if they are part of an initializing cohort, or if their parents were born before 
	 DiscardBefore.  */
	
	/* things having to do with migration */
	double ***MigPr;  /* array of probs that an individual from a population will emigrate out of it.  subscripted by time sex and age (from 1 to MA) */		
	
	double ****Dest;  /* the destination probability matrix, subscripted by time, sex, age (1 to MA) and  
	 population of destination.  For example, if the matrix for 1-year-olds looked like:
	 
	 0		.1		.9
	 .2		0		.8
	 .5		.5		0
	 
	 It means that 10% of the migrating 1-year-olds in population 0 will go to population 1, and 90% to pop 2.
	 And so forth.  I will constrain things in the input so that Dest[k][k] = 0 always.  (i.e. migrants don't have
	 any chance of actually showing up where they started in the year they migrate.    
	 
	 */		
	
	
	double *SneakyMProb;  /* array of probabilities.  Each of these is the probability that an individual from a different locale
	 might sneak a mating with a female in this locale.  */
	
	/* here, if I want to coerce spawning in groups of SpawnGroupSize females and SpawnGroupSize males, I set this variable.
	 By default it should be set to -1.  This should be useful for simulating full-factorial mating in groups of size SpawnGroupSize
	 males and females.  */
	int SpawnGroupSize;
	
	/*  To take all the spawners in a single year and divvy them into SpawnGroupCnt groups of males and SpawnGroupCnt groups of females that are somewhat
	 randomly sized (drawn from a multinomial dsn) and are only allowed to mate amongst members of the same spawner group, then this is option to use.  By default 
	 it should be set to -1.  This should be useful for paritioning a year's spawners into SpawnGroupCnt day buckets.  */
	int SpawnGroupCnt; 
	
	
	/* When memory gets allocated to the pointers for the next cohort, this is the factor by which the number of spaces is 
	 greater than the expected number.  (basically it determines mspaces and fspaces in MakeBabies).  By default it will be 
	 set to 10.  But should you be doing lots of generations and not wanting to eat up so much memory, you can set it smaller 
	 (like 2 for the simulations Oystein was doing).  Note that ultimately I should just realloc all those things after the 
	 offspring are created (right at the end of MakeBabies() to reclaim that space. */
	int OffspringAllocOverload;
	
	
} PopPars;

typedef struct { /* a struct to hold data and variables related to sampling of individuals
 from the population.  */
	int NumLoc;  /* number of indepedently-segregating loci to sample */
	
	int *NumAlleles;
	double **AlleFreqs;
	
	/* the following pertain to the --newborn-samp option for sampling only the newborns */
	int SampleNewbornsOnly;
	int *NumF; /* the number of females from each year to randomly sample */
	int *NumM; /* the number of males from each year to randomly sample */  
	int nc;  /* number of different cohorts to sample from */
	int *ts;  /* the years to sample from */
	
	
	
	/* the following pertain to the --gtyps-for-all option */
	int SampleAllFromCohorts;
	int SA_lo;  /* low and high values of years to sample from */
	int SA_hi;
	
	
	/* the following pertain to the --gtyp-a-proportion type options */
	double *pre_male_samp;
	double *post_male_samp;
	double *pre_fem_samp;
	double *post_fem_samp;
	
	double dur_repro_fem_samp;
	double dur_repro_male_samp;
	
	/*  and here are the arrays that give 1's for the times when the sampling is to be 
	 done and 0's otherwise */
	int *pre_male_samp_years;
	int *pre_fem_samp_years;
	int *post_male_samp_years;
	int *post_fem_samp_years;
	int *dur_male_samp_years;
	int *dur_fem_samp_years;
	
	double LethalityProb;
	
	
} SampPars;

typedef struct {
	char CohortSizeMode[1000]; /* can be either "const" or "var" , for the time being */
	int *CohortSizes; /* if using the var option, this array holds the pop sizes */
	int ConstSize;  /* if using the const option, this is our cohort size */
} PopSizePars;

typedef struct indiv_ptr {  /* a structure to hold info about individuals */
	
	char ID[50];
	struct indiv_ptr *ma;
	struct indiv_ptr *pa;
	struct indiv_ptr *OldestOffsp;
	int dead;
	int t;  /* time at which the individual was born */
	int idx;  /* the index of the individual (they are subscripted starting from 0) */
	
	int birth_pop;  /* the index of the population in which it was born.  It gets set in the FillIndiv function.  */
	
	double ron;  /* relative offspring number */
	int *MaFounder;  /* an array of the indices of the founders of the maternal gene at each sampled locus (i.e. this is subscripted by locus) */
	int *PaFounder;  /* same as above, but for the paternal gene */
	
	int N_PreK;  /* number of times this individual was sampled in a pre-killing sampling episode */
	int N_PostK;  /* number of times this individual was sampled in a post-killing sampling episode */
	int N_DurRepro;  /* number of times this individual was sampled during reproduction */
	
	int PreSpace;  /* the amount of space currently allocated to the PreK_ts or the PostK_ts of the DurRepro_ts */
	int PostSpace;
	int DurReproSpace;
	
	int *PreK_ts; /* the years in which the individual was sampled in a pre-killing sampling episode */
	int *PostK_ts; /* the years in which the individual was sampled in a pre-killing sampling episode */
	int *DurRepro_ts; /* the years in which the individual was sampled while in the act of reproduction */
	
	int **SampLocales;
	
	int ReportMe;  /* this takes the value 1 if this indiv struct is the instance of an individual which should be
	 used for reporting genotype data */
	
	int IsReproducing;
	
	double sticky_birth_rate;  /* the gamma rv the individual will get at birth (of it's a male) to implement 
								variance in repro sucess that is maintained across different years by the indiv. */
	
} indiv;


/* this is a structure that holds all the information relevant to a single population
 (or locale) */
typedef struct {
	PopPars P;
	indiv **Males, **Females;
	int *NM, *NF;  /* for initial numbers of males and females */
	int *NMR, *NFR;  /* for the number of remaining individuals after some of them die...*/
	int MA;
	int *MSpace,*FSpace;  /* amount of memory allocated to the array of males and females each year */
	
	SampPars S;
	PopSizePars Pop;
	
	int T; /* the number of years to run on this one */
	
	int HasLocFile;
	char LocFileName[10000];
} Locale;


/* FUNCTION PROTOTYPES */
int GetALocale(Locale **L, int which, int *argc, char **argv[]);
void PrintMigrantSummary(int ****M, int ***P, int NP, int MA);
void AYearInTheLife(Locale *L, int t, int bp, Locale **LA, int NumPops);
void CopyIndiv(indiv *From, indiv *To);
void Migrate_O_Rama(Locale **LA, int NumPops, int t);
double ***AllocAndInitMigPr(int T, int MA);
double ****AllocAndInitMigDest(int T, int MA, int NumPops, int CurrentPop);
int PopSizeNonZero(int t, int MA, int *NFR, int *NMR);
void PreDeathSample(int t, int MA, SampPars S, int *NM, int *NF, indiv **Males, indiv **Females, int locale, int *NMR, int *NFR );
void PostDeathSample(int t, int MA, SampPars S, int *NM, int *NF, indiv **Males, indiv **Females, int locale, int *NMR, int *NFR );
void StandardAnnualDieOff(int t, int MA, int *NM, int *NF, indiv **Males, indiv **Females, int *NMR, int *NFR, PopPars P);
void FillIndiv(indiv *baby, indiv *mom, indiv *dad, int *counter, int year, int num, char sex, int birthplace, PopPars P);
void PrintAndDiscardPedigree(indiv *self, int year, int DiscardBefore, int TossAll, int MaxAge);
void InitCohorts(int MA, int *NF, int *NM, int T, indiv ***M, indiv ***F, PopPars P, int birthplace, int *MSpace, int *FSpace);
void KillOff(indiv *I, int N, int *NR, double s, int t);
void MakeBabies(indiv **F, indiv **M, int t, int MA, PopPars P, int N, int *NFR, int *NF, int *NM, int *NMR, SampPars S, int birthplace, int *MSpace, int *FSpace, Locale **LA, int NumPops);
int NumOffspringRV(double mu, double cv, enum dsn_type Type);
void SibshipPapas(indiv **Dads, indiv **PossibleDads, int n, int N, double *Wts, double S, double r);
void AncestryOfSampledGenes(int **indics, int n, indiv *ind, int *nf,  int *ma, int *pa, int ReInit, int *AncsIdx, int *AncsGenes, int tf, int **FounderLocs, int *FLM);
void OldAncestryOfSampledGenes(int **indics, int n, indiv *ind, int *nf,  int *ma, int *pa, int ReInit);
void PrintProbDistOfPapas(double r, int maxn);
void ReadLocFile(SampPars *S, char *FileName);

int ProcessOptions(int *output_argc, char **output_argv[], PopPars *P,int *T, int **NM, int **NF, 
				   int **NMR, int **NFR, SampPars *S, char *locfilename, PopSizePars *Pop, int *N, int **MSpace, int **FSpace, int *HasLocFile);

void MarkAsSampled(indiv *I, int N, double prob, enum samtype IsPostKilling, int t, int locale, double LethalityProb, int *NR);
void SampleMarkedAtLocusl(/* indiv **M, indiv **F */ int NumPops, Locale **LA, /*int *NM, int *NF,*/ int LastT, /* int no_ibd_t, */ int l, int L, int *nf, int **fl, int *flm, int n);
void PrintAllelesOfSampledInds(indiv **M, indiv **F, int *NM, int *NF, int L, int start, int stop, int **FounderTypes, int **FounderLocations);
int **DrawFounderTypes(int n, int *na, Locale **LA, int *nf, int **FL);
void MarkIndAsSampled(indiv *I, enum samtype IsPostKilling, int t, int locale);
void PrintCensusSizes(int t, int MA, int *NFR, int *NMR, const char *string, int age_start, int whichpop);
void MarkSampledNewborns(SampPars S, int *NM, int *NF, indiv **Males, indiv **Females, int whichpop);
int ReturnCohortSize(PopSizePars *Pop,int t);
char *dsn_type2string(enum dsn_type d);
void PrintPars(PopPars P, PopSizePars Pop, SampPars S, int T);
int DispersionParError(double d, char *whichoption);
void PrintBillHillPars(PopPars P, PopSizePars Pop);
void NEW_SampleMarkedAtLocusl(/* indiv **M, indiv **F */ int NumPops, Locale **LA, /*int *NM, int *NF,*/ int LastT, /* int no_ibd_t, */ int l, int L, int *nf, int **fl, int *flm, int n);
void MultiPop_InitCohorts(Locale **LA, int NumPops, int T);

/* SOME GLOBAL VARIABLES */
int gID;  /* for keeping track of individual ID numbers */
int gRunTests = 0; /* for a verbosity setting that spews extra information that can be used to 
 test the variance in family size */
int gNumPops;
int gNumPops_Flag;
int gCurrentPopNum = 0;   /* sort of a hack.  This is here to let us know which population we are processing on the commmand line */
int ***gPotMig, ****gMig;  /* for counting overall the number of migrants and "potential migrants"  gPotMig is subscripted by [MALE/FEMALE][POP][AGE] 
 and gMig is subscripted by [MALE/FEMALE][POP][AGE][DESTINATION] */


int main(int argc, char **argv) 
{
	int i,l,t;
	int T=0;  /* the number of new cohorts to produce */
	int more_pops=1;
	
	/* for tracing genes once the pedigree is established */
	int *NumFounding;  /* number of founding genes of the sample, subscripted by locale and then by locus */
	int **FounderTypes, **FounderLocations=NULL, *FoundLocMem=NULL;  /* FounderLocations is for keeping track of which locale each founder is from, and FoundLocMem records how much
	 space has been allocated to each locus' FounderLocations array */
	FILE *sfile;
	long seed1, seed2;
	
	int bp, NumPops = 0;
	Locale **LocaleArray = (Locale **)ECA_CALLOC(1000,sizeof(Locale *));
	Locale *L;
	int which = 0;
	int CheckMA;
	int WeAreGenotyping = 0;
	int OverallNumLoc, *OverallNumAlle;
	
	/* initialize the counter for indiv ID's  and other global variables */
	gID = 0;
	gNumPops = 0;
	gNumPops_Flag = 0;
	
	/* work through the command line, collecting information about each population,
	 and stop when there are no more --new-pop options.  Note that the argc and the 
	 argv get modified each time there is a --new-pop. */
	for(which=0;more_pops;which++) {
		more_pops = 0;
		NumPops++;  /* each time we go through here we get another population */
		more_pops = GetALocale(LocaleArray,which,&argc, &argv);
		
		/* increment a global counter that tells us what population we are working on
		 inside of GetALocale (a bit of a hack, but it works)  */
		gCurrentPopNum++;
		printf("DONE_COLLECTING_INFORMATION_FOR_POPULATION : %d\n",gCurrentPopNum-1);
		
		/* get the LocusData if necessary */
		L = LocaleArray[which];
		/*if(L->S.nc > 0 || L->S.pre_male_samp != NULL || L->S.post_male_samp != NULL || L->S.pre_fem_samp != NULL || L->S.post_fem_samp != NULL ) { */
		if(L->HasLocFile==1) {
			ReadLocFile(&(L->S),L->LocFileName);
		}
	}
	
	if(NumPops > 1  && gNumPops != NumPops)  {
		fprintf(stderr,"Error! Argument of --num-pops option is %d but number of populations entered using the --new-pop option is %d.  Exiting...\n\n",gNumPops,NumPops);
		exit(1);
	}
	gNumPops = NumPops;
	
	
	/* allocate memory to the global migration counters */
	if(NumPops>1) {
		gMig = (int ****)ECA_CALLOC(2,sizeof(int ***));
		gPotMig = (int ***)ECA_CALLOC(2,sizeof(int **));
		for(i=0;i<2;i++)  {
			gMig[i] = (int ***)ECA_CALLOC(NumPops,sizeof(int **));
			gPotMig[i] = (int **)ECA_CALLOC(NumPops,sizeof(int *));
			for(l=0;l<NumPops;l++)  {
				gMig[i][l] = (int **)ECA_CALLOC(LocaleArray[0]->P.MaxAge+1,sizeof(int *));
				gPotMig[i][l] = (int *)ECA_CALLOC(LocaleArray[0]->P.MaxAge+1,sizeof(int));
				for(t=0;t<=LocaleArray[0]->P.MaxAge;t++) {
					gMig[i][l][t] = (int *)ECA_CALLOC(NumPops,sizeof(int));
				}
			}
		}
	}
	
	/* spew out how many populations there are */
	printf("NUMBER_OF_POPULATIONS : %d\n",NumPops);
	
	/* now print out all the information about each of the populations in series */
	for(which=0;which<NumPops;which++) {
		printf("SUMMARY_OF_SETTINGS_FOR_POPULATION : %d\n",which);
		L = LocaleArray[which];
		/* set T to the max of the L->T's */
		if(L->T > T) {
			T = L->T;
		}
		
		/* print out the demographic parameters */
		PrintPars(L->P,L->Pop,L->S,T);
		
		if(gRunTests) {  /* only print these out if running some tests */
			PrintBillHillPars(L->P,L->Pop);
		}
		/* here we have a section to print out all the parameter values.  It is still incomplete */
		PrintProbDistOfPapas(L->P.MateSpec, 20);
	}
	
	
	/* I OUGHT TO CHECK TO MAKE SURE THE NUMBER OF LOCI AND NUMBER OF ALLELES ARE THE SAME FOR ALL POPS,
	 AND GIVE AN ERROR IF THE LOCUS AND ALLELE NUMBERS MISMATCH. ALSO, SET NumLoc to the global number of loci */
	OverallNumLoc = LocaleArray[0]->S.NumLoc;
	OverallNumAlle = (int *)ECA_CALLOC(OverallNumLoc, sizeof(int));
	for(l=0;l<OverallNumLoc;l++)  {
		OverallNumAlle[l] = LocaleArray[0]->S.NumAlleles[l];
	}
	for(bp=1;bp<NumPops;bp++)  {
		L=LocaleArray[bp];
		if(OverallNumLoc != L->S.NumLoc) {
			fprintf(stderr,"Error!  Numbers of loci mismatching in different locales.  %d in Locale %d, and %d in previous locales.  Must be equal.  Exiting...\n",
					L->S.NumLoc,bp,OverallNumLoc);
			exit(1);
		}
		for(l=0;l<OverallNumLoc;l++)  {
			if(OverallNumAlle[l] != LocaleArray[bp]->S.NumAlleles[l]) {
				fprintf(stderr,"Error!  Allele Number mismatch Locus %d : Locale %d has %d alleles, previous locales had %d.  Exiting...\n",
						l,bp,LocaleArray[bp]->S.NumAlleles[l],OverallNumAlle[l]);
				exit(1);
			}
		}
		
	}
	
	
	/* initialize the RNG */
	SeedFromFile("spip_seeds");
	
	
	/* cycle over the pops and initialize the cohorts in each, AND check to make
	 sure that the max age is the same in all of them. */
	for(bp=0;bp<NumPops;bp++)  {
		L = LocaleArray[bp];
		L->MA = L->P.MaxAge;  /* a convenient  shorthand */
		if(bp==0)
			CheckMA = L->MA;
		else {
			if(CheckMA != L->MA)  {
				fprintf(stderr,"Error! MaxAge of different populations must be equal.\n");
				fprintf(stderr,"\tMaxAge of population %d is %d, but MaxAge of previous populations is %d\n",bp,L->MA,CheckMA);
				exit(1);
			}
		}
		
		
		/*InitCohorts(L->MA,L->NF,L->NM,T,&(L->Males),&(L->Females), L->P,bp, L->MSpace,L->FSpace); */
	}
	
	
	/* initialize all the cohorts in all the populations */
	printf("INITIALIZING_COHORTS_FOR_ALL_POPULATIONS\n");
	MultiPop_InitCohorts(LocaleArray,NumPops,T);
	
	/* down here report how long the simulation will be done for */
	printf("VERBIAGE : Simulation will be run for : %d : years as that was the largest T for all the populations\n",T);
	
	for(t=CheckMA;t<CheckMA+T;t++)  {  /* cycle over all the new cohort years from MA to MA+T-1 */
		
		/* The first thing each year is a migration step.  Note that the youngest individuals migrating
		 are one-year-olds */
		Migrate_O_Rama(LocaleArray,NumPops,t);
		
		for(bp=0;bp<NumPops;bp++)  {  /* cycle over all the locales */
			L = LocaleArray[bp];
			AYearInTheLife(L, t, bp,LocaleArray,NumPops);
		} /* closes loop over the different populations */
		
		/* HERE IS WHERE I WILL HAVE THE MIGRATION STEP! THIS'LL BE BEFORE THE DEATH EPISODE, I GUESS */
		
	}  /* done with simulating the pedigree of the population */
	
	
	/* OK, now we can sprinkle genes around */
	t--;  /* now t is the absolute time of the last cohort born */
	
	
	/**  This is where individuals get marked if they are to be "cohort-sampled"  THIS HAS BEEN DEPRECATED **/
	/*for(bp=0;bp<NumPops;bp++)  {   cycle over all the locales */
	/*	L = LocaleArray[bp];
	 if(L->S.SampleNewbornsOnly==1) {
	 MarkSampledNewborns(L->S,L->NM,L->NF,L->Males,L->Females,bp);
	 }
	 } */
	
	
	PrintMigrantSummary(gMig, gPotMig, NumPops, LocaleArray[0]->P.MaxAge);
	
	
	/* First, we have to determine if any of the populations are getting sampled, and, if they are,
	 we have to cycle over ALL of the individuals that ever existed, because, even those who
	 were born in a different locale may have been sampled after migrating to another locale */
	for(bp=0;bp<NumPops;bp++)  {  /* cycle over all the locales */
		L = LocaleArray[bp];
		if(L->S.pre_male_samp != NULL || L->S.post_male_samp != NULL || L->S.pre_fem_samp != NULL || 
		   L->S.post_fem_samp != NULL   /*|| L->S.SampleNewbornsOnly==1 */ || L->S.dur_repro_fem_samp > 0.0  || L->S.dur_repro_male_samp > 0.0) { 
			WeAreGenotyping++;
		}
	}
	
	
	/* if we are genotyping, then we have to allocate memory to NumFounding for all pops, first, then trace ancestry */
	if(WeAreGenotyping>0) {
		
		/* Allocate memory */
		NumFounding = (int *)ECA_CALLOC(OverallNumLoc, sizeof(int));  
		/* Allocate to the first fold of FounderLocations */
		FounderLocations = (int **)ECA_CALLOC(OverallNumLoc, sizeof(int*));
		FoundLocMem = (int *)ECA_CALLOC(OverallNumLoc, sizeof(int));
		
		/* cycle over the loci and trace the ancestry */
		for(l=0;l<L->S.NumLoc;l++) {
			NEW_SampleMarkedAtLocusl(NumPops, LocaleArray,t,l,OverallNumLoc, &(NumFounding[l]),&(FounderLocations[l]),&(FoundLocMem[l]), gID);
			printf("PROGRESS :  Done tracing ancestry at locus %d\n",l+1);
		}
		
	}
	
	/* then print the allelic ancestry (and when they were sampled */
	for(bp=0;bp<NumPops;bp++)  {  /* cycle over all the locales */
		L = LocaleArray[bp];
		PrintAllelesOfSampledInds(L->Males, L->Females, L->NM, L->NF, OverallNumLoc, L->P.DiscardBefore, t, NULL,NULL);
	}
	
	/* then print the locale origin of the founder alleles (and when they were sampled */
	for(bp=0;bp<NumPops;bp++)  {  /* cycle over all the locales */
		L = LocaleArray[bp];
		PrintAllelesOfSampledInds(L->Males, L->Females, L->NM, L->NF, OverallNumLoc, L->P.DiscardBefore, t, NULL,FounderLocations);
	}
	
	/* then determine the types of all those founder gene copies */
	FounderTypes = DrawFounderTypes(OverallNumLoc, OverallNumAlle, LocaleArray, NumFounding, FounderLocations);
	
	/* then print out the genotypes by allelic type */
	for(bp=0;bp<NumPops;bp++)  {  /* cycle over all the locales */
		L = LocaleArray[bp];
		PrintAllelesOfSampledInds(L->Males, L->Females, L->NM, L->NF, L->S.NumLoc, L->P.DiscardBefore, t, FounderTypes,NULL);
	}
	
	
	/* put the seeds out */
	SeedToFile("spip_seeds");
	
	return(0);
}



/* print a summary report about migration that occurred during the whole simulation */
void PrintMigrantSummary(int ****M, int ***P, int NP, int MA)
{
	int sex,age,from,i;
	double tot;
	
	if(NP==1) return;
	
	for(sex=0;sex<2;sex++) {
		for(from=0;from<NP;from++)  {
			for(age=1;age<=MA;age++)  {
				for(tot=0.0,i=0;i<NP;i++)  {
					tot += M[sex][from][age][i];
				}
				printf("MIGRANT_SUMMARY : ");
				if(sex==0) printf("MALE ");
				else printf("FEM ");
				printf(" %d  :  %d  : %d :  %.0f : %f : ",from,age,P[sex][from][age], tot, tot/P[sex][from][age]);
				for(i=0;i<NP;i++) {
					if(tot>0.0) {
						printf("%f ",M[sex][from][age][i]/tot);
					}
					else {
						printf("0.0  ");
					}
				}
				printf("\n");
			}
		}
	}
}








void CopyIndiv(indiv *From, indiv *To)
{
	sprintf(To->ID, "%s", From->ID);
	To->ma = From->ma;
	To->pa = From->pa;
	To->dead = From->dead;
	To->t = From->t;  
	To->idx = From->idx;
	
	To->birth_pop = From->birth_pop;
	
	To->ron = From->ron;
	To->MaFounder = From->MaFounder;
	To->PaFounder = From->MaFounder;
	
	To->N_PreK = From->N_PreK;
	To->N_PostK = From->N_PostK;
	To->N_DurRepro = From->N_DurRepro;
	
	To->PreSpace = From->PreSpace;
	To->PostSpace = From->PostSpace;
	To->DurReproSpace = From->DurReproSpace;
	
	To->PreK_ts = From->PreK_ts;
	To->PostK_ts = From->PostK_ts;
	To->DurRepro_ts = From->DurRepro_ts;
	
	To->SampLocales = From->SampLocales;
}









/* a function that does all the migrants in and out of the array of populations.  The basic
 way that this operates is that it goes through the populations in order.  In each population
 it cycles over all the individuals and if the individual is not dead, it chooses it to be a 
 migrant according to its age-specific out-migration rate.   Then, if chosen to outmigrate, its
 destination is chosen according to its in-migration rate.  It gets added to the end of that populations
 array of individuals, while it becomes "dead" in the current population.  The NM and NF numbers (number
 of males and females in each population) don't get changed until the very end, so it is not possible
 for someone to migrate out of a population to a new one, then migrate out of that one as well.  
 */
void Migrate_O_Rama(Locale **LA, int NumPops, int t) 
{
	int i,at,sex,j,new_idx;
	int age;
	int MA = LA[0]->P.MaxAge;
	int ***AddNums = (int ***)ECA_CALLOC(2,sizeof(int**));
	double rando;
	int to;
	int this_time;
	
	
	if(NumPops==1) 
		return;
	
	
	/* allocate some memory */
	AddNums[MALE] = (int **)ECA_CALLOC(NumPops, sizeof(int *));
	AddNums[FEMALE] = (int **)ECA_CALLOC(NumPops, sizeof(int *));
	
	for(i=0;i<NumPops;i++)  {
		AddNums[MALE][i] = (int *)ECA_CALLOC(MA+1,sizeof(int*));
		AddNums[FEMALE][i] = (int *)ECA_CALLOC(MA+1,sizeof(int*));
	}
	
	
	for(age=1;age<=MA;age++)  {  /* cycle over the ages */
		at = t - age;  /* absolute time subscript for an individual of age "age" at time t */ 
		
		
		for(i=0;i<NumPops;i++)  {  /* cycle over the populations individuals will migrate out of */
			
			
			/* cycle over the males */
			this_time=0;
			for(j=0;j<LA[i]->NM[at];j++)  {
				
				if(LA[i]->Males[at][j].dead==0) {
					gPotMig[MALE][i][age]++;
					
					rando = (double)ranf();
					if(rando < LA[i]->P.MigPr[t][MALE][age]) {
						
						/* choose the destination */
						to = IntFromProbsRV(LA[i]->P.Dest[t][MALE][age],0,NumPops);
						
						/* increment AddNums for the destination */
						AddNums[MALE][to][age]++;
						
						/* check to make sure there is space there, and realloc if not. */
						while (LA[to]->MSpace[at] < LA[to]->NM[at] + AddNums[MALE][to][age]) {
							LA[to]->Males[at] = (indiv *)realloc( (void *)LA[to]->Males[at], sizeof(indiv) * (LA[to]->MSpace[at] + REALLOC_BLOCK) );
							LA[to]->MSpace[at] += REALLOC_BLOCK;
						}
						
						/* copy the migrant to the end of his destination population */
						new_idx = LA[to]->NM[at] + AddNums[MALE][to][age] - 1;
						CopyIndiv(&(LA[i]->Males[at][j]), &(LA[to]->Males[at][ new_idx  ]) );
						
						/* kill the migrant in his current population */
						LA[i]->Males[at][j].dead = 1;
						
						/* decrement the number of remaining males */
						LA[i]->NMR[at]--;
						
						/* print out the migration event information */
						if(this_time==0) {
							printf("MALE_MIGRATION :  %d  : %d  :  %d  :  %s %d -> %d",t,age,i,LA[i]->Males[at][j].ID, i,to);
						}
						else {
							printf(" | %s %d -> %d ",LA[i]->Males[at][j].ID, i,to);
						}
						this_time++;
						gMig[MALE][i][age][to]++;
					}
				}
			}
			if(this_time>0)
				printf("\n");
			
			this_time = 0;
			/* cycle over the females */
			for(j=0;j<LA[i]->NF[at];j++)  {
				if(LA[i]->Females[at][j].dead==0) {
					gPotMig[FEMALE][i][age]++;
					rando = (double)ranf();
					if(rando < LA[i]->P.MigPr[t][FEMALE][age]) {
						
						/* choose the destination */
						to = IntFromProbsRV(LA[i]->P.Dest[t][FEMALE][age],0,NumPops);
						
						/* increment AddNums for the destination */
						AddNums[FEMALE][to][age]++;
						
						/* check to make sure there is space there, and realloc if not.  NOT YET IMPLEMENTED!!! */
						while (LA[to]->FSpace[at] < LA[to]->NF[at] + AddNums[FEMALE][to][age]) {
							LA[to]->Females[at] = (indiv *)realloc( (void *)LA[to]->Females[at], sizeof(indiv) * (LA[to]->FSpace[at] + REALLOC_BLOCK) );
							LA[to]->FSpace[at] += REALLOC_BLOCK;
						}
						
						/* copy the migrant to the end of her destination population */
						new_idx = LA[to]->NF[at] + AddNums[FEMALE][to][age] - 1;
						CopyIndiv(&(LA[i]->Females[at][j]), &(LA[to]->Females[at][ new_idx  ]) );
						
						/* kill the migrant in her current population */
						LA[i]->Females[at][j].dead = 1;
						
						/* decrement the number of remaining females */
						LA[i]->NFR[at]--;
						
						/* print out the migration event information */
						/* print out the migration event information */
						if(this_time==0) {
							printf("FEM_MIGRATION :  %d  : %d  :  %d  :  %s %d -> %d",t,age,i,LA[i]->Females[at][j].ID, i,to);
						}
						else {
							printf(" | %s %d -> %d ",LA[i]->Females[at][j].ID, i,to);
						}
						this_time++;
						gMig[FEMALE][i][age][to]++;
					}
				}
			}
			if(this_time>0)
				printf("\n");
			
			
		}  /* close the cycle over pops */
		
		/* down here we update the NM and NF and NMR and NFR for this age group in each of the populations */
		for(i=0;i<NumPops;i++)  {
			LA[i]->NM[at] += AddNums[MALE][i][age];
			LA[i]->NF[at] += AddNums[FEMALE][i][age];
			LA[i]->NMR[at] += AddNums[MALE][i][age];
			LA[i]->NFR[at] += AddNums[FEMALE][i][age];
		}
		
	}  /* close the cycle over ages */	
	
	
	
	/* free the memory */
	for(sex=0;sex<2;sex++)  {
		for(i=0;i<NumPops;i++)  {
			free(AddNums[sex][i]);
		}
		free(AddNums[sex]);
	}
	free(AddNums);
	
}




/* a function to retrieve a bunch of information about the demographic
 parameters in a locale.  It exits when it hits the --end-locale option
 and will check to make sure that it has everything that it needs for
 the current population
 */

int GetALocale(Locale **M, int which, int *argc, char **argv[])
{
	Locale *L;
	int new_pop_flag = 0;
	
	M[which] = (Locale *)malloc(sizeof(Locale));
	
	L = M[which];
	
	L->HasLocFile=0;
	
	new_pop_flag = ProcessOptions(argc,argv,&(L->P),&(L->T), &(L->NM), &(L->NF),&(L->NMR), &(L->NFR), &(L->S), (L->LocFileName),&(L->Pop),&(L->Pop.ConstSize),&(L->MSpace),&(L->FSpace),
								  &(L->HasLocFile));
	
	/* here, if no migration parameters got set, we will put them to the defaults */
	if(L->P.MigPr==NULL) {
		L->P.MigPr = AllocAndInitMigPr(L->T, L->P.MaxAge);
	}
	if(L->P.Dest==NULL) {
		L->P.Dest = AllocAndInitMigDest(L->T, L->P.MaxAge, gNumPops, which);
	}
	
	return(new_pop_flag);
}





double ***AllocAndInitMigPr(int T, int MA)
{
	int t,sex;
	
	double ***m=(double ***)ECA_CALLOC(T+MA,sizeof(double **));
	
	for(t=0;t<T+MA;t++) {
		m[t] = (double **)ECA_CALLOC(2,sizeof(double *));
		for(sex=0;sex<2;sex++) {
			m[t][sex] = (double *)ECA_CALLOC(MA+1,sizeof(double));
		}
	}
	return(m);
} 

double ****AllocAndInitMigDest(int T, int MA, int NumPops, int CurrentPop)
{
	int t,sex,age,i;
	
	double ****m=(double ****)ECA_CALLOC(T+MA,sizeof(double ***));
	
	for(t=0;t<T+MA;t++) {
		m[t] = (double ***)ECA_CALLOC(2,sizeof(double **));
		for(sex=0;sex<2;sex++) {
			m[t][sex] = (double **)ECA_CALLOC(MA+1,sizeof(double *));
			for(age=0;age<=MA;age++)  {
				m[t][sex][age] = (double *)ECA_CALLOC(NumPops,sizeof(double));
				for(i=0;i<NumPops;i++)  {
					if(i != CurrentPop)
						m[t][sex][age][i] = 1.0/(NumPops-1.0);
				}
			}
		}
	}
	return(m);
} 







/* high-level function in which everything happens in a 
 single population in a year */
void AYearInTheLife(Locale *L, int t, int bp, Locale **LA, int NumPops) 
{
	int N;
	
	/*  pre-death random sampling of juveniles and adults */
	PreDeathSample(t, L->MA, L->S, L->NM, L->NF, L->Males, L->Females,bp, L->NMR, L->NFR );
	
	/*  print out the pre-kill census for males and females */
	PrintCensusSizes(t,L->MA,L->NFR,L->NMR,"PREKILL_CENSUS_",1,bp);
	
	/* survival from the last year to the current year's reproducing group */
	StandardAnnualDieOff(t, L->MA, L->NM, L->NF, L->Males, L->Females, L->NMR, L->NFR, L->P);
	
	/* post-death random sampling of juveniles and adults */
	PostDeathSample(t, L->MA, L->S, L->NM, L->NF, L->Males, L->Females, bp, L->NMR, L->NFR );
	
	/*  print out the post-kill census for males and females */
	PrintCensusSizes(t,L->MA,L->NFR,L->NMR,"POSTKILL_CENSUS_",1,bp);
	
	/* set the cohort size for this time */
	N = ReturnCohortSize(&(L->Pop),t);
	
	/* Now we make offspring from the remaining females and males */
	if(PopSizeNonZero(t,L->MA,L->NFR,L->NMR)) {
		MakeBabies(L->Females, L->Males, t, L->MA, L->P, N, L->NFR, L->NF, L->NM, L->NMR, L->S, bp, L->MSpace, L->FSpace,LA,NumPops);
	}
	
	/* print the post-reproduction census sizes for males and females */
	PrintCensusSizes(t,L->MA,L->NFR,L->NMR,"POST_REPROD_CENSUS_",0,bp);
	
}










/* this is a simple function that initializes some necessary information
 when and individual born.  the variable sex should get 'M' or 'F' depending on
 whether the individual is a male or a female. */
void FillIndiv(indiv *baby, indiv *mom, indiv *dad, int *counter, int year, int num, char sex, int birthplace, PopPars P)
{
	sprintf(baby->ID,"%c%d_%d_%d",sex,year,birthplace,num);
	baby->ma = mom;
	baby->pa = dad;
	baby->OldestOffsp = NULL;
	baby->dead = 0;
	baby->idx = (*counter) = (*counter) + 1;
	baby->t = year;
	baby->birth_pop = birthplace;
	baby->MaFounder = NULL;
	baby->PaFounder = NULL;
	baby->PreK_ts = NULL;
	baby->PostK_ts = NULL;
	baby->DurRepro_ts = NULL;
	baby->N_PreK = 0;
	baby->N_PostK = 0;
	baby->N_DurRepro = 0;
	baby->PreSpace = 0;
	baby->PostSpace = 0;
	baby->DurReproSpace = 0;
	
	baby->SampLocales = (int **)ECA_CALLOC(3,sizeof(int*));
	
	baby->ReportMe = 0;
	
	/* finally, down here, if it is a male, we give him the sticky-repro-var.  Simulate it from a 
	 gamma distribution with shape parameter alpha, and we may as well use 1 for the scale parameter */
	if(P.sticky_m_alpha>0 && sex=='M') {
		baby->sticky_birth_rate = (double)gengam(1.0f ,(float)P.sticky_m_alpha);
	}
	else {
		baby->sticky_birth_rate = -1.0;
	}
}


/* these are high-level functions that operate on single populations. Mostly these are 
 collections of functions that have been lumped together so that things are cleaner in 
 main(). */
void PreDeathSample(int t, int MA, SampPars S, int *NM, int *NF, indiv **Males, indiv **Females, int locale, int *NMR, int *NFR)
{
	int at,age;
	for(age=1;age<=MA;age++)  {
		at = t - age;  /* absolute time subscript for an individual of age "age" at time t */ 
		if(S.pre_male_samp != NULL && S.pre_male_samp_years[t]) MarkAsSampled(Males[at],NM[at],S.pre_male_samp[age],BEF_KILL,t,locale, S.LethalityProb, &(NMR[at]) );
		if(S.pre_fem_samp != NULL && S.pre_fem_samp_years[t]) MarkAsSampled(Females[at],NF[at],S.pre_fem_samp[age],BEF_KILL,t,locale, S.LethalityProb, &(NFR[at]) );
	}
}
void PostDeathSample(int t, int MA, SampPars S, int *NM, int *NF, indiv **Males, indiv **Females, int locale, int *NMR, int *NFR )
{
	int at,age;
	for(age=1;age<=MA;age++)  {
		at = t - age;  /* time subscript for an individual of age "age" at time t */ 
		if(S.post_male_samp != NULL && S.post_male_samp_years[t]) MarkAsSampled(Males[at],NM[at],S.post_male_samp[age],AFT_KILL,t,locale, S.LethalityProb, &(NMR[at]));
		if(S.post_fem_samp != NULL && S.post_fem_samp_years[t]) MarkAsSampled(Females[at],NF[at],S.post_fem_samp[age],AFT_KILL,t,locale, S.LethalityProb, &(NFR[at]));
	}
}
void StandardAnnualDieOff(int t, int MA, int *NM, int *NF, indiv **Males, indiv **Females, int *NMR, int *NFR, PopPars P)
{
	int age, at;
	for(age=1;age<=MA;age++)  {
		at = t - age;  /* this is the absolute time subscript for an individual who could be a parent at age "age" of the current cohort at t */ 
		KillOff(Males[at],NM[at],&(NMR)[at],P.malesurv[age],t);
		KillOff(Females[at],NF[at],&(NFR)[at],P.femsurv[age],t);
	}
	/* down here need to add a step to kill off individuals that will be age MA+1 in next time step.
	 since age will be MA+1 at the end of the above loop we should be able to do this, so long as at>0 */
	at = t - age;  /* this is the absolute time subscript for an individual who could be a parent at age "age" of the current cohort at t */ 
	if(at>=0) {
		printf("ABOUT_TO_KILL_REMAINING_OLD_MALES\n");
		KillOff(Males[at],NM[at],&(NMR)[at],0.0,t);
		printf("ABOUT_TO_KILL_REMAINING_OLD_FEMALES\n");
		KillOff(Females[at],NF[at],&(NFR)[at],0.0,t);
	}
	
}



char *dsn_type2string(enum dsn_type d) {
	switch(d) {
		case(NEG_BINOM):
			return("NegativeBinomial");
		case(POISSON):
			return("POISSON");
		case(BINARY):
			return("BINARY");
		default:
			return("Unknown");
	}	
}

void PrintPars(PopPars P, PopSizePars Pop, SampPars S, int T)
{
	int i,sex,t;
	
	printf("SETUP : MaxAge : %d\n", P.MaxAge);
	printf("SETUP : FemaleSurvivalRates : ");
	for(i=1;i<=P.MaxAge;i++) {
		printf("%.5f ",P.femsurv[i]);
	} 
	printf("\n");
	
	printf("SETUP : MaleSurvivalRates : ");
	for(i=1;i<=P.MaxAge;i++) {
		printf("%.5f ",P.malesurv[i]);
	}
	printf("\n");
	
	printf("SETUP : FemaleAgeSpecificFecundities : ");
	for(i=1;i<=P.MaxAge;i++) {
		printf("%.5f ",P.Fasrf[i]);
	}
	printf("\n");
	printf("SETUP : MaleAgeSpecificFecundities : ");
	for(i=1;i<=P.MaxAge;i++) {
		printf("%.5f ",P.Masrp[i]);
	}
	printf("\n");
	
	printf("SETUP : FemaleAgeSpecificProbabilitiesOfReproduction : ");
	for(i=1;i<=P.MaxAge;i++) {
		printf("%.5f ",P.FPR[i]);
	}
	printf("\n");
	printf("SETUP : MaleAgeSpecificProbabilitiesOfReproduction : ");
	for(i=1;i<=P.MaxAge;i++) {
		printf("%.5f ",P.MPR[i]);
	}
	printf("\n");
	
	printf("SETUP : FemaleAgeSpecificProbabilitiesOfPostReproductiveDeath : ");
	for(i=1;i<=P.MaxAge;i++) {
		printf("%.5f ",P.FPDPR[i]);
	}
	printf("\n");
	printf("SETUP : MaleAgeSpecificProbabilitiesOfPostReproductiveDeath : ");
	for(i=1;i<=P.MaxAge;i++) {
		printf("%.5f ",P.MPDPR[i]);
	}
	printf("\n");
	
	
	printf("SETUP : FemaleRatioOfMeanToVarianceOfOffspringNumber : %f\n", P.Fcv);
	printf("SETUP : MaleRatioOfMeanToVarianceOfOffspringNumber : %f\n", P.MCVrs);
	printf("SETUP : MateFidelityParameter : %f\n", P.MateSpec);
	
	printf("SETUP : SexRatio : %f\n", P.SexRatio);
	
	printf("SETUP : OffspringDistribution : %s\n", dsn_type2string(P.OffsDsn));
	
	printf("SETUP : CohortSizesRandomOrFixed : ");
	if(P.CohortSizesRandom==1) printf("Random\n");
	else printf("Fixed\n");

	if(P.MigPr==NULL) {
		printf("SETUP : MigPars_OutAndIn :  NO_MIGRATION_PARAMETERS_SET_UP \n");
	}
	else for(t=0;t<P.MaxAge+T;t++) {
		for(sex=0;sex<2;sex++)  {
			for(i=1;i<=P.MaxAge;i++)  { int j;
				printf("SETUP : MigPars_OutAndIn : time= %3d  : sex= %s  :  age= %3d  :  ",t,(sex==0 ? "MALE  " : "FEMALE" ),i);
				printf(" outrate= %.5f   :  inrates= NOT REPORTED AT THE MOMENT",P.MigPr[t][sex][i]);
				
				for(j=0;j<gNumPops;j++)  {
					; //printf(" %.5f ",P.Dest[t][sex][i][j]);
				}
				printf("\n");
			}
		}
	}
	
	printf("SETUP : DiscardBefore : %d\n", P.DiscardBefore);
	
	printf("SETUP : DiscardParentsOrAll : ");
	if(P.TossAll==1) printf("All\n");
	else printf("Parents\n");
	
	
	printf("SETUP : CohortSizeModel : ");
	if(strcmp(Pop.CohortSizeMode,"const")==0) {
		printf("Constant\n");
		printf("SETUP : ConstantCohortSize : %d\n", Pop.ConstSize);
	}
	else if(strcmp(Pop.CohortSizeMode,"var")==0) {
		printf("Variable\n");
		printf("SETUP : VariableCohortSizes : ");
		for(i=P.MaxAge;i<P.MaxAge+T;i++)  {
			printf("%d ",Pop.CohortSizes[i]);
		}
		printf("\n");
	}
	printf("SETUP : NumberOfYears : %d\n",T);	
	
	
	if(S.pre_male_samp != NULL) {
		printf("SETUP : PreDeathMaleSamplingProportions :  ");
		for(i=1;i<=P.MaxAge;i++) {
			printf(" %f",S.pre_male_samp[i]);
		}
		printf("\n");
		printf("SETUP : PreDeathMaleSamplingYears :  ");
		for(i=0;i<P.MaxAge+T;i++) {
			if(S.pre_male_samp_years[i]==1) printf(" %d",i);
		}
		printf("\n");
	}
	
	if(S.pre_fem_samp != NULL) {
		printf("SETUP : PreDeathFemSamplingProportions :  ");
		for(i=1;i<=P.MaxAge;i++) {
			printf(" %f",S.pre_fem_samp[i]);
		}
		printf("\n");
		printf("SETUP : PreDeathFemSamplingYears :  ");
		for(i=0;i<P.MaxAge+T;i++) {
			if(S.pre_fem_samp_years[i]==1) printf(" %d",i);
		}
		printf("\n");
	}
	if(S.post_male_samp != NULL) {
		printf("SETUP : PostDeathMaleSamplingProportions :  ");
		for(i=1;i<P.MaxAge;i++) {
			printf(" %f",S.post_male_samp[i]);
		}
		printf("\n");
		printf("SETUP : PostDeathMaleSamplingYears :  ");
		for(i=0;i<P.MaxAge+T;i++) {
			if(S.post_male_samp_years[i]==1) printf(" %d",i);
		}
		printf("\n");
	}
	
	if(S.post_fem_samp != NULL) {
		printf("SETUP : PostDeathFemSamplingProportions :  ");
		for(i=1;i<P.MaxAge;i++) {
			printf(" %f",S.post_fem_samp[i]);
		}
		printf("\n");
		printf("SETUP : PostDeathFemSamplingYears :  ");
		for(i=0;i<P.MaxAge+T;i++) {
			if(S.post_fem_samp_years[i]==1) printf(" %d",i);
		}
		printf("\n");
	}
}



void PrintBillHillPars(PopPars P, PopSizePars Pop)
{
	int i,j;
	double x,y,s;
	double *thetam,*thetaf,*exf,*exm, ex2;
	int MA = P.MaxAge;
	double *pm = P.MPR, *pf = P.FPR, *fm = P.Masrp, *ff = P.Fasrf;
	double numer,denom;
	double LM,LF;
	double varm,varf,mlifeout,flifeout;
	
	thetam = (double *)ECA_CALLOC(MA+1, sizeof(double));
	thetaf = (double *)ECA_CALLOC(MA+1, sizeof(double));
	exf = (double *)ECA_CALLOC(MA+1, sizeof(double));
	exm = (double *)ECA_CALLOC(MA+1, sizeof(double));
	
	/* first print the line that says what years we are dealing with */
	printf("SUMMARY : StableAgeDistAges : ");
	for(i=0;i<=P.MaxAge;i++)  {
		printf("%d ",i);
	}
	printf("\n");
	
	/* now, down here we also want to print some summaries of things */
	y = 1.0; /* this starts at 1.0 because we have to account for all the newborns that are 
	 in the population */
	x = 1.0;
	thetaf[0] = 1.0;
	for(i=1;i<=P.MaxAge;i++)  {
		x *= P.femsurv[i];
		thetaf[i] = x;
		y += x;
	}
	printf("SUMMARY : FemaleStableAgeDist : ");
	for(i=0;i<=P.MaxAge;i++)  {
		thetaf[i] /= y;
		printf("%f  ",thetaf[i]);
	}
	printf("\n");
	
	/* now, we do the same for males */
	y = 1.0; /* this starts at 1.0 because we have to account for all the newborns that are 
	 in the population */
	x = 1.0;
	thetam[0] = 1.0;
	for(i=1;i<=P.MaxAge;i++)  {
		x *= P.malesurv[i];
		thetam[i] = x;
		y += x;
	}
	printf("SUMMARY : MaleStableAgeDist : ");
	for(i=0;i<=P.MaxAge;i++)  {
		thetam[i] /= y;
		printf("%f  ",thetam[i]);
	}
	printf("\n");
	
	
	/* OK, now we are going to compute Lm and Lf, since the prob of reproducing at age 0 is zero, we
	 just cycle from 1 to MA inclusive  */
	for(denom=0.0,numer=0.0,i=1;i<=MA;i++)  {
		numer += i * thetaf[i] * pf[i] * ff[i];
		denom += thetaf[i] * pf[i] * ff[i];
	}
	LF = numer/denom;
	printf("SUMMARY : FemaleGenerationInterval : %f\n",LF);
	for(denom=0.0,numer=0.0,i=1;i<=MA;i++)  {
		numer += i * thetam[i] * pm[i] * fm[i];
		denom += thetam[i] * pm[i] * fm[i];
	}
	LM = numer/denom;
	printf("SUMMARY : MaleGenerationInterval : %f\n",LM);
	printf("SUMMARY : BothSexesGenerationInterval : %f\n",(LM+LF)/2.0);
	
	
	
	/* and now we want to compute the means and variances of lifetime 
	 offspring output.  Computing these is somewhat long and nasty. */
	
	/* first step is to compute the expected number of offspring for each
	 male or female of a given age, given that they are actually reproducing 
	 that year */
	/* go for females first */
	for(denom=0.0,i=1;i<=MA;i++)  {
		denom += thetaf[i] * pf[i] * ff[i];
	}
	printf("SUMMARY : ExpectedNumoffsForBreedingFemales : ");
	for(i=1;i<=MA;i++)  {
		exf[i] = (1.0/(1.0-P.SexRatio)) * (ff[i]/denom) * thetaf[0];
		printf("%f  ", exf[i] );
	}
	printf("\n");
	/* then do the males */
	for(denom=0.0,i=1;i<=MA;i++)  {
		denom += thetam[i] * pm[i] * fm[i];
	}
	printf("SUMMARY : ExpectedNumoffsForBreedingMales : ");
	for(i=1;i<=MA;i++)  {
		exm[i] = (1.0/P.SexRatio) * (fm[i]/denom) * thetam[0];
		printf("%f  ", exm[i] );
	}
	printf("\n");
	
	/* now we will compute the lifetime-expected number of offspring per male and female. */
	/* note that numer and denom are not actually a numerator and denominator here, and note that
	 s is the prob of surviving to age i.*/
	for(numer=0.0,s=1.0,i=1;i<=MA;i++)  {
		s *= P.femsurv[i];
		numer += s * pf[i] * exf[i];
	}
	flifeout = numer;
	printf("SUMMARY : ExpectedFemLifetimeOffspringOutput : %f\n",flifeout);
	for(numer=0.0,s=1.0,i=1;i<=MA;i++)  {
		s *= P.malesurv[i];
		numer += s * pm[i] * exm[i];
	}
	mlifeout = numer;
	printf("SUMMARY : ExpectedMaleLifetimeOffspringOutput : %f\n",mlifeout);
	
	
	/* finally, we will compute the variance in lifetime family size */
	for(numer=-4.0,s=1.0,i=1;i<=MA;i++)  {
		s *= P.malesurv[i];
		ex2 = (exm[i] / P.MCVrs) + exm[i]*exm[i];
 		numer += s * pm[i] * ex2;
 		/* now we sum over the cross-terms, too */
 		for(j=1;j<i;j++)  {
 			numer += 2.0 * exm[i] * exm[j] * pm[i] * pm[j] * s;
 		}
	}
	varm = numer;
	printf("SUMMARY : MaleVarianceInLifetimeReproduction : %f\n",varm);
	
	
	for(numer=-4.0,s=1.0,i=1;i<=MA;i++)  {
		s *= P.femsurv[i];
		ex2 = (exf[i] / P.Fcv) + exf[i]*exf[i];
 		numer += s * pf[i] * ex2;
 		/* now we sum over the cross-terms, too */
 		for(j=1;j<i;j++)  {
 			numer += 2.0 * exf[i] * exf[j] * pf[i] * pf[j] * s;
 		}
	}
	varf = numer;
	printf("SUMMARY : FemaleVarianceInLifetimeReproduction : %f\n",varf);
	
	/* then at the end, assuming an equal sex ratio, calculate Hill's Ne */
	printf("SUMMARY : HillsNeAssumingEqualSexRatio : %f\n",8.0 * ((LM+LF)/2.0) * Pop.ConstSize / (varf + varm + 4.0));
	printf("SUMMARY : HillsAnnualNeAssumingEqualSexRatio : %f\n",8.0 * ((LM+LF)/2.0) * ((LM+LF)/2.0) * Pop.ConstSize / (varf + varm + 4.0));
	
	
	free(thetam);
	free(thetaf);
}


/*  
 This is a function that prints out the appropriate pedigree lines.
 If year is less than DiscardBefore value + MaxAge, then 
 1) it prints out the unmodified pedigree on a "DISCARDED_PDGEE" line. 
 2) It changes the ma and pa fields to NULL according to TossAll
 If TossAll is 1 then self's ma and pa fields get NULL.
 If TossAll is 0 then self's ma and pa fields get NULL only if it ma and pa
 pointers, respectively, point to either NULL, or to an individual
 with a time field (t) less than DiscardBefore.  This way of doing
 things is a little more realistic and reasonable, I believe, because,
 otherwise, none of the individuals in the first MaxAge years are related.   
 
 If year is greater than or equal to DiscardBefore, then it also prints out the 
 MODIFIED (i.e. after the above changes) information on "PEDIGREE" lines.
 */
void PrintAndDiscardPedigree(indiv *self, int year, int DiscardBefore, int TossAll, int MaxAge)  
{
	char mom[100], dad[100];
	
	/* store the strings for the parents */
	if(self->ma==NULL) 
		sprintf(mom,"0");
	else 
		sprintf(mom,"%s",self->ma->ID);
	if(self->pa==NULL) 
		sprintf(dad,"0");
	else 
		sprintf(dad,"%s",self->pa->ID);
	
	
	/* print the discarded pedigree line and make the modifications, if year is early enough */
	if(year < DiscardBefore + MaxAge) {
		printf("DISCARDED_PDGEE : %d : %d :  %s    %s    %s\n",year,self->birth_pop,self->ID,dad,mom);
		
		if(TossAll==1) {
			self->ma=NULL;
			self->pa=NULL;
		}
		else {
			/* rub out ma if she is too old */
			if(self->ma != NULL) 
				if(self->ma->t < DiscardBefore)
					self->ma=NULL;
			/* rub out dad if he is too old, too */
			if(self->pa != NULL) 
				if(self->pa->t < DiscardBefore)
					self->pa=NULL;
			
			
		}
	}
	
	/* then, if year is >= DiscardBefore, we print out the normal PEDIGREE lines */
	if(year >= DiscardBefore) {
		/* store the strings for the parents */
		if(self->ma==NULL) 
			sprintf(mom,"0");
		else 
			sprintf(mom,"%s",self->ma->ID);
		if(self->pa==NULL) 
			sprintf(dad,"0");
		else 
			sprintf(dad,"%s",self->pa->ID);
		
		/* then print the PEDIGREE line based on the newly modified information */
		printf("PEDIGREE : %d : %d :  %s    %s    %s\n",year,self->birth_pop,self->ID,dad,mom);
		
	}
	
	/*	printf("CurrentValueOfgID = %d\n",gID);  */		
}


/* read data out of the named locus parameters file */
void ReadLocFile(SampPars *S, char *FileName)
{
	int j,k;
	double sum;
	FILE *in;
	
	/* now get all the locus pars */
	if( (in=fopen(FileName,"r"))==NULL ) {
		fprintf(stderr,"\n\nCouldn't open file \"%s\" for locus information.\nExiting...\n\n",FileName);
		exit(1);
	}
	
	printf("ALLEFREQS : opened file \"%s\" to get genetic data\n",FileName);
	while(eat_comments(in,'&')) ;
	/* get the number of loci */
	fscanf(in," %d",&(S->NumLoc) );
	
	printf("ALLEFREQS : number of loci in file %s is %d\n",FileName,S->NumLoc);
	
	while(eat_comments(in,'&')) ;
	/* then get all the rest of the locus parameters */
	S->NumAlleles = (int *)ECA_CALLOC(S->NumLoc, sizeof(int));
	
	while(eat_comments(in,'&')) ;
	S->AlleFreqs = (double **)ECA_CALLOC(S->NumLoc, sizeof(double *));
	for(j=0;j<S->NumLoc;j++)  {
		while(eat_comments(in,'&')) ;
		fscanf(in," %d",&(S->NumAlleles[j]));
		S->AlleFreqs[j] = (double *)ECA_CALLOC(S->NumAlleles[j], sizeof(double));
		for(sum=0.0,k=0;k<S->NumAlleles[j];k++)  {
			while(eat_comments(in,'&')) ;
			fscanf(in," %lf",&(S->AlleFreqs[j][k]));
			sum += S->AlleFreqs[j][k];
		}
		printf("ALLEFREQS : Locus %d : %d Alleles : ", j+1,S->NumAlleles[j]);
		for(k=0;k<S->NumAlleles[j];k++)  {
			S->AlleFreqs[j][k] /= sum;
			printf(" %f ",S->AlleFreqs[j][k]);
		}
		printf("\n");
		
	}
	fclose(in);
}

int DispersionParError(double d, char *whichoption) 
{
	if(d<=0.0)  {
		fprintf(stderr,"Error! The argument to %s cannot be less than or equal to zero!\n",whichoption);
		return(1);
	}
	if(d>1.0) {
		fprintf(stderr,"Error! The argument to %s cannot be greater than one!  That would represent underdispersion relative to the Poisson distribution, and has not been implemented in spip.\n",whichoption);
		return(1);
	}
	
	return(0);
	
}



/*  
 This goes through and with probability "prob" we "genetically sample" each individual
 that is still living.  At this stage, "genetically sampling" just means that we mark it
 as having been sampled, (we also say which year it was sampled, and whether it was
 before or after the killing spree) so we will include it later in tracing the ancestry of loci, etc.
 
 N is the number of individuals (alive or dead) that are in the array I. 
 
 IsPostKilling = BEF_KILL means sampling was before the death episode.  and ==1 means it was
 after the death episode.   
 
 t is the year in which this sampling is taking place.
 
 NR is for the number of males remaining.
 */
void MarkAsSampled(indiv *I, int N, double prob, enum samtype IsPostKilling, int t, int locale, double LethalityProb, int *NR)
{
	int i;
	
	for(i=0;i<N;i++)  {
		if(I[i].dead==0 && ranf() < prob) {
			if(IsPostKilling==AFT_KILL) {
				if(I[i].PostSpace <= I[i].N_PostK)  {  /* realloc space if necessary */
					I[i].PostK_ts = (int *)realloc( (void *)I[i].PostK_ts, sizeof(int) * (I[i].PostSpace + REALLOC_BLOCK) );
					I[i].SampLocales[IsPostKilling] = (int *)realloc( (void *)I[i].SampLocales[IsPostKilling], sizeof(int) * (I[i].PostSpace + REALLOC_BLOCK) );
					I[i].PostSpace += REALLOC_BLOCK;
				}
				I[i].SampLocales[IsPostKilling][I[i].N_PostK] = locale;
				I[i].PostK_ts[ I[i].N_PostK++ ] = t;
			}
			else if(IsPostKilling==BEF_KILL) {
				if(I[i].PreSpace <= I[i].N_PreK)  {  /* realloc space if necessary */
					I[i].PreK_ts = (int *)realloc( (void *)I[i].PreK_ts, sizeof(int) * (I[i].PreSpace + REALLOC_BLOCK) );
					I[i].SampLocales[IsPostKilling] = (int *)realloc( (void *)I[i].SampLocales[IsPostKilling], sizeof(int) * (I[i].PreSpace + REALLOC_BLOCK) );
					I[i].PreSpace += REALLOC_BLOCK;
				}
				I[i].SampLocales[IsPostKilling][I[i].N_PreK] = locale;
				I[i].PreK_ts[ I[i].N_PreK++ ] = t;
			}
			else if(IsPostKilling==WHILE_REPRO) {
				if(I[i].DurReproSpace <= I[i].N_DurRepro)  {  /* realloc space if necessary */
					I[i].DurRepro_ts = (int *)realloc( (void *)I[i].DurRepro_ts, sizeof(int) * (I[i].DurReproSpace + REALLOC_BLOCK) );
					I[i].SampLocales[IsPostKilling] = (int *)realloc( (void *)I[i].SampLocales[IsPostKilling], sizeof(int) * (I[i].DurReproSpace + REALLOC_BLOCK) );
					I[i].DurReproSpace += REALLOC_BLOCK;
				}
				I[i].SampLocales[IsPostKilling][I[i].N_DurRepro] = locale;
				I[i].DurRepro_ts[ I[i].N_DurRepro++ ] = t;
			}
			else {
				fprintf(stderr,"IsPostKilling not equal to 0 or 1 or 2.  Equal to %d\nExiting\n",IsPostKilling);
				exit(1);
			}
			/* down here we enforce the lethality of the sampling that got done here */
			if(ranf() < LethalityProb) {
				I[i].dead = 1;
				(*NR)--;
				printf("KILLING:  Just Killed individual : %s  at time %d at age %d as a sampling lethality!\n",I[i].ID,t-1,t-I[i].t-1);
			}
		}
	}
}

/* a bit of a kluge.  This was written after MarkAsSampled.  This function
 operates on a pointer to an individual, and it just marks it as sampled with probability
 1, regardless of whether it is dead or not */
void MarkIndAsSampled(indiv *I, enum samtype IsPostKilling, int t, int locale)
{
	if(IsPostKilling==AFT_KILL) {
		if(I->PostSpace <= I->N_PostK)  {  /* realloc space if necessary */
			I->PostK_ts = (int *)realloc( (void *)(I->PostK_ts), sizeof(int) * (I->PostSpace + REALLOC_BLOCK) );
			I->SampLocales[IsPostKilling] = (int *)realloc( (void *)I->SampLocales[IsPostKilling], sizeof(int) * (I->PostSpace + REALLOC_BLOCK) );
			I->PostSpace += REALLOC_BLOCK;
		}
		I->SampLocales[IsPostKilling][I->N_PostK] = locale;
		I->PostK_ts[ I->N_PostK++ ] = t;
	}
	else if(IsPostKilling==BEF_KILL) {
		if(I->PreSpace <= I->N_PreK)  {  /* realloc space if necessary */
			I->PreK_ts = (int *)realloc( (void *)(I->PreK_ts), sizeof(int) * (I->PreSpace + REALLOC_BLOCK) );
			I->SampLocales[IsPostKilling] = (int *)realloc( (void *)I->SampLocales[IsPostKilling], sizeof(int) * (I->PreSpace + REALLOC_BLOCK) );
			I->PreSpace += REALLOC_BLOCK;
		}
		I->SampLocales[IsPostKilling][I->N_PreK] = locale;
		I->PreK_ts[ I->N_PreK++ ] = t;
	}
	else if(IsPostKilling==WHILE_REPRO) {
		if(I->DurReproSpace <= I->N_DurRepro)  {  /* realloc space if necessary */
			I->DurRepro_ts = (int *)realloc( (void *)I->DurRepro_ts, sizeof(int) * (I->DurReproSpace + REALLOC_BLOCK) );
			I->SampLocales[IsPostKilling] = (int *)realloc( (void *)I->SampLocales[IsPostKilling], sizeof(int) * (I->DurReproSpace + REALLOC_BLOCK) );
			I->DurReproSpace += REALLOC_BLOCK;
		}
		I->SampLocales[IsPostKilling][I->N_DurRepro] = locale;
		I->DurRepro_ts[ I->N_DurRepro++ ] = t;
	}
	else {
		fprintf(stderr,"IsPostKilling not equal to 0 or 1 or 2.  Equal to %d\nExiting\n",IsPostKilling);
		exit(1);
	}
}

/*
 This cycles through all the individuals in M and in the array F.  For each of those
 individuals that has been genetically sampled (N_PreK>0 || N_PostK>0) we trace its
 ancestry at locus l.
 
 NumPops is the number of different locales
 
 LA is an array of pointers to each locale
 
 NM and NF are the numbers of females in the different cohorts.
 
 LastT is the time at which the last simulated individuals were born
 
 no_ibd_t is the time at which all gene copies born into individuals at that
 time or before are considered to be non-identical by descent 
 
 l is the index of the current locus.
 
 L is the total number of loci, which we include so that we can allocate memory as we go.
 
 nf is an output parameter which returns the number of founding genes at the locus.
 
 fl is an array which will holds the locales from which each founder gene comes from, since it gets modified within this
 function, it is actually the address of the pointer that points to an array of ints.
 
 flm is an output parameter that holds the amount of memory allocated to fl---memory gets realloc'd to fl in this function
 
 
 n is the total number of individuals in the pedigree
 
 After this function has run, each individual who was marked to be sampled will 
 have entries in his or her MaFounder and PaFounder fields that give the index of
 the founding gene copies there. 
 
 THIS IS THE MULTI-POP VERSION OF THIS FUNCTION, WHICH MEANS IT HAS TO CYCLE OVER ALL INDIVIDUALS
 IN ALL POPULATIONS FOR EACH LOCUS.  SO, THE INPUT IS THE ARRAY OF LOCALES.
 
 */
void SampleMarkedAtLocusl(/* indiv **M, indiv **F */ int NumPops, Locale **LA, /*int *NM, int *NF,*/ int LastT, /* int no_ibd_t, */ int l, int L, int *nf, int **fl, int *flm, int n)
{
	int t,i,bp;
	int ReInit=1;
	int ma, pa;
	int **indics;
	int *AncsIdx, *AncsGenes;  /* to be used repeatedly by AncestryOfSampledGenes()  */
	int *NM, *NF; 
	indiv **M,**F;
	int no_ibd_t;
	
	/* allocate space for indics.  Don't forget to free it at the end of the function */
	indics = (int **)ECA_CALLOC(n+1, sizeof(int *));
	for(i=0;i<=n;i++)  {
		indics[i] = (int *)ECA_CALLOC(2,sizeof(int));
	}
	AncsIdx = (int *)ECA_CALLOC(n+1, sizeof(int *));
	AncsGenes = (int *)ECA_CALLOC(n+1, sizeof(int *));
	
	/* cycle over locales to record the earliest --discard-before time.  This we record as
	 no_ibd_t, and we use it later */
	no_ibd_t = LA[0]->P.DiscardBefore;
	for(bp=1;bp<NumPops;bp++)  {
		if(LA[bp]->P.DiscardBefore < no_ibd_t) {
			no_ibd_t = LA[bp]->P.DiscardBefore;
		}
	}
	
	
	/* HERE I SHOULD CYCLE OVER ALL LOCALES */
	for(bp=0;bp<NumPops;bp++)  {
		
		/* cycle over time from discard-before forward */  /* THIS IBD TIME CAN BE SPECIFIC TO EACH POPULATION */
		for(t=LA[bp]->P.DiscardBefore;t<=LastT;t++) {
			
			/* define some convenient shorthand variables */
			NF = LA[bp]->NF;
			NM = LA[bp]->NM;
			F = LA[bp]->Females;
			M = LA[bp]->Males;
			
			/* then cycle over every female born in the locale */
			for(i=0;i<NF[t];i++) {
				
				/* then if this individual was marked for sampling at some point... */
				if(F[t][i].N_PreK + F[t][i].N_PostK + F[t][i].N_DurRepro > 0) {
					/* when we call this one, we reset the "indics" when we do the first marked individual.
					 that is what the ReInit=1 takes care of the first time.  Then it gets set to 0. */
					/*printf("Entering AncOfSampGenes for female t=%d i=%d with ReInit=%d\n", t,i,ReInit); */
					AncestryOfSampledGenes(indics, gID, &(F[t][i]),nf,&ma,&pa, ReInit, AncsIdx, AncsGenes, no_ibd_t, fl, flm );  /* THIS WILL HAVE TO BE CHANGED SO THAT NO_IBD_T is the minimum --discard-before TIME (it is only used as a 
					 limit for resetting indics if need be,
					 we have to pass FounderLocations into it. */
					ReInit = 0;
					
					/* now allocate space to hold the founder gene indexes of the individual */
					if(F[t][i].MaFounder==NULL)  {
						F[t][i].MaFounder = (int *)ECA_CALLOC(L,sizeof(int));
					}
					if(F[t][i].PaFounder==NULL) {
						F[t][i].PaFounder = (int *)ECA_CALLOC(L,sizeof(int));
					}
					/* then assign the ma and pa values to this individual's ma and pa founder indexes for this locus */
					F[t][i].PaFounder[l] = pa;
					F[t][i].MaFounder[l] = ma;
				}
			}
			/* then cycle over every male */
			for(i=0;i<NM[t];i++) {
				
				/* then if this individual was marked for sampling at some point... */
				if(M[t][i].N_PreK + M[t][i].N_PostK +  M[t][i].N_DurRepro > 0) {
					/* when we call this one, we reset the "indics" when we do the first marked individual.
					 that is what the ReInit=1 takes care of the first time.  Then it gets set to 0. */
					/* printf("Entering AncOfSampGenes for male t=%d i=%d with ReInit=%d\n", t,i,ReInit);  */
					AncestryOfSampledGenes(indics, gID, &(M[t][i]),nf,&ma,&pa, ReInit, AncsIdx, AncsGenes, no_ibd_t, fl, flm);
					ReInit = 0;
					
					/* now allocate space to hold the founder gene indexes of the individual */
					if(M[t][i].MaFounder==NULL)  {
						M[t][i].MaFounder = (int *)ECA_CALLOC(L,sizeof(int));
					}
					if(M[t][i].PaFounder==NULL) {
						M[t][i].PaFounder = (int *)ECA_CALLOC(L,sizeof(int));
					}
					/* then assign the ma and pa values to this individual's ma and pa founder indexes for this locus */
					M[t][i].PaFounder[l] = pa;
					M[t][i].MaFounder[l] = ma;
				}
			}
		} /* closes loop over times in locale */
	} /* closes loop over locales */
	
	/* at the end, free up indics */
	for(i=0;i<=n;i++)  {
		free(indics[i]);
	}
	free(indics);
	free(AncsIdx);
	free(AncsGenes);
}



/*
 
 THIS ONE HAS BEEN MODIFIED TO DEAL WITH THE FACT THAT WITH MIGRATION THERE
 WILL BE CURRENT AND "GHOST" INDIVIDUALS.  TO DEAL WITH THAT, AT EACH AGE, WE
 GO THROUGH AND MAKE A LIST OF POINTERS TO INDIVIDUALS, THAT WILL BE SUBSCRIPTED
 VIA THEIR GLOBAL ID'S.  WE GO THROUGH AND ASSIGN POINTERS TO THE INDIVIDUALS, AND
 IF WE GET TO ANOTHER INSTANCE OF AN INDIVIDUAL (FOR INSTANCE, ONE THAT MIGRATED) THEN,
 IF IT HAS BEEN GENETICALLY SAMPLED MORE TIMES, THEN WE ASSIGN ITS ADDRESS TO THE POINTER.
 
 THEN WE CYCLE OVER THE POINTERS TO DETERMINE IF WE SHOULD TRACE THEIR ANCESTRY.  AT THIS
 PHASE WE ASSIGN A "USE THIS ONE" TAG TO THOSE INDIVIDUALS SO THAT THE PRINT_ALLELES_OF_SAMPLED_INDS
 FUNCTION KNOW WHOM TO OPERATE ON.  
 
 This cycles through all the individuals in M and in the array F.  For each of those
 individuals that has been genetically sampled (N_PreK>0 || N_PostK>0) we trace its
 ancestry at locus l.
 
 NumPops is the number of different locales
 
 LA is an array of pointers to each locale
 
 NM and NF are the numbers of females in the different cohorts.
 
 LastT is the time at which the last simulated individuals were born
 
 no_ibd_t is the time at which all gene copies born into individuals at that
 time or before are considered to be non-identical by descent 
 
 l is the index of the current locus.
 
 L is the total number of loci, which we include so that we can allocate memory as we go.
 
 nf is an output parameter which returns the number of founding genes at the locus.
 
 fl is an array which will holds the locales from which each founder gene comes from, since it gets modified within this
 function, it is actually the address of the pointer that points to an array of ints.
 
 flm is an output parameter that holds the amount of memory allocated to fl---memory gets realloc'd to fl in this function
 
 
 n is the total number of individuals in the pedigree
 
 After this function has run, each individual who was marked to be sampled will 
 have entries in his or her MaFounder and PaFounder fields that give the index of
 the founding gene copies there. 
 
 THIS IS THE MULTI-POP VERSION OF THIS FUNCTION, WHICH MEANS IT HAS TO CYCLE OVER ALL INDIVIDUALS
 IN ALL POPULATIONS FOR EACH LOCUS.  SO, THE INPUT IS THE ARRAY OF LOCALES.
 
 */
void NEW_SampleMarkedAtLocusl(/* indiv **M, indiv **F */ int NumPops, Locale **LA, /*int *NM, int *NF,*/ int LastT, /* int no_ibd_t, */ int l, int L, int *nf, int **fl, int *flm, int n)
{
	int t,i,bp;
	int ReInit=1;
	int ma, pa;
	int **indics;
	int *AncsIdx, *AncsGenes;  /* to be used repeatedly by AncestryOfSampledGenes()  */
	int *NM, *NF; 
	indiv **M,**F;
	int no_ibd_t;
	indiv **BigArray;
	int minidx,maxidx,idx;
	
	/* allocate space for indics.  Don't forget to free it at the end of the function */
	indics = (int **)ECA_CALLOC(n+1, sizeof(int *));
	for(i=0;i<=n;i++)  {
		indics[i] = (int *)ECA_CALLOC(2,sizeof(int));
	}
	AncsIdx = (int *)ECA_CALLOC(n+1, sizeof(int *));
	AncsGenes = (int *)ECA_CALLOC(n+1, sizeof(int *));
	
	/* allocate space to the pointers to inds */
	BigArray = (indiv **)ECA_CALLOC(gID+1,sizeof(indiv *));
	
	
	/* cycle over locales to record the earliest --discard-before time.  This we record as
	 no_ibd_t, and we use it later */
	no_ibd_t = LA[0]->P.DiscardBefore;
	for(bp=1;bp<NumPops;bp++)  {
		if(LA[bp]->P.DiscardBefore < no_ibd_t) {
			no_ibd_t = LA[bp]->P.DiscardBefore;
		}
	}
	
	/* cycle over times */
	for(t=no_ibd_t;t<=LastT;t++) {
		minidx = gID + 100;
		maxidx = 0;
		for(bp=0;bp<NumPops;bp++)  { /* cycle over pops to make a linear array of pointers */
			/* define some convenient shorthand variables */
			NF = LA[bp]->NF;
			NM = LA[bp]->NM;
			F = LA[bp]->Females;
			M = LA[bp]->Males;
			
			/* get the pointers assigned to females */
			for(i=0;i<NF[t];i++) {
				idx = F[t][i].idx;
				if(minidx > idx) minidx = idx;
				if(maxidx < idx) maxidx = idx;
				if(BigArray[idx]==NULL)  { /* if this is the first time we have seen this individual */
					BigArray[idx]=&(F[t][i]);
					/*printf("BigArray[idx=%d] gets address of F[t=%d][i=%d] at bp=%d\n",idx,t,i,bp);*/
				}
				/* otherwise, if this one has been sampled more times than the earlier instance, use it */
				else if(F[t][i].N_PreK + F[t][i].N_PostK + F[t][i].N_DurRepro > 
						BigArray[idx]->N_PreK + BigArray[idx]->N_PostK + BigArray[idx]->N_DurRepro) {
					BigArray[idx]=&(F[t][i]);
				}
			}
			/* get the pointers assigned to males */
			for(i=0;i<NM[t];i++) {
				idx = M[t][i].idx;
				if(minidx > idx) minidx = idx;
				if(maxidx < idx) maxidx = idx;
				if(BigArray[idx]==NULL)  { /* if this is the first time we have seen this individual */
					BigArray[idx]=&(M[t][i]);
					/*printf("BigArray[idx=%d] gets address of M[t=%d][i=%d] at bp=%d\n",idx,t,i,bp);*/
				}
				/* otherwise, if this one has been sampled more times than the earlier instance, use it */
				else if(M[t][i].N_PreK +M[t][i].N_PostK +M[t][i].N_DurRepro > 
						BigArray[idx]->N_PreK + BigArray[idx]->N_PostK + BigArray[idx]->N_DurRepro) {
					BigArray[idx]=&(M[t][i]);
				}
			}
		}  /* done with getting the pointers set up */
		
		
		/* cycle over the linear array of pointers */
		for(idx=minidx;idx<=maxidx;idx++)   {
			
			/* then if this individual was marked for sampling at some point... */
			if(BigArray[idx]->N_PreK + BigArray[idx]->N_PostK + BigArray[idx]->N_DurRepro > 0) {
				/* when we call this one, we reset the "indics" when we do the first marked individual.
				 that is what the ReInit=1 takes care of the first time.  Then it gets set to 0. */
				/*printf("Entering AncOfSampGenes for female t=%d i=%d with ReInit=%d\n", t,i,ReInit); */
				AncestryOfSampledGenes(indics, gID,BigArray[idx],nf,&ma,&pa, ReInit, AncsIdx, AncsGenes, no_ibd_t, fl, flm );  /* THIS WILL HAVE TO BE CHANGED SO THAT NO_IBD_T is the minimum --discard-before TIME (it is only used as a 
				 limit for resetting indics if need be,
				 we have to pass FounderLocations into it. */
				ReInit = 0;
				
				/* now allocate space to hold the founder gene indexes of the individual */
				if(BigArray[idx]->MaFounder==NULL)  {
					BigArray[idx]->MaFounder = (int *)ECA_CALLOC(L,sizeof(int));
				}
				if(BigArray[idx]->PaFounder==NULL) {
					BigArray[idx]->PaFounder = (int *)ECA_CALLOC(L,sizeof(int));
				}
				/* then assign the ma and pa values to this individual's ma and pa founder indexes for this locus */
				BigArray[idx]->PaFounder[l] = pa;
				BigArray[idx]->MaFounder[l] = ma;
				
				BigArray[idx]->ReportMe = 1;
			}
			
		}  /* closes loop over the linear array of pointers */
	} /* closes loop over times */
	
	
	/* at the end, free up indics and other memory */
	for(i=0;i<=n;i++)  {
		free(indics[i]);
	}
	free(indics);
	free(AncsIdx);
	free(AncsGenes);
	free(BigArray);
}





/* 
 this cycles over all the individuals in Males and Females from time start to stop, inclusive.  For those that were sampled (have a 
 positive number in the N_PreK or N_PostK) we print out the alleles that their
 maternal and paternal gene copies descended from. 
 
 FounderTypes is an array giving the allelic type of the distinctly numbered founding gene copies.
 
 If FounderTypes is NULL, then it just prints out that distinct number of each founding gene
 copy that that gene is descended from.  Otherwise, it prints out the allelic type.  
 
 L is number of loci.
 
 */
void PrintAllelesOfSampledInds(indiv **M, indiv **F, int *NM, int *NF, int L, int start, int stop, int **FounderTypes, int **FounderLocations)
{
	int i,l,t,k,pre,post,dur;
	
	if(FounderTypes != NULL && FounderLocations != NULL)  {
		fprintf(stderr,"Error!  Somehow PrintAllelesOfSampledInds() was called with both FounderTypes and FounderLocations being NULL.  Exiting...\n");
		exit(1);
	}
	
	for(t=start;t<=stop;t++)  {
		for(i=0;i<NF[t];i++) {
			pre = F[t][i].N_PreK;
			post = F[t][i].N_PostK;
			dur = F[t][i].N_DurRepro;
			if(F[t][i].ReportMe==1 && pre + post + dur > 0) {
				if(FounderTypes==NULL && FounderLocations==NULL) {
					printf("GenotypesByFounderAlleles : ");
				}
				else if(FounderTypes != NULL && FounderLocations==NULL) {
					printf("GenotypesByAllelicType : ");
				}
				else if(FounderTypes == NULL && FounderLocations != NULL) {
					printf("GenotypesByFounderLocales : ");
				}
				printf("%s : ",F[t][i].ID);
				for(k=0;k<pre;k++)  {
					printf("%d ",F[t][i].PreK_ts[k]);
				}
				printf(" : ");
				for(k=0;k<pre;k++)  {
					printf("%d ",F[t][i].SampLocales[BEF_KILL][k]);
				}
				printf(" : ");
				for(k=0;k<post;k++)  {
					printf("%d ",F[t][i].PostK_ts[k]);
				}
				printf(" : ");
				for(k=0;k<post;k++)  {
					printf("%d ",F[t][i].SampLocales[AFT_KILL][k]);
				}
				printf(" : ");
				for(k=0;k<dur;k++)  {
					printf("%d ",F[t][i].DurRepro_ts[k]);
				}
				printf(" : ");
				for(k=0;k<dur;k++)  {
					printf("%d ",F[t][i].SampLocales[WHILE_REPRO][k]);
				}
				printf(" : ");
				for(l=0;l<L;l++)  {
					if(FounderTypes==NULL && FounderLocations==NULL) {
						printf("%d/%d\t",F[t][i].MaFounder[l],F[t][i].PaFounder[l]);
					}
					else if(FounderTypes != NULL && FounderLocations == NULL)  {
						/* note that the the founding alleles are named 1 to NumFounders, so we have to use one minus that
						 to subscript the FounderTypes */
						printf("%d/%d\t",FounderTypes[l][ F[t][i].MaFounder[l]-1 ]+1,
							   FounderTypes[l][ F[t][i].PaFounder[l]-1 ]+1);
					}
					else if(FounderTypes == NULL && FounderLocations != NULL) {
						printf("%d/%d\t",FounderLocations[l][ F[t][i].MaFounder[l]-1 ],
							   FounderLocations[l][ F[t][i].PaFounder[l]-1 ]);
					}
				}
				printf("\n");
			}
		}
		/* then do the males */
		for(i=0;i<NM[t];i++) {
			pre = M[t][i].N_PreK;
			post = M[t][i].N_PostK;
			dur = M[t][i].N_DurRepro;
			if(M[t][i].ReportMe==1 && pre + post + dur > 0) {
				if(FounderTypes==NULL && FounderLocations==NULL) {
					printf("GenotypesByFounderAlleles : ");
				}
				else if(FounderTypes != NULL && FounderLocations==NULL) {
					printf("GenotypesByAllelicType : ");
				}
				else if(FounderTypes == NULL && FounderLocations != NULL) {
					printf("GenotypesByFounderLocales : ");
				}
				printf("%s : ",M[t][i].ID);
				for(k=0;k<pre;k++)  {
					printf("%d ",M[t][i].PreK_ts[k]);
				}
				printf(" : ");
				for(k=0;k<pre;k++)  {
					printf("%d ",M[t][i].SampLocales[BEF_KILL][k]);
				}
				printf(" : ");
				for(k=0;k<post;k++)  {
					printf("%d ",M[t][i].PostK_ts[k]);
				}
				printf(" : ");
				for(k=0;k<post;k++)  {
					printf("%d ",M[t][i].SampLocales[AFT_KILL][k]);
				}
				printf(" : ");
				for(k=0;k<dur;k++)  {
					printf("%d ",M[t][i].DurRepro_ts[k]);
				}
				printf(" : ");
				for(k=0;k<dur;k++)  {
					printf("%d ",M[t][i].SampLocales[WHILE_REPRO][k]);
				}
				printf(" : ");
				for(l=0;l<L;l++)  {
					if(FounderTypes==NULL && FounderLocations==NULL) {
						printf("%d/%d\t",M[t][i].MaFounder[l],M[t][i].PaFounder[l]);
					}
					else if(FounderTypes != NULL && FounderLocations==NULL)   {
						printf("%d/%d\t",FounderTypes[l][ M[t][i].MaFounder[l]-1 ]+1,
							   FounderTypes[l][ M[t][i].PaFounder[l]-1 ]+1);
					}
					else if(FounderTypes == NULL && FounderLocations != NULL) {
						printf("%d/%d\t",FounderLocations[l][ M[t][i].MaFounder[l]-1 ],
							   FounderLocations[l][ M[t][i].PaFounder[l]-1 ]);
					}
				}
				printf("\n");
			}
		}
	}
	
}


/* 
 Given the number of founding gene copies at each locus, and the allele frequencies at 
 each of those loci, this allocates memory to, and returns, an array of the types of
 of each of those founders.
 
 n -- number of loci
 
 na -- number of alleles at each locus 
 
 LA -- array of locales used to get to the allele freqs at the n loci
 
 nf -- number of founding gene copies in the pedigree at each locus 
 
 FL -- the locations of the founder alleless  
 */
int **DrawFounderTypes(int n, int *na, Locale **LA, int *nf, int **FL)
{
	int i,k;
	int **T;
	int bp;
	
	/* allocate memory straight off */
	T = (int **)ECA_CALLOC(n, sizeof(int *));
	for(i=0;i<n;i++)  {
		/* allocate memory for each locus */
		T[i] = (int *)ECA_CALLOC(nf[i],sizeof(int));
		
		printf("FounderTypes : Number_of_Founders_at_Locus : %d : %d\n",i+1,nf[i]);
		
		
		printf("FounderTypes : Founders_locales_at_locus : %d : ",i+1);
		for(k=0;k<nf[i];k++)  {
			printf("%d ",FL[i][k]);
		}
		printf("\n");
		
		
		printf("FounderTypes : Allelic_types_assigned_to_founders_at_locus : %d : ",i+1);
		/* then assign the types to each founder gene copy */
		for(k=0;k<nf[i];k++)  {
			bp = FL[i][k];  /* record the founders population of origin in a convenient variable */
			
			T[i][k] = IntFromProbsRV(LA[bp]->S.AlleFreqs[i], 0, na[i]);
			printf("%d ",T[i][k]+1);
		}
		printf("\n");
		
		
	}
	
	/* that's it!  now we just return T */
	return(T);
	
}






void PrintCensusSizes(int t, int MA, int *NFR, int *NMR, const char *string, int age_start, int whichpop) 
{
	int at,age;
	
	printf("%sAGES :           %d  :  %d :", string,t,whichpop);
	for(age=age_start;age<=MA;age++) {
		printf("\t%d",age);
	}
	printf("\n");
	printf("%sCOUNTS : MALES  : %d  :  %d :", string,t,whichpop);
	for(age=age_start;age<=MA;age++) {
		at = t - age;
		printf("\t%d",NMR[at]);
	}
	printf("\n");
	printf("%sCOUNTS : FEM    : %d  :  %d :",string, t,whichpop);
	for(age=age_start;age<=MA;age++) {
		at = t - age;
		printf("\t%d",NFR[at]);
	}
	printf("\n");
}

/* This just returns what the cohort size should be for the current year */
int ReturnCohortSize(PopSizePars *Pop,int t) 
{
	if(strcmp(Pop->CohortSizeMode,"const")==0) {
		return(Pop->ConstSize);
	}
	else if(strcmp(Pop->CohortSizeMode,"var")==0) {
		return(Pop->CohortSizes[t]);
	}
	else {
		fprintf(stderr,"Error! Unrecognized CohortSizeModel in function ReturnCohortSize().\nExiting\n");
		exit(1);
	}
	return(-999999);
}


/*
 This returns a 1 if there are more than zero potentially reproducing
 males and females at time t.  If there are zero, then it returns a zero and it also
 prints out a ZEROPOP message.
 */
int PopSizeNonZero(int t, int MA, int *NFR, int *NMR)
{
	int age,at,nm=0,nf=0;
	
	for(age=1;age<=MA;age++)  {
		at = t - age;
		nf += NFR[at];
		nm += NMR[at];   /* here we also take advantage of the loop to count up the number of 
		 males that can be contributing */
	}
	if(nm>0 && nf>0) 
		return(1);
	else {
		if(nm==0) {
			printf("ZEROPOPSIZE : %d : Zero reproducing Males\n",t);
		}
		if(nf==0) {
			printf("ZEROPOPSIZE : %d : Zero reproducing Females\n",t);
		}
	}
	return(0);
}





void MarkSampledNewborns(SampPars S, int *NM, int *NF, indiv **Males, indiv **Females, int whichpop)
{	
	int i,j;
	long **SSM,**SSF;
	
	/* randomly draw the subscripts of all those juveniles that will be sampled
	 by just making a random permutation of all of the subscripts, and using only
	 the first n (where n is the number of individuals you will be sampling  */
	SSM = (long **)ECA_CALLOC(S.nc,sizeof(long *));
	SSF = (long **)ECA_CALLOC(S.nc,sizeof(long *));
	for(i=0;i<S.nc;i++)  {
		SSM[i] = (long *)ECA_CALLOC(NM[ S.ts[i] ], sizeof(long));
		SSF[i] = (long *)ECA_CALLOC(NF[ S.ts[i] ], sizeof(long));
		
		/* fill up the arrays with all the possible numbers */
		for(j=0;j<NM[ S.ts[i] ]; j++)  {
			SSM[i][j] = j;
		}
		for(j=0;j<NF[ S.ts[i] ]; j++)  {
			SSF[i][j] = j;
		}
		/* then permute all those */
		genprm(SSM[i],NM[ S.ts[i] ]);
		genprm(SSF[i],NF[ S.ts[i] ]);
		
		/* then print out who we are going to be sampling, and whether or not we got all of them, and while we are doing 
		 that we will allocate memory to the MaFounder and PaFounder fields for each of them.  */
		printf("NEWBORN_SAMPLING : %d : FEMALES : AskedFor : %d : Available : %d :  ",S.ts[i],S.NumF[i], NF[ S.ts[i] ]);
		for(j=0;j<ECA_MIN(NF[ S.ts[i] ], S.NumF[i]);j++)  {
			MarkIndAsSampled( &(Females[S.ts[i]][ SSF[i][j]]), 0, S.ts[i], whichpop);
			printf("%s ",Females[S.ts[i]][ SSF[i][j]].ID);
			
		}
		printf("\n");
		printf("NEWBORN_SAMPLING : %d : MALES : AskedFor : %d : Available : %d :  ",S.ts[i],S.NumM[i], NM[ S.ts[i] ]);
		for(j=0;j<ECA_MIN(NM[ S.ts[i] ], S.NumM[i]);j++)  {
			MarkIndAsSampled( &(Males[S.ts[i]][ SSM[i][j]]), 0, S.ts[i], whichpop);
			printf("%s ",Males[S.ts[i]][ SSM[i][j]].ID);
		}
		printf("\n");
		
	}
}	


/*
 Initialize the male and female cohorts in M and F.  
 MA is MaxAge
 NF is number of females for cohorts born from 0 to MA (but will be used for later gens too)
 NM is number of males for cohorts born from  0 to MA  (but will be used for later gens too)
 T is the number of generations forward to allow space for 
 */
void InitCohorts(int MA, int *NF, int *NM, int T, indiv ***M, indiv ***F, PopPars P, int birthplace, int *MSpace, int *FSpace)
{
	int i,j;
	
	/* allocate to all the years we will need */
	*M = (indiv **)ECA_CALLOC(T+MA+1,sizeof(indiv *));
	*F = (indiv **)ECA_CALLOC(T+MA+1,sizeof(indiv *));
	
	/* cycle over the MA year classes and initialize */
	for(i=0;i<MA;i++)  {
		(*M)[i] = (indiv *)ECA_CALLOC(NM[i],sizeof(indiv));
		/* record how much memory is allocated to that row of males */
		MSpace[i] = NM[i];
		
		for(j=0;j<NM[i];j++)  {
			FillIndiv( &((*M)[i][j]), NULL, NULL, &gID, i, j,'M',birthplace,P);
			PrintAndDiscardPedigree( &((*M)[i][j]), i, P.DiscardBefore,P.TossAll, P.MaxAge);
		}
		
		(*F)[i] = (indiv *)ECA_CALLOC(NF[i],sizeof(indiv));
		/* record how much memory is allocated to that row of females */
		FSpace[i] = NF[i];
		for(j=0;j<NF[i];j++)  {
			FillIndiv(&((*F)[i][j]), NULL, NULL, &gID, i, j,'F',birthplace,P);
			PrintAndDiscardPedigree( &((*F)[i][j]), i,P.DiscardBefore,P.TossAll, P.MaxAge);
		}
	}
	
}


/*
 THIS IS A MULTIPOP VERSION THAT KEEPS THE gID's IN GOOD ORDER, WHICH MAKES
 THINGS MUCH BETTER WHEN IT COMES TIME TO SIMULATE GENETIC DATA ON THE 
 PEDIGREE.
 
 Initialize the male and female cohorts in M and F.  
 MA is MaxAge
 NF is number of females for cohorts born from 0 to MA (but will be used for later gens too)
 NM is number of males for cohorts born from  0 to MA  (but will be used for later gens too)
 T is the number of generations forward to allow space for 
 */
void MultiPop_InitCohorts(Locale **LA, int NumPops, int T)           
{
	int i,j,bp;
	int MA;
	
	for(bp=0;bp<NumPops;bp++)  {
		MA = LA[bp]->MA;
		/* allocate to all the years we will need */
		LA[bp]->Males = (indiv **)ECA_CALLOC(T+MA+1,sizeof(indiv *));
		LA[bp]->Females = (indiv **)ECA_CALLOC(T+MA+1,sizeof(indiv *));
	}
	/* cycle over the MA year classes and initialize */
	for(i=0;i<MA;i++)  {
		for(bp=0;bp<NumPops;bp++) {
			
			LA[bp]->Males[i] = (indiv *)ECA_CALLOC(LA[bp]->NM[i],sizeof(indiv));
			/* record how much memory is allocated to that row of males */
			LA[bp]->MSpace[i] = LA[bp]->NM[i];
			
			for(j=0;j<LA[bp]->NM[i];j++)  {
				FillIndiv( &(LA[bp]->Males[i][j]), NULL, NULL, &gID, i, j,'M',bp,LA[bp]->P);
				PrintAndDiscardPedigree( &(LA[bp]->Males[i][j]), i, LA[bp]->P.DiscardBefore,LA[bp]->P.TossAll, LA[bp]->P.MaxAge);
			}
			
			LA[bp]->Females[i] = (indiv *)ECA_CALLOC(LA[bp]->NF[i],sizeof(indiv));
			/* record how much memory is allocated to that row of females */
			LA[bp]->FSpace[i] = LA[bp]->NF[i];
			for(j=0;j<LA[bp]->NF[i];j++)  {
				FillIndiv(&(LA[bp]->Females[i][j]), NULL, NULL, &gID, i, j,'F',bp,LA[bp]->P);
				PrintAndDiscardPedigree( &(LA[bp]->Females[i][j]), i,LA[bp]->P.DiscardBefore,LA[bp]->P.TossAll, LA[bp]->P.MaxAge);
			}
		} /* close loop over pops */
		
	}  /* close loop over times */
	
}


/* given an array I of N indivs, go through and kill each living one
 off with probability s.  For each one you kill, decrement the number remaining
 NR.  We pass t in so we can print an informative message aout when they die.
 */
void KillOff(indiv *I, int N, int *NR, double s, int t)
{
	int i;
	
	for(i=0;i<N;i++)  {
		if(I[i].dead == 0) {
			if( ranf() > s) {
				I[i].dead = 1;
				(*NR)--;
								/* note that the age and time of death here used to reflect dying at the end of the previous time
								   step, rather than dying at the beginning of the current time step. But I have changed that
								   to be the actual age in the actual time step. */
								printf("KILLING:  Just Killed individual : %s  at time %d at age %d\n",I[i].ID,t,t-I[i].t);
				
			}
		}
	}
}


/* 
 if SizeNotCount!=0:
	given that you have N individuals that are going to be put into N/n + (N%n>0) groups of n 
	individuals (with the last one being of size N%n, if that quantity is greater
	than 0), randomly permute the individuals into these groups which are numbered 0, 1, 2,...
 if SizeNotCount==0:
	take the N individuals and put them into n groups according to a multinomial dsn with N trials
	and n cells each with prob 1.0/n.  Then permute them.  return the indexes of the cells of each of the
	n individuals.
 
	MaxGroup is an output variable that tells us the max number of individuals in any spawning group (to be used for
	memory allocation purposes).  NumGroup tells us how many groups there are.
 */ 
int *ReturnLinearSpawnGroupIds(int N, int n, int SizeNotCount, int *MaxGroup, int *NumGroup)
{
	int i,j,k;
	int *ret=(int *)calloc(N, sizeof(int));
	int *x;
	double *p;
	int mx=0;
	
	if(SizeNotCount) {
		for(i=0;i<N;i++)  {
			ret[i] = i/n;
		}
		mx=n;
		*NumGroup = N/n + (N%n>0);
	}
	else {
		x = (int *)calloc(n,sizeof(int));
		p = (double *)calloc(n, sizeof(double));
		for(i=0;i<n;i++) p[i]=1.0/(double)n;
		MultinomialRV(N,p,n,x);
		for(i=0,k=0;i<n;i++)  {
			if(x[i]>mx) {
				mx=x[i];
			}
			for(j=0;j<x[i];j++) {
				ret[k++]=i;
			}
		}
		free(x);
		free(p);
		*NumGroup = n;
	}
	
	genprmINT(ret,N);
	*MaxGroup=mx;
	return(ret);
}


/* make offspring for the cohort at time t.  The expected number
 to make is N.  The newborns go into M[t] and F[t].  THis is where the reproduction-mediated 
 death takes place, too.  */
void MakeBabies(indiv **F, indiv **M, int t, int MA, PopPars P, int N, int *NFR, int *NF, int *NM, int *NMR, SampPars S, int birthplace, int *MSpace, int *FSpace, Locale **LA, int NumPops)
{
	int age,at,i,m,f,x,j,n;
	double wt=0.0;
	double mu[600];
	int fspaces = (int)(P.OffspringAllocOverload*N*(1.0-P.SexRatio)+1);  /* amount of space allocated to new females */
	int mspaces = (int)(P.OffspringAllocOverload*N*P.SexRatio+1);
	indiv **MalePointers, **tempFathers;
	double *MaleProbs;
	int NumMaleParents=0, NumMalesContributing=0;
	int **FemNumOffs; /* array subscripted by [age][idx] for the number of offspring for each female in a gCohortSizesRandom==0 situation */
	double *FemDirPars;  /* array for the dirichlet (or multinomial) pars for females */
	int *FemNumOffsLin;  /* array of number of offspring for each female in linear format (so it can be output from
	 the multinomial of cmd functions */
	int *FemSpawnGroupsLin;  /* linear array telling which spawner group each reproducing female belongs to */
	int TotNumFemales=0;  /* to count the total number of contributing females */
	double dirsum=0.0;  /* for normalizing some pars */
	double MaleProbSum=0.0;
	int *FemDiedPostReprod, *MaleDiedPostReprod, NumLocalMales,NumFemsReproducing;
	int **SpawnGroup2D;  /* array subscripted by [age][idx] for the spawner group of each female in the SpawnGroupSize or SpawnGroupCnt situation */
	int *MaleSpawnGroupsLin; /* linear array telling us which spawning group the corresponding males belong to */
	indiv ***MalePointersSG;  /* indexed by spawner group and then by male within that spawner group */
	double **MaleProbsSG;  /* indexed by spawner group and then male within that spawner group */
	int *NumMalesContributingSG; /* indexed by spawner group */
	int MaxMalesInSpawnerGroup;  /* will store the largest number of males in a spawner group if we are doing SpawnerGroupCnt*/
	int NumMaleSpawnGroups; /* if we are doing SpawnerGroupSize, then this will tell us how many groups that makes */
	
	
	/* 	printf("\nVERBIAGE :   Entering MakeBabies : t=%d     N=%d\n",t,N);  */
	
	/* allocate lots of memory to F[t] and M[t].  Make it ten times the expected number, 
	 and maybe we can reallocate later.   THIS IS HORRIBLY INEFFICIENT ON MEMORY.  I SHOULD
	 FIX THIS AT SOME POINT IN THE FUTURE */
	M[t] = (indiv *)ECA_CALLOC(mspaces, sizeof(indiv));
	F[t] = (indiv *)ECA_CALLOC(fspaces, sizeof(indiv));
	MSpace[t] = mspaces;
	FSpace[t] = fspaces;
	/* initialize the number of males and females at t */
	NF[t] = 0;
	NM[t] = 0;
	NFR[t]=0;
	NMR[t]=0;
	
	/* allocate space to count how many individuals of different ages crump due to reproduction-related death */
	FemDiedPostReprod = (int *)ECA_CALLOC(MA+1,sizeof(int));
	MaleDiedPostReprod = (int *)ECA_CALLOC(MA+1,sizeof(int));
	
	/* first count up the expected number of offspring if each female contributes, on average,
	 a number of offspring equal to P.Fasrf[age], and also weighted by the prob that one of these
	 individuals will reproduce at all. */
	for(age=1;age<=MA;age++)  {
		at = t - age;
		wt += NFR[at] * P.Fasrf[age] * P.FPR[age];
		NumMaleParents += NMR[at];   /* here we also take advantage of the loop to count up the number of 
		 males that can be contributing */
	}
	NumLocalMales = NumMaleParents;
	
	/* IMPLEMENTING THE SNEAKY MALE THING */
	/*  If there is the possibility of sneaky males, we count them up too */
	if(P.SneakyMProb!=NULL)  {
		for(i=0;i<NumPops;i++) {
			if(i != birthplace && P.SneakyMProb[i]>0.0) {
				for(age=1;age<=MA;age++)  {
					at = t - age;
					NumMaleParents += LA[i]->NMR[at];
				}
			}
		}
	}
	
	/* allocate memory for the array of possible males and for the assignment pointers to
	 male parents  */
	MaleProbs = (double *)ECA_CALLOC( NumMaleParents, sizeof(double));
	MalePointers = (indiv **)ECA_CALLOC( NumMaleParents, sizeof(indiv *));
	tempFathers = (indiv **)ECA_CALLOC( N, sizeof(indiv *));
	
	/* then determine the mean number of offspring per female for each age class */
	mu[0] = 0.0;
	wt = (double)N/wt;
	/* printf("VERBIAGE :  Mean Num Offspring per Female  (NumFemalesRemaining) and [ProbOfReproducing] of Different Ages  : "); */
	for(age=1;age<=MA;age++)  {
		at = t - age;
		mu[age] = wt *  P.Fasrf[age];
		
		/* and count the total number of females (dead or alive) while we are at it */
		TotNumFemales += NF[at];
	}
	printf("\n");
	
	
	
	
	/* here, in order to accommodate the SpawnGroupStuff, we need to decide now who is going to be reproducing
	 and who isn't and count up the total number of those indivs.  We store the information in the indiv_struct. */
	NumFemsReproducing=0;
	for(age=1;age<=MA;age++)  {
		at = t - age;
		for(i=0;i<NF[at];i++)  {
			if(F[at][i].dead == 0  && ranf() < P.FPR[age] ) { 
				F[at][i].IsReproducing = 1;
				NumFemsReproducing++;
			}
			else {
				F[at][i].IsReproducing = 0;
			}
		}
	}
	
	
	
	/* if we are doing spawning groups, then make the list for the females and translate it into an array
	 subscripted by [age][i].  */
	if(P.SpawnGroupSize!=-1) { int dummy;
		FemSpawnGroupsLin = ReturnLinearSpawnGroupIds(NumFemsReproducing,P.SpawnGroupSize,1,&dummy,&dummy);
	}
	if(P.SpawnGroupCnt!=-1) {int dummy;
		FemSpawnGroupsLin = ReturnLinearSpawnGroupIds(NumFemsReproducing,P.SpawnGroupCnt,0,&dummy,&dummy);
	}
	if(P.SpawnGroupCnt!=-1 || P.SpawnGroupSize!=-1) {
		SpawnGroup2D = (int **)ECA_CALLOC(MA+1,sizeof(int *));
		for(j=0,age=1;age<=MA;age++)  {
			at = t - age;
			SpawnGroup2D[age] = (int *)ECA_CALLOC(NF[at], sizeof(int));
			for(i=0;i<NF[at];i++) {
				if(F[at][i].IsReproducing==1) {
					SpawnGroup2D[age][i] = FemSpawnGroupsLin[j++];
					printf("SPAWNING_GROUPS:  %d : %s : %d  \n",t,F[at][i].ID,SpawnGroup2D[age][i]);
				}
				else {
					SpawnGroup2D[age][i] = -1;
				}
				//printf("TESTING_SPAWN_GROUP_FEM : %s   isRepro= %d  SpawnGroup= %d\n",F[at][i].ID,F[at][i].IsReproducing,SpawnGroup2D[age][i]);
			}
		}
	}
	
	/********************************************************************************/	
	
	/* and here, if gCohortSizesRandom==0 we are doing fixed cohort sizes, so we are doing either 
	 poisson or negative binomial, but conditional on the total number of offspring in the cohort,
	 so we need to do it as binomial or beta-binomial sampling.  What we do here is make an array
	 of how many offspring each female has. */
	if(P.CohortSizesRandom==0) {
		/* first deal with allocating the necessary memory */
		FemDirPars = (double *)ECA_CALLOC(TotNumFemales, sizeof(double));
		FemNumOffsLin = (int *)ECA_CALLOC(TotNumFemales, sizeof(int));
		FemNumOffs = (int **)ECA_CALLOC(MA+1,sizeof(int *));
		/*  MODDDD: ADD ALLOCATION TO FEMSPAWNGROUPS HERE */
		for(age=1;age<=MA;age++)  {
			at = t - age;
			FemNumOffs[age] = (int *)ECA_CALLOC(NF[at], sizeof(int));
		}
		
		/* then assign to each female the weight she will get.  ultimately these
		 will turn out to be dirichlet or multinomial parameters */
		j=0;
		for(age=1;age<=MA;age++)  {
			at = t - age;
			for(i=0;i<NF[at];i++)  {
				if(F[at][i].IsReproducing == 1 ) {   /* if they are not dead yet and they are chosen to reproduce this year */
					/* put the mean number of offspring that they should produce into the FemDirPars array */
					FemDirPars[j] = mu[age];
					dirsum += FemDirPars[j];
					/* sample them while reproducing, if indicated */
					if(ranf() < S.dur_repro_fem_samp && S.dur_fem_samp_years[t]) {
						MarkIndAsSampled(&(F[at][i]), WHILE_REPRO, t,birthplace);
					}
					
					/* kill them off after mating if indicated */
					if(ranf() < P.FPDPR[age]) {
						F[at][i].dead = 1;
						FemDiedPostReprod[age]++;
					}
					
					if(gRunTests==1) {printf("TESTING : FemMeanNumOffs : %d  %s   %f\n",t,F[at][i].ID,mu[age]);  }
				}
				else {  /* otherwise, put 0 into that array */
					FemDirPars[j] = 0.0;
				}
				j++;
			}
		}
		
		
		
		
		/* and now we are ready to make those parameters work for us.  There are two cases: (1) Poisson (P.Fcv==1.0) and 
		 Negative Binomial (P.Fcv<1).  Which are implemented here as Multinomial and CMD, respectively */
		if(N>0) {
			/* RIGHT HERE I SHOULD IMPLEMENT A SAMPLING WITHOUT REPLACEMENT TO DO THE BINARY CASE...LATER */
			if(P.Fcv==1.0) {
				/* first, make the dir pars multinomial pars by making them sum to one */
				for(i=0;i<j;i++)  {
					FemDirPars[i] /= dirsum;
				}
				/* then draw a multinomial RV */  
				D_MultinomialRV(N,FemDirPars,j,FemNumOffsLin);
			}
			else if(P.Fcv<1.0 && P.Fcv>0.0) {
				/* first, make the dir pars reflect the FCVrs variable  */
				for(i=0;i<j;i++)  {
					FemDirPars[i] *= P.Fcv / (1.0 - P.Fcv);
				}
				/* then draw a CMD RV */  
				CompoundDirMultRV(FemDirPars,N,j,FemNumOffsLin);
			}
			else {
				for(i=0;i<j;i++)  {
					FemDirPars[i] = 0.0;
				}
				printf("ERROR:  P.Fcv less than zero or greater than 1.0.  Exiting...\n\n");
				exit(1);
			}
			/* and here we translate the FemNumOffsLin back to the 2-D array by age and idx, and it will be ready for use later */
			j=0;
			for(age=1;age<=MA;age++)  {
				at = t - age;
				for(i=0;i<NF[at];i++)  {
					FemNumOffs[age][i] = FemNumOffsLin[j];
					j++;	
				}
			}
			
		}
		else {  /* if the cohort size is zero the do nothing, really */
			printf("ZERO_COHORT : %d : No offspring produced, but individuals may try to reproduce and hence may die or be sampled at that time\n",t);
		}
	}
	/*******************************************************************/	
	
	
	
	
	
	
	
	
	
	
	/* compute the weights for the males.  These will be simply their normalized relative age-specific reproductive "prowesses".
	 If we have some extra-bimomial varaince in male reproductive success, these weights will get modified by using a weighted
	 version of them as a Dirichlet dsn parameter to make their weights a Dirichlet r.v.  Note that we only include non-dead
	 individuals who have been chosen to reproduce in this array of weights.   */
	for(j=0,age=1;age<=MA;age++)  {
		at = t - age;
		for(i=0;i<NM[at];i++) { 
			if(M[at][i].dead == 0  && ranf() < P.MPR[age]) {
				MalePointers[j] = &(M[at][i]);
				/* here we require some way to use the ron variable if we had stickiness across years of
				 male reproductive prowess, but I'll deal with that later */
				MaleProbs[j] = P.Masrp[age];
				
				/* this is the simple step to incorporate the extra variance in male reproductive success which is a maintained across years
				 by the males */
				if(P.sticky_m_alpha>0) {
					MaleProbs[j] *= MalePointers[j]->sticky_birth_rate;
				}
				
				MaleProbSum += MaleProbs[j];		
				j++;
				/* sample them while reproducing, if indicated */
				if(ranf() < S.dur_repro_male_samp && S.dur_male_samp_years[t]) {
					MarkIndAsSampled(&(M[at][i]), WHILE_REPRO, t,birthplace);
				}
				if(ranf() < P.MPDPR[age]) {
					M[at][i].dead = 1;
					MaleDiedPostReprod[age]++;
				}
			}	
		}
	}
	/*  If there is the possibility of sneaky males, we factor them in too*/
	if(P.SneakyMProb!=NULL)  {
		for(n=0;n<NumPops;n++) {
			if(n != birthplace && P.SneakyMProb[n]>0.0) {
				for(age=1;age<=MA;age++)  {float temprando;
					at = t - age;
					for(i=0;i<LA[n]->NM[at];i++) { 
						temprando = ranf();
						if(LA[n]->Males[at][i].dead == 0  && temprando < (LA[n]->P.MPR[age] * P.SneakyMProb[n]) ) {
							printf("SNEAKY_MALE_IN_ARRAY : %d : %d : %s : %d : %f\n",t,birthplace,LA[n]->Males[at][i].ID,n,temprando);
							MalePointers[j] = &(LA[n]->Males[at][i]);
							/* here we require some way to use the ron variable if we had stickiness across years of
							 male reproductive prowess, but I'll deal with that later */
							MaleProbs[j] = P.Masrp[age];
							MaleProbSum += MaleProbs[j];		
							j++;
							/* sample them while reproducing, if indicated */
							if(ranf() < S.dur_repro_male_samp && S.dur_male_samp_years[t]) {
								MarkIndAsSampled(&(LA[n]->Males[at][i]), WHILE_REPRO, t,birthplace);
							}
							/* but don't kill them off */
							/*if(ranf() < P.MPDPR[age]) {
							 M[at][i].dead = 1;
							 MaleDiedPostReprod[age]++;
							 }*/
						}	
					}  /* closes loop over the males */
					
				}  /* closes loop over ages */
			}
		}
	}
	
	
	
	NumMalesContributing = j;
	/* now we do the normalizing.  If the male reproductive variance ratio is 1 we just are going to
	 be doing multinomial sampling, so we normalize MaleProbs to actually be probabilities */
	if(P.MCVrs == 1.0)  {
		for(j=0;j<NumMalesContributing;j++)  {
			MaleProbs[j] /= MaleProbSum;
			
			/* here we print out some stuff to test for the right distribution of male offspring number */
			if(gRunTests==1) {printf("TESTING : MaleMeanNumOffs : %d  %s   %f\n",t,MalePointers[j]->ID,MaleProbs[j] * N);  }
		}
	}
	else if(P.MCVrs < 1.0 && P.MCVrs > 0.0) {  double *TempMaleDirPars;
		/* in this case, the weight they get is their expected number of offspring times
		 the theta/(1-theta) factor that makes it an appropriately overdispersed 
		 Dirichlet parameter, then we make a probability vector from that */
		TempMaleDirPars = (double *)ECA_CALLOC(NumMalesContributing, sizeof(double));
		for(j=0;j<NumMalesContributing;j++)  {  
			TempMaleDirPars[j] =  MaleProbs[j] * (double)N * (P.MCVrs  / (1.0 - P.MCVrs)) / MaleProbSum;
			/* here we print out some stuff to test for the right distribution of male offspring number */
			if(gRunTests==1) {printf("TESTING : MaleMeanNumOffs : %d  %s   %f\n",t,MalePointers[j]->ID,MaleProbs[j]/MaleProbSum * N);  }
		}
		if(N>0) { DirichletRV(TempMaleDirPars,NumMalesContributing,MaleProbs);}
		free(TempMaleDirPars);
	}
	else {
		printf("ERROR: Male variance ratio of reproductive success >1 or <= 0.   == %f\n\n",P.MCVrs);
		exit(1);
	}
	
	
	/* now that everything has been computed for these males, we can put them into separate spawning groups if that is indicated */
	if(P.SpawnGroupSize!=-1) {
		MaleSpawnGroupsLin = ReturnLinearSpawnGroupIds(NumMalesContributing,P.SpawnGroupSize,1,&MaxMalesInSpawnerGroup,&NumMaleSpawnGroups);
	}
	if(P.SpawnGroupCnt!=-1) {
		MaleSpawnGroupsLin = ReturnLinearSpawnGroupIds(NumMalesContributing,P.SpawnGroupCnt,0,&MaxMalesInSpawnerGroup,&NumMaleSpawnGroups);
	}
	if(P.SpawnGroupCnt!=-1 || P.SpawnGroupSize!=-1) {
		/* allocate a bunch of memory to hold each spawner group */
		MalePointersSG = (indiv ***)calloc(NumMaleSpawnGroups,sizeof(indiv **));
		MaleProbsSG = (double **)calloc(NumMaleSpawnGroups,sizeof(double *));
		NumMalesContributingSG = (int *)calloc(NumMaleSpawnGroups,sizeof(int));
		for(i=0;i<NumMaleSpawnGroups;i++)  {
			MalePointersSG[i] = (indiv **)calloc(MaxMalesInSpawnerGroup,sizeof(indiv *));
			MaleProbsSG[i] = (double *)calloc(MaxMalesInSpawnerGroup,sizeof(double));
		}
		
		/* then cycle over the males in MalePointers and according to MaleSpawnGroupsLin point the corresponding elements
		 in MalePointersSG to them and put the right values in MaleProbsSG, and count how many are in each spawn group */
		for(i=0;i<NumMalesContributing;i++)  { int g; int k;
			g = MaleSpawnGroupsLin[i];
			k = NumMalesContributingSG[g]++;
			MalePointersSG[g][k] = MalePointers[i];
			MaleProbsSG[g][k] = MaleProbs[i];
			
			/* while we are at it, we will print out the SpawnGroup information */
			printf("SPAWNING_GROUPS:  %d : %s : %d  \n",t,MalePointers[i]->ID,MaleSpawnGroupsLin[i]);
		}
		
		/* now let's put some code in for testing */
		/*for(i=0;i<NumMaleSpawnGroups;i++)  {
			printf("MALE_SPAWN_GROUP: %d : %d : %d :",t,i,NumMalesContributingSG[i]);
			for(j=0;j<NumMalesContributingSG[i];j++) {
				printf(" %s ",MalePointersSG[i][j]->ID);
			}
			printf("\n");
		}  */
		
		/* and here we normalize the MaleProbsSG's and then make them cumulative probs */
		for(i=0;i<NumMaleSpawnGroups;i++)  { double normo;
		  normo = 0.0;
			for(j=0;j<NumMalesContributingSG[i];j++) {
				normo += MaleProbsSG[i][j];
			}
			for(j=0;j<NumMalesContributingSG[i];j++) {
				MaleProbsSG[i][j]/=normo;
			}
			for(j=1;j<NumMalesContributingSG[i];j++) {
				MaleProbsSG[i][j] += MaleProbsSG[i][j-1];
			}
		}
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/* now, make MaleProbs a cumulative array so that it works with the binary search routine I use */
	for(j=1;j<NumMalesContributing;j++)  {
		MaleProbs[j] += MaleProbs[j-1];
	}
	
	
	/* HEY!  HERE IS A COOL AND QUICK AND EASY WAY THAT I COULD IMPLEMENT EXACT NUMBER OF OFFSPRING (AT LEAST
	 FOR BINOMIAL AND NEGATIVE BINOMIAL OFFSPRING NUMBER DISTRIBUTIONS).  I COULD CREATE AN ARRAY INDEXED BY
	 AGE AND BY INDEX OF FEMALE HERE THAT INCLUDED THE NUMBER OF OFFSPRING FOR EACH FEMALE (THIS COULD BE
	 DONE BY MULTINOMIAL OR COMPOUND-MULTINOMIAL SAMPLING), THEN THE ELEMENTS OF THAT ARRAY COULD BE USED 
	 TO SET x DOWN BELOW INSTEAD OF USING NumOffspringRV() FOR THAT.  WE JUST HAVE TO PUT A CONDITIONAL
	 IN THERE TELLING WHETHER TO USE "ABOLUTE" OR "APPROXIMATE" COHORT SIZES.  I now have that implemented! */  
	/* then cycle over the females and make them have offspring */
	if(NumMalesContributing==0) {
		printf("ZEROPOP : %d : Individuals alive, but apparently no males were chosen to reproduce\n",t); 
	}
	else for(age=1;age<=MA;age++)  {
		at = t - age;
		for(i=0;i<NF[at];i++) {  int HasNoKids; /* cycle over all females */
			
			/* do some stuff for the Reproductive Inhibition */
			if(P.ReproInhib > 0 && F[at][i].OldestOffsp != NULL &&  F[at][i].OldestOffsp->dead==0   && t - F[at][i].OldestOffsp->t <= P.ReproInhib) {
				HasNoKids=0;
			}
			else {
				HasNoKids=1;
			}
			if(P.CohortSizesRandom==0) {  /* if we are doing the fixed cohort sizes, then we have already chosen and we just assign it */
				x = FemNumOffs[age][i];
			}
			else if(F[at][i].IsReproducing==1 && HasNoKids) {   /* if they are not dead yet and if they are chosen to reproduce this year */
				x = NumOffspringRV(mu[age],P.Fcv, P.OffsDsn);  /* choose number of offspring */
				/* sample them while reproducing, if indicated */
				if(ranf() < S.dur_repro_fem_samp &&  S.dur_fem_samp_years[t]) {
					MarkIndAsSampled(&(F[at][i]), WHILE_REPRO, t,birthplace);
				}
				/* rub them out, afterward, if appropriate */
				if(ranf() < P.FPDPR[age]) {
					F[at][i].dead = 1;
					FemDiedPostReprod[age]++;
				}
				
				if(gRunTests==1) {printf("TESTING : FemMeanNumOffs : %d  %s   %f\n",t,F[at][i].ID,mu[age]);  }
				
			}
			else {
				x = 0;  /* otherwise, they get no offspring */
			}
			
			/* then assign number of different sexes to the offspring */
			if(x>0) {  /* we have to make sure there are fathers here for them, too */
				m = (int)ignbin((long)x,(float)P.SexRatio);
				f = x - m;
				/* now, we get a list of fathers and store them in tempFathers */
				/* RIGHT HERE IS WHERE I CAN PUT IN SOME CODE TO FORCE MALES TO BE MEMBERS OF A PARTICULAR SPAWNING GROUP */
				/* if we NumMalesContributingSG[female's spawn group] is greater than 0 we are good to go, otherwise we just
				 set m=f=0 and go on to the next female */
				if(P.SpawnGroupCnt!=-1 || P.SpawnGroupSize!=-1) { int g;
					g = SpawnGroup2D[age][i];
					if(g<NumMaleSpawnGroups && NumMalesContributingSG[g]>0) {
						SibshipPapas(tempFathers, MalePointersSG[g], x, NumMalesContributingSG[g], 
									 MaleProbsSG[g], MaleProbsSG[g][NumMalesContributingSG[g]-1], P.MateSpec);
					}
					else {
						m=f=0;
					}
				}
				else {
					SibshipPapas(tempFathers, MalePointers, x, NumMalesContributing, 
							 MaleProbs, MaleProbs[NumMalesContributing-1], P.MateSpec);
				}
			}
			else {
				m=f=0;
			}
			
			/* now use x as in index variable, so initialize it.  */
			x = 0;
			
			
			/* then actually put the females into this year's array */
			if(NF[t]+f > fspaces) {  /* check to make sure there is room for them */
				printf("\n\nToo many females at time %d.  fspaces = %d\n\nExiting\n\n",t,fspaces);
				exit(1);
			}
			for(j=NF[t];j<NF[t]+f;j++,x++)  {
				
				FillIndiv( &(F[t][j]), &(F[at][i]), tempFathers[x], &gID, t, j,'F',birthplace,P);
				F[at][i].OldestOffsp = &(F[t][j]);
				PrintAndDiscardPedigree( &(F[t][j]), t, P.DiscardBefore,P.TossAll, P.MaxAge);
			}
			NF[t] += f;
			NFR[t]+= f;
			
			/* then put the males into this year's array */
			if(NM[t]+m > mspaces) {  /* check to make sure there is room for them */
				printf("\n\nToo many males at time %d.  mspaces = %d\n\nExiting\n\n",t,mspaces);
				exit(1);
			}
			for(j=NM[t];j<NM[t]+m;j++,x++)  {
				FillIndiv( &(M[t][j]), &(F[at][i]), tempFathers[x], &gID, t, j,'M',birthplace,P);
				F[at][i].OldestOffsp = &(M[t][j]);
				PrintAndDiscardPedigree( &(M[t][j]), t, P.DiscardBefore,P.TossAll, P.MaxAge);
			}
			NM[t] += m;
			NMR[t] += m;
			
			
		} /* closes the loop over i */
	}  /* closes the loop over age */
	
	
	/* and now print out how many died post-reproduction */
	/*	printf("POSTREPRO_DIED : AGES :  %d  : ",t);
	 for(age=1;age<=MA;age++)  {
	 printf("%d  ",age);
	 }
	 printf("\n");
	 printf("POSTREPRO_DIED : MALES :  %d  : ",t);
	 for(age=1;age<=MA;age++)  {
	 printf("%d  ",MaleDiedPostReprod[age]);
	 }
	 printf("\n");
	 printf("POSTREPRO_DIED : FEM :  %d  : ",t);
	 for(age=1;age<=MA;age++)  {
	 printf("%d  ",MaleDiedPostReprod[age]);
	 }
	 printf("\n");*/
	
	
	/* down here we need to decrement the counters that tell how many individuals are remaining---this is
	 because we need to account for any post-reprod death that they may have suffered */
	for(age=1;age<=MA;age++)  {
		at = t - age;
		NMR[at] -= MaleDiedPostReprod[age];
		NFR[at] -= FemDiedPostReprod[age];
	}
	
	
	
	/*	printf("VERBIAGE : Done with this generation\n"); */
	
	
	
	/* free up the memory we have allocated */
	free(MaleProbs);
	free(MalePointers);
	free(tempFathers);
	if(P.CohortSizesRandom==0) {
		free(FemDirPars);
		free(FemNumOffsLin);
		for(age=1;age<=MA;age++)  {
			at = t - age;
			free(FemNumOffs[age]);
		}
		free(FemNumOffs);
	}
	if(P.SpawnGroupSize!=-1 || P.SpawnGroupCnt!=-1) {
		free(FemSpawnGroupsLin);
		for(age=1;age<=MA;age++)  {
			free(SpawnGroup2D[age]);
		}
		free(SpawnGroup2D);
		for(i=0;i<NumMaleSpawnGroups;i++)  {
			free(MalePointersSG[i]);
			free(MaleProbsSG[i]);
		}
		free(MalePointersSG);
		free(MaleProbsSG);
		free(NumMalesContributingSG);
	}
	
	free(FemDiedPostReprod);
	free(MaleDiedPostReprod);
	
	/*	printf("VERBIAGE :  Exiting MakeBabies()\n"); */
	
}


/* generate number of offspring of females (or weights for males) given mean number of offspring,
 the coefficient of variation, and the Type of distribution.  Make that an Enum. */
int NumOffspringRV(double mu, double cv, enum dsn_type Type)
{
	if(mu==0.0) return (0);
	if(cv==1.0 && Type==NEG_BINOM)
		Type = POISSON;
	
	
	switch(Type) {
			
		case(NEG_BINOM):
			/* return a neg binom rv as gamma-mixed poisson --- this parameterization
			 is easier to understand in this case than the ranlib parameterization (which 
			 doesn't allow both variables to be doubles anyway */
			return(NegBinomAB_RV( mu*cv/(1.0-cv), (1.0-cv)/cv )); 
			break;
		case(BINARY):
			if(ranf() < mu) return(1);
			else return(0);
			break;
		default:
		case(POISSON):
			return(ignpoi((float)mu));
	}
}

/* 
 return in the array of pointers to indiv's Dads (which must have memory allocated to it outside of this function)
 a list of pointers to the fathers for n members of a full sibship.  All the potential N fathers
 have weight in the array Wts, the sum of which is S.  The pointers to the fathers themselves are
 in the array PossibleDads which is indexed the same way as Wts is.  The mate fidelity is dealt with by a 
 sort of Dirichlet Process / Chinese Restaurant process sort of thing having parameter r:
 
 With probability
 r/(i-1+r)
 
 the father of the i-th child in the sibship is drawn from the array of PossibleDads according to the
 probabilities in Wts.   With probability (i-1)/r, the father of child i is the father of one of the 
 i-1 children with subscript less than i.  Pretty straightforward.  
 
 r == infinity corresponds to totally random (no mate fidelity), but we denote that here by 
 passing any negative value of r
 
 */
void SibshipPapas(indiv **Dads, indiv **PossibleDads, int n, int N, double *Wts, double S, double r)
{
	int i,idx,j;
	double rando,k;
	
	if(r < 0.0) { /* each father is independently assigned according to Wts */
		for(i=0;i<n;i++)  {
			rando = (double)ranf() * S;
			idx = IntFromArray_BinarySearch(rando, Wts, 0, N-1);
			Dads[i] = PossibleDads[idx];
			
		}
	}
	else if(r==0.0) {  /* this is total mate fidelity --- i.e. the father of the first child is the father of every child */
		rando = (double)ranf() * S;
		idx = IntFromArray_BinarySearch(rando, Wts, 0, N-1);
		for(i=0;i<n;i++)  {
			Dads[i] = PossibleDads[idx];
		}
	}	
	else { /* in this case we do the Dirichlet Process */
		/* choose the father of the first child */
		rando = (double)ranf() * S;
		idx = IntFromArray_BinarySearch(rando, Wts, 0, N-1);
		Dads[0] = PossibleDads[idx];
		/* then do the rest */
		for(j=1;j<n;j++)  {
			k = (double)j + 1.0;
			if(ranf() < r / ( k - 1.0 + r) ) {  /* here we draw a father at random from all possible fathers */
				rando = (double)ranf() * S;
				idx = IntFromArray_BinarySearch(rando, Wts, 0, N-1);
				Dads[j] = PossibleDads[idx];
			}
			else {  /* otherwise we choose a child at random and assign his/her father to the current father */
				idx = UniformRV(0,j-1);
				Dads[j] = Dads[idx];
			}
		}
	}
}




/*
 this traces a random genealogy back for the two gene copies in ind, until it either coalesces with
 the lineage of one of the nf founding genes currently known, or it reaches a founder (the last gene
 copy, going back in time, before tf, or before a null parent) in which case that lineage becomes a lineage
 of the ++nf-th known founder.
 
 So, all gene copies present at time tf or earlier are assumed to be non-identical by descent.
 
 n is the total number of individuals in the pedigree.
 
 indics is a useful variable for keeping track of the lineages.  indics[idx][0] is the paternal
 gene of idx and indics[idx][1] is the maternal. The values are initialized to -2.  They get a -1 if 
 a lineage goes through them, and that gets changed to the index of an allele, eventually. 
 
 ma and pa are output variables.  Their values get set to the subscript of the founder gene ancestral
 to ma and the founder gene ancestral to pa, respectively.
 
 indics must have memory allocated to it outside of this function.  If you want to initialize all the 
 values in it to -2, then ReInit  should be 1.    
 
 tf is the time before which everyone is discarded.  
 
 AncsIdx and AncsGenes are arrays for holding information about who the ancestors were.
 They should have memory allocated outside this function.  You may as well allocate enough
 space to hold all individuals, but only a tiny fraction of that will be used. These arrays
 do not need to be initialized before calling AncestryOfSampledGenes() 
 
 FounderLocs is an array of the locales in which each founder gene was found. since it gets modified
 (new memory gets realloced to it)  within this
 function it is actually the address of the pointer to the array fo ints.
 
 *FLM is an output parameter telling how much memory has been allocated to FounderLocs
 
 */
void AncestryOfSampledGenes(int **indics, int n, indiv *ind, int *nf,  int *ma, int *pa, int ReInit, int *AncsIdx, int *AncsGenes, int tf, int **FounderLocs, int *FLM)
{	
	int i,j;
	indiv *temp, *old_temp;
	int route;
	int allele_idx;
	int done;
	int NumAncs;  /* for a more efficient accounting of who the ancestry gets traced through */
	
	
	if(ReInit==1) {  /* if null, then allocate to and initialize it */
		/*		printf("VERBIAGE :  initializing indics to -2 in AncestryOfSampledGenes\n"); */
		for(i=tf;i<=n;i++)  {  /* no need to worry about times before tf */
			indics[i][0] = -2;
			indics[i][1] = -2; 
		}
	}
	
	/* cycle over doing first the paternal then the maternal route out of the individual */
	for(i=0;i<2;i++)  {
		
		NumAncs = 0;  /* initialize this to count the number of ancestors for the maternal or paternal gene */
		indics[ind->idx][i] = -1;
		AncsIdx[NumAncs] = ind->idx;
		AncsGenes[NumAncs++] = i;
		
		if(i==0)
			temp = ind->pa;
		else 
			temp = ind->ma;
		
		/* have to define this in case the parents of this individual are both null */	
		old_temp = ind;
		
		done=0;
		while(done==0) {
			
			if(temp==NULL)  {  /* in this case we have reached a new founding gene
			 and so we record the lineage as such */
				/* printf("DEBUGGING: inside temp==NULL.\n"); */
				*nf = *nf + 1;
				allele_idx = *nf;
				
				/* and add the founder's location */
				if(*nf > *FLM) {  /* realloc more memory if necessary */
					(*FounderLocs) = (int *)realloc( (void *)(*FounderLocs), sizeof(int) * (*FLM + REALLOC_BLOCK) );
					*FLM = *FLM + REALLOC_BLOCK;
				} 
				(*FounderLocs)[*nf - 1] = old_temp->birth_pop;
				
				for(j=0;j<NumAncs;j++)  {
					if(indics[ AncsIdx[j] ][ AncsGenes[j] ] != -1)  {
						printf("WARNING : indics of a putative ancestor is %d not -1 with j = %d\n", indics[ AncsIdx[j] ][ AncsGenes[j] ] , j);
					}
					else {
						indics[ AncsIdx[j] ][ AncsGenes[j] ] = allele_idx;
					}
				}
				if(i==0)
					*pa = allele_idx;
				else
					*ma = allele_idx;
				
				done = 1;
			}
			else {
				/* printf("DEBUGGING: about to simulate route\n"); */
				route = UniformRV(0,1);
				if(indics[temp->idx][route] >= 0) {  /* we have hit an existing lineage */
 					allele_idx = indics[temp->idx][route];
 					for(j=0;j<NumAncs;j++)  {
						if(indics[ AncsIdx[j] ][ AncsGenes[j] ] != -1)  {
							printf("WARNING : indics of a putative ancestor is %d not -1\n", indics[ AncsIdx[j] ][ AncsGenes[j] ] );
						}
						else {
							indics[ AncsIdx[j] ][ AncsGenes[j] ] = allele_idx;
						}
					}
					if(i==0)
						*pa = allele_idx;
					else
						*ma = allele_idx;
					done = 1;
				}
				else if(indics[temp->idx][route] == -1) {
					printf("\n\nBad News!  We have a -1 where it is not expected in AncestryOfSampledGenes()\n\nExiting\n\n");
					exit(1);
				}
				else {  /* otherwise, record that this, as yet unnamed, lineage has been through here
				 and proceed to pa or ma depending on the value of route */
					/* printf("DEBUGGING: Moving on from temp->ID=%s which has temp->idx= %d \n",temp->ID,temp->idx); 
					printf("DEBUGGING: And that indiv has ma pointer %x and pa pointer %x\n",temp->ma, temp->pa); */
					indics[temp->idx][route] = -1;
					AncsIdx[NumAncs] = temp->idx;
					AncsGenes[NumAncs++] = route;
					
					old_temp=temp;
					
					
					if(route == 0)
						temp = temp->pa;
					else
						temp = temp->ma;
				}
			}
			
		}
	}
	
}





/* this is just a silly little function that will print the probability of different numbers of fathers 
 (actually different numbers of "dipping into the pools of fathers") contributing to a sibship.
 
 I will put this on the stdout stream from main. 
 
 r is the Dirichlet process parameter (like theta in the IAM).
 */
void PrintProbDistOfPapas(double r, int maxn)
{
	double **probs = (double **)ECA_CALLOC(maxn+1, sizeof(double *));
	int n,i;
	
	probs[1] = (double *)ECA_CALLOC(2,sizeof(double));
	
	
	/* initialize */
	probs[1][1] = 1.0;
	probs[1][0] = 0.0;
	
	printf("MATE_SPECIFITY : Printing matrix of probability of number of \"new\" males paternal to sibship of given size.  This is the number of \"dippings\" into the pool of males\n");
	printf("MATE_SPECIFITY : 1 :\t1.0000\n");
	
	if(r>=0.0)  {
		for(n=2;n<=maxn;n++)  {
			probs[n] = (double *)ECA_CALLOC(n+1,sizeof(double));
			printf("MATE_SPECIFITY : %d : \t",n);
			for(i=1;i<=n;i++)  {
				probs[n][i] = probs[n-1][i] * (n-1.0)/(r+n-1.0) + probs[n-1][i-1] * r/(r+n-1.0);
				printf("%.4f\t",probs[n][i]);
			}
			printf("\n");
		}
	}
	else {
		for(n=2;n<=maxn;n++)  { double temp;
			printf("MATE_SPECIFITY : %d : \t",n);
			for(i=1;i<=n;i++)  {
				if(i==n) temp = 1.0;
				else temp = 0.0;
				printf("%.4f\t",temp);
			}
			printf("\n");
		}
	}
}



#define DOING_GENETIC_SAMPLING NewbornSampFlag + GenosForAll_Flag + GtypAPropFemPre_Flag + GtypAPropMalePre_Flag + \
GtypAPropFemPost_Flag + GtypAPropMalePost_Flag + GtypAPropMaleDur_Flag + GtypAPropFemDur_Flag

#define REQ_FOR_SAMPLING ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ALREADY_HAS(TFlag, -T or --number-of-years) &&  \
ALREADY_HAS(DiscardAll_Flag, --discard-all) 

/* ProcessOptions(&P,&T, &NM, &NF,&NMR, &NFR, &S) */
/* This returns 1 if it hit a new_pop, and zero otherwise */
int ProcessOptions(int *output_argc, char **output_argv[], PopPars *P,int *T, int **NM, int **NF, 
				   int **NMR, int **NFR, SampPars *S, char *locfilename, PopSizePars *Pop, int *N, int **MSpace, int **FSpace, int *HasLocFile)
{
	int j;
	int argc = *output_argc;
	char **argv = *(output_argv);
	int CommFilesInserted = 0;
	/* here are some flags for the options */
	int MaxAgeFlag = 0,
	SurvivalRatesFlag = 0,
	FemSurvRatesFlag = 0,
	MaleSurvRatesFlag = 0,
	FasrfFlag = 0,
	MasrpFlag = 0,
	FPRFlag = 0,
	MPRFlag = 0,
	RCSFlag = 0,
	FCSFlag = 0,
	FCVFlag = 0,
	OffsDsnFlag = 0,
	MCVFlag = 0,
	StickyM_Flag = 0,
	MateSpecFlag = 0,
	YearVarFlag = 0,
	SexRatioFlag = 0,
	TFlag = 0,
	InitMaleFlag = 0,
	InitFemaleFlag = 0,
	NewbornSampFlag = 0,
	LocFileFlag = 0,
	Cohort_Size_Flag = 0,
	GenosForAll_Flag = 0,
	GtypAPropFemPre_Flag = 0,
	GtypAPropMalePre_Flag = 0,
	GtypAPropFemDur_Flag = 0,
	GtypAPropMaleDur_Flag = 0,		
	GtypAPropFemPost_Flag = 0,
	GtypAPropMalePost_Flag = 0,
	LethalSampling_Flag = 0,
	DiscardAll_Flag = 0,
	DiscardParents_Flag = 0,
	FPDPR_Flag = 0,
	NewPop_Flag = 0,
	MaleProbMigIn_Flag = 0,
	MaleProbMigOut_Flag = 0,
	FemProbMigIn_Flag = 0,
	FemProbMigOut_Flag = 0,
	MPDPR_Flag = 0,
	SneakyMProb_Flag = 0,
	ReproInhibFlag = 0,
	SpawnGroupCnt_Flag = 0,
	SpawnGroupSize_Flag = 0,
	alloc_extra_f = 0;
	
	DECLARE_ECA_OPT_VARS;
	
	SET_OPT_WIDTH(26);
	SET_ARG_WIDTH(25);
	SET_VERSION("\nspip_m -- a program for simulating pedigrees within populations exchanging migrants\n\nVERSION: \
				1.0 Beta\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 3 November 2005 \nCOPYRIGHT: U.S. Federal Government Work\n\n")
	SET_VERSION_HISTORY("\n\nIN PROGRESS\n")
	SET_DESCRIBE("spip_m  ---  a program for simulating pedigrees within populations exchanging migrants")
	
	
	/* here I set some defaults for non-required options */
	P->OffsDsn = NEG_BINOM;
	P->Fcv = 1.0;
	P->ReproInhib = 0;
	P->MCVrs = 1.0;
	P->sticky_m_alpha = -1.0;  /* if it is > 0 it is a flag to be using it */ 
	P->SexRatio = .5;
	P->MateSpec = -1.0;
	P->DiscardBefore = 0;
	P->TossAll = 0;
	P->CohortSizesRandom=1;
	P->MigPr = NULL;
	P->Dest = NULL;
	P->SneakyMProb = NULL;
	P->SpawnGroupSize = -1;
	P->SpawnGroupCnt = -1;
	S->nc = 0;  /* don't need to sample any individuals if you don't want to */
	S->NumLoc = 0;
	P->OffspringAllocOverload = 10;
	
	/* the following NULL's serve as flags that the proportional random sampling is 
	 to be implemented */
	S->pre_male_samp = NULL;
	S->post_male_samp = NULL;
	S->pre_fem_samp = NULL;
	S->post_fem_samp = NULL;
	S->pre_male_samp_years = NULL;
	S->pre_fem_samp_years = NULL;
	S->post_male_samp_years = NULL;
	S->post_fem_samp_years = NULL;
	S->dur_male_samp_years = NULL;
	S->dur_fem_samp_years = NULL;
	S->dur_repro_fem_samp = 0.0;
	S->dur_repro_male_samp = 0.0;
	S->LethalityProb = 0.0;
	
	BEGIN_OPT_LOOP
	
	/* Here, if there was a command-file that got inserted into the command line
	 (that would have happened from some code in BEGIN_OPT_LOOP),
	 argc and argv will have been modified, and we will want that modification to 
	 be reflected in output_argc and output argv.  So, we make the assignment here */
	/* To get this to work right, of course, we only want to make these assignments if 
	 new command files have been inserted this time through.  The if statement takes 
	 care of that. */
	if(__OptCommFileFlag != CommFilesInserted) {
		*output_argc = argc;
		*output_argv = argv;
		CommFilesInserted = __OptCommFileFlag;
	}
	
	
	
	/*		if(__OptCommentLevel > 0) {
	 printf("\n   ****  options for dealing with multiple-population simulations ****\n\n");
	 } */
	OPEN_SUBSET(Options for dealing with multiple population simulations,
				Multi-Pop Options,
				options for dealing with multiple-population simulations);
	
	if(OPTION(NewPop_Flag,
			  ,
			  new-pop,
			  ,
			  signifies completion of options for the current population,
			  use this option to signify that you are done inputting options
			  for the current population.  After this option is issued the program
			  will check to make sure that all necessary options were given to 
			  describe the population fully and the program will prepare itself
			  to read in options for the next population.  DO NOT include a --new-pop
			  on the command line after the options for the final population or before
			  the options for the first population!))  {
		/* here we just advance the locale counter and the argv pointer, etc. but we will
		 do that later */
		printf("NEWPOP :  Received the --new-pop option at element %d of the currently modified version of the command line.\n",__Opt_i);
		
		/* modify the command line to just include the part past where it got to */
		*output_argc = argc - __Opt_i;
		*output_argv = &((*output_argv)[__Opt_i]);
		(*output_argv)[0] = argv[0];  /* put the program name back in there */
		printf("REMAINING_COMMAND_LINE_LOOKS_LIKE output_argc %d : ",*output_argc);
		for(j=0;j<*output_argc;j++)  {
			printf("%s ",(*output_argv)[j]);
		}
		printf("\n");
		
		/* reducing argc ensures that the next time it looks for missing options */
		argc = __Opt_i+1;		
		
	}
	
	if(OPTION (gNumPops_Flag,
			   ,
			   num-pops,
			   J,
			   number of populations in the simulation,
			   This option must be given if you want the populations defined using the 
			   --new-pop option to exchange migrants.  J must be the number of populations
			   in the simulation---i.e. one plus the number of times the --new-pop option will be given.
			   This option may be used only once and must be given before the --fem-prob-mig-in or --male-prob-mig-in options.
			   
			   )) {
		if(gNumPops_Flag > 1) {
			fprintf(stderr,"Error! Attempting to issue the --num-pops option twice.  Not allowed!\n");
			GET_INT;
			OPT_ERROR;
		}
		else if(ARGS_EQ(1)) {
			gNumPops = GET_INT;
		}
	}
	
	CLOSE_SUBSET;
	
	OPEN_SUBSET(Age-structure and reproduction parameters,
				Age Structure/Reproduction,
				Age-structure and reproduction parameters);
	
	if(REQUIRED_OPTION(MaxAgeFlag,A,max-age,J , the maximum age the organism can reach , this sets the maximum age that any individual \
					   in the population can reach.  All organisms are programmed to have died with probability one before they exceed\
					   this age. Required.) ) {
		if(ARGS_EQ(1)) {
			P->MaxAge = GET_INT;
		}
		/* take some time to allocate memory and set some default values for the prob of reproduction-induced death */
		P->FPDPR = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
		P->MPDPR = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
		
	}
	if(COND_REQ_OPTION(SurvivalRatesFlag,s,survival-probs,R1~R2:R_MA , age specific probabilities of survival  ,  \
					   age-specific probability of survival.  R1 is the probabiilty of of surviving from age
					   0 to age 1.  R2 is the probability of surviving from age 1 to age 2.  There should be 
					   MaxAge of these R arguments.  The last of them is the probability of surviving from age 
					   MaxAge-1 to age MaxAge. This option is required unless sex-specific survival rates are set
					   using the fem-surv-probs and the male-surv-probs options.,
					   FemSurvRatesFlag + MaleSurvRatesFlag == 0,
					   when both fem-surv-probs and male-surv-probs are not given
					   ) ) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			P->femsurv = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
			P->malesurv = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				P->femsurv[j] = GET_DUB;
				P->malesurv[j] = P->femsurv[j];
			}	
		}
	}
	if(COND_REQ_OPTION(FemSurvRatesFlag,
					   ,
					   fem-surv-probs,
					   R1~R2:R_MA ,
					   age-specific survival probs for females,
					   Age-specific probability of survival for females.  R1 is the probability of of surviving from age
					   0 to age 1.  R2 is the probability of surviving from age 1 to age 2.  There should be 
					   MaxAge of these R arguments.  The last of them is the probability of surviving from age 
					   MaxAge-1 to age MaxAge.  Both fem-surv-probs and male-surv-probs are required unless the survival
					   probabilities for males and females are the same.  In that case then male and female survival rates
					   may both be set by the arguements to the single option survival-probs.  , 
					   SurvivalRatesFlag==0, when the survival-probs is not given
					   )) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			P->femsurv = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				P->femsurv[j] = GET_DUB;
			}	
		}		
	}
	if(COND_REQ_OPTION(MaleSurvRatesFlag,
					   ,
					   male-surv-probs,
					   R1~R2:R_MA,
					   age-specific survival probs for males,
					   Age-specific probability of survival for males.  R1 is the probability of of surviving from age
					   0 to age 1.  R2 is the probability of surviving from age 1 to age 2.  There should be 
					   MaxAge of these R arguments.  The last of them is the probability of surviving from age 
					   MaxAge-1 to age MaxAge.  Either both --fem-surv-probs and --male-surv-probs are required; or just
					   --survival-probs (-s).  , 
					   SurvivalRatesFlag==0, when the survival-probs is not given
					   )) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			P->malesurv = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				P->malesurv[j] = GET_DUB;
			}	
		}		
	}
	if(REQUIRED_OPTION(FasrfFlag,f,fem-asrf,R1~R2:R_MA, female age specific relative fecundities  ,  \
					   age-specific relative fecundities.  R1 is the relative fecundity of a one-year old female.  R2 
					   is the same for a 2 year-old female.  There should be 
					   MaxAge of these R arguments.  The last one is the relative fecundity of MaxAge-year-old females.
					   It is assumed that 0-year olds do not reproduce.  These numbers need not sum to one.  All that 
					   matters is their _relative_ sizes. Required. ) ) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			P->Fasrf = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				P->Fasrf[j] = GET_DUB;
			}	
		}
	}
	if(REQUIRED_OPTION(FPRFlag,,fem-prob-repro,R1~R2:R_MA, female age specific probabilities of reproducing ,  \
					   R1 is the absolute probability that a one-year old female will reproduce.  R2 
					   is the same for a 2 year-old female.  There should be 
					   MaxAge of these R arguments.  The last one is the probability that a MaxAge-year-old female 
					   reproduces.
					   It is assumed that 0-year olds do not reproduce.  These numbers need not sum to one.  They must
					   be between zero and 1.
					   They are not relative values.  Required.) ) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			P->FPR = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				P->FPR[j] = GET_DUB;
				if(P->FPR[j] > 1.0 || P->FPR[j] < 0.0) {
					fprintf(stderr,"Error! Argument %d to option --fem-prob-repro is %f which is >1 or <0 and not a valid probability.\n",j,P->FPR[j]);
					OPT_ERROR;
				}
			}	
		}
	}
	if(OPTION(ReproInhibFlag,
			  ,
			  repro-inhib,
			  J,
			  age of offspring below which reproduction is inhibited,
			  If a female has an offspring which is less than or equal to J then she will not reproduce.  This option
			  is used to deal with organisms like whales in which a female will not produce another calf if she is currently
			  nursing one.  This option may only be used with J > 0 only when the binary variant of the --offsp-dsn option has been invoked.  Note
			  that if the calf dies before the age of J then it no longer inhibits reproduction. ) ) {
		if(ARGS_EQ(1)) {
			P->ReproInhib = GET_INT;
		}
		
	}
	if(REQUIRED_OPTION(MasrpFlag,m,male-asrp,R1~R2:R_MA, male age specific reproductive potential  ,  \
					   age-specific relative male reproductive potential.  R1 is the relative reproductive potential of 
					   a one-year old male.  R2 
					   is the same for a 2 year-old male.  There should be 
					   MaxAge of these R arguments.  The last one is the relative reproductive potential of MaxAge-year-old males.
					   It is assumed that 0-year olds do not reproduce.  These numbers need not sum to one.  All that 
					   matters is their _relative_ sizes. Required. ) ) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			P->Masrp = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				P->Masrp[j] = GET_DUB;
			}	
		}
	}
	if(REQUIRED_OPTION(MPRFlag,,male-prob-repro,R1~R2:R_MA, male age specific probabilities of reproducing ,  \
					   R1 is the absolute probability that a one-year old male will reproduce.  R2 
					   is the same for a 2 year-old male.  There should be 
					   MaxAge of these R arguments.  The last one is the probability that a MaxAge-year-old male 
					   reproduces.
					   It is assumed that 0-year olds do not reproduce.  These numbers need not sum to one. They must
					   be between zero and 1.
					   They are not relative values.  Required.) ) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			P->MPR = (double *)ECA_CALLOC(P->MaxAge+1,sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				P->MPR[j] = GET_DUB;
				if(P->MPR[j] > 1.0 || P->MPR[j] < 0.0) {
					fprintf(stderr,"Error! Argument %d to option --male-prob-repro is %f which is >1 or <0 and not a valid probability.\n",j,P->MPR[j]);
					OPT_ERROR;
				}
			}	
		}
	}
	if( OPTION(OffsDsnFlag,,offsp-dsn ,C{pois|negbin|binary},distribution of offspring number , this sets the
			   distribution of the number of offspring per female.  C may be either pois or negbin
			   or binary
			   to make the distribution either Poisson or Negative Binomial or Binary.  A Binary distribution
			   means that each reproductive female will have either one of zero offspring.  The binary option
			   may only be used if --fixed-cohort-size has not been issued.  If 
			   it is negative binomial then the parameters of that dsn are set according to the value for 
			   fem-rep-disp-par.  Optional.  Default is negbin.) )  { char tempstr[100];
		if(ARGS_EQ(1)) {
			GET_STR(tempstr);
			if(strcmp(tempstr,"negbin")==0) {
				P->OffsDsn = NEG_BINOM;
			}
			else if(strcmp(tempstr,"pois")==0) {
				P->OffsDsn = POISSON;
			}
			else if(strcmp(tempstr,"binary")==0) {
				P->OffsDsn = BINARY;
			}
			else {
				fprintf(stderr, "Error! Unknown parameter %s to option offsp-dsn\n", tempstr);
				__OptError=1;
			}
		}		
	}
	if( CLASHABLE_OPTION(FCVFlag,,fem-rep-disp-par,R,female dispersion of reprod. success , this variable
						 controls the variance in female reproductive success.  It is the ratio of the mean offspring number
						 to the variance in offspring number per female.  The argument must be greater than 0 and less 
						 than or equal to 1.  Optional.  Default value is 1.,
						 (P->OffsDsn == POISSON && P->Fcv != 1.0) , is set to something other than 1.0 and offsp-dsn has been 
						 given the pois argument.
						 )  )  {
		if(ARGS_EQ(1)) {
			P->Fcv = GET_DUB;
			if(DispersionParError(P->Fcv, "--fem-rep-disp-par") ) {
				OPT_ERROR;
			} 
		}		
	}
	if( CLASHABLE_OPTION(MCVFlag,,male-rep-disp-par,R,male dispersion of reprod. success , this variable
			   controls the variance in male reproductive success.  It is the ratio of the mean offspring number
			   to the variance in offspring number per male. The argument must be greater than 0 and less 
			   than or equal to 1. Optional.  Default value is 1.,
			   StickyM_Flag,
				is issued and the male-sticky-rep-var has also been given.)  )  {
		if(ARGS_EQ(1)) {
			P->MCVrs = GET_DUB;
			if(DispersionParError(P->MCVrs, "--male-rep-disp-par") ) {
				OPT_ERROR;
			} 
		}		
	}
	if( CLASHABLE_OPTION(StickyM_Flag,
						 ,
						 male-sticky-rep-var,
						 R,
						 male variance in reproductive success maintained over years , 
						 this variable controls the variance in male reproductive success.  Unlike the male-rep-disp-par option\054
						 however\054 the male retains his mating rate across different years.  Specifically\054 R is the shape
						 paramater of the gamma distribution from which every male obtains (upon birth) a mating success rate.  His
						 overall rate in any given year is then then product of this mating success rate and his age-specific relative 
						 reproductive prowess\054 given his age.  These rates are normalized to sum to one and thus become probabilities
						 of fathering children.  It is not really possible to say how this relates to the ratio of the mean to the variance
						 in offspring number because the mean is not necessarily known.  But it can be figured out.  This option cannot be 
						 used along with the male-rep-disp-par option.  WARNING: THIS IS SOMEWHAT EXPERIMENTAL AT THE MOMENT AND HAS NOT BEEN
						 TESTED WITH LOTS OF OTHER OPTIONS LIKE SNEAKY-MALES\054 ETC.,
						 (MCVFlag) , is issued and the male-rep-disp-par has also been given.
						 )  )  {
		if(ARGS_EQ(1)) {
			P->sticky_m_alpha = GET_DUB;
			if( P->sticky_m_alpha <=0  ) {
				fprintf(stderr,"The argument of the male-sticky-rep-var option must be greater than 0!  Not %f. \n",P->sticky_m_alpha);
				OPT_ERROR;
			} 
		}		
	}
	if( OPTION(MateSpecFlag,,mate-fidelity,R,mate fidelity parameter , 
			   controls how faithful females are to mating with a single male.  If [R] is negative then the father for 
			   each successive child of a female is chosen randomly.  If [R] is zero then there will only be a single
			   father for all members of a sibship.  (i.e. each female mates with only one male). As [R] gets larger the fidelity
			   decreases.  [R] turns out to be the parameter in the urn scheme of Hoppe. Optional.  Default value is -1---
			   no mate fidelity.)  )  {
		if(ARGS_EQ(1)) {
			P->MateSpec = GET_DUB;
		}		
	}
	if(OPTION(FPDPR_Flag,
			  ,
			  fem-postrep-die,
			  R1~R2:R_MA,
			  additional prob of death for a female that mates,
			  R1 is the probability that a female will die after mating. Whether or not an individual mates
			  is determined by the corresponding argument to the --fem-prob-repro option.  This is an extra source
			  of mortality above and beyond that specified by the survival probabilities.  R2 is the prob that 
			  a 2-year-old female dies after mating.  There should be MaxAge such values.    This option is useful 
			  for semelparous organisms that and others that have an elevated
			  probability of death in years when they mate. Optional.  The default is
			  0 for all age groups.)) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			for(j=1;j<=P->MaxAge;j++)  {
				P->FPDPR[j] = GET_DUB;
			}	
		}
	}
	if(OPTION(MPDPR_Flag,
			  ,
			  male-postrep-die,
			  R1~R2:R_MA,
			  additional prob of death for a male that mates,
			  R1 is the probability that a male will die after mating. Whether or not an individual mates
			  is determined by the corresponding argument to the --male-prob-repro option.  This is an extra source
			  of mortality above and beyond that specified by the survival probabilities.  R2 is the prob that 
			  a 2-year-old male dies after mating.  There should be MaxAge such values.    This option is useful 
			  for semelparous organisms that and others that have an elevated
			  probability of death in years when they mate. Optional.  The default is
			  0 for all age groups.)) {
		if( ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge)) {
			for(j=1;j<=P->MaxAge;j++)  {
				P->MPDPR[j] = GET_DUB;
			}	
		}
	}
	if( OPTION(SexRatioFlag,,sex-ratio ,R, probability an offspring is male,
			   R sets the probability that an individual will be born as a male.  R should be between 0 and 1.
			   Optional.  Default value is .5.
			   )  )  {
		if(ARGS_EQ(1)) {
			P->SexRatio = GET_DUB;
		}		
	}
	if( OPTION(YearVarFlag , ,year-var-switch , , not currently implemented ,)  )  {
		;		
	}
	
	
	/***************  SPAWNER GROUP SORTS OF THINGS **************/
	OPEN_SUBSET(Spawning group restriction parameters,
				Spawning Group,
				Parameters for restricting spawning to be among certain clusters of individuals in each year);
	
	if(OPTION(
			SpawnGroupSize_Flag,
			,
			spawn-grp-size,
			J,
			Set spawning group sizes of J females and J males each year,
			When this is invoked\054 the spawners each years are randomly partitioned into groups of size J males and J females.  These groups are given simple indexes
				0\054 1\054... as many as needed.  Females in a group with index X may only mate with males from a group with index X.  If there are more males than
				females\054 some of those males may be in spawning groups that do not get the chance to mate because there is not corresponding female spawning group.  
				The same is true of some females when there are more females than males.  This can be used to simulate full-factorial mating of small groups of individuals
				at salmon hatcheries. This option conflicts with the spawn-grp-cnt option.)) {
		if(SpawnGroupCnt_Flag) {
			fprintf(stderr,"Error!  The --spawn-grp-size option cannot be issued if --spawn-grp-cnt has already been given within a population!");
			OPT_ERROR;
		}
		if(ARGS_EQ(1)) {
			P->SpawnGroupSize = GET_INT;
		}
	}
	
	if(OPTION(
						SpawnGroupCnt_Flag,
						,
						spawn-grp-cnt,
						J,
						Randomly partition each year spawners into J randomly sized spawner groups,
						When this is invoked\054 the male and female spawners each years are randomly partitioned into J groups with sizes drawn from a multinomial distribution
						with cell probabilities equal to 1/J.  These groups are given simple indexes
						0\054 1\054...\054 J.  Females in a group with index X may only mate with males from a group with index X.  This option should be useful for simulating
						the occurrence of day-buckets at salmon hatcheries.  This option conflicts with the spawn-grp-size option. )) {
		if(SpawnGroupSize_Flag) {
			fprintf(stderr,"Error!  The --spawn-grp-cnt option cannot be issued if the --spawn-grp-size has already been given within a population!");
			OPT_ERROR;
		}
		if(ARGS_EQ(1)) {
			P->SpawnGroupCnt = GET_INT;
		}
	}
	
	CLOSE_SUBSET;
	
	/**************** MIGRATION RATE OPTIONS ********************/
	CLOSE_SUBSET;
	OPEN_SUBSET(Migration rate and related parameters,
				Migration,
				Parameters for setting migration rates into and out of demes etc.);
	
	if(MULT_USE_OPTION(
					   MaleProbMigOut_Flag,
					   ,
					   male-prob-mig-out,
					   G1 G2 R,
					   time-and-age-specific prob. that a male migrates out of a population,
					   R is the probability in each year given in range G1 that a male whose age is given in the
					   range G2 will migrate out of the population he is currently resident in.  This is 
					   population-specific.  Must appear after the --max-age or -A option for a population.
					   This option may be used multiple times for a single population.  If given more than once for the same age
					   group the program will use the last-given value.  In other words later uses of this
					   option can overwrite earlier uses of it within the same population.  Note that the default
					   migration probability is zero.  BIG NOTE: no migration is going to be done at any time point which is
					   less than the max age because that corresponds to the user-given initial counts of organisms.  If G1 includes
					   any times less than MaxAge\054 they will be discarded and you will get a WARNING from function StringToIntArray. ,
					   999999999) ) { enum sex_enum sexnum = MALE;
		if(ALREADY_HAS(MaxAgeFlag,-A or --max-age) &&  ALREADY_HAS(TFlag, -T or --number-of-years  ) && ARGS_EQ(3) ) { char timerange[10000]; 
			char agerange[10000]; 
			int *timearray; 
			int *agearray; 
			double temp_pr;
			int t;
			GET_STR(timerange);
			GET_STR(agerange);
			temp_pr = GET_DUB;
			timearray = StringToIntArray(timerange, P->MaxAge, P->MaxAge + *T - 1, 1);
			agearray = StringToIntArray(agerange, 0, P->MaxAge, 1);
			if(P->MigPr==NULL) {
				P->MigPr = AllocAndInitMigPr(*T, P->MaxAge);
			}
			for(t=P->MaxAge;t<*T+P->MaxAge;t++) {
				if(timearray[t]==1) {
					for(j=1;j<=P->MaxAge;j++)  {
						if(agearray[j]==1) P->MigPr[t][sexnum][j] = temp_pr;
					}
				}
			}
			free(timearray);
			free(agearray);
		}
		
	}
	
	if(MULT_USE_OPTION(
					   MaleProbMigIn_Flag,
					   ,
					   male-prob-mig-in,
					   G1 G2 R0~R1:R_NP ,
					   time-and-age-specific prob of a male migrating to a population,
					   Rk is the probability in the year given in range G1 that a male whose age is given in the
					   range G2 will migrate to population k given that it has migrated out of the population
					   it is currently resident in.  Populations are indexed from 0 up to the number of populations minus one. 
					   After G1 and G2 there should be one probability for 
					   each population.  So the total number of arguments---including G1 and G2---should be
					   the number of populations plus 2. This option is 
					   population-specific.  The option cannot appear before the --max-age or -A option for a population and it
					   also cannot appear before the --num-pops option is given nor before the -T or --number-of-years option.
					   This option may be used multiple times for a single population.  If given more than once for the same age
					   group the program will use the last-given value.  In other words later uses of this
					   option can overwrite earlier uses of it within the same population.  NOTE:  the probability of migrating
					   to the population of origin will always be set to zero regardless of the choice of the user on the 
					   command line.  And the remaining entries will be scaled so as to sum to 1.  If this option is not given
					   or if it is not  given for some age groups the default behavior is to have equal probability of migration
					   into any of the populations. BIG NOTE: no migration is going to be done at any time point which is
					   less than the max age because that corresponds to the user-given initial counts of organisms.  If G1 includes
					   any times less than MaxAge\054 they will be discarded and you will get a WARNING from function StringToIntArray.,
					   999999999) ) {  enum sex_enum sexnum = MALE;
		if(ALREADY_HAS(MaxAgeFlag,-A or --max-age) && ALREADY_HAS(TFlag,-T or --number-of-years) && ALREADY_HAS(gNumPops_Flag, --num-pops) ) {
			if(ARGS_EQ(gNumPops + 2) ) { 
				int k; 
				char timerange[10000]; 
				char agerange[10000]; 
				int *timearray; 
				int *agearray; 
				int t;
				double tempprobs[9999]; 
				double sum=0.0;
				GET_STR(timerange);
				GET_STR(agerange);
				for(j=0;j<gNumPops;j++)  { 
					tempprobs[j] = GET_DUB;
					if(j != gCurrentPopNum) sum += tempprobs[j];
				}
				timearray = StringToIntArray(timerange, P->MaxAge, P->MaxAge + *T - 1, 1);
				agearray = StringToIntArray(agerange, 0, P->MaxAge, 1);
				if(P->Dest==NULL) {
					P->Dest = AllocAndInitMigDest(*T, P->MaxAge, gNumPops, gCurrentPopNum);
				}
				for(t=P->MaxAge;t<*T+P->MaxAge;t++) {
					if(timearray[t]==1) {
						for(k=1;k<=P->MaxAge;k++) {
							if(agearray[k]==1) {
								for(j=0;j<gNumPops;j++)  {
									if(j == gCurrentPopNum) {
										P->Dest[t][sexnum][k][j] = 0.0;
									}
									else {
										P->Dest[t][sexnum][k][j] = tempprobs[j]/sum;
									}
								}
							}
						}
					}
				}
			}/* closes if ARGS_EQ */
			else {
				printf("Warning! Mismatch between argument of --num-pops and the number of arguments to --male-prob-mig-in\n");
			}
		}  
		
		
	}
	
	
	if(MULT_USE_OPTION(
					   FemProbMigOut_Flag,
					   ,
					   fem-prob-mig-out,
					   G1 G2 R,
					   time-and-age-specific prob. that a female migrates out of a population,
					   R is the probability in each year given in range G1 that a female whose age is given in the
					   range G2 will migrate out of the population she is currently resident in.  This is 
					   population-specific.  Must appear after the --max-age or -A option for a population.
					   This option may be used multiple times for a single population.  If given more than once for the same age
					   group the program will use the last-given value.  In other words later uses of this
					   option can overwrite earlier uses of it within the same population.  Note that the default
					   migration probability is zero. BIG NOTE: no migration is going to be done at any time point which is
					   less than the max age because that corresponds to the user-given initial counts of organisms.  If G1 includes
					   any times less than MaxAge\054 they will be discarded and you will get a WARNING from function StringToIntArray.,
					   999999999) ) { enum sex_enum sexnum = FEMALE;
		if(ALREADY_HAS(MaxAgeFlag,-A or --max-age) &&  ALREADY_HAS(TFlag, -T or --number-of-years  ) && ARGS_EQ(3) ) { char timerange[10000]; 
			char agerange[10000]; 
			int *timearray; 
			int *agearray; 
			double temp_pr;
			int t;
			GET_STR(timerange);
			GET_STR(agerange);
			temp_pr = GET_DUB;
			timearray = StringToIntArray(timerange, P->MaxAge, P->MaxAge + *T - 1, 1);
			agearray = StringToIntArray(agerange, 0, P->MaxAge, 1);
			if(P->MigPr==NULL) {
				P->MigPr = AllocAndInitMigPr(*T, P->MaxAge);
			}
			for(t=P->MaxAge;t<*T+P->MaxAge;t++) {
				if(timearray[t]==1) {
					for(j=1;j<=P->MaxAge;j++)  {
						if(agearray[j]==1) P->MigPr[t][sexnum][j] = temp_pr;
					}
				}
			}
			free(timearray);
			free(agearray);
		}
		
	}
	
	if(MULT_USE_OPTION(
					   FemProbMigIn_Flag,
					   ,
					   fem-prob-mig-in,
					   G1 G2 R0~R1:R_NP-1,
					   time-and-age-specific prob of a female migrating to a population,
					   Rk is the probability in the year given in range G1 that a female whose age is given in the range
					   G2 will migrate to population k given that it has migrated out of the population
					   it is currently resident in. Populations are indexed from 0 up to the number of populations minus one.
					   After G1 and G2 there should be one probability for 
					   each population.  So the total number of arguments---including G1 and G2---should be
					   the number of populations plus 2. This option is 
					   population-specific.  The option cannot appear before the --max-age or -A option for a population and it
					   also cannot appear before the --num-pops option is given nor before the -T or --number-of-years option.
					   This option may be used multiple times for a single population.  If given more than once for the same age
					   group the program will use the last-given value.  In other words later uses of this
					   option can overwrite earlier uses of it within the same population.  NOTE:  the probability of migrating
					   to the population of origin will always be set to zero regardless of the choice of the user on the 
					   command line.  And the remaining entries will be scaled so as to sum to 1.  If this option is not given
					   or if it is not  given for some age groups the default behavior is to have equal probability of migration
					   into any of the populations. BIG NOTE: no migration is going to be done at any time point which is
					   less than the max age because that corresponds to the user-given initial counts of organisms.  If G1 includes
					   any times less than MaxAge\054 they will be discarded and you will get a WARNING from function StringToIntArray.,
					   999999999) ) {  enum sex_enum sexnum = FEMALE;
		if(ALREADY_HAS(MaxAgeFlag,-A or --max-age) && ALREADY_HAS(TFlag,-T or --number-of-years) && ALREADY_HAS(gNumPops_Flag, --num-pops) ) {
			if(ARGS_EQ(gNumPops + 2) ) { 
				int k; 
				char timerange[10000]; 
				char agerange[10000]; 
				int *timearray; 
				int *agearray; 
				int t;
				double tempprobs[9999]; 
				double sum=0.0;
				GET_STR(timerange);
				GET_STR(agerange);
				for(j=0;j<gNumPops;j++)  { 
					tempprobs[j] = GET_DUB;
					if(j != gCurrentPopNum) sum += tempprobs[j];
				}
				timearray = StringToIntArray(timerange, P->MaxAge, P->MaxAge + *T - 1, 1);
				agearray = StringToIntArray(agerange, 0, P->MaxAge, 1);
				if(P->Dest==NULL) {
					P->Dest = AllocAndInitMigDest(*T, P->MaxAge, gNumPops, gCurrentPopNum);
				}
				for(t=P->MaxAge;t<*T+P->MaxAge;t++) {
					if(timearray[t]==1) {
						for(k=1;k<=P->MaxAge;k++) {
							if(agearray[k]==1) {
								for(j=0;j<gNumPops;j++)  {
									if(j == gCurrentPopNum) {
										P->Dest[t][sexnum][k][j] = 0.0;
									}
									else {
										P->Dest[t][sexnum][k][j] = tempprobs[j]/sum;
									}
								}
							}
						}
					}
				}
			}/* closes if ARGS_EQ */
			else {
				printf("Warning! Mismatch between argument of --num-pops and the number of arguments to --fem-prob-mig-in\n");
			}
		}  		
	}
	if(OPTION(SneakyMProb_Flag,
			  ,
			  sneaky-males,
			  R0~R1:R_NP-1 ,
			  rates of sneaky-male matings,
			  Rk is the probability that a male from locale k will sire offspring with a female from
			  a population that is not k.  Populations are indexed from 0 up to the number of populations minus one. 
			  MUST EXPLAIN MORE LATER.  Number of args should
			  be the same as the number of pops given by the --num-pops flag.  Note that as implemented we do not
			  impose any extra probability of dying on these non-native males mates.  Also note that the male age specific
			  relative reproductive prowess is taken from the population he is sneaking into because those values are relative
			  measures and cannot be compared easily between populations.  ) ) {
		if(ALREADY_HAS(gNumPops_Flag, --num-pops)) {
			if(ARGS_EQ(gNumPops)) { int k; double burn;
				P->SneakyMProb = (double *)ECA_CALLOC(gNumPops, sizeof(double));
				for(k=0;k<gNumPops;k++)  {
					if(k==gCurrentPopNum) {
						P->SneakyMProb[k] = 0.0;
						burn = GET_DUB;	
					}
					else 
						P->SneakyMProb[k] = GET_DUB;
					
					if(P->SneakyMProb[k]<0.0 || P->SneakyMProb[k] > 1.0) {
						fprintf(stderr,"Error!  Argument %d of --sneaky-males is %f, which is less than zero or greater than 1\n",k+1,P->SneakyMProb[k]);
						OPT_ERROR;
					}
				}
			}
		}
		
	}
	
	/************** INITIAL CONDITIONS AND DURATION ************/
	CLOSE_SUBSET;
	OPEN_SUBSET(Initial-condition and duration-of-simulation parameters,
				Starting Conditions,
				Parameters for setting the starting state of the population and controlling the number of generations to simulate.);
	
	
	if(REQUIRED_OPTION(TFlag,T,number-of-years,J,number of years to simulate,
					   J sets the number of years to carry out in the simulations.  The years for which  individuals 
					   "exist" in the simulations are from year 0 to year T + MaxAge - 1 inclusive.  The first simulated offspring cohort is born at time
					   MaxAge. ) ) {
		if(ARGS_EQ(1) ) {
			*T = GET_INT;
		}	
	}
	if(REQUIRED_OPTION(InitMaleFlag,,initial-males,J0~J1:J_MA-1, initial number of males of different ages,
					   Initial number of males of different ages.  J0 is the number of newborn males at time MaxAge; J1 is the number of 
					   one-year-old males still alive at time MaxAge and so forth.  There should be MaxAge such Js because the first thing
					   to happen is an episode of death going from time MaxAge to time MaxAge+1; thus any MaxAge-year-olds would die before reproducing
					   anyway.  The last J is the number of MaxAge-1-year-old males at the outset (i.e. at time MaxAge).
					   ) ) {
		if( ALREADY_HAS(TFlag, -T or --number-of-years) && ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge) ) {
			/* allocate necessary memory */
			*NM = (int *)ECA_CALLOC(P->MaxAge+1+(*T),sizeof(int));
			*NMR = (int *)ECA_CALLOC(P->MaxAge+1+(*T),sizeof(int));
			*MSpace = (int *)ECA_CALLOC(P->MaxAge+1+(*T),sizeof(int));
			for(j=P->MaxAge-1;j>=0;j--) {
				(*NM)[j] = GET_INT;
				(*NMR)[j] = (*NM)[j];
			}
		}	
	}
	if(REQUIRED_OPTION(InitFemaleFlag,,initial-females,J0~J1:J_MA-1, initial number of females of different ages,
					   Initial number of females of different ages.  J0 is the number of newborn females at time MaxAge; J1 is the number of 
					   one-year-old females still alive at time MaxAge and so forth.  There should be MaxAge such Js because the first thing
					   to happen is an episode of death going from time MaxAge to time MaxAge+1; thus any MaxAge-year-olds would die before reproducing
					   anyway.  The last J is the number of MaxAge-1-year-old females at the outset (i.e. at time MaxAge).
					   ) ) {
		if( ALREADY_HAS(TFlag, -T or --number-of-years) && ALREADY_HAS(MaxAgeFlag, -A or --max-age) &&  ARGS_EQ(P->MaxAge) ) {
			/* allocate necessary memory */
			*NF = (int *)ECA_CALLOC(P->MaxAge+1+(*T),sizeof(int));
			*NFR = (int *)ECA_CALLOC(P->MaxAge+1+(*T),sizeof(int));
			*FSpace = (int *)ECA_CALLOC(P->MaxAge+1+(*T),sizeof(int));
			for(j=P->MaxAge-1;j>=0;j--) {
				(*NF)[j] = GET_INT;
				(*NFR)[j] = (*NF)[j];
			}
		}	
	}
	
	
	/********** POPULATION SIZE STUFF *******************/
	CLOSE_SUBSET;
	OPEN_SUBSET(Population Size Parameters and Option,
				PopSize,
				Parameters for setting the size of the population and features thereof);
	
	if(CLASHABLE_OPTION(RCSFlag , ,random-cohort-size, ,makes the cohort sizes random variables,
						random-cohort-sizes is the default.  This option takes
						no arguments.  Optional,
						FCSFlag,
						is used in conjunction with fixed-cohort-size) ) {
		P->CohortSizesRandom=1;		
	}
	if(CLASHABLE_OPTION(FCSFlag , ,fixed-cohort-size, ,makes the cohort sizes fixed quantities,
						random-cohort-sizes is the default.  This option takes
						no arguments.  Optional,
						RCSFlag,
						is used in conjunction with random-cohort-size) ) {
		P->CohortSizesRandom=0;		
	}
	if( REQUIRED_OPTION(Cohort_Size_Flag,
						,
						cohort-size,
						C{const[J]|var[J1~J2:JT$80x5$]}, 
						cohort size determined by model C with pars ..., 
						This option is used to specify how the cohort size is determined.  The 
						different models for determining cohort size are specified by the choice string C.
						Each model might take a different set of parameters. Listed below are the choices
						C for each model followed by the set of extra arguments to the option with a description
						of how each model works.\n\n\tconst  J   ---   this makes the cohort size (or at least
																								   the _expected_ cohort size (if using the --random-cohort-size option)) J every
						year regardless of how large the parental population is.\n\n\tvar J1 J2...JT  ---  this makes the 
						expected cohort size each year from MaxAge to T+MaxAge-1 equal to J1 J2...JT respectively. Hence there
						must be T such arguments J following the var option (where T is the number of years that the simulation will run---the argument to the
																			 -T or the --number-of-years option.)
						)  )  {
		GET_STR(Pop->CohortSizeMode);
		if(strcmp(Pop->CohortSizeMode,"const")==0) {
			if(ARGS_EQ(1)) {
				*N = GET_INT;
			}		
			else {
				fprintf(stderr,"\tThe \"const\" mode specifier should be followed by a single integer\n");
				OPT_ERROR;
			}
		}
		else if(strcmp(Pop->CohortSizeMode,"var")==0) {
			if(ALREADY_HAS(MaxAgeFlag, -A or --max-age) && ALREADY_HAS(TFlag, -T or --number-of-years )) {
				if(ARGS_EQ( *T ) ) {
					Pop->CohortSizes = (int *)ECA_CALLOC(*T + P->MaxAge, sizeof(int));
					for(j=P->MaxAge; j< (*T)+P->MaxAge; j++) {
						Pop->CohortSizes[j] = GET_INT;
					}
				}
				else {
					printf("when using the \"var\" specifier.\n");
				}
			}
		}
		else {
			fprintf(stderr,"Error! Unrecognized argument \"%s\" to option --cohort-size\n",Pop->CohortSizeMode);
			OPT_ERROR;
		}
	}
	
	
	
	
	
	/************* GENETIC SAMPLING STUFF ********************/
	CLOSE_SUBSET;
	OPEN_SUBSET(Genetic Sampling Parameters,
				Genetics,
				Parameters controlling initial allele frequencies and the pattern of genetic sampling);
	
	/* Allow a warm-up period, so that you can scrub the previously simulated pedigree */
	if (CLASHABLE_OPTION(DiscardAll_Flag,
						 , 
						 discard-all, 
						 J , 
						 make all pedigree members born b/t J and J+MaxAge-1 the founders, 
						 This option renders all individuals born between the years of J and J+MaxAge-1 (inclusive)
						 founders of the pedigree.  This means that all those individuals will be regarded as having
						 unknown parents.  (i.e. their mother and father fields in the PEDIGREE will take the value
											0.)  The pedigree that would have been simulated without the discard-all option gets printed 
						 on the DISCARDED_PDGEE lines of output.  When spip simulates genetic data all the individuals
						 born between J and J+MaxAge-1 (inclusive) are considered the founders. Optional.  The
						 default is --discard-before 0. , 
						 DiscardParents_Flag,
						 is used at same time as --discard-par
						 ) ) {
		if( ALREADY_HAS(TFlag, -T or --number-of-years) && ALREADY_HAS(MaxAgeFlag,-A or --max-age) && ARGS_EQ(1) ) {
			P->DiscardBefore = GET_INT;
			P->TossAll = 1;
			/* Test to see if it is an appropriate value */
			if(P->DiscardBefore > (*T + P->MaxAge + 1) ) {
				fprintf(stderr,"Error! Argument of --discard-all = %d is greater than arg of -T (or --number-of-years) + MaxAge - 1 = %d\n",
						P->DiscardBefore,*T + P->MaxAge - 1);
				OPT_ERROR;
			}
			if(P->DiscardBefore < 0) {
				fprintf(stderr,"Error! Argument of --discard-all = %d is less than zero\n",P->DiscardBefore);
				OPT_ERROR;
			}
		}
	}
	
	/*		if (CLASHABLE_OPTION(DiscardParents_Flag,
	 , 
	 discard-parents, 
	 [J] ,
	 pedigree genes from parents born before J are founders, 
	 This option discards in the pedigree any parents who were born before time J.  If a parent (let us say
	 a father) born before
	 time J produces an offspring then that child will be considered to have an unkown father.  (i.e. in the 
	 pedigree the child will have a 0 for its father but not necessarily for its mother.) The pedigree 
	 that would have been simulated without the discard-parents option gets printed 
	 on the DISCARDED_PDGEE lines of output.  There will be discarded pedigree outputs until the cohort
	 born at time J + MaxAge because after that time it is not possible for a child to have a parent that
	 was born before time J.  When spip simulates genetic data the paternal allele founder is the child
	 if the father was discarded.  This option is not terribly useful.  Unfortunately the programs of the 
	 MORGAN package cannot deal with pedigree with only one parent having a 0.  It is possible to get a similar
	 effect by using the --discard-all option and running the simulation a little longer.  , 
	 DiscardAll_Flag,
	 is used at same time as --discard-all
	 ) ) {
	 if( ALREADY_HAS(TFlag, -T or --number-of-years) && ALREADY_HAS(MaxAgeFlag,-A or --max-age) && ARGS_EQ(1) ) {
	 P->DiscardBefore = GET_INT;
	 P->TossAll = 0;  */
	/* Test to see if it is an appropriate value */
	/*				if(P->DiscardBefore > *T) {
	 fprintf(stderr,"Error! Argument of --discard-parents = %d is greater than arg of -T (or --number-of-years) = %d\n",P->DiscardBefore,*T);
	 OPT_ERROR;
	 }
	 if(P->DiscardBefore < 0) {
	 fprintf(stderr,"Error! Argument of --discard-parents = %d is less than zero\n",P->DiscardBefore,*T);
	 OPT_ERROR;
	 }
	 else {
	 printf("DISCARD_BEFORE : %d\n",P->DiscardBefore);
	 }
	 }
	 }
	 */
	
	if(CLASHABLE_OPTION(GtypAPropFemPre_Flag, 
						,
						gtyp-ppn-fem-pre,
						G R1~R2:R_MA,
						genotype a proportion of females randomly prior to death,
						this causes the program to randomly sample living females
						in the years specified by range G before the episode of death and survival.   The probability 
						that a living female is sampled is age-dependent.  R1 is the probability 
						that a female born in the previous year (but before any death episode) is 
						sampled.  Note that this is equivalent to sampling newborns before any of
						them have the chance to die.  R2 is the probability that a two-year
						old is sampled.  There should be MaxAge  such arguments. The last one is
						the probability that a MaxAge-year old is sampled before the death episode
						each year.  Years specified in G which are before the year specified in 
						--discard-all or which are after T + MaxAge - 1 will be ignored. ,
						NewbornSampFlag,
						is used in conjunction with --newborn-samp) ) {
		if(REQ_FOR_SAMPLING && ARGS_EQ(P->MaxAge+1) )	{ char therange[10000];
			GET_STR(therange);
			S->pre_fem_samp_years = StringToIntArray(therange, P->DiscardBefore, *T + P->MaxAge-1, 1);
			if(S->pre_fem_samp==NULL) S->pre_fem_samp = (double *)ECA_CALLOC(P->MaxAge+1, sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				S->pre_fem_samp[j] = GET_DUB;
			}
		}	
	}
	if(CLASHABLE_OPTION(GtypAPropMalePre_Flag, 
						,
						gtyp-ppn-male-pre,
						G R1~R2:R_MA,
						genotype a proportion of males randomly prior to death,
						this causes the program to randomly sample living males
						in the years specified by range G before the episode of death and survival.   The probability 
						that a living male is sampled is age-dependent.  R1 is the probability 
						that a newborn male of the previous year is sampled before it has a chance
						to die.  R2 is the probability that a two-year
						old is sampled.  There should be MaxAge such arguments. The last one is
						the probability that a MaxAge-year old is sampled before the death episode
						each year. Years specified in G which are before the year specified in 
						--discard-all or which are after T + MaxAge - 1 will be ignored. ,
						NewbornSampFlag,
						is used in conjunction with --newborn-samp) ) {
		if(REQ_FOR_SAMPLING && ARGS_EQ(P->MaxAge+1))	{ char therange[10000];
			GET_STR(therange);
			S->pre_male_samp_years = StringToIntArray(therange, P->DiscardBefore, *T + P->MaxAge-1, 1);
			if(S->pre_male_samp==NULL) S->pre_male_samp = (double *)ECA_CALLOC(P->MaxAge+1, sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				S->pre_male_samp[j] = GET_DUB;
			}
		}	
	}
	if(CLASHABLE_OPTION(GtypAPropFemPost_Flag, 
						,
						gtyp-ppn-fem-post,
						G R1~R2:R_MA,
						genotype a proportion of females randomly after death episode,
						this causes the program to randomly sample living females
						in the years specified by range G after the episode of death and survival.   The probability 
						that a living female is sampled is age-dependent.  R1 is the probability 
						that a one-year-old female is sampled.  R2 is the probability that a two-year
						old is sampled.  There should be MaxAge-1 such arguments. The last one is
						the probability that a MaxAge-1 year-old is sampled after the death episode
						each year.  Years specified in G which are before the year specified in 
						--discard-all or which are after T + MaxAge - 1 will be ignored. ,
						NewbornSampFlag,
						is used in conjunction with --newborn-samp) ) {
		if(REQ_FOR_SAMPLING && ARGS_EQ(P->MaxAge+1) )	{ char therange[10000];
			GET_STR(therange);
			S->post_fem_samp_years = StringToIntArray(therange, P->DiscardBefore, *T + P->MaxAge-1, 1);
			if(S->post_fem_samp==NULL) S->post_fem_samp = (double *)ECA_CALLOC(P->MaxAge+1, sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				S->post_fem_samp[j] = GET_DUB;
			}
		}	
	}
	if(CLASHABLE_OPTION(GtypAPropMalePost_Flag, 
						,
						gtyp-ppn-male-post,
						G R1~R2:R_MA,
						genotype a proportion of males randomly after death episode,
						this causes the program to randomly sample living males
						in the years specified by range G after the episode of death and survival.   The probability 
						that a living male is sampled is age-dependent.  R1 is the probability 
						that a one-year-old male is sampled.  R2 is the probability that a two-year
						old is sampled.  There should be MaxAge-1 such arguments. The last one is
						the probability that a MaxAge-1 year-old is sampled after the death episode
						each year.  Years specified in G which are before the year specified in 
						--discard-all or which are after T + MaxAge - 1 will be ignored. ,
						NewbornSampFlag,
						is used in conjunction with --newborn-samp) ) {
		if(REQ_FOR_SAMPLING && ARGS_EQ(P->MaxAge+1) )	{  char therange[10000];
			GET_STR(therange);
			S->post_male_samp_years = StringToIntArray(therange, P->DiscardBefore, *T + P->MaxAge-1, 1);
			if(S->post_male_samp==NULL) S->post_male_samp = (double *)ECA_CALLOC(P->MaxAge+1, sizeof(double));
			for(j=1;j<=P->MaxAge;j++)  {
				S->post_male_samp[j] = GET_DUB;
			}
		}	
	}
	if(CLASHABLE_OPTION(
						GtypAPropFemDur_Flag,
						,
						gtyp-ppn-fem-dur,
						G R,
						genotype a proportion R of reproducing females,
						This causes the program to randomly sample (with probability R) females
						who are reproducing in the years specified by range G.  This is particularly useful for simulating species which
						are easily sampled when they congregate to spawn or reproduce. Years specified in G which are before the year specified in 
						--discard-all or which are after T + MaxAge - 1 will be ignored. , 
						NewbornSampFlag, 
						is used in conjunction with --newborn-samp ) ) {
		if(REQ_FOR_SAMPLING && ARGS_EQ(2)) { char therange[10000];
			GET_STR(therange);
			S->dur_fem_samp_years = StringToIntArray(therange, P->DiscardBefore, *T + P->MaxAge-1, 1);
			S->dur_repro_fem_samp = GET_DUB;
		}	
	}
	if(CLASHABLE_OPTION(
						GtypAPropMaleDur_Flag,
						,
						gtyp-ppn-male-dur,
						G R,
						genotype a proportion R of reproducing males,
						This causes the program to randomly sample (with probability R) males
						who are reproducing in the years specified by range G.  This is particularly useful for simulating species which
						are easily sampled when they congregate to spawn or reproduce. Years specified in G which are before the year specified in 
						--discard-all or which are after T + MaxAge - 1 will be ignored. , 
						NewbornSampFlag, 
						is used in conjunction with --newborn-samp ) ) {
		if(REQ_FOR_SAMPLING && ARGS_EQ(2)) { char therange[10000];
			GET_STR(therange);
			S->dur_male_samp_years = StringToIntArray(therange, P->DiscardBefore, *T + P->MaxAge-1, 1);
			S->dur_repro_male_samp = GET_DUB;
		}	
	}
	if(CLASHABLE_OPTION(GenosForAll_Flag, 
						,
						gtyps-for-all,
						G,
						draw genotypes for all individuals born at times in range G,
						this causes genotypes to be simulated for all individuals born into 
						cohorts at the times specified in range G.  G must be a list of nonnegative
						numbers in increasing order separated by commas and dashes as described above. Note that 
						as currently implemented this option overwrites
						the  --gtyp-ppn-fem-pre and --gtyp-ppn-male-pre options.  Also incompatible with the --newborn-samp
						option.,
						NewbornSampFlag,
						is used in conjunction with --newborn-samp) ) {
		if(REQ_FOR_SAMPLING && ARGS_EQ(1) )	{char therange[10000];
			GET_STR(therange);
			S->pre_male_samp_years = StringToIntArray(therange, P->DiscardBefore, *T + P->MaxAge-1, 1);
			S->pre_fem_samp_years = StringToIntArray(therange, P->DiscardBefore, *T + P->MaxAge-1, 1);
			/* this is equivalent to setting the prob that a zero-year old is sampled before the death episode to 1.0 */
			if(S->pre_male_samp==NULL) S->pre_male_samp = (double *)ECA_CALLOC(P->MaxAge+1, sizeof(double));
			if(S->pre_fem_samp==NULL) S->pre_fem_samp = (double *)ECA_CALLOC(P->MaxAge+1, sizeof(double));
			S->pre_fem_samp[1] = 1.0;
			S->pre_male_samp[1] = 1.0;
		}	
	}
	
	if(OPTION(
			LethalSampling_Flag,
			,
			lethal-sampling,
			R,
			makes the -pre and -post sampling lethal with probability R,
			This is somewhat experimental.  The purpose of this option is to kill off with probability R individuals that
			  have been sampled.  I believe this should apply only to those sampled from the --gtyp...-pre and --gtyp...-post
			  options and not to those sampled during reproduction not those sampled as newborns nor sampled
			  using the --gtyps-for-all option (but I am not sure of that at this point).  Use with extreme care.  This is sort of 
			for the developer only at this point.)) {
		if(ARGS_EQ(1)) {
			S->LethalityProb = GET_DUB;
		}
			
	}
	/*
	 if( OPTION(NewbornSampFlag,,newborn-samp,[J1] [J2] ... ,the years in which to sample newborns, The
	 years in which to sample newborns.  All the Js should be less than or equal to T+MaxAge-1.  The
	 Js should also be greater than or equal to the value given in the --discard-all option.  Optional. Default is none. )) {
	 S->nc = COUNT_ARGS;
	 if(ALREADY_HAS(MaxAgeFlag, -A or --max-age)) { int bail=0;
	 if(S->nc == 0) {
	 fprintf(stderr,"Error! --newborn-samp option given with no arguments.\n");
	 OPT_ERROR;
	 }
	 S->SampleNewbornsOnly = 1;
	 S->ts = (int *)ECA_CALLOC(S->nc,sizeof(int));
	 for(j=0;j<S->nc;j++) {
	 S->ts[j] = GET_INT;
	 if(S->ts[j] > P->MaxAge + (*T) - 1) {
	 fprintf(stderr,"Error! Requesting newborn sampling at t=%d which is greater than MaxAge+T-1=%d\n",S->ts[j], P->MaxAge + (*T) - 1 );
	 bail=1;
	 }
	 }
	 if(bail) {
	 OPT_ERROR;
	 }
	 }
	 }
	 if (COND_REQ_OPTION(MaleSampFlag,, newborn-males-samp, [J1] [J2] ... ,how many males sampled from sampled cohorts, , NewbornSampFlag && !FemaleSampFlag  , --newborn-samp   and --females-sampled is not used) ) {
	 if(ALREADY_HAS(NewbornSampFlag, newborn-samp) && ARGS_EQ(S->nc) ) {
	 S->NumM = (int *)ECA_CALLOC(S->nc,sizeof(int));
	 for(j=0;j<S->nc;j++) {
	 S->NumM[j] = GET_INT;
	 }
	 }
	 }
	 if (COND_REQ_OPTION(FemaleSampFlag,
	 , 
	 newborn-fem-samp, 
	 [J1] [J2] ... ,
	 how many females sampled from sampled cohorts, 
	 , NewbornSampFlag && !MaleSampFlag 
	 , --newborn-samp and --males-sampled is not used
	 ) ) {
	 if(ALREADY_HAS(NewbornSampFlag, newborn-samp) && ARGS_EQ(S->nc) ) {
	 S->NumF = (int *)ECA_CALLOC(S->nc,sizeof(int));
	 for(j=0;j<S->nc;j++) {
	 S->NumF[j] = GET_INT;
	 }
	 }
	 }
	 */
	if (COND_REQ_OPTION(LocFileFlag,, locus-file, F ,F=pathname to file with locus information,
						F=pathname to file with locus information.  The format of this allele-frequency information
						file is simple.  It includes only integers (or comments enclosed by a pair of & characters).
						The first integer is the number of loci in the file.  Then for each locus you give the
						number of alleles at that locus followed by the allele frequencies of each allele at the locus.
						All these things should be separated by whitespace.  No punctuation!  If the allele frequencies
						at a locus do not sum to one they will be rescaled so that the do sum to one.  If the number of
						allele frequencies listed does not match the number of alleles given for the locus then the program
						may fail without warning or give otherwise unexpected results.  The allele frequencies in this 
						file are assumed to be the allele frequencies in all years before the year given by the
						year-of-no-ibd option.  , 
						DOING_GENETIC_SAMPLING  , 
						when using --newborn-samp or any other genetic sampling options ) ) {
		if(  ARGS_EQ(1) ) {
			GET_STR(locfilename);
			*HasLocFile = 1;
		}
	}
	
	CLOSE_SUBSET;
	
	OPEN_SUBSET(Program Implementation Options,
				Program Implementation,
				Program Implementation Options);
	if(OPTION(
			  alloc_extra_f,
			  ,
			  alloc-extra,
			  J,
			  Integer factor by which to increase memory allocation.,
			  Integer factor by which to increase memory allocation to offspring arrays.  When memory gets allocated to the
			  pointers for the next cohort\054 this is the factor by which the number of spaces is 
			  greater than the expected number.  Basically it determines mspaces and fspaces in function MakeBabies.  By default it will be 
			  set to 10.  But should you be doing lots of generations and not wanting to eat up so much memory\054 you can set it smaller 
			  (like 2 for the simulations Oystein was doing).  Note that ultimately I should just realloc all those things after the 
			  offspring are created (right at the end of MakeBabies() to reclaim that space). But\054 until I do that and debug it
			  we have to suffer along with this.
			) ) {
		if(ARGS_EQ(1)) {
			P->OffspringAllocOverload = GET_INT;
		}
		
	}	
	
	CLOSE_SUBSET;
	
	
	/* here are a couple of little things to check each time to make sure that we haven't gotten any clashing
	 options (that are too complex to catch with the CLASHING_OPTION macro) */
	if( __OptSeekUnusedRequiredOpts == 1 && (
											 (SurvivalRatesFlag + FemSurvRatesFlag) > 1 ||  (SurvivalRatesFlag + MaleSurvRatesFlag) > 1 ) ) {
		fprintf(stderr,"Error! You cannot issue either the fem-surv-probs or the male-surv-probs option along with the survival-probs option.\n");
		OPT_ERROR
	} 
	
	END_OPT_LOOP
	
	/* a little bit more error checking down here */
	if(P->CohortSizesRandom==0 && P->OffsDsn==BINARY)  {
		fprintf(stderr,"Error! The \"binary\" argument to the --offsp-dsn option cannot be used if the --fixed-cohort-size option is in effect.\nExiting\n");
		exit(1);
	}
	if(	P->ReproInhib > 0 && P->OffsDsn != BINARY)  {
		fprintf(stderr,"Error! The --repro-inhib option cannot be used with [J] > 0 if the \"binary\" argument to the --offsp-dsn is not in effect.\nExiting\n");
		exit(1);
	}
	if(	P->ReproInhib > 0 && P->CohortSizesRandom==0)  {
		fprintf(stderr,"Error! The --repro-inhib option cannot be used with [J] > 0 if the --fixed-cohort-size option is in effect..\nExiting\n");
		exit(1);
	}
	return(NewPop_Flag);
}






















