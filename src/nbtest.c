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


int main(void) {
	double mu = 0.005817;
	double Fcv = 0.7;
	int i;
	int RV;

	for(i=0;i<=500000000; i++) {
		RV = NegBinomAB_RV(mu, Fcv);
	 	printf("%d\n", RV);
	}

	return(0);
	
}