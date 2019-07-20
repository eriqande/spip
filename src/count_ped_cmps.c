/*************************************************************************
 count_ped_cmps.c

 This program is implemented in MORGAN_V2.4.1, with parameter statements
 added for kin.

 For each pedigree component:
   Lists the number of the component and the ID's of the individuals in that
   component.
   
                                                Jan 2004  ECA

***********************************************************************/

#define PNAM "kin"
#include "pedcomp.h"

/**
 * @defgroup kin kin
 * @ingroup PedComp
 * @{
 */

int   main (int argc, char *argv[])
{
   int   c, i, j, nfou, wid;
   double x, g4;
   
   struct _Indlist *eca_temp;

   Ind  *pr;

   FILE *pedF;

   atexit (print_tail);
   PCleanUp = clean_check;
   atexit (clean_morg);
   print_head (*argv);
   if (!proc_pars (argc, argv, PNAM))   /* process input parameter statements */
      exit (EXIT_FAILURE);

   if (!proc_kin ())                    /* check parameter statements;        */
      exit (EXIT_FAILURE);              /* prepare controls                   */

   if (!proc_fils ())                   /* set up pedigree file parameters    */
      exit (EXIT_FAILURE);

   if (!chk_pedfil ())                  /* reads & checks pedigree file for   */
      give_up (ERR_qit);                /* errors                             */

   if (!mak_opedfil ())                 /* write new pedigree file, if there  */
      give_up (ERR_qit);                /* are changes                        */

   if (!(pedF = open_ped ()))           /* reopen pedigree file               */
      give_up (ERR_qit);

   alloc_ped ();
   input_ped (pedF, 0, NULL, NULL);     /*input pedigree                      */
   fclose (pedF);

   nfou = count_founders ();
   printf (" Number of founders is %d\n\n", nfou);
   printf (" Number of components:  %d\n\n", n_compnts);

   printf("GIVING LISTING OF COMPONENTS AND INDIVIDUALS WITHIN EACH\n\n");
   for(i=1;i<=n_compnts;i++)  {
   	  printf("CPC: Component %d\nCPC:",i);
   	  eca_temp = cmplist[i].people;
   	  while (eca_temp != NULL) {
   	  	printf(" %s",eca_temp->current->self);
   	  	eca_temp = eca_temp->next;
   	  } 
   	  printf("\n");
   }
	
   exit (EXIT_SUCCESS);
}


/** @} */
