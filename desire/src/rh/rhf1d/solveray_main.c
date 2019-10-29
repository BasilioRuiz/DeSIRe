/* ------- file: --------------------- solveray_main.c -------------- */

/* 15/05/19: Created by Esperanza Paez (epm@iac.es)
 * Adaptation to the DeSIRe project - Interoperability C/Fortran.
 * Conversion of main() as solveray().
 */

#include <stdio.h>


/* --- Function prototypes --                          -------------- */

int solveray(int argc, char *argv[], int *nspect, int *nsize,
             double *lambda, double *I, double *Q, double *U, double *V);


/* ------- begin --------------------- solveray_main.c -------------- */

int main( int argc, char *argv[] )
{
   int nspect, nsize = 0;
   return(solveray(argc, argv, &nspect, &nsize,
                   NULL, NULL, NULL, NULL, NULL));
}

/* ------- end ----------------------- solveray_main.c -------------- */
