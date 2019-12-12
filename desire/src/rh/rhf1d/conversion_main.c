/* ------- file: --------------------- conversion_main.c ------------ */

/* 19/06/19: Created by Esperanza Paez (epm@iac.es)
 * Adaptation to the DeSIRe project - Interoperability C/Fortran.
 * Conversion of main() as conversion().
 */

#include <stdio.h>


/* --- Function prototypes --                          -------------- */

int conversion(int argc, char *argv[],
               int *nspace, int *nsize, double *tau, double *mass,
               double *z, double *T, double *ne, double *vel, double *vturb);


/* ------- begin --------------------- conversion_main.c ------------ */

int main( int argc, char *argv[] )
{
   int nspace, nsize = 0;
   return(conversion(argc, argv, &nspace, &nsize,
                     NULL, NULL, NULL, NULL, NULL, NULL, NULL));
}

/* ------- end ----------------------- conversion_main.c ------------ */
