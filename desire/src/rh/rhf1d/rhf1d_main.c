/* ------- file: --------------------- rhf1d_main.c ----------------- */

/* 15/05/19: Created by Esperanza Paez (epm@iac.es)
 * Adaptation to the DeSIRe project - Interoperability C/Fortran.
 * Conversion of main() as rhf1d().
 */

#include <stdio.h>


/* --- Function prototypes --                          -------------- */

int rhf1d(int argc, char *argv[],
          int *nspace, int *nlevel, int *nsize, double *coefs);


/* ------- begin --------------------- rhf1d_main.c ----------------- */

int main( int argc, char *argv[] )
{
   int nspace, nlevel, nsize = 0;
   return(rhf1d(argc, argv, &nspace, &nlevel, &nsize, NULL));
}

/* ------- end ----------------------- rhf1d_main.c ----------------- */
