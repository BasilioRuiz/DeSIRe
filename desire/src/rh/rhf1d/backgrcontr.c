/* ------- file: -------------------------- backgrcontr.c -----------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr 22 09:45:42 2009 --

       --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "background.h"
#include "inputs.h"
#include "error.h"

#define CONTR_INPUT_FILE "contribute.input"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

// 31/07/19 epm: Global variables definition is now in "rh_glob.c".
extern Atmosphere atmos;
extern Geometry geometry;
extern char messageStr[];


/* ------- begin -------------------------- backgrcontr.c ----------- */

int main(int argc, char *argv[])
{
  register int n;

  int     Nlambda, result;
  double *lambda;
  FILE   *fp;

  /* --- Read input data and initialize --             -------------- */

  setOptions(argc, argv);
  SetFPEtraps();

  readInput();
  MULTIatmos(&atmos, &geometry);
  readAtomicModels();
  readMolecularModels();

  if ((fp = fopen(CONTR_INPUT_FILE, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", CONTR_INPUT_FILE);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }

  result = fscanf(fp, "%d", &Nlambda);
  lambda = (double *) malloc(Nlambda * sizeof(double));
  for (n = 0;  n < Nlambda;  n++)
    result = fscanf(fp, "%lf", &lambda[n]);
  fclose(fp);

  backgrOpac(Nlambda, lambda);

  free(lambda);
}

/* ------- end ---------------------------- backgrcontr.c ----------- */
