/* ------- file: -------------------------- readb_xdr.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan 19 16:51:09 2012 --

       --------------------------                      ----------RH-- */

/* --- Reads the atmospheric magnetic field.
       XDR (external data representation) version. --  -------------- */


/* --- Read the absolute value of the magnetic field, and the angles
       gamma (the angle of the field with the normal direction to the
       atmosphere) and chi (the angle of the field with the positive
       x-axis, measured counter-clockwise).

 Note: Magnetic field strength should be given in Tesla (1 T = 10^4 Gauss) 
       Angles should be given in radians.

 Modifications:

       - 11/11/20 epm:
         We supply the atmospheric magnetic field from SIR.

       --                                              -------------- */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "desire.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern InputData input;
extern char messageStr[];

// 11/11/20 epm: Structure to load the atmospheric magnetic field from SIR.
extern SIRatmos siratmos;


/* ------- begin -------------------------- readB.c ----------------- */

bool_t readB(Atmosphere *atmos)
{
  const char routineName[] = "readB";

  bool_t  result = TRUE;
  size_t  recordsize;
  off_t   offset;
  FILE   *fp_stokes;
  XDR     xdrs;
  int     i, nitems;

  // if (!strcmp(input.Stokes_input, "none")) return FALSE;

  if (siratmos.nB == 0) return FALSE;
  if (siratmos.nB != atmos->Nspace)
  {
    sprintf(messageStr, "Wrong number of depths from SIR: %d", siratmos.nB);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  // if ((fp_stokes = fopen(input.Stokes_input, "r")) == NULL) {
  //   sprintf(messageStr, "Unable to open inputfile %s", input.Stokes_input);
  //   Error(ERROR_LEVEL_2, routineName, messageStr);
  // }

  atmos->B = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->gamma_B = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->chi_B   = (double *) malloc(atmos->Nspace * sizeof(double));

  // if (input.xdr_endian) {
  //   xdrstdio_create(&xdrs, fp_stokes, XDR_DECODE);
  //
  //   result &= xdr_vector(&xdrs, (char *) atmos->B, atmos->Nspace,
  //                        sizeof(double), (xdrproc_t) xdr_double);
  //   result &= xdr_vector(&xdrs, (char *) atmos->gamma_B, atmos->Nspace,
  //                        sizeof(double), (xdrproc_t) xdr_double);
  //   result &= xdr_vector(&xdrs, (char *) atmos->chi_B, atmos->Nspace,
  //                        sizeof(double), (xdrproc_t) xdr_double);
  //   xdr_destroy(&xdrs);
  // } else {
  //   recordsize = atmos->Nspace;
  //
  //   for (result = TRUE, i = 0; i < recordsize; i++) {
  //     nitems = fscanf(fp_stokes, "%lf %lf %lf",
  //                     &atmos->B[i], &atmos->gamma_B[i], &atmos->chi_B[i]);
  //     if (nitems != 3) {
  //       result = FALSE;
  //       break;
  //     }
  //   }
  // }
  //
  // if (!result) {
  //   sprintf(messageStr, "Unable to read from input file %s",
  //           input.Stokes_input);
  //   Error(ERROR_LEVEL_1, routineName, messageStr);
  // }
  //
  // fclose(fp_stokes);

  for (i = 0; i < atmos->Nspace; i++)
  {
    atmos->B[i] = siratmos.B[i];
    atmos->gamma_B[i] = siratmos.gamma[i];
    atmos->chi_B[i] = siratmos.phi[i];
  }

  return result;
}
/* ------- end ---------------------------- readB.c ----------------- */
