/* ------- file: -------------------------- writespect_xdr.c --------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 11 12:18:36 2021 --

       --------------------------                      ----------RH-- */

/* --- Writes spectroscopic data to output file.
       XDR (external data representation) version.

       Modifications:

       - 04/04/20 epm:
         Disk access is avoided as much as possible.
         J_FILE is saved in memory rather than disk.
         The keyword LIMIT_MEMORY has to be set to FALSE.
         Be aware that 'rhf1d' and 'solveray' need the keyword
         STARTING_J to be set to NEW_J as J_FILE is missing.

       --                                              -------------- */

#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern enum Topology topology;

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeSpectrum.c --------- */

void writeSpectrum(Spectrum *spectrum)
{
  const char routineName[] = "writeSpectrum";
  register int nspect, nact;

  bool_t     result = TRUE;
  int        Nintensity;
  double    *lambda_air, vacuum_to_air_limit = VACUUM_TO_AIR_LIMIT;
  FILE      *fp_spectrum;
  XDR        xdrs;
  ActiveSet *as;

  if (strcmp(input.spectrum_output, "none") != 0)
  {
    if ((fp_spectrum = fopen(input.spectrum_output, "w")) == NULL) {
      sprintf(messageStr, "Unable to open output file %s",
              input.spectrum_output);
      Error(ERROR_LEVEL_1, routineName, messageStr);
      return;
    }
    xdrstdio_create(&xdrs, fp_spectrum, XDR_ENCODE);

    result &= xdr_int(&xdrs, &spectrum->Nspect);

    if (spectrum->vacuum_to_air) {
      lambda_air = (double *) malloc(spectrum->Nspect * sizeof(double));
      vacuum_to_air(spectrum->Nspect, spectrum->lambda, lambda_air);
      result &= xdr_vector(&xdrs, (char *) lambda_air, spectrum->Nspect,
                           sizeof(double), (xdrproc_t) xdr_double);
      free(lambda_air);
    } else
      result &= xdr_vector(&xdrs, (char *) spectrum->lambda, spectrum->Nspect,
                           sizeof(double), (xdrproc_t) xdr_double);

    switch (topology) {
    case ONE_D_PLANE:
      Nintensity = atmos.Nrays * spectrum->Nspect;
      break;
    case TWO_D_PLANE:
      Nintensity = atmos.N[0] * atmos.Nrays * spectrum->Nspect;
      break;
    case THREE_D_PLANE:
      Nintensity = atmos.N[0]*atmos.N[1] * atmos.Nrays * spectrum->Nspect;
      break;
    case SPHERICAL_SYMMETRIC:
      Nintensity = atmos.Nrays * spectrum->Nspect;
      break;
    default:
      Nintensity = 0;
      break;
    }
    result &= xdr_vector(&xdrs, (char *) spectrum->I[0], Nintensity, 
                         sizeof(double), (xdrproc_t) xdr_double);

    result &= xdr_bool(&xdrs, &spectrum->vacuum_to_air);
    result &= xdr_double(&xdrs, &vacuum_to_air_limit);

    if (atmos.Stokes || input.backgr_pol) {
      result &= xdr_vector(&xdrs, (char *) spectrum->Stokes_Q[0], Nintensity, 
                           sizeof(double), (xdrproc_t) xdr_double);
      result &= xdr_vector(&xdrs, (char *) spectrum->Stokes_U[0], Nintensity, 
                           sizeof(double), (xdrproc_t) xdr_double);
      result &= xdr_vector(&xdrs, (char *) spectrum->Stokes_V[0], Nintensity, 
                           sizeof(double), (xdrproc_t) xdr_double);
    }

    if (!result) {
      sprintf(messageStr, "Unable to write proper amount to output file %s",
              input.spectrum_output);
      Error(ERROR_LEVEL_1, routineName, messageStr);
    }
    xdr_destroy(&xdrs);
    fclose(fp_spectrum);
  }

  /* --- Write angle-averaged mean intensity to file -- ------------- */

  if (spectrum->updateJ && !input.limit_memory) {
    // 04/04/20 epm: From now on, we save J in memory rather than disk.
    // if ((spectrum->fd_J =
    //      open(input.JFile, O_WRONLY | O_CREAT, PERMISSIONS)) == -1) {
    //   sprintf(messageStr, "Unable to open output file %s", input.JFile);
    //   Error(ERROR_LEVEL_1, routineName, messageStr);
    //   return;
    // }
    for (nspect = 0;  nspect < spectrum->Nspect;  nspect++)
      writeJlambda(nspect, spectrum->J[nspect]);

    // close(spectrum->fd_J);

    /* --- Write the anisotropy J^2_0 in the z-direction -- --------- */

    if (input.backgr_pol) {
      if ((spectrum->fd_J20 =
           open(J20_DOT_OUT, O_WRONLY | O_CREAT, PERMISSIONS)) == -1) {
        sprintf(messageStr, "Unable to open output file %s", J20_DOT_OUT);
        Error(ERROR_LEVEL_1, routineName, messageStr);
        return;
      }
      for (nspect = 0;  nspect < spectrum->Nspect;  nspect++)
        writeJ20lambda(nspect, spectrum->J20[nspect]);

      close(spectrum->fd_J20);
    }
  }

  /* --- Clean up memory allocations --                -------------- */

  free(spectrum->lambda);

  for (nspect = 0; nspect < spectrum->Nspect; nspect++) {
    as = &spectrum->as[nspect];

    for (nact = 0; nact < atmos.Nactiveatom; nact++)
    {
      if (as->Nactiveatomrt[nact] > 0)
        {
          free(as->art[nact]);
          free(as->lower_levels[nact]);
          free(as->upper_levels[nact]);
      }
    }
    /* --- Reallocate space for the molecular transition arrays -- -- */

    for (nact = 0; nact < atmos.Nactivemol; nact++)
    {
      if (as->Nactivemolrt[nact] > 0)
      {
        free(as->mrt[nact]);
      }
    }
  }
}
/* ------- end ---------------------------- writeSpectrum.c --------- */
