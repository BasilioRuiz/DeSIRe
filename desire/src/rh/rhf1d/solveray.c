/* ------- file: -------------------------- solveray.c --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May  1 09:57:08 2020 --

       Updates (epm): Adaptation to the DeSIRe project.
                      Interoperability C-Fortran.
                      Conversion of main() as a function for Fortran.

       --------------------------                      ----------RH-- */

/* --- Solves radiative transfer for given atmosphere and model atom
       along a ray with arbitrary \mu_z, assuming the atom's population
       numbers and angle-averaged radiation field is given.

       Expects input file "ray.input" containing two lines of the form
         muz
         Nsource  wave_index1  ....   wave_indexNsource

       Modifications:

       - 04/04/20 epm:
         Disk access is avoided as much as possible.
         BACKGROUND_FILE is saved in memory rather than disk.
         Two keywords have to be set: OLD_BACKGROUND=FALSE
         (input.old_background=FALSE in "background.c") and,
         VMACRO_TRESH=0.0 (atmos->moving=TRUE in "multiatmos.c").

         We don't need "spectrum_%4.2f" anymore.

         Cosine of heliocentric angle as input argument to solveray().

       --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"


#define COMMENT_CHAR        "#"

#define RAY_INPUT_FILE      "ray.input"
#define ASCII_SPECTRUM_FILE "spectrum_%4.2f.asc"


/* --- Function prototypes --                          -------------- */

void closefiles(void);

int solveray(int argc, char *argv[], int *nspect, int *nsize, double *lambda,
             double *I, double *Q, double *U, double *V, double *xmu);


/* --- Global variables --                             -------------- */

// 31/07/19 epm: Global variables definition is now in "rh_glob.c".
extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- solveray.c -------------- */

// 29/05/19 epm: Function to be called from Fortran.
int solveray_( int *nspect, int *nsize, double *lambda, double *stokes_I,
               double *stokes_Q, double *stokes_U, double *stokes_V,
               double *xmu )
{
  static char *argv[] = {"solveray"};
  int argc = sizeof(argv) / sizeof(argv[0]);
  return(solveray(argc, argv, nspect, nsize, lambda,
                  stokes_I, stokes_Q, stokes_U, stokes_V, xmu));
}

int solveray( int argc, char *argv[],
              int *nspect, int *nsize, double *lambda, double *stokes_I,
              double *stokes_Q, double *stokes_U, double *stokes_V,
              double *xmu )
{
  register int n, k, la;

  char    rayFileName[14], inputLine[MAX_LINE_SIZE], ascFilename[18];
  bool_t  result, exit_on_EOF, to_obs, initialize, crosscoupling,
          analyze_output, equilibria_only;
  int     Nsource, Nread, Nrequired, checkPoint, *wave_index = NULL;
  int     ier;
  double  muz, *S, *chi, *J;
  FILE   *fp_out, *fp_ray, *fp_stokes, *fp_out_asc;
  XDR     xdrs;
  ActiveSet *as;

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);

  // 18/06/19 epm: La siguiente llamada habilita el lanzamiento de floating
  // point exceptions (por defecto 1/0 lanza una excepcion pero 1.0/0.0 no).
  // SetFPEtraps() llama a feenableexcept(int excepts) -OS dependent- cuyo
  // manual en linux dice que 'excepts' puede ser:
  // FE_DIVBYZERO  division by zero
  // FE_INEXACT    inexact result: rounding was necessary to store the result
  // FE_INVALID    domain error (e.g. 0.0/0.0, sqrt(-1))
  // FE_OVERFLOW   operation too large to be representable
  // FE_UNDERFLOW  operation subnormal with a loss of precision
  // FE_ALL_EXCEPT bitwise OR of all supported floating-point exceptions
  // Examples:
  // 0.0/0.0 = nan       => exceptions raised: FE_INVALID
  // 1.0/0.0 = inf       => exceptions raised: FE_DIVBYZERO
  // 1.0/10.0 = 0.100000 => exceptions raised: FE_INEXACT
  // sqrt(-1) = -nan     => exceptions raised: FE_INVALID
  // DBL_MAX*2.0 = inf   => exceptions raised: FE_INEXACT & FE_OVERFLOW
  // nextafter(DBL_MIN/pow(2.0,52),0.0) = 0.0 => FE_INEXACT & FE_UNDERFLOW
  //
  // El caso es que al convertir 'solveray' en rutina, las excepciones
  // se habilitan tambien en el codigo Fortran y este no esta preparado para
  // ello. Es decir, el codigo de SIR admite que se produzcan operaciones
  // invalidas sin lanzar ningun aviso porque simplemente la solucion en
  // cuestion no convergera. Para evitar que este codigo Fortran aborte con
  // una floating point exception, debemos dejarlas inhabilitadas en C.
  //
  // Trap floating point exceptions on various machines.
  // SetFPEtraps();

  /* --- Read input data and initialize --             -------------- */

  readInput();
  spectrum.updateJ = FALSE;

  /* --- Read input data for atmosphere --             -------------- */

  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry);

  /* --- Read direction cosine for ray --              -------------- */

  // 04/04/20 epm: We don't need this file anymore since the cosine of
  // heliocentric angle is now an input argument.
  // if ((fp_ray = fopen(RAY_INPUT_FILE, "r")) == NULL) {
  //   sprintf(messageStr, "Unable to open inputfile %s", RAY_INPUT_FILE);
  //   Error(ERROR_LEVEL_2, argv[0], messageStr);
  // }
  // getLine(fp_ray, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  // Nread = sscanf(inputLine, "%lf", &muz);
  // checkNread(Nread, Nrequired=1, argv[0], checkPoint=1);
  muz = *xmu;

  if (muz <= 0.0  ||  muz > 1.0) {
    sprintf(messageStr,
            "Value of muz = %f does not lie in interval <0.0, 1.0]", muz);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }

  if ((atmos.Stokes && input.StokesMode == FIELD_FREE) ||
      input.backgr_pol) {

    /* --- Want formal solution to be polarized in these cases -- --- */

    input.StokesMode = FULL_STOKES;
  }

  /* --- redefine geometry for just this one ray --    -------------- */

  atmos.Nrays = geometry.Nrays = 1;
  geometry.muz[0] = muz;
  geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
  geometry.muy[0] = 0.0;
  geometry.wmu[0] = 1.0;
  if (atmos.Stokes) Bproject();

  input.startJ = OLD_J;

  readAtomicModels();
  readMolecularModels();
  SortLambda();

  getBoundary(&geometry);

  /* --- Open file with background opacities --        -------------- */

  // 04/04/20 epm: The keyword VMACRO_TRESH should be 0.0 (atmos.moving=TRUE).
  // Email from Han suggesting how to avoid problems if we save
  // "background.dat|ray" in memory rather than disk: "I think the best way
  // would be to assume that we have velocities (atmos.moving=TRUE) and/or
  // magnetic field in an inversion, so we can avoid readBRS()".
  if (atmos.moving || atmos.Stokes) {
    
    /* --- Case of angle-dependent background opacities -- ---------- */
      
    strcpy(input.background_File, input.background_ray_File);
    Background(analyze_output=FALSE, equilibria_only=FALSE);
  } else {
    Background(analyze_output=FALSE, equilibria_only=TRUE);

    if ((atmos.fd_background =
         open(input.background_File, O_RDONLY, 0)) == -1) {
      sprintf(messageStr, "Unable to open inputfile %s",
              input.background_File);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
    }
    readBRS();
  }
  convertScales(&atmos, &geometry);

  getProfiles();
  initSolution();
  initScatter();

  getCPU(1, TIME_POLL, "Total initialize");

  /* --- Solve radiative transfer equations --         -------------- */

  solveSpectrum(FALSE, FALSE);

  /* --- Write emergent spectrum to output file --     -------------- */

  // 04/04/20 epm: We don't need this file anymore.
  // sprintf(rayFileName, "spectrum_%4.2f", muz);
  // if ((fp_out = fopen(rayFileName, "w")) == NULL) {
  //   sprintf(messageStr, "Unable to open output file %s", rayFileName);
  //   Error(ERROR_LEVEL_2, argv[0], messageStr);
  // }
  // xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);
  //
  // result = xdr_double(&xdrs, &muz);
  // result = xdr_vector(&xdrs, (char *) spectrum.I[0], spectrum.Nspect,
  //                     sizeof(double), (xdrproc_t) xdr_double);

  /* --- Write ASCII table for special applications -- -------------- */

  // if (!input.xdr_endian) {
  //   sprintf(ascFilename, ASCII_SPECTRUM_FILE, muz);
  //   if ((fp_out_asc = fopen(ascFilename, "w" )) == NULL) {
  //     sprintf(messageStr, "Unable to open output file %s", ascFilename);
  //     Error(ERROR_LEVEL_2, argv[0], messageStr);
  //   }
  //   fprintf(fp_out_asc, "%d\n", spectrum.Nspect);
  //
  //   if (atmos.Stokes || input.backgr_pol) {
  //     for (la = 0;  la < spectrum.Nspect;  la++) {
  //       fprintf(fp_out_asc,
  //               "%19.13lf %17.10lg %17.10lg %17.10lg %17.10lg\n",
  //               spectrum.lambda[la],
  //               spectrum.I[0][la], spectrum.Stokes_Q[0][la],
  //               spectrum.Stokes_U[0][la], spectrum.Stokes_V[0][la]);
  //     }
  //   } else {
  //     for (la = 0;  la < spectrum.Nspect;  la++) {
  //       fprintf(fp_out_asc,
  //               "%19.13lf %17.10lg %17.10lg %17.10lg %17.10lg\n",
  //               spectrum.lambda[la],
  //               spectrum.I[0][la], 0.0, 0.0, 0.0);
  //     }
  //   }
  //   fclose(fp_out_asc);
  // }

  /* --- 10/06/19 epm: Save the spectrum in the array arguments ----- */

  *nspect = spectrum.Nspect;
  if (*nsize >= *nspect)
  {
    if (atmos.Stokes || input.backgr_pol)
    {
      for (la = 0; la < *nspect; la++)
      {
        lambda[la]   = spectrum.lambda[la];
        stokes_I[la] = spectrum.I[0][la];
        stokes_Q[la] = spectrum.Stokes_Q[0][la];
        stokes_U[la] = spectrum.Stokes_U[0][la];
        stokes_V[la] = spectrum.Stokes_V[0][la];
      }
    }
    else
    {
      for (la = 0; la < *nspect; la++)
      {
        lambda[la]   = spectrum.lambda[la];
        stokes_I[la] = spectrum.I[0][la];
        stokes_Q[la] = 0.0;
        stokes_U[la] = 0.0;
        stokes_V[la] = 0.0;
      }
    }
    ier = 0;
  }
  else if (*nsize == 0)   // we don't want the data
  {
    ier = 0;
  }
  else   // unsufficient memory for the arrays
  {
    ier = -1;
  }

  /* --- Read wavelength indices for which chi and S are to be
         written out for the specified direction --    -------------- */

  // 04/04/20 RH Manual: Nsource is the number of source functions and
  // opacity wavelengths that are stored with the solution. Typically,
  // this is set to zero.
  // Nread = fscanf(fp_ray, "%d", &Nsource);
  // checkNread(Nread, 1, argv[0], checkPoint=2);

  Nsource = 0;

  // if (Nsource > 0) {
  //   wave_index = (int *) malloc(Nsource * sizeof(int));
  //   Nread = 0;
  //   while (fscanf(fp_ray, "%d", &wave_index[Nread]) != EOF) Nread++;
  //   checkNread(Nread, Nsource, argv[0], checkPoint=3);
  //   fclose(fp_ray);
  //
  //   chi = (double *) malloc(atmos.Nspace * sizeof(double));
  //   if (atmos.Stokes)
  //     S = (double *) malloc(4 * atmos.Nspace * sizeof(double));
  //   else
  //     S = (double *) malloc(atmos.Nspace * sizeof(double));
  // }
  // result = xdr_int(&xdrs, &Nsource);

  /* --- Go through the list of wavelengths --         -------------- */

  // if (Nsource > 0  &&  input.limit_memory)
  //   J = (double *) malloc(atmos.Nspace * sizeof(double));
  //
  // for (n = 0;  n < Nsource;  n++) {
  //   if (wave_index[n] < 0  ||  wave_index[n] >= spectrum.Nspect) {
  //     sprintf(messageStr, "Illegal value of wave_index[n]: %4d\n"
  //             " Value has to be between 0 and %4d",
  //             wave_index[n], spectrum.Nspect);
  //     Error(ERROR_LEVEL_2, argv[0], messageStr);
  //     continue;
  //   }
  //   sprintf(messageStr, " Processing n = %4d, lambda = %9.3f [nm]\n\n",
  //           wave_index[n], spectrum.lambda[wave_index[n]]);
  //   Error(MESSAGE, NULL, messageStr);
  //
  //   as = &spectrum.as[wave_index[n]];
  //   alloc_as(wave_index[n], crosscoupling=FALSE);
  //   Opacity(wave_index[n], 0, to_obs=TRUE, initialize=TRUE);
  //   readBackground(wave_index[n], 0, to_obs=TRUE);
  //
  //   if (input.limit_memory) {
  //     readJlambda(wave_index[n], J);
  //   } else
  //     J = spectrum.J[wave_index[n]];
  //
  //   /* --- Add the continuum opacity and emissivity -- ----------- */
  //
  //   for (k = 0;  k < atmos.Nspace;  k++) {
  //     chi[k] = as->chi[k] + as->chi_c[k];
  //     S[k]   = (as->eta[k] + as->eta_c[k] + as->sca_c[k]*J[k]) / chi[k];
  //   }
  //   result = xdr_int(&xdrs, &wave_index[n]);
  //   result = xdr_vector(&xdrs, (char *) chi, atmos.Nspace,
  //                       sizeof(double), (xdrproc_t) xdr_double);
  //   result = xdr_vector(&xdrs, (char *) S, atmos.Nspace,
  //                       sizeof(double), (xdrproc_t) xdr_double);
  //
  //   free_as(wave_index[n], crosscoupling=FALSE);
  // }

  /* --- If magnetic fields are present or background is polarized -- */

  // if (atmos.Stokes || input.backgr_pol) {
  //   result = xdr_vector(&xdrs, (char *) spectrum.Stokes_Q[0],
  //                       spectrum.Nspect, sizeof(double),
  //                       (xdrproc_t) xdr_double);
  //   result = xdr_vector(&xdrs, (char *) spectrum.Stokes_U[0],
  //                       spectrum.Nspect, sizeof(double),
  //                       (xdrproc_t) xdr_double);
  //   result = xdr_vector(&xdrs, (char *) spectrum.Stokes_V[0],
  //                       spectrum.Nspect, sizeof(double),
  //                       (xdrproc_t) xdr_double);
  // }

  // if (Nsource > 0  &&  input.limit_memory)
  //   free(J);

  // xdr_destroy(&xdrs);
  // fclose(fp_out);
  printTotalCPU();

  // 09/09/19 epm: Just before finishing, close open writing descriptors.
  closefiles();

  // 02/02/21 epm: C memory leak detection library.
  // PrintTotalAllocatedMemory();
  // PrintMemoryReservedByCMemDbgLibrary();

  return(ier);
}
/* ------- end ---------------------------- solveray.c -------------- */
