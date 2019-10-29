/* ------- file: -------------------------- conversion.c ------------

       19/02/19 brc : Extract of "rhf1d.c" just to convert models.

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 24 16:40:14 2011 --

       First change:  Tue Jun 18 2019 by Esperanza Paez (epm@iac.es)
                      Adaptation to the DeSIRe project.
                      Interoperability C-Fortran.
                      Conversion of main() as a function for Fortran.

       --------------------------                      ----------RH-- */


#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "inputs.h"


/* --- Function prototypes --                          -------------- */

void closefiles(void);

int conversion(int argc, char *argv[],
               int *nspace, int *nsize, double *tau, double *mass,
               double *z, double *T, double *ne, double *vel, double *vturb);


/* --- Global variables --                             -------------- */

// 31/07/19 epm: Global variables definition is now in "rh_glob.c".
extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;


/* ------- begin -------------------------- conversion.c ------------ */

// Function to be called from Fortran.
int conversion_( int *nspace, int *nsize, double *tau, double *mass,
                 double *z, double *T, double *ne, double *vel, double *vturb )
{
  static char *argv[] = {"conversion", NULL};
  return(conversion(1, argv, nspace, nsize, tau, mass, z, T, ne, vel, vturb));
}

int conversion( int argc, char *argv[],
                int *nspace, int *nsize, double *tau, double *mass,
                double *z, double *T, double *ne, double *vel, double *vturb )
{
  bool_t write_analyze_output, equilibria_only;
  FILE  *f;
  int    i;
  int    ier;

  Atom *atom;
  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */

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
  // El caso es que al convertir 'conversion' en rutina, las excepciones
  // se habilitan tambien en el codigo Fortran y este no esta preparado para
  // ello. Es decir, el codigo de DeSIRe admite que se produzcan operaciones
  // invalidas sin lanzar ningun aviso porque simplemente la solucion en
  // cuestion no convergera. Para evitar que el codigo Fortran aborte con una
  // floating point exception, debemos dejarlas inhabilitadas en C.
  //
  // Trap floating point exceptions on various machines.
  // SetFPEtraps();

  readInput();
  spectrum.updateJ = TRUE;

  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry);
  if (atmos.Stokes) Bproject();

  readAtomicModels();
  readMolecularModels();
  SortLambda();

  getBoundary(&geometry);

  Background(write_analyze_output=TRUE, equilibria_only=FALSE);
  convertScales(&atmos, &geometry);

  getCPU(1, TIME_POLL, "Total Initialize");

  // 19/06/19 epm: We don't need this anymore.
  //
  // printf("Writing to output file = %s\n", "atmos.txt.out");
  // if ((f = fopen("atmos.txt.out", "w")) == NULL)
  // {
  //   printf("Error opening file %s!!!\n", "atmos.txt.out");
  //   exit(1);
  // }
  // fprintf(f,"%d\n", atmos.Nspace);
  // for (i = 0; i < atmos.Nspace; i++)
  // {
  //   fprintf(f, "%e %e %e %e %e %e %e\n",geometry.tau_ref[i],geometry.cmass[i], geometry.height[i],atmos.T[i],atmos.ne[i], geometry.vel[i], atmos.vturb[i]);
  // }
  // fclose(f);

  /* --- 19/06/19 epm: Save the model in the array arguments. ------- */

  *nspace = atmos.Nspace;
  if (*nsize >= *nspace)
  {
    for (i = 0; i < *nspace; i++)
    {
      tau[i]   = geometry.tau_ref[i];
      mass[i]  = geometry.cmass[i];
      z[i]     = geometry.height[i];
      T[i]     = atmos.T[i];
      ne[i]    = atmos.ne[i];
      vel[i]   = geometry.vel[i];
      vturb[i] = atmos.vturb[i];
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

  printTotalCPU();

  // 09/09/19 epm: Just before finishing, close open writing descriptors.
  closefiles();

  return(ier);
}
/* ------- end ---------------------------- conversion.c ------------ */
