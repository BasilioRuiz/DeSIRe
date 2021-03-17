/* ------- file: -------------------------- rhf1d.c -----------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 11 12:03:51 2021

       Updates (epm): Adaptation to the DeSIRe project.
                      Interoperability C-Fortran.
                      Conversion of main() as a function for Fortran.

       --------------------------                      ----------RH-- */

/* --- Main routine of 1D plane-parallel radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

       Formal solution is performed with Feautrier difference scheme
       in static atmospheres, and with piecewise quadratic integration
       in moving atmospheres.

       --                                              -------------- */

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


/* --- Function prototypes --                          -------------- */

void closefiles(void);

int rhf1d(int argc, char *argv[],
          int *nspace, int *nlevel, int *nsize, double *coefs);


/* --- Global variables --                             -------------- */

// 31/07/19 epm: Global variables definition is now in "rh_glob.c".
extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;

// 02/07/19 epm: Boolean variables to control new RH executions.
extern bool_t new_background;   // for background.c :: Background()
extern bool_t new_hminus_ff;    // for hydrogen.c :: Hminus_ff()
extern bool_t new_h2minus_ff;   // for hydrogen.c :: H2minus_ff()
extern bool_t new_h2plus_ff;    // for hydrogen.c :: H2plus_ff()
extern bool_t new_passive_bb;   // for metal.c :: passive_bb()


/* ------- begin -------------------------- rhf1d() ----------------- */

// 23/05/19 epm: Function to be called from Fortran.
int rhf1d_( int *nspace, int *nlevel, int *nsize, double *coefs )
{
  static char *argv[] = {"rhf1d"};
  int argc = sizeof(argv) / sizeof(argv[0]);
  return(rhf1d(argc, argv, nspace, nlevel, nsize, coefs));
}

int rhf1d( int argc, char *argv[],
           int *nspace, int *nlevel, int *nsize, double *coefs )
{
  bool_t write_analyze_output, equilibria_only;
  double deltaJ;
  int    niter, nact;
  int    index, i, k;
  int    ier;

  Atom *atom;
  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */

  // 02/07/19 epm: New RH execution.
  new_background = TRUE;
  new_hminus_ff  = TRUE;
  new_h2minus_ff = TRUE;
  new_h2plus_ff  = TRUE;
  new_passive_bb = TRUE;

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
  // El caso es que al convertir 'rhf1d' en rutina, las excepciones
  // se habilitan tambien en el codigo Fortran y este no esta preparado para
  // ello. Es decir, el codigo de SIR admite que se produzcan operaciones
  // invalidas sin lanzar ningun aviso porque simplemente la solucion en
  // cuestion no convergera. Para evitar que este codigo Fortran aborte con
  // una floating point exception, debemos dejarlas inhabilitadas en C.
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

  getProfiles();
  initSolution();
  initScatter();

  getCPU(1, TIME_POLL, "Total Initialize");

  /* --- Solve radiative transfer for active ingredients -- --------- */

  Iterate(input.NmaxIter, input.iterLimit);

  /* --- If appropriate final solution(s) should be polarized -- ---- */

  adjustStokesMode();
  niter = 0;
  while (niter < input.NmaxScatter) {
    deltaJ = solveSpectrum(FALSE, FALSE);
    if (!input.backgr_pol && deltaJ <= input.iterLimit) break;
    niter++;
  }

  /* --- Write output files --                         -------------- */

  // 04/04/20 epm: Where possible, avoid writing on the hard disk.

  if (atmos.hydrostatic) {
    geometry.scale = COLUMN_MASS;
    convertScales(&atmos, &geometry);
  }
  getCPU(1, TIME_START, NULL);

  // writeInput();              // we don't want "input.out"
  writeAtmos(&atmos);           // avoided setting ATMOS_OUTPUT = none
  writeGeometry(&geometry);     // avoided setting GEOMETRY_OUTPUT = none
  // writeFlux(FLUX_DOT_OUT);   // we don't want "flux.out"

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    // writeAtom(atom);         // we don't want "atom.%s.out"
    writePopulations(atom);     // we save populations in memory
    writeRadRate(atom);         // avoided setting RADRATE_OUTPUT = none
    writeCollisionRate(atom);   // avoided setting COLLRATE_OUTPUT = none
    writeDamping(atom);         // avoided setting DAMPING_OUTPUT = none
  }
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    writeMolPops(molecule);     // at the moment we don't use molecules
  }

  writeSpectrum(&spectrum);     // avoided setting SPECTRUM_OUTPUT = none
  writeOpacity();               // avoided setting OPACITY_OUTPUT = none

  getCPU(1, TIME_POLL, "Write output");

  /* --- 06/06/19 epm: Save the populations in the array argument --- */

  *nspace = atmos.Nspace;
  *nlevel = atom->Nlevel;
  if (*nsize >= (*nspace) * (*nlevel))
  {
    for (index = 0, k = 0;  k < *nspace;  k++)
    {
      for (i = 0;  i < *nlevel;  i++)
      {
        coefs[index++] = atom->n[i][k] / atom->nstar[i][k];
      }
    }
    ier = 0;
  }
  else if (*nsize == 0)   // we don't want the data
  {
    ier = 0;
  }
  else   // unsufficient memory for 'coefs'
  {
    ier = -1;
  }

  printTotalCPU();

  // 09/09/19 epm: Just before finishing, close open writing descriptors.
  closefiles();
  // 04/04/20 epm: Print a new line at the end of rhf1d.
  Error(MESSAGE, NULL, "\n");

  return(ier);
}
/* ------- end ---------------------------- rhf1d() ----------------- */

/* ------- begin -------------------------- load_hpop_rh_() --------- */

// 11/11/20 epm: Function to be called from Fortran.
void load_hpop_rh_( int *ntau, float nh[*ntau][6] )
{
  int i, j, k;
  //FILE *pepe = fopen("pepe","w");

  // Pass the H populations to SIR (Fortran). Be aware that Fortran and
  // C have different conventions for storing multi-dimensional arrays.
  // Fortran stores multi-dimensional arrays in column-major order while
  // C arrays are stored in row-major order.

  // Unlike SIR, RH's atmosphere model runs from upper layers to lower
  // ones, so we need to turn around the array.
  for (k = *ntau-1, i = 0;  k >= 0;  k--, i++)
  {
    for (j = 0;  j < 6;  j++)
    {
       nh[i][j] = atmos.H->n[j][k]*1e-6;  // convert to CGS system
    }
    //fprintf(pepe, "%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
    //        nh[i][0],nh[i][1],nh[i][2],nh[i][3],nh[i][4],nh[i][5]);
  }
  //fclose(pepe);
}
/* ------- end ---------------------------- load_hpop_rh_() --------- */
