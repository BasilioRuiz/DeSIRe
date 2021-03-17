/* ------- file: -------------------------- barklem.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Nov 5, 2014, 15:24 --

       --------------------------                      ----------RH-- */

/* --- Routines to deal with Collisional broadening as formulated by
       Barklem et al.

  See: - Anstee & O'Mara 1995, MNRAS 276, 859-866 (s-p, p-s)
       - Barklem & O'Mara 1997, MNRAS 290, 102 (p-d, d-p)
       - Barklem, O'Mara & Ross 1998, MNRAS 296, 1057-1060 (d-f, f-d)
       - Barklem, O'Mara 1998, MNRAS 300, 863-871

  Modifications:

       - 05/11/16 JdlCR:
         The tables from Barklem do not work with ionized species
         like Ca II or Mg II. Added the possibility to provide the
         cross-sections in the atom file.

       - 21/06/19 epm:
         The directory holding Barklem files is now an input from
         "keyword.input".

       - 31/07/19 epm:
         We supply collisional broadening for some lines from SIR.

       --                                              -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"
#include "desire.h"


// 21/06/19 epm: The directory is now an input from "keyword.input".
//#define BARKLEM_SP_DATA   "../src/rh/Atoms/Barklem_spdata.dat"
#define BARKLEM_SP_DATA     "Barklem_spdata.dat"
#define BARKLEM_SP_NS       21
#define BARKLEM_SP_NP       18
#define BARKLEM_SP_NEFF1    1.0
#define BARKLEM_SP_NEFF2    1.3

//#define BARKLEM_PD_DATA   "../src/rh/Atoms/Barklem_pddata.dat"
#define BARKLEM_PD_DATA     "Barklem_pddata.dat"
#define BARKLEM_PD_NP       18
#define BARKLEM_PD_ND       18
#define BARKLEM_PD_NEFF1    1.3
#define BARKLEM_PD_NEFF2    2.3

//#define BARKLEM_DF_DATA   "../src/rh/Atoms/Barklem_dfdata.dat"
#define BARKLEM_DF_DATA     "Barklem_dfdata.dat"
#define BARKLEM_DF_ND       18
#define BARKLEM_DF_NF       18
#define BARKLEM_DF_NEFF1    2.3
#define BARKLEM_DF_NEFF2    3.3

#define BARKLEM_DELTA_NEFF  0.1

#define BOHR_RADIUS_SQUARED 2.8002e-17


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];

// 31/07/19 epm: New struct to load information about SIR lines.
extern SIRBarklem sirbarklem;


/* ------- begin -------------------------- readBarklemTable.c ------ */

bool_t readBarklemTable(enum Barklemtype type, Barklemstruct *bs)
{
  register int n, i, j;
  const char routineName[] = "readBarklemTable";

  char    filename[MAX_LINE_SIZE], inputLine[MAX_LINE_SIZE], *charptr;
  int     nread;
  double  neff1_0, neff2_0;
  FILE   *fp_Barklem;

  switch (type) {
  case SP:
  //strcpy(filename, BARKLEM_SP_DATA);
    sprintf(filename, "%s/%s", input.barklem_dir, BARKLEM_SP_DATA);
    bs->N1 = BARKLEM_SP_NS;
    bs->N2 = BARKLEM_SP_NP;

    neff1_0 = BARKLEM_SP_NEFF1;
    neff2_0 = BARKLEM_SP_NEFF2;
    break;

  case PD:
  //strcpy(filename, BARKLEM_PD_DATA);
    sprintf(filename, "%s/%s", input.barklem_dir, BARKLEM_PD_DATA);
    bs->N1 = BARKLEM_PD_NP;
    bs->N2 = BARKLEM_PD_ND;

    neff1_0 = BARKLEM_PD_NEFF1;
    neff2_0 = BARKLEM_PD_NEFF2;
    break;

  case DF:
  //strcpy(filename, BARKLEM_DF_DATA);
    sprintf(filename, "%s/%s", input.barklem_dir, BARKLEM_DF_DATA);
    bs->N1 = BARKLEM_DF_ND;
    bs->N2 = BARKLEM_DF_NF;

    neff1_0 = BARKLEM_DF_NEFF1;
    neff2_0 = BARKLEM_DF_NEFF2;
    break;
  }

  if ((fp_Barklem = fopen(filename, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", filename);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }

  bs->neff1 = (double *) malloc(bs->N1 * sizeof(double));
  for (n = 0;  n < bs->N1;  n++)
    bs->neff1[n] = neff1_0 + n * BARKLEM_DELTA_NEFF;

  bs->neff2 = (double *) malloc(bs->N2 * sizeof(double));
  for (n = 0;  n < bs->N2;  n++)
    bs->neff2[n] = neff2_0 + n * BARKLEM_DELTA_NEFF;

  bs->cross = matrix_double(bs->N1, bs->N2);
  bs->alpha = matrix_double(bs->N1, bs->N2);

  for (n = 0;  n < 3;  n++)
    charptr = fgets(inputLine, MAX_LINE_SIZE, fp_Barklem);

  for (i = 0;  i < bs->N1;  i++)
    for (j = 0;  j < bs->N2;  j++) {
      nread = fscanf(fp_Barklem, "%lf", &bs->cross[i][j]);
  }
  for (n = 0;  n < 2;  n++)
    charptr = fgets(inputLine, MAX_LINE_SIZE, fp_Barklem);

  for (i = 0;  i < bs->N1;  i++)
    for (j = 0;  j < bs->N2;  j++) {
      nread = fscanf(fp_Barklem, "%lf", &bs->alpha[i][j]);
  }

  fclose(fp_Barklem);
  return TRUE;
}
/* ------- end ---------------------------- readBarklemTable.c------- */

/* ------- begin -------------------------- getBarklemcross.c ------- */

// 02/10/19 epm: getBarklemcross() receives the wavelength (last argument).

bool_t getBarklemcross(Barklemstruct *bs, RLK_Line *rlk, double lambda_air)
{
  const char routineName[] = "getBarklemcross";

  int index, k;
  double diff;
  double Z, neff1, neff2, findex1, findex2, reducedmass, meanvelocity,
         crossmean, E_Rydberg, deltaEi, deltaEj;
  Element *element;

  element = &atmos.elements[rlk->pt_index - 1];

  /* --- 30/09/19 epm:
   *     We supply collisional broadening for some lines from SIR --- */

  for (k = 0; k < sirbarklem.nlines; k++)
  {
    diff = fabs(sirbarklem.wave[k] * 100 - lambda_air * 1000);
    if (diff < 2.0)   // 20 miliAngstrom
    {
      if (sirbarklem.alpha[k] > 1e-20 && sirbarklem.sigma[k] > 1e-20)
      {
        rlk->alpha = sirbarklem.alpha[k];
        rlk->cross = sirbarklem.sigma[k]/BOHR_RADIUS_SQUARED;
        break;
      }
      else
      {
        sprintf(messageStr,
                "%s %.3lf: Taking UNSOLD collisional broadening",
                element->ID, sirbarklem.wave[k]);
        Error(WARNING, routineName, messageStr);
        return(FALSE);
      }
    }
  }
  if (k == sirbarklem.nlines)   // if the line is not a SIR line
  {
     return(FALSE);   // use UNSOLD
  }

  // 30/09/19 epm: We don't want to interpole Barklem table from RH anymore.
  //
  // /* --- ABO tabulations are valid only for neutral atoms -- -- */
  //
  // if (rlk->stage > 0)
  //   return FALSE;
  //
  // if ((deltaEi = element->ionpot[rlk->stage] - rlk->Ei) <= 0.0)
  //   return FALSE;
  // if ((deltaEj = element->ionpot[rlk->stage] - rlk->Ej) <= 0.0)
  //   return FALSE;
  //
  // Z = (double) (rlk->stage + 1);
  // E_Rydberg = E_RYDBERG / (1.0 + M_ELECTRON / (element->weight * AMU));
  // neff1 = Z * sqrt(E_Rydberg / deltaEi);
  // neff2 = Z * sqrt(E_Rydberg / deltaEj);
  //
  // if (rlk->Li > rlk->Lj) SWAPDOUBLE(neff1, neff2);
  //
  // if (neff1 < bs->neff1[0] || neff1 > bs->neff1[bs->N1-1])
  //   return FALSE;
  // Locate(bs->N1, bs->neff1, neff1, &index);
  // findex1 =
  //   (double) index + (neff1 - bs->neff1[index]) / BARKLEM_DELTA_NEFF;
  //
  // if (neff2 < bs->neff2[0] || neff2 > bs->neff2[bs->N2-1])
  //   return FALSE;
  // Locate(bs->N2, bs->neff2, neff2, &index);
  // findex2 =
  //   (double) index + (neff2 - bs->neff2[index]) / BARKLEM_DELTA_NEFF;
  //
  // /* --- Find interpolation in table --             -------------- */
  //
  // rlk->cross = cubeconvol(bs->N2, bs->N1,
  //                         bs->cross[0], findex2, findex1);
  // rlk->alpha = cubeconvol(bs->N2, bs->N1,
  //                         bs->alpha[0], findex2, findex1);

  reducedmass  = AMU / (1.0/atmos.H->weight + 1.0/element->weight);
  meanvelocity = sqrt(8.0 * KBOLTZMANN / (PI * reducedmass));
  crossmean    = SQ(RBOHR) * pow(meanvelocity / 1.0E4, -rlk->alpha);

  rlk->cross *= 2.0 * pow(4.0/PI, rlk->alpha/2.0) *
    exp(gammln((4.0 - rlk->alpha)/2.0)) * meanvelocity * crossmean;

  rlk->vdwaals = BARKLEM;
  return TRUE;
}
/* ------- end ---------------------------- getBarklemcross.c ------- */

/* ------- begin -------------------------- getBarklemactivecross.c - */

bool_t getBarklemactivecross(AtomicLine *line)
{
  const char routineName[] = "getBarklemactivecross";

  bool_t determined = TRUE, useBarklem = FALSE;
  int index, Ll, Lu, nq, i, j, ic, k;
  double Sl, Su, Jl, Ju;
  double Z, neff1, neff2, findex1, findex2, reducedmass, meanvelocity,
         crossmean, E_Rydberg, deltaEi, deltaEj;
  Atom *atom;
  Barklemstruct bs;

  atom = line->atom;
  j = line->j;
  i = line->i;

  /* --- 31/07/19 epm:
   *     We supply collisional broadening for some lines from SIR --- */

  for (k = 0; k < sirbarklem.nlines; k++)
  {
    // 31/07/19 epm: In NLTE use SIR values (if not null) written
    // in SIR atomic file when Barklem is on in the RH atomic file.

    if (sirbarklem.anumber[k] == (atom->periodic_table+1) &&
        sirbarklem.stage[k]   == (atom->stage[i]+1)       &&
        sirbarklem.low[k]     == i                        &&
        sirbarklem.up[k]      == j                         )
    {
      if (sirbarklem.alpha[k] > 1e-20 && sirbarklem.sigma[k] > 1e-20)
      {
        line->cvdWaals[1] = sirbarklem.alpha[k];
        line->cvdWaals[0] = sirbarklem.sigma[k]/BOHR_RADIUS_SQUARED;
        break;
      }
      else
      {
        sprintf(messageStr,
                "%s_%d (low %d, up %d): Taking UNSOLD collisional broadening",
                atom->ID, atom->stage[i]+1, i, j);
        Error(WARNING, routineName, messageStr);
        return(FALSE);
      }
    }
  }

  /* ---
     JdlCR: Interpolate the tables only if sigma is smaller than 20
     otherwise assume that we are already giving Barklem cross-sections
     --- */

  if (line->cvdWaals[0] < 20.0)
  {

    /* --- ABO tabulations are only valid for neutral atoms  -- ------- */

    if (atom->stage[i] > 0)
      return FALSE;

    /* --- Get the quantum numbers for orbital angular momentum -- ---- */

    determined &= determinate(atom->label[i], atom->g[i],
                              &nq, &Sl, &Ll, &Jl);
    determined &= determinate(atom->label[j], atom->g[j],
                              &nq, &Su, &Lu, &Ju);

    /* --- See if one of the Barklem cases applies --    -------------- */

    if (determined) {
      if ((Ll == S_ORBIT && Lu == P_ORBIT) ||
          (Ll == P_ORBIT && Lu == S_ORBIT)) {
          useBarklem = readBarklemTable(SP, &bs);
      } else if ((Ll == P_ORBIT && Lu == D_ORBIT) ||
                 (Ll == D_ORBIT && Lu == P_ORBIT)) {
          useBarklem = readBarklemTable(PD, &bs);
      } else if ((Ll == D_ORBIT && Lu == F_ORBIT) ||
                 (Ll == F_ORBIT && Lu == D_ORBIT)) {
          useBarklem = readBarklemTable(DF, &bs);
      }
    }
    if (!determined || !useBarklem) return FALSE;

    /* --- Determine the index of the appropriate continuum level -- -- */

    Z = atom->stage[j] + 1;
    for (ic = j + 1;  atom->stage[ic] < atom->stage[j]+1;  ic++);

    deltaEi   = atom->E[ic] - atom->E[i];
    deltaEj   = atom->E[ic] - atom->E[j];
    E_Rydberg = E_RYDBERG / (1.0 + M_ELECTRON / (atom->weight * AMU));

    neff1 = Z * sqrt(E_Rydberg / deltaEi);
    neff2 = Z * sqrt(E_Rydberg / deltaEj);

    if (Ll > Lu) SWAPDOUBLE(neff1, neff2);

    /* --- Interpolate according to effective principal quantum number  */

    if (neff1 < bs.neff1[0] || neff1 > bs.neff1[bs.N1-1])
      return FALSE;
    Locate(bs.N1, bs.neff1, neff1, &index);
    findex1 =
      (double) index + (neff1 - bs.neff1[index]) / BARKLEM_DELTA_NEFF;

    if (neff2 < bs.neff2[0] || neff2 > bs.neff2[bs.N2-1])
      return FALSE;
    Locate(bs.N2, bs.neff2, neff2, &index);
    findex2 =
      (double) index + (neff2 - bs.neff2[index]) / BARKLEM_DELTA_NEFF;

    /* --- Find interpolation in table --                -------------- */

    line->cvdWaals[0] = cubeconvol(bs.N2, bs.N1,
                                   bs.cross[0], findex2, findex1);
    line->cvdWaals[1] = cubeconvol(bs.N2, bs.N1,
                                   bs.alpha[0], findex2, findex1);
  }

  reducedmass  = AMU / (1.0/atmos.atoms[0].weight + 1.0/atom->weight);
  meanvelocity = sqrt(8.0 * KBOLTZMANN / (PI * reducedmass));
  crossmean    = SQ(RBOHR) * pow(meanvelocity / 1.0E4, -line->cvdWaals[1]);

  line->cvdWaals[0] *= 2.0 * pow(4.0/PI, line->cvdWaals[1]/2.0) *
    exp(gammln((4.0 - line->cvdWaals[1])/2.0)) * meanvelocity * crossmean;

  /* --- Use UNSOLD for the contribution of Helium atoms -- ---------- */

  line->cvdWaals[2] = 1.0;
  line->cvdWaals[3] = 0.0;

  return TRUE;
}
/* ------- end ---------------------------- getBarklemactivecross.c -- */
