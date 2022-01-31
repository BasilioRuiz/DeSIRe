//
// File          _______ : rh_sir.c
// Description   _______ : Functions to pass information between SIR and RH.
// Project       _______ : DeSIRe
// Creation date _______ : 31/07/19
// Author        _______ : epm@iac.es
//

/* 
 http://indico.ictp.it/event/a13229/session/2/contribution/11/material/0/0.pdf
 *
 * CALLING C FROM F77
 * Need to make C function look like Fortran 77.
 * Append underscore (except on AIX, HP-UX).
 * Call by reference conventions.
 * Best only used for "subroutine" constructs (cf. MPI) as passing return
 * value of functions varies a lot:
 *    void add_abs_(int *v1,int *v2,int *res){*res = abs(*v1)+abs(*v2);}
 * Strings are tricky (no terminal 0, length added).
 * Arrays are always passed as "flat" 1d arrays by providing a pointer to the
 * first array element.
 *
 * CALLING F77 FROM C
 * Difficult for anything but Fortran 77 style calls since Fortran 90+
 * features need extra info:
 *    - shaped arrays, optional parameters, modules.
 * Arrays need to be "flat". C-style multi-dimensional arrays are lists of
 * pointers to individual pieces of storage, which may not be consecutive
 * => use 1d and compute position.
 */


#include <stdio.h>
#include <strings.h>

#include "desire.h"

extern SIRabun     sirabun;
extern SIRatmos    siratmos;
extern SIRBarklem  sirbarklem;
extern SIRcpl      sircpl[2];
extern SIRflags    sirflags;
extern SIRkeywords sirkeywords;
extern SIRKurucz   sirkurucz;
extern SIRwave     sirwave;
extern RHdc       *rhdc;
extern RHhpop      rhhpop;


//____________________________________________________________________________


/*------------------- PASSING INFORMATION FROM SIR TO RH -------------------*/

//____________________________________________________________________________
//
//  Method: sirabundances_()
//
/** Pass abundance values from SIR to RH.
 *  @param  nelem number of elements.
 *  @param  abun array with the abundances.
 */
//____________________________________________________________________________

void sirabundances_( int *nelem, float *abun )
{
   sirabun.nelem = *nelem;
   sirabun.abun  = abun;
}


//____________________________________________________________________________
//
//  Method: siratmos_()
//
/** Pass atmospheric models from SIR to RH.
 *  @param  idascii array with the ASCII codes for the ID model string.
 *  @param  nidascii size of 'idascii'.
 *  @param  scale scale identifier ('M', 'T' or 'H').
 *  @param  logg gravity decimal logarithm.
 *  @param  ntau number of optical depths.
 *  @param  dscale array with depth scales.
 *  @param  T array with the temperatures.
 *  @param  elecdens array with the electronic densities (cm-3).
 *  @param  vz array with the line of sight velocity (km/s).
 *  @param  vmic array with the microturbulence velocity (km/s).
 *  @param  nB number of elements in B.
 *  @param  B array with the magnetic field strength (Teslas).
 *  @param  gamma array with magnetic field inclination (radians).
 *  @param  phi array with magnetic field azimuth (radians).
 */
//____________________________________________________________________________

void siratmos_( int *idascii, int *nidascii, char *scale,
                double *logg, int *ntau, double *dscale,
                float *T, float *elecdens, float *vz, float *vmic,
                int *nB, float *B, float *gamma, float *phi )
{
   int i, n;

   n = sizeof(siratmos.id);
   for (i = 0; i < *nidascii && i < (n-1); i++)
   {
      siratmos.id[i] = (char) idascii[i];
   }
   siratmos.id[i]    = '\0';
   siratmos.scale    = *scale;
   siratmos.logg     = *logg;
   siratmos.ntau     = *ntau;
   siratmos.dscale   = dscale;
   siratmos.T        = T;
   siratmos.elecdens = elecdens;
   siratmos.vz       = vz;
   siratmos.vmic     = vmic;
   siratmos.nB       = *nB;
   siratmos.B        = B;
   siratmos.gamma    = gamma;
   siratmos.phi      = phi;
}


//____________________________________________________________________________
//
//  Method: siratmosnh_()
//
/** Pass hydrogen populations from SIR to RH.
 *  @param  nh1 population of the fundamental level.
 *  @param  nh2 population of the first excited level.
 *  @param  nh3 population of the second excited level.
 *  @param  nh4 population of the third excited level.
 *  @param  nh5 population of the fourth excited level.
 *  @param  nh6 number of protons per cm3.
 */
//____________________________________________________________________________

void siratmosnh_( float *nh1, float *nh2, float *nh3,
                  float *nh4, float *nh5, float *nh6 )
{
   siratmos.nh1  = nh1;
   siratmos.nh2  = nh2;
   siratmos.nh3  = nh3;
   siratmos.nh4  = nh4;
   siratmos.nh5  = nh5;
   siratmos.nh6  = nh6;
}


//____________________________________________________________________________
//
//  Method: sirbarklem_()
//
/** Pass Barklem data from SIR to RH.
 *  @param  nlines number of lines.
 *  @param  anumber array with the atomic numbers.
 *  @param  stage array with the ionization stages.
 *  @param  low array with the low levels.
 *  @param  up array with the up levels.
 *  @param  alpha array with the alpha values.
 *  @param  sigma array with the sigma values.
 *  @param  wave array with the air wavelength values.
 */
//____________________________________________________________________________

void sirbarklem_( int *nlines, int *anumber, int *stage, int *low, int *up,
                  float *alpha, float *sigma, float *wave )
{
   sirbarklem.nlines  = *nlines;
   sirbarklem.anumber = anumber;
   sirbarklem.stage   = stage;
   sirbarklem.low     = low;
   sirbarklem.up      = up;
   sirbarklem.alpha   = alpha;
   sirbarklem.sigma   = sigma;
   sirbarklem.wave    = wave;
}


//____________________________________________________________________________
//
//  Method: sircpl_()
//
/** Pass coupling data (quantum numbers) of the lower and upper level from
 *  SIR to RH. The quantum number depends on the coupling scheme:
 *  - qn1 = J if LS | J  if JK | J  if JJ
 *  - qn2 = L if LS | L1 if JK | j1 if JJ
 *  - qn3 = S if LS | S1 if JK | l1 if JJ
 *  - qn4 = 0 if LS | J1 if JK | j2 if JJ
 *  - qn5 = 0 if LS | l  if JK | l2 if JJ
 *  - qn6 = 0 if LS | K  if JK | 0  if JJ
 *  @param  nlines number of LTE lines.
 *  @param  cpl_low_up array with the coupling scheme for each line.
 *  @param  qn1_low/up array with the quantum number 1 for each line.
 *  @param  qn2_low/up array with the quantum number 2 for each line.
 *  @param  qn3_low/up array with the quantum number 3 for each line.
 *  @param  qn4_low/up array with the quantum number 4 for each line.
 *  @param  qn5_low/up array with the quantum number 5 for each line.
 *  @param  qn6_low/up array with the quantum number 6 for each line.
 */
//____________________________________________________________________________

void sircpl_( int *nlines, int *cpl_low, int *cpl_up,
              float *qn1_low, float *qn1_up, float *qn2_low, float *qn2_up,
              float *qn3_low, float *qn3_up, float *qn4_low, float *qn4_up,
              float *qn5_low, float *qn5_up, float *qn6_low, float *qn6_up )
{
   sircpl[0].nlines = *nlines;
   sircpl[0].cpl    = cpl_low;
   sircpl[0].qn1    = qn1_low;
   sircpl[0].qn2    = qn2_low;
   sircpl[0].qn3    = qn3_low;
   sircpl[0].qn4    = qn4_low;
   sircpl[0].qn5    = qn5_low;
   sircpl[0].qn6    = qn6_low;
   sircpl[1].nlines = *nlines;
   sircpl[1].cpl    = cpl_up;
   sircpl[1].qn1    = qn1_up;
   sircpl[1].qn2    = qn2_up;
   sircpl[1].qn3    = qn3_up;
   sircpl[1].qn4    = qn4_up;
   sircpl[1].qn5    = qn5_up;
   sircpl[1].qn6    = qn6_up;
}


//____________________________________________________________________________
//
//  Method: sirflags_()
//
/** Pass some command line flags from SIR to RH.
 *  @param  flagv verbose flag.
 */
//____________________________________________________________________________

void sirflags_( int *flagv )
{
   sirflags.flagv = *flagv;
}


//____________________________________________________________________________
//
//  Method: sirkeywords_()
//
/** Pass some keywords from SIR to RH.
 *  @param  nrays number of rays.
 *  @param  startj starting solution for J.
 */
//____________________________________________________________________________

void sirkeywords_( int *nrays, int *startj )
{
   sirkeywords.nrays  = *nrays;
   sirkeywords.startj = *startj;
}


//____________________________________________________________________________
//
//  Method: sirkurucz_()
//
/** Pass Kurucz data from SIR to RH.
 *  @param  nlines number of LTE lines.
 *  @param  lambda array with the wavelength in air (nm).
 *  @param  gf array with the log gf.
 *  @param  elemcode array with the element code.
 *  @param  ei array with the first energy level (cm-1).
 *  @param  ej array with the second energy level (cm-1).
 *  @param  grad array with the log of radiative damping constant.
 *  @param  gstark array with the log of stark damping const/electron number.
 *  @param  gvdwaals array with the log of van der Waals damping.
 * 
 *  @note The parameters contain a list in the format described on the Kurucz
 *        web page http://kurucz.harvard.edu/linelists.html
 */
//____________________________________________________________________________

void sirkurucz_( int *nlines,
                 double *lambda, double *gf, double *elemcode,
                 double *ei, double *ej,
                 double *grad, double *gstark, double *gvdwaals )
{
   sirkurucz.nlines   = *nlines;
   sirkurucz.lambda   = lambda;
   sirkurucz.gf       = gf;
   sirkurucz.elemcode = elemcode;
   sirkurucz.Ei       = ei;
   sirkurucz.Ej       = ej;
   sirkurucz.Grad     = grad;
   sirkurucz.GStark   = gstark;
   sirkurucz.GvdWaals = gvdwaals;
}


//____________________________________________________________________________
//
//  Method: sirwave_()
//
/** Pass the wavetable from SIR to RH.
 *  @param  nwave number of wavelengths.
 *  @param  wavetable array with the wavetable.
 */
//____________________________________________________________________________

void sirwave_( int *nwave, double *wavetable )
{
   sirwave.nwave     = *nwave;
   sirwave.wavetable = wavetable;
}


//____________________________________________________________________________


/*------------------- PASSING INFORMATION FROM RH TO SIR -------------------*/

//____________________________________________________________________________
//
//  Method: rhdeparture_()
//
/** Pass the departure coefficients for the atom argument from RH to SIR.
 *  @param  atomi array with the ASCII codes for the ID of the atom.
 *  @param  ilow lower level index (0 for the fundamental level).
 *  @param  iup upper level index.
 *  @param  nsize dimension of the array departure.
 *  @param  departure array with the coefficients.
 *  @return error code (0 if successful).
 */
//____________________________________________________________________________

int rhdeparture_( int *atomi, int *ilow, int *iup, int *nsize,
                  double *departure )
{
   char id[3];
   int  ier = 0;
   int  i, k, jtau, index;
   int  nactiveatom, nspace;

   for (i = 0; i < 2; i++) id[i] = (char) atomi[i];
   id[2] = '\0';

   if (rhdc != NULL)
   {
      nactiveatom = rhdc[0].nactiveatom;
      nspace      = rhdc[0].nspace;

      for (i = 0; i < nactiveatom; i++)
      {
         if (strcasecmp(id, rhdc[i].atomid) == 0)
         {
            if (*nsize >= nspace * 2)  // only two levels: lower & upper
            {
               for (index = 0, k = 0;  k < nspace;  k++)
               {
                  jtau = rhdc[i].nlevels * k;  // beggining of the tau segment
                  departure[index++] = rhdc[i].coef[*ilow + jtau];
                  departure[index++] = rhdc[i].coef[*iup + jtau];
               }
            }
            else   // unsufficient memory for 'departure'
            {
               ier = -1;
            }
            break;
         }
      }
      if (i == nactiveatom)   // atom not found
      {
         ier = -2;
      }
   }
   else
   {
      ier = -2;
   }

   return(ier);
}


//____________________________________________________________________________
//
//  Method: rhhpop_()
//
/** Pass the H populations from RH to SIR.
 *  @param  ntau number of optical depth.
 *  @param  nh array with the populations.
 */
//____________________________________________________________________________

void rhhpop_( int *ntau, float *nh )
{
   register int i;

   if (*ntau == rhhpop.nspace)
   {
      for (i = 0;  i < (*ntau) * 6;  i++)
      {
         nh[i] = rhhpop.pop[i];
      }
   }
}


//____________________________________________________________________________
