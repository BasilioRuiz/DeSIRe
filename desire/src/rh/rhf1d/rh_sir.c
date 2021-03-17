//
// File          _______ : rh_sir.c
// Description   _______ : Functions to pass information from SIR to RH.
// Project       _______ : DeSIRe
// Creation date _______ : 31/07/19
// Author        _______ : epm@iac.es
//


#include "desire.h"

extern SIRabun     sirabun;
extern SIRatmos    siratmos;
extern SIRBarklem  sirbarklem;
extern SIRflags    sirflags;
extern SIRkeywords sirkeywords;
extern SIRKurucz   sirkurucz;
extern SIRwave     sirwave;


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
 *  @param  ntau number of depths.
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
      siratmos.id[i] = idascii[i];
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
/** Pass SIR lines with Barklem to RH.
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
   // We supply collisional broadening for some lines from SIR.
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
//  Method: sirflags_()
//
/** Pass some command line flags from SIR to RH.
 *  @param  flagv verbose flag.
 */
//____________________________________________________________________________

void sirflags_( int *flagv )
{
   // Verbose flag.
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
 *  @param  nlines number of lines.
 *  @param  field1 array with the wavelength in air (nm).
 *  @param  field2 array with the log gf.
 *  @param  field3 array with the element code.
 *  @param  field4 array with the first energy level (cm-1).
 *  @param  field5 array with the J for first level.
 *  @param  field7 array with the second energy level (cm-1).
 *  @param  field8 array with the J for second level.
 *  @param  field10 array with the log of radiative damping constant.
 *  @param  field11 array with the log of stark damping const/electron number.
 *  @param  field12 array with the log of van der Waals damping.
 *  @param  field6a array with the spin for the first level.
 *  @param  field6b array with the orbital angular moment for the first level.
 *  @param  field9a array with the spin for the second level.
 *  @param  field9b array with the orbital angular moment for second level.
 * 
 *  @note The parameters contain a list in the format described on the Kurucz
 *        web page http://kurucz.harvard.edu/linelists.html
 */
//____________________________________________________________________________

void sirkurucz_( int *nlines,
                 double *field1, double *field2,  double *field3,
                 double *field4, double *field5,  double *field7,
                 double *field8, double *field10, double *field11,
                 double *field12,
                 int *field6a, unsigned char *field6b,
                 int *field9a, unsigned char *field9b )
{
   // We supply Kurucz data from SIR for LTE lines.
   sirkurucz.nlines   = *nlines;
   sirkurucz.lambda   = field1;
   sirkurucz.gf       = field2;
   sirkurucz.elemcode = field3;
   sirkurucz.Ei       = field4;
   sirkurucz.Ji       = field5;
   sirkurucz.spini    = field6a;
   sirkurucz.orbiti   = field6b;
   sirkurucz.Ej       = field7;
   sirkurucz.Jj       = field8;
   sirkurucz.spinj    = field9a;
   sirkurucz.orbitj   = field9b;
   sirkurucz.Grad     = field10;
   sirkurucz.GStark   = field11;
   sirkurucz.GvdWaals = field12;
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
