//
// File          _______ : rh_desire.c
// Description   _______ : Functions to pass information from DeSIRe to RH.
// Project       _______ : DeSIRe
// Creation date _______ : 31/07/19
// Author        _______ : epm@iac.es
//


#include "desire.h"

//extern DesireConfig desireconfig;
extern DesireLines desirelines;


//____________________________________________________________________________
//
//  Method: desireconfig_()
//
/** Pass DeSIRe configuration to RH.
 *  @param  flagv command line verbose flag.
 *  @param  nray number of rays.
 */
//____________________________________________________________________________

void desireconfig_( int *flagv, int *nray )
{
   // DeSIRe command line.
   //desireconfig.flagv = *flagv;

   // We supply values for some RH keywords from DeSIRe.
   //desireconfig.nrays = *nrays;
}


//____________________________________________________________________________
//
//  Method: desirelines_()
//
/** Pass DeSIRe lines to RH.
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

void desirelines_( int *nlines, int *anumber, int *stage, int *low, int *up,
                   float *alpha, float *sigma, float *wave )
{
   // We supply collisional broadening for some lines from DeSIRe.
   desirelines.nlines = *nlines;
   desirelines.anumber = anumber;
   desirelines.stage = stage;
   desirelines.low = low;
   desirelines.up = up;
   desirelines.alpha = alpha;
   desirelines.sigma = sigma;
   desirelines.wave = wave;
}


//____________________________________________________________________________
