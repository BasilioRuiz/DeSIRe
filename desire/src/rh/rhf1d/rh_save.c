//
// File          _______ : rh_save.c
// Description   _______ : Functions to save RH information.
// Project       _______ : DeSIRe
// Creation date _______ : 06/06/21
// Author        _______ : epm@iac.es
//


#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "desire.h"

extern Atmosphere atmos;

extern RHdc  *rhdc;
extern RHhpop rhhpop;


//____________________________________________________________________________
//
//  Method: save_rhdc()
//
/** Save RH departure coefficients.
 */
//____________________________________________________________________________

void save_rhdc( void )
{
   Atom  *atom;
   int    n, i, j, k;
   size_t npop;

   // I am using realloc(NULL,...) to not free this memory in case of
   // calling FreeMemoryLeakNoRealloc().
   if (rhdc == NULL)
   {
      rhdc = realloc(NULL, atmos.Nactiveatom * sizeof(RHdc));
      for (n = 0; n < atmos.Nactiveatom; n++) rhdc[n].coef = NULL;
   }
   for (n = 0; n < atmos.Nactiveatom; n++)
   {
      atom = atmos.activeatoms[n];
      npop = atom->Nlevel * atmos.Nspace * sizeof(double);

      if (rhdc != NULL && rhdc[n].coef == NULL)
      {
         strcpy(rhdc[n].atomid, atom->ID);
         rhdc[n].nlevels = atom->Nlevel;
         rhdc[n].nspace = atmos.Nspace;
         rhdc[n].nactiveatom = atmos.Nactiveatom;
         rhdc[n].coef = realloc(NULL, npop);
      }
      if (rhdc == NULL || rhdc[n].coef == NULL)
      {
         Error(ERROR_LEVEL_1, "save_rhdc",
               "Memory problem saving departure coefficients");
      }
      if (atom->nstar[0] != NULL)
      {
         for (i = 0, k = 0; k < atmos.Nspace; k++)
            for (j = 0; j < atom->Nlevel; j++)
               rhdc[n].coef[i++] = atom->n[j][k]/atom->nstar[j][k];
      }
      else
      {
         for (i = 0, k = 0; k < atmos.Nspace; k++)
            for (j = 0; j < atom->Nlevel; j++)
               rhdc[n].coef[i++] = 1.0;
      }
   }
}


//____________________________________________________________________________
//
//  Method: save_rhhpop()
//
/** Save RH H populations.
 */
//____________________________________________________________________________

void save_rhhpop( void )
{
   int i, j, k;

   // I am using realloc(NULL,...) to not free this memory in case of
   // calling FreeMemoryLeakNoRealloc().
   if (rhhpop.pop == NULL)
   {
      rhhpop.nspace = atmos.Nspace;
      rhhpop.pop = realloc(NULL, atmos.Nspace * 6 * sizeof(float));
   }
   if (rhhpop.pop == NULL)
   {
      Error(ERROR_LEVEL_1, "save_rhhpop",
            "Memory problem saving H populations");
   }
   // Unlike SIR, RH's atmosphere model runs from upper layers to lower
   // ones, so we need to turn around the array.
   for (i = 0, k = atmos.Nspace-1;  k >= 0;  k--)
   {
      for (j = 0; j < 6; j++)
      {
         rhhpop.pop[i++] = atmos.H->n[j][k]*1e-6;  // convert to CGS system
      }
   }
}


//____________________________________________________________________________
