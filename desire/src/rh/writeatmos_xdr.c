/* ------- file: -------------------------- writeatmos_xdr.c --------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May  7 16:09:49 2021 --

       --------------------------                      ----------RH-- */

/* --- Writes atmospheric data to output file.
       XDR (external data representation) version. --  -------------- */

 
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- xdr_counted_string.c ---- */

bool_t xdr_counted_string(XDR *xdrs, char **p)
{
  bool_t output = (xdrs->x_op == XDR_ENCODE) ? TRUE : FALSE; 
  short  length;

  /* --- Reads or writes string of length length from or to xdr stream
         pointed to by xdrs.

    See: IDL User's Guide 17-34
         --                                            -------------- */

  if (output) length = strlen(*p);

  if (!xdr_short(xdrs, &length))  return FALSE;
  if (!output) {
    *p = (char *) malloc((length + 1) * sizeof(char));
    (*p)[length] = '\0';
  }
  return (length ? xdr_string(xdrs, p, length) : TRUE);
}
/* ------- end ---------------------------- xdr_counted_string.c ---- */

/* ------- begin -------------------------- writeAtmos.c ------------ */

void writeAtmos(Atmosphere *atmos)
{
  const char routineName[] = "writeAtmos";
  register int n;

  bool_t result = TRUE;
  int    Nspace = atmos->Nspace;
  char  *elemID, *atmosID;
  FILE  *fp_out;
  XDR    xdrs;

  if (!strcmp(input.atmos_output, "none")) return;

  if ((fp_out = fopen(input.atmos_output, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", input.atmos_output);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  result &= xdr_int(&xdrs, &atmos->H->Nlevel);
  result &= xdr_int(&xdrs, &atmos->Nelem);
  result &= xdr_bool(&xdrs, &atmos->moving);

  result &= xdr_vector(&xdrs, (char *) atmos->T, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(&xdrs, (char *) atmos->ne, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(&xdrs, (char *) atmos->vturb, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(&xdrs, (char *) atmos->H->n[0],
		       Nspace * atmos->H->Nlevel,
		       sizeof(double), (xdrproc_t) xdr_double);

  atmosID = atmos->ID;
  result &= xdr_counted_string(&xdrs, &atmosID);

  for (n = 0;  n < atmos->Nelem;  n++) {
    elemID = atmos->elements[n].ID;
    result &= xdr_counted_string(&xdrs, &elemID);
    result &= xdr_double(&xdrs, &atmos->elements[n].weight);
    result &= xdr_double(&xdrs, &atmos->elements[n].abund);
  }
  /* --- In case magnetic fields were present --       -------------- */

  if (atmos->Stokes) {
    result &= xdr_bool(&xdrs, &atmos->Stokes);

    result &= xdr_vector(&xdrs, (char *) atmos->B, Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(&xdrs, (char *) atmos->gamma_B, Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(&xdrs, (char *) atmos->chi_B, Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
  }    

  if (!result) {
    sprintf(messageStr, "Unable to write proper amount to output file %s",
	    input.atmos_output);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeAtmos.c ------------ */


/* ------- begin -------------------------- freeAtmos.c ------------- */

void freeAtmos(Atmosphere *atmos)
{
  register int n;
  Element *element;

  /* --- Free memory associated with atmospheric structure -- ------- */

  free(atmos->N);      // 06/06/21 epm: add this free()
  free(atmos->T);
  free(atmos->ne);
  free(atmos->vturb);  // 06/06/21 epm: add this free()
  free(atmos->nHtot);
  free(atmos->nHmin);

  if (atmos->Stokes) {
    free(atmos->B);
    free(atmos->gamma_B);
    free(atmos->chi_B);

    freeMatrix((void **) atmos->cos_gamma);
    freeMatrix((void **) atmos->cos_2chi);
    freeMatrix((void **) atmos->sin_2chi);
  }

  for (n = 0;  n < atmos->Natom;  n++) {
    if (atmos->atoms[n].active  ||  
        atmos->hydrostatic  ||
        input.solve_ne == ITERATION) 
      freeAtom(&atmos->atoms[n]);
  }
  free(atmos->activeatoms);
  
  for (n = 0;  n < atmos->Nelem;  n++) {
    element = &atmos->elements[n];
    if (element->Nmolecule > 0)
      free(element->mol_index);
  }
  free(atmos->elements);

  for (n = 0;  n < atmos->Nmolecule;  n++) {
    if (atmos->molecules[n].active  ||  
        atmos->hydrostatic  ||
        input.solve_ne == ITERATION) 
      freeMolecule(&atmos->molecules[n]);
  }
  free(atmos->activemols);

  free(atmos->backgrrecno);
  free(atmos->backgrflags);

  // 06/06/21 epm: Free also this array (allocated en hydrogen.c line 110).
  if (!atmos->H_LTE)
    if (!atmos->H->active) freeMatrix((void **) atmos->H->n);

  // 06/06/21 epm: atmos.nH is freed in hydrogen.c line 151.
  // freeMatrix((void **) atmos.nH);
}
/* ------- end ---------------------------- freeAtmos.c ------------- */
