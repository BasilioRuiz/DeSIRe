/* ------- file: -------------------------- writen_xdr.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan 26 14:50:56 2012 --

       --------------------------                      ----------RH-- */

/* --- Routines for reading and writing populations from and to file.
       XDR (external data representation) version.

       Modifications:

       - 04/04/20 epm:
         Disk access is avoided as much as possible.
         Populations are saved in memory rather than disk.
         'popsinFile' and 'popsoutFile' are ignored.

       --                                              -------------- */

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

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];

static void **pmem = NULL;


/* ------- begin -------------------------- xdr_populations.c ------- */

bool_t xdr_populations(XDR *xdrs, char *atmosID, int Nlevel, int Nspace,
		       double *n, double *nstar)
{
  const char routineName[] = "xdr_populations";

  char  *ID;
  bool_t result = TRUE;
  int    Npop = Nlevel * Nspace, Nl, Ns;

  /* --- The actual reading/writing routine. Upon input the values
         for atmosID, Nlevel and Nspace in the file are checked against
         the values supplied in the call. --           -------------- */

  if (xdrs->x_op == XDR_ENCODE) {
    result &= xdr_counted_string(xdrs, &atmosID);
    result &= xdr_int(xdrs, &Nlevel);
    result &= xdr_int(xdrs, &Nspace);
  } else {
    result &= xdr_counted_string(xdrs, &ID);
    // 04/04/20 epm: Do not compare these two variables because, if we have
    // OLD_POPULATIONS, 'ID' is delayed respect to 'atmosID' but it's right.
    // if (!strstr(ID, atmosID)) {
    //   sprintf(messageStr, "Populations were computed with different"
    //                       " atmosphere (%s) than the current one", ID);
    //   Error(WARNING, routineName, messageStr);
    // }
    free(ID);

    result &= xdr_int(xdrs, &Nl);
    result &= xdr_int(xdrs, &Ns);
    if ((Nl != Nlevel)  ||  (Ns != Nspace)) return FALSE;
  }
  /* --- Exit if true populations do not exist --      -------------- */

  if (n != NULL) {
    result &= xdr_vector(xdrs, (char *) n, Npop,
			 sizeof(double), (xdrproc_t) xdr_double);
  } else
    return FALSE;

  /* --- We can live without the LTE values since they can always be
         constructed with routine LTEpops --           -------------- */

  if (nstar != NULL)
    result &= xdr_vector(xdrs, (char *) nstar, Npop,
			 sizeof(double), (xdrproc_t) xdr_double);

  return result;
}
/* ------- end ---------------------------- xdr_populations.c ------- */

/* ------- begin -------------------------- writePopulations.c ------ */

void writePopulations( Atom *atom )
{
   const char routineName[] = "writePopulations";

   char  *atmosid    = (char *) atmos.ID;
   int    nlevel     = atom->Nlevel;
   int    nspace     = atmos.Nspace;
   off_t  offset     = 0;
   size_t npop       = atom->Nlevel * atmos.Nspace * sizeof(double);
   size_t recordsize = sizeof(atmos.ID) + sizeof(atom->Nlevel) +
                       sizeof(atmos.Nspace) + npop + npop;
   int    i;

   // 04/04/20 epm: Save the atom population in memory rather than disk.
   // I am using realloc(NULL,...) to not free this memory in case of
   // calling FreeMemoryLeakNoRealloc().
   if (pmem == NULL)
   {
     pmem = (void **) realloc(NULL, atmos.Nactiveatom * sizeof(void *));
     for (i = 0; i < atmos.Nactiveatom; i++) pmem[i] = NULL;
   }
   if (pmem != NULL && pmem[atom->activeindex] == NULL)
   {
      pmem[atom->activeindex] = realloc(NULL, recordsize);
   }
   if (pmem == NULL || pmem[atom->activeindex] == NULL)
   {
      sprintf(messageStr, "Memory problem writing populations for atom %s",
              atom->ID);
      Error(ERROR_LEVEL_1, routineName, messageStr);
   }
   memcpy(pmem[atom->activeindex] + offset, atmosid, sizeof(atmos.ID));
   offset += sizeof(atmos.ID);
   memcpy(pmem[atom->activeindex] + offset, &nlevel, sizeof(atom->Nlevel));
   offset += sizeof(atom->Nlevel);
   memcpy(pmem[atom->activeindex] + offset, &nspace, sizeof(atmos.Nspace));
   offset += sizeof(atmos.Nspace);
   if (atom->n[0] != NULL)
   {
      memcpy(pmem[atom->activeindex] + offset, atom->n[0], npop);
      offset += npop;
   }
   else   // exit if true populations do not exist
   {
      Error(ERROR_LEVEL_1, routineName,
            "Error writing populations: true populations do not exist");
   }
   if (atom->nstar[0] != NULL)   // but we can live without LTE values
   {
      memcpy(pmem[atom->activeindex] + offset, atom->nstar[0], npop);
      offset += npop;
   }
}

void writePopulations_OLD(Atom *atom)
{
  const char routineName[] = "writePopulations";

  FILE *fp_out;
  XDR   xdrs;

  // /* BRC and DOS playing */
  FILE *f;
  register int i, k;
  char filetxt[strlen(atom->popsoutFile) + 4 + 1];
  strcpy(filetxt, atom->popsoutFile);
  strcat(filetxt,".txt");

  if ((fp_out = fopen(atom->popsoutFile, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", atom->popsoutFile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  if (!xdr_populations(&xdrs, atmos.ID, atom->Nlevel, atmos.Nspace,
                       atom->n[0], atom->nstar[0])) {
    sprintf(messageStr, "Unable to write to output file %s",
            atom->popsoutFile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_out);

  // /* BRC and DOS playing */
  if ((f = fopen(filetxt, "w")) == NULL)
  {
    sprintf(messageStr, "Error opening file %s", filetxt);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  for (k = 0;  k < atmos.Nspace;  k++)
  {
    for (i = 0;  i < atom->Nlevel;  i++)
    {
      fprintf(f, "%d %d %lf\n", i, k, atom->n[i][k] / atom->nstar[i][k]);
    }
  }
  fclose(f);
}
/* ------- end ---------------------------- writePopulations.c ------ */

/* ------- begin -------------------------- readPopulations.c ------- */

void readPopulations( Atom *atom )
{
   const char routineName[] = "readPopulations";

   char   atmosid[sizeof(atmos.ID)];
   int    nlevel;
   int    nspace;
   off_t  offset     = 0;
   size_t npop       = atom->Nlevel * atmos.Nspace * sizeof(double);
   size_t recordsize = sizeof(atmos.ID) + sizeof(atom->Nlevel) +
                       sizeof(atmos.Nspace) + npop + npop;

   if (pmem != NULL && pmem[atom->activeindex] != NULL)
   {
      memcpy(atmosid, pmem[atom->activeindex] + offset, sizeof(atmos.ID));
      offset += sizeof(atmos.ID);
      memcpy(&nlevel, pmem[atom->activeindex] + offset, sizeof(atom->Nlevel));
      offset += sizeof(atom->Nlevel);
      memcpy(&nspace, pmem[atom->activeindex] + offset, sizeof(atmos.Nspace));
      offset += sizeof(atmos.Nspace);
      if (atom->n[0] != NULL)
      {
         memcpy(atom->n[0], pmem[atom->activeindex] + offset, npop);
         offset += npop;
      }
      else   // exit if true populations do not exist
      {
         Error(ERROR_LEVEL_2, routineName,
               "Error reading populations: true populations do not exist");
      }
      if (atom->nstar[0] != NULL)   // but we can live without LTE values
      {
         memcpy(atom->nstar[0], pmem[atom->activeindex] + offset, npop);
         offset += npop;
      }

      // 04/04/20 epm: Do not compare these two variables because, if we have
      // OLD_POPULATIONS, 'atmosid' is delayed respect to 'ID' but it's right.
      // if (memcmp(atmosid, atmos.ID, strlen(atmos.ID)) != 0)
      // {
      //    sprintf(messageStr, "Populations were computed with different"
      //                        " atmosphere (%s) than the current one", id);
      //    Error(WARNING, routineName, messageStr);
      // }
      if (nlevel != atom->Nlevel)
      {
         sprintf(messageStr, "Populations were computed with different"
                             " nlevel (%d) than the current one", nlevel);
         Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      if (nspace != atmos.Nspace)
      {
         sprintf(messageStr, "Populations were computed with different"
                             " nspace (%d) than the current one", nspace);
         Error(ERROR_LEVEL_2, routineName, messageStr);
      }
   }
   else
   {
      sprintf(messageStr, "Pointer null reading populations for"
              " atom %s (offset %lld, size %zu)",
              atom->ID, (long long) offset, recordsize);
      Error(ERROR_LEVEL_2, routineName, messageStr);
   }
}

void readPopulations_OLD(Atom *atom)
{
  const char routineName[] = "readPopulations";

  FILE *fp_in;
  XDR   xdrs;

  /* --- Read populations from file.

   Note: readPopulations could read the true populations and not
         the LTE populations. To this effect it should pass a NULL
         pointer to xdr_populations as its last argument.
         --                                            -------------- */

  if ((fp_in = fopen(atom->popsinFile, "r")) == NULL) {
    // 31/07/19 epm: Write the message with a suggestion.
    sprintf(messageStr,
      "Unable to open input file %s\n"
      " Are chemical symbols in the last column of %s capitalized?",
      atom->popsinFile, input.atoms_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  xdrstdio_create(&xdrs, fp_in, XDR_DECODE);

  if (!xdr_populations(&xdrs, atmos.ID, atom->Nlevel, atmos.Nspace,
		       atom->n[0], atom->nstar[0])) {
    sprintf(messageStr, "Unable to read from input file %s",
	    atom->popsinFile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_in);
}
/* ------- end ---------------------------- readPopulations.c ------- */

/* ------- begin -------------------------- writeMolPops.c ---------- */

void writeMolPops(struct Molecule *molecule)
{
  const char routineName[] = "writeMolPops";

  FILE *fp_out;
  XDR   xdrs;

  /* --- Write molecular (vibration) populations to file. -- -------- */


  if ((fp_out = fopen(molecule->popsFile, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s",
	    molecule->popsFile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  if (!xdr_populations(&xdrs, atmos.ID, molecule->Nv, atmos.Nspace,
		       molecule->nv[0], molecule->nvstar[0])) {
    sprintf(messageStr, "Unable to write to output file %s",
	    molecule->popsFile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeMolPops.c ---------- */

/* ------- begin -------------------------- readMolPops.c ----------- */

void readMolPops(struct Molecule *molecule)
{
  const char routineName[] = "readMolPops";

  FILE *fp_in;
  XDR   xdrs;

  /* --- Read molecular (vibration) populations from file.

   Note: readMolPops only reads the true populations and not
         the LTE populations. To this effect it passes a NULL pointer
         to xdr_populations as its last argument.
         --                                            -------------- */

  if ((fp_in = fopen(molecule->popsFile, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s",
	    molecule->popsFile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  xdrstdio_create(&xdrs, fp_in, XDR_DECODE);

  if (!xdr_populations(&xdrs, atmos.ID, molecule->Nv, atmos.Nspace,
		       molecule->nv[0], NULL)) {
    sprintf(messageStr, "Unable to read from input file %s",
	    molecule->popsFile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_in);
}
/* ------- end ---------------------------- readPopulations.c ------- */
