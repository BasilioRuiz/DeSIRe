/* ------- file: -------------------------- desire.h ----------------

       Version:       DeSIRe v2.10
       Author:        epm@iac.es
       Last modified: November 11 2020

       --------------------------                      ----------RH-- */

#ifndef __DESIRE_H__
#define __DESIRE_H__


// 10/10/19 epm: Structure to pass Barklem data from SIR to RH.
typedef struct
{
   int    nlines;
   int   *anumber, *stage, *low, *up;
   float *alpha, *sigma, *wave;
}
SIRBarklem;

// 04/04/20 epm: Structure to pass some command line flags from SIR to RH.
typedef struct
{
   int flagv;
}
SIRflags;

// 07/07/20 epm: Structure to pass Kurucz data from SIR to RH.
typedef struct
{
   int            nlines;
   int           *spini, *spinj;
   double        *lambda, *gf, *elemcode;
   double        *Ei, *Ji, *Ej, *Jj, *Grad, *GStark, *GvdWaals;
   unsigned char *orbiti, *orbitj;
}
SIRKurucz;

// 08/08/20 epm: Structure to pass abundance values from SIR to RH.
typedef struct
{
   int    nelem;
   float *abun;
   char  *id[99];
}
SIRabun;

// 10/10/20 epm: Structure to pass the wavetable from SIR to RH.
typedef struct
{
   int     nwave;
   double *wavetable;
}
SIRwave;

// 10/10/20 epm: Structure to pass some keywords from SIR to RH.
typedef struct
{
   int nrays;
   int startj;
}
SIRkeywords;

// 11/11/20 epm: Structure to pass atmospheric models from SIR to RH.
typedef struct
{
   char    id[100];
   char    scale;
   int     ntau;
   int     nB;
   double  logg;
   double *dscale;
   float  *T, *elecdens, *vz, *vmic;
   float  *B, *gamma, *phi;
   float  *nh1, *nh2, *nh3, *nh4, *nh5, *nh6;
}
SIRatmos;


#endif /* !__DESIRE_H__ */

/* ------- end ---------------------------- desire.h ---------------- */
