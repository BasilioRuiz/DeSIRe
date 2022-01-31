/* ------- file: -------------------------- desire.h ----------------

       Version:       DeSIRe v3.10
       Author:        epm@iac.es
       Last modified: October 10 2021

       --------------------------                      ----------RH-- */

#ifndef __DESIRE_H__
#define __DESIRE_H__


/*------------------- PASSING INFORMATION FROM SIR TO RH -------------------*/

// 08/08/20 epm: Structure to pass abundance values from SIR to RH.
typedef struct
{
   int    nelem;    // number of elements
   float *abun;     // array with the abundances
   char  *id[99];   // element (initialized in "rh_glob.c")
}
SIRabun;

// 11/11/20 epm: Structure to pass atmospheric models from SIR to RH.
typedef struct
{
   // Atmospheric models (written in "siratmos_.c").
   char    id[100];    // ID model string
   char    scale;      // scale identifier ('M', 'T' or 'H')
   int     ntau;       // number of optical depths
   int     nB;         // number of elements in B
   double  logg;       // gravity decimal logarithm
   double *dscale;     // array with depth scales
   float  *T;          // array with the temperatures
   float  *elecdens;   // array with the electronic densities (cm-3)
   float  *vz;         // array with the line of sight velocity (km/s)
   float  *vmic;       // array with the microturbulence velocity (km/s)
   float  *B;          // array with the magnetic field strength (Teslas)
   float  *gamma;      // array with magnetic field inclination (radians)
   float  *phi;        // array with magnetic field azimuth (radians)
   // Hydrogen populations (written in "siratmosnh_.c").
   float  *nh1;        // population of the fundamental level
   float  *nh2;        // population of the first excited level
   float  *nh3;        // population of the second excited level
   float  *nh4;        // population of the third excited level
   float  *nh5;        // population of the fourth excited level
   float  *nh6;        // number of protons per cm3
}
SIRatmos;

// 10/10/19 epm: Structure to pass Barklem data from SIR to RH.
typedef struct
{
   int    nlines;    // number of lines
   int   *anumber;   // array with the atomic numbers
   int   *stage;     // array with the ionization stages
   int   *low;       // array with the low levels
   int   *up;        // array with the up levels
   float *alpha;     // array with the alpha values
   float *sigma;     // array with the sigma values
   float *wave;      // array with the air wavelength values
}
SIRBarklem;

// 04/04/21 epm: Structure to pass coupling data from SIR to RH.
typedef struct
{
   // J : total angular momentum
   // L : angular (orbital) momentum of the atom
   // S : spin of the atom
   // L1: angular momentum of the parent
   // S1: parent's spin
   // J1: total angular momentum of the parent
   // l : electron's angular momentum
   // K : orbital angular momentum of the coupling
   // j1: first component of total angular momentum
   // l1: first component of (orbital) angular momentum 
   // j2: second component of total angular momentum
   // l2: second component of (orbital) angular momentum 
   int    nlines;   // number of LTE lines
   int   *cpl;      // coupling scheme (0=LS, 1=JK, 2=JJ)
   float *qn1;      // J if LS | J  if JK | J  if JJ
   float *qn2;      // L if LS | L1 if JK | j1 if JJ
   float *qn3;      // S if LS | S1 if JK | l1 if JJ
   float *qn4;      // 0 if LS | J1 if JK | j2 if JJ
   float *qn5;      // 0 if LS | l  if JK | l2 if JJ
   float *qn6;      // 0 if LS | K  if JK | 0  if JJ
}
SIRcpl;

// 04/04/20 epm: Structure to pass some command line flags from SIR to RH.
typedef struct
{
   int flagv;   // verbose flag
}
SIRflags;

// 10/10/20 epm: Structure to pass some keywords from SIR to RH.
typedef struct
{
   int nrays;    // number of rays
   int startj;   // starting solution for J
}
SIRkeywords;

// 07/07/20 epm: Structure to pass Kurucz data from SIR to RH.
typedef struct
{
   int     nlines;     // number of LTE lines
   double *lambda;     // array with the wavelength in air (nm)
   double *gf;         // array with the log gf
   double *elemcode;   // array with the element code
   double *Ei;         // array with the first energy level (cm-1)
   double *Ej;         // array with the second energy level (cm-1)
   double *Grad;       // array with the log of radiative damping
   double *GStark;     // array with the log of stark damping
   double *GvdWaals;   // array with the log of van der Waals damping
}
SIRKurucz;

// 10/10/20 epm: Structure to pass the wavetable from SIR to RH.
typedef struct
{
   int     nwave;       // number of wavelengths
   double *wavetable;   // array with the wavetable
}
SIRwave;


/*------------------- PASSING INFORMATION FROM RH TO SIR -------------------*/

// 06/06/21 epm: Structure to pass departure coefficients from RH to SIR.
typedef struct
{
   char    atomid[3];     // ID of the atom
   int     nlevels;       // number of levels
   int     nspace;        // number of optical depths
   int     nactiveatom;   // number of active atoms
   double *coef;          // departure coefficients
}
RHdc;

// 06/06/21 epm: Structure to pass H populations from RH to SIR.
typedef struct
{
   int    nspace;   // number of optical depths
   float *pop;      // array with the H populations
}
RHhpop;


#endif /* !__DESIRE_H__ */

/* ------- end ---------------------------- desire.h ---------------- */
