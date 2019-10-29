/* ------- file: -------------------------- desire.h ----------------

       Version:       DeSIRe v1.4
       Author:        epm@iac.es
       Last modified: July 31 2019

       --------------------------                      ----------RH-- */

#ifndef __DESIRE_H__
#define __DESIRE_H__


typedef struct
{
  int    nlines;
  int   *anumber, *stage, *low, *up;
  float *alpha, *sigma, *wave;
}
DesireLines;


#endif /* !__DESIRE_H__ */

/* ------- end ---------------------------- desire.h ---------------- */
