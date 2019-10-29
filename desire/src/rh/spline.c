/* ------- file: -------------------------- spline.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 16 14:57:56 1999   --

       --------------------------                      ----------RH-- */

/* --- Cubic spline interpolation routines. --         -------------- */

 
#include <stdio.h>
#include <stdlib.h>

#include "rh.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

// 01/07/19 epm: Rename these global variables appending "2"
// to avoid the same name that is used in "expspline.c".
static bool_t  ascend2;
static int     Ntable2;
static double *xtable2, xmin2, xmax2, *M2 = NULL;
static double *ytable2;


/* ------- begin -------------------------- splineCoef.c ------------ */

void splineCoef(int N, double *x, double *y)
{
  register int j;
  static double *u = NULL;

  double  p, *q, hj, hj1, D, D1, mu;

  ascend2 = (x[1] > x[0]) ? TRUE : FALSE;
  xmin2 = (ascend2) ? x[0] : x[N-1];
  xmax2 = (ascend2) ? x[N-1] : x[0];

  q = M2 = (double *) realloc(M2, N * sizeof(double));
  u = (double *) realloc(u, N * sizeof(double));
  hj = x[1] - x[0];
  D  = (y[1] - y[0]) / hj;

  q[0] = u[0] = 0.0;
  for (j = 1;   j < N-1;  j++) {
    hj1 = x[j+1] - x[j];
    mu  = hj / (hj + hj1);
    D1  = (y[j+1] - y[j]) / hj1;

    p = mu*q[j-1] + 2;
    q[j] = (mu - 1) / p;
    u[j] = ((D1 - D) * 6/(hj + hj1) - mu*u[j-1]) / p;

    hj = hj1;  D = D1;
  }

  M2[N - 1] = 0.0;
  for (j = N-2;  j >= 0;  j--) {
    M2[j] = q[j]*M2[j+1] + u[j];
  }
  Ntable2 = N;
  xtable2 = x;  ytable2 = y;
}
/* ------- end ---------------------------- splineCoef.c ------------ */

/* ---------------------------------------- splineEval.c ------------ */

void splineEval(int N, double *x, double *y, bool_t hunt)
{
  register int n;

  int    j = 0;
  double hj, fx, fx1;

  for (n = 0;  n < N;  n++) {
    if (x[n] <= xmin2)
      y[n] = (ascend2) ? ytable2[0] : ytable2[Ntable2-1];
    else if (x[n] >= xmax2)
      y[n] = (ascend2) ? ytable2[Ntable2-1] : ytable2[0];
    else {
      if (hunt) 
	Hunt(Ntable2, xtable2, x[n], &j);
      else
	Locate(Ntable2, xtable2, x[n], &j);

      hj  = xtable2[j+1] - xtable2[j];
      fx  = (x[n] - xtable2[j]) / hj;
      fx1 = 1 - fx;

      y[n] = fx1*ytable2[j] + fx*ytable2[j+1] +
	(fx1*(SQ(fx1) - 1) * M2[j] + fx*(SQ(fx) - 1) * M2[j+1]) * SQ(hj)/6.0;
    }
  }
}
/* ------- end ---------------------------- splineEval.c ------------ */
