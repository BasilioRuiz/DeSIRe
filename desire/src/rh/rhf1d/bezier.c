/* ------- file: -------------------------- bezier.c ----------------

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   Modifications:
           2017-03-12, JdlCR: Created!

       Last modified: Wed May 23 09:05:19 2018 --

       --------------------------                      ----------RH-- */


#include <math.h>
#include <string.h>    // memcpy, memset
#include <x86intrin.h> // Intrinsic SSE instructions

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "bezier.h"


/* --- Macros --                                       -------------- */

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];


/* --- Identity matrix --                              -------------- */

static const double ident[4][4] =
  {{1.0, 0.0, 0.0, 0.0},
   {0.0, 1.0, 0.0, 0.0},
   {0.0, 0.0, 1.0, 0.0},
   {0.0, 0.0, 0.0, 1.0}};


/* ------- begin -------------------------- cent_deriv.c ------------ */

inline double cent_deriv(double dsup,double dsdn, 
			 double chiup,double chic, double chidn)
{
  /* --- Derivative Fritsch & Butland (1984) --        -------------- */

  double fim1, fi, alpha, wprime;
  
  fim1=(chic-chiup)/dsup;
  fi=(chidn-chic)/dsdn;

  if (fim1*fi > 0) {
    alpha = 0.333333333333333333333333 * ( 1.0 + dsdn/(dsdn+dsup) );
    wprime = (fim1*fi) / ( (1.0-alpha) * fim1 + alpha*fi );
  } else {
    wprime=0.0;
  }
  return wprime;
}
/* ------- end ---------------------------- cent_deriv.c ------------ */


/* ------- begin -------------------------- cent_deriv_mat.c -------- */

inline void cent_deriv_mat(double wprime[4][4], double dsup, double dsdn,
			   double chiup[4][4], double chic[4][4],
			   double chidn[4][4])
{
  register int i,j;
  
  for(j=0;j<4;j++)
    for(i=0;i<4;i++)
      wprime[j][i] = cent_deriv(dsup, dsdn, chiup[j][i], chic[j][i],
				chidn[j][i]);
}
/* ------- end ---------------------------- cent_deriv_mat.c -------- */


/* ------- begin -------------------------- cent_deriv_vec.c -------- */

inline void cent_deriv_vec(double wprime[4], double dsup, double dsdn,
		    double chiup[4], double chic[4], double chidn[4])
{
  register int i;
  
  for(i=0;i<4;i++)
    wprime[i] = cent_deriv(dsup, dsdn, chiup[i], chic[i], chidn[i]);
  
}
/* ------- end ---------------------------- cent_deriv_vec.c -------- */


/* ------- begin -------------------------- m4m.c ------------------- */


inline void m4m(double a[4][4], double b[4][4], double c[4][4])
{

  /* --- Matrix multiplication --                      -------------- */

  register int i, j, k;
  memset(&c[0][0],0,sizeof(double)*16);
  for(j = 0; j<4; j++)
    for(i = 0; i<4; i++)
      for(k = 0; k<4; k++)
	c[j][i] += a[k][i]*b[j][k]; 
}
/* ------- end ---------------------------- m4m.c ------------------- */


/* ------- begin -------------------------- m4v.c ------------------- */

/* --- Matrix/vector multiplication.
       We use matrix as float to be able to use Intel's
       matrix inversion --                         ------------------ */

inline void m4v(float a[4][4], double b[4], double c[4])
{
  register int k, i;
  memset(&c[0],0,sizeof(double)*4);
  for(i = 0; i<4; i++)
    for(k = 0; k<4; k++)
      c[i] += ((double)a[i][k]) * b[k];
}
/* ------- end ---------------------------- m4v.c ------------------- */


/* ------- begin -------------------------- Svec.c ------------------ */

inline void Svec(int k, double **S, double *Sf)
{
  /* --- Extracts the Source vector at depth-points k -- ------------ */

  Sf[0] = S[0][k], Sf[1] = S[1][k], Sf[2] = S[2][k], Sf[3] = S[3][k];
}
/* ------- end ---------------------------- Svec.c ------------------ */


/* ------- begin -------------------------- SIMD_MatInv.c ----------- */

void SIMD_MatInv(float* src)
{
  /* --- 

     Very fast in-place 4x4 Matrix inversion using SIMD instrutions
     Only works with 32-bits floats. It uses Cramer's rule.
     
     Provided by Intel

     Requires SSE instructions but all x86 machines since 
     Pentium III have them.
     
     --                                            ------------------ */
  
  __m128 minor0, minor1, minor2, minor3;
  __m128 row0, row1, row2, row3;
  __m128 det, tmp1;
  
  // -----------------------------------------------
  tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src)), (__m64*)(src+ 4));
  row1 = _mm_loadh_pi(_mm_loadl_pi(row1, (__m64*)(src+8)), (__m64*)(src+12));
  row0 = _mm_shuffle_ps(tmp1, row1, 0x88);
  row1 = _mm_shuffle_ps(row1, tmp1, 0xDD);
  tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src+ 2)), (__m64*)(src+ 6));
  row3 = _mm_loadh_pi(_mm_loadl_pi(row3, (__m64*)(src+10)), (__m64*)(src+14));
  row2 = _mm_shuffle_ps(tmp1, row3, 0x88);
  row3 = _mm_shuffle_ps(row3, tmp1, 0xDD);
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row2, row3);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor0 = _mm_mul_ps(row1, tmp1);
  minor1 = _mm_mul_ps(row0, tmp1);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
  minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
  minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row1, row2);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
  minor3 = _mm_mul_ps(row0, tmp1);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
  minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
  minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  row2 = _mm_shuffle_ps(row2, row2, 0x4E);
  minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
  minor2 = _mm_mul_ps(row0, tmp1);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
  minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
  minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row0, row1);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
  minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
  minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row0, row3);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
  minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
  minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));
  // -----------------------------------------------
  tmp1 = _mm_mul_ps(row0, row2);
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
  minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
  minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));
  tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
  minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
  minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);
  // -----------------------------------------------
  det = _mm_mul_ps(row0, minor0);
  det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
  det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
  tmp1 = _mm_rcp_ss(det);
  det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1),
		   _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
  det = _mm_shuffle_ps(det, det, 0x00);
  minor0 = _mm_mul_ps(det, minor0);
  _mm_storel_pi((__m64*)(src), minor0);
  _mm_storeh_pi((__m64*)(src+2), minor0);
  minor1 = _mm_mul_ps(det, minor1);
  _mm_storel_pi((__m64*)(src+4), minor1);
  _mm_storeh_pi((__m64*)(src+6), minor1);
  minor2 = _mm_mul_ps(det, minor2);
  _mm_storel_pi((__m64*)(src+ 8), minor2);
  _mm_storeh_pi((__m64*)(src+10), minor2);
  minor3 = _mm_mul_ps(det, minor3);
  _mm_storel_pi((__m64*)(src+12), minor3);
  _mm_storeh_pi((__m64*)(src+14), minor3);
}
/* ------- end ---------------------------- SIMD_MatInv.c ----------- */


/* ------- begin -------------------------- Bezier3_coeffs ---------- */

inline void Bezier3_coeffs(double dt, double *alpha, double *beta,
		    double *gamma, double *theta, double *eps)
{
  /* --- Integration coeffs. for cubic Bezier interpolants
         Use Taylor expansion if dtau is small --  ------------------ */
  
  double dt2 = dt*dt, dt3 = dt2 * dt,  dt4;
    
  if(dt >= 5.e-2){

    *eps = exp(-dt);

    *alpha = (-6.0 + 6.0 * dt - 3.0 * dt2 + dt3 + 6.0 * eps[0])        / dt3;
    dt3 = 1.0/dt3;
    *beta  = (6.0 + (-6.0 - dt * (6.0 + dt * (3.0 + dt))) * eps[0])    * dt3;
    *gamma = 3.0 * (6.0 + (-4.0 + dt)*dt - 2.0 * (3.0 + dt) * eps[0])  * dt3;
    *theta = 3.0 * ( eps[0] * (6.0 + dt2 + 4.0 * dt) + 2.0 * dt - 6.0) * dt3;
  } else{
    dt4 = dt2*dt2;
    *eps = 1.0 - dt + 0.5 * dt2 - dt3 / 6.0 + dt4 / 24.0;

    *alpha = 0.25 * dt - 0.05 * dt2 + dt3 / 120.0 - dt4 / 840.0;
    *beta  = 0.25 * dt - 0.20 * dt2 + dt3 / 12.0  - dt4 / 42.0; 
    *gamma = 0.25 * dt - 0.10 * dt2 + dt3 * 0.025 - dt4 / 210.0; 
    *theta = 0.25 * dt - 0.15 * dt2 + dt3 * 0.05  - dt4 / 84.0; 
  }
}
/* ------- end ---------------------------- Bezier3_coeffs ---------- */


/* ------- begin -------------------------- PiecewiseStokesBezier3 -- */

void PiecewiseStokesBezier3(int nspect, int mu, bool_t to_obs,
			    double *chi, double **S, double **I,
			    double *Psi)
{
  /* --- Cubic DELO-Bezier solver for polarized light
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

         Reference(s):
         J. de la Cruz Rodriguez & N. Piskunov (2013)
         --                                        ------------------ */
  
  const char routineName[] = "PiecewiseStokesBezier3";
  register int k, n, m, i, j;
  
  int    Ndep = geometry.Ndep, k_start, k_end, dk;
  double dtau_uw, dtau_dw = 0.0, c1, c2, w[3], dsdn2, dchi_dn,
         I_upw[4], Bnu[2];
  double dchi_up,dchi_c,dt03;
  double dsup,dsdn,dt,eps=0,alpha=0,beta=0,gamma=0,theta=0;
  double Ku[4][4], K0[4][4], Kd[4][4], dKu[4][4], dK0[4][4];
  double Su[4], S0[4], Sd[4], dSu[4], dS0[4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4], V1[4];
  double imu = 1.0 / geometry.muz[mu];
  float Md[4][4];
  double *z = geometry.height;

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = 0.5 * imu * (chi[k_start] + chi[k_start+dk]) *
    fabs(z[k_start] - z[k_start+dk]);
  
  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      for (n = 0;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_upw[0] = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case IRRADIATED:
      I_upw[0] = geometry.Ibottom[nspect][mu];
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      for (n = 0;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case IRRADIATED:
      I_upw[0] = geometry.Itop[nspect][mu];
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    default:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[TOP]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }

  for (n = 0;  n < 4;  n++) I[n][k_start] = I_upw[n];
  if (Psi) Psi[k_start] = 0.0;
  
  k=k_start+dk;
  dsup = fabs(z[k] - z[k-dk]) * imu;
  dsdn = fabs(z[k+dk] - z[k]) * imu;
  dchi_up= (chi[k] - chi[k-dk])/dsup;

  
  /* ---  dchi/ds at central point--               ------------------ */
  
  dchi_c = cent_deriv(dsup,dsdn,chi[k-dk],chi[k],chi[k+dk]);
  
  
  /* --- Upwind path_length (BEzier3 integration) -- ---------------- */

  c2 = max(chi[k]    - (dsup/3.0) * dchi_c,  0.0);
  c1 = max(chi[k-dk] + (dsup/3.0) * dchi_up, 0.0);
  
  dtau_uw = 0.25 * dsup * (chi[k] + chi[k-dk] + c1 + c2);

  
  /* --- Ku, K0 and dKu, dSu -                     ------------------ */
  
  StokesK(nspect, k_start,    chi[k_start],    Ku);
  StokesK(nspect, k_start+dk, chi[k_start+dk], K0);

  Svec(k_start,    S, Su);
  Svec(k_start+dk, S, S0);

  /* --- Assume side derivative in the first interval -- ------------ */
  
  for(n = 0;  n < 4;  n++){
    dSu[n] = (S0[n] - Su[n]) / dtau_uw;
    
    for(m = 0;  m < 4;  m++)
      dKu[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
  }
  
  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end;  k += dk) {      
      
    /* --- dchi/ds at downwind point --                -------------- */
      
    dsdn = fabs(z[k+dk] - z[k]) * imu;
      
    if(abs(k - k_end) > 1){
      dsdn2 = fabs(z[k+2*dk] - z[k+dk]) * imu;
      dchi_dn = cent_deriv(dsdn, dsdn2, chi[k], chi[k+dk], chi[k+2*dk]);       
    } else
      dchi_dn = (chi[k+dk] - chi[k])/dsdn;      
      
    /* --- Make sure that c1 and c2 don't do below zero -- ---------- */
      
    c2 = max(chi[k]    + (dsdn/3.0) * dchi_c , 0.0);
    c1 = max(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);
          
    /* --- Bezier3 integrated dtau --              ------------------ */
      
    dtau_dw = 0.25 * dsdn * (chi[k] + chi[k+dk] + c1 + c2);
    dt = dtau_uw, dt03 = dt / 3.0;
  
    /* --- Bezier3 coeffs. --                      ------------------ */
      
    Bezier3_coeffs(dt, &alpha, &beta, &gamma, &theta, &eps);
   
    /* --- Diagonal operator --                    ------------------ */
      
    if(Psi) Psi[k] = alpha + gamma;
   
    /* ---- get algebra in place --                ------------------ */
      
    StokesK(nspect, k+dk, chi[k+dk], Kd);
    Svec(k+dk, S, Sd);

    cent_deriv_mat(dK0, dtau_uw, dtau_dw, Ku, K0, Kd);
    cent_deriv_vec(dS0, dtau_uw, dtau_dw, Su, S0, Sd);

    m4m(Ku, Ku, Ma); // Ku # Ku
    m4m(K0, K0, A ); // K0 # K0

    for(j = 0;  j < 4;  j++){
      for(i = 0;  i < 4;  i++){
	Md[j][i] = ident[j][i] + alpha * K0[j][i] - gamma *
	  -(dt03 * (A[j][i] + dK0[j][i] + K0[j][i]) + K0[j][i]);
	  
	Ma[j][i] = eps * ident[j][i] - beta * Ku[j][i] + theta *
	  (dt03 * (Ma[j][i] + dKu[j][i] + Ku[j][i]) - Ku[j][i]);
	  
	Mb[j][i] = beta * ident[j][i] + theta * (ident[j][i] -
						 dt03 * Ku[j][i]);
	Mc[j][i] = alpha* ident[j][i] + gamma * (ident[j][i] +
						 dt03 * K0[j][i]);
      }
    }
      
    /* --- Here I am doing Ma*stk + Mb * Su + Mc * S0 + 
           (gam * dS0 - theta * dSu) * dtau / 3.0 to compute the 
           right-hand term
           --                                      ------------------ */
    
    memset(V0, 0, 4*sizeof(double));
    
    for(i = 0;  i < 4;  i++){
      for(j = 0;  j < 4;  j++){
	V0[i] += Ma[i][j] * I[j][k-dk] + Mb[i][j] * Su[j] +
	  Mc[i][j] * S0[j];
      }
      V0[i] += dt03 * (gamma * dS0[i] - theta * dSu[i]);
    }
    /* --- Solve linear system to get the intensity -- -------------- */
      
    SIMD_MatInv(Md[0]);   // Invert Md
    m4v(Md, V0, V1);      // Multiply Md^-1 * V0

    for(i=0;i<4;i++) I[i][k] = V1[i];
      
      
    /* --- Shift values for next depth --          ------------------ */
      
    memcpy(Su,   S0, 4*sizeof(double));
    memcpy(S0,   Sd, 4*sizeof(double));
    memcpy(dSu, dS0, 4*sizeof(double));
      
    memcpy(Ku[0],   K0[0], 16*sizeof(double));
    memcpy(K0[0],   Kd[0], 16*sizeof(double));
    memcpy(dKu[0], dK0[0], 16*sizeof(double));
      
    dtau_uw = dtau_dw;
    dsup    = dsdn;
    dchi_up = dchi_c;
    dchi_c  = dchi_dn;      
  }
      
  /* --- Linear integration in the last interval -- ----------------- */
  
  k = k_end;
  dtau_uw = 0.5*imu * (chi[k] + chi[k-dk]) *
    fabs(geometry.height[k] - geometry.height[k-dk]);
  w3(dtau_uw, w);

  /* --- dSu is defined negative in Han's implementation ------------ */
  
  for (n = 0;  n < 4;  n++)
    V0[n] = w[0]*S[n][k] + w[1] * -dSu[n];
  
  if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
      
  for (n = 0;  n < 4;  n++) {
    for (m = 0;  m < 4;  m++) {
      A[n][m]  = -w[1]/dtau_uw * Ku[n][m];
      Md[n][m] = (w[0] - w[1]/dtau_uw) * K0[n][m];
    }
    A[n][n]  = 1.0 - w[0];
    Md[n][n] = 1.0;
  }
      
  for (n = 0;  n < 4;  n++) 
    for (m = 0;  m < 4;  m++) 
      V0[n] += A[n][m] * I[m][k-dk];

  /* --- Solve linear system --                    ------------------ */
  
  SIMD_MatInv(Md[0]); // Invert Md
  m4v(Md,V0,V1);      // Multiply Md^-1 * V0
  
  for (n = 0;  n < 4;  n++) I[n][k] = V1[n];
  }
/* ------- end ------------------------- PiecewiseStokesBezier3 -- -- */


/* ------- begin ----------------------- Piecewise_Bezier3 -- ------- */

void Piecewise_Bezier3(int nspect, int mu, bool_t to_obs,
		       double *chi, double *S, double *I, double *Psi)
{
  
  /* --- Cubic Bezier solver for unpolarized light
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

         Reference:
         J. de la Cruz Rodriguez & N. Piskunov (2013), Auer (2003)

         Comments: 
           JdlCR: We only check that the control points of the opacity
	          and source function are always above zero to avoid having
	          a negative interpolant.
         --                                            -------------- */
  
  register int k;
  const char routineName[] = "Piecewise_Bezier3";

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, c1, c2, w[3],
         zmu, Bnu[2];
  double dsup,dsdn,dt,dt03,eps=0,alpha=0,beta=0,gamma=0,theta=0;
  double dS_up,dS_c,dchi_up,dchi_c,dchi_dn,dsdn2;

  zmu = 1.0 / geometry.muz[mu];

  /* --- Distinguish between rays going from BOTTOM to TOP
         (to_obs == TRUE), and vice versa --           -------------- */

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  
  dtau_uw = 0.5 * zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);
  dS_uw = (S[k_start] - S[k_start+dk]) / dtau_uw;

  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      I_upw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_upw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      break;
    case IRRADIATED:
      I_upw = geometry.Ibottom[nspect][mu];
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      I_upw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[0], spectrum.lambda[nspect], Bnu);
      I_upw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
      break;
    case IRRADIATED:
      I_upw = geometry.Itop[nspect][mu];
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
  
  I[k_start] = I_upw;
  if (Psi) Psi[k_start] = 0.0;
  
  /* --- Set variables for first iteration to allow simple 
         shift for all next iterations --              -------------- */

  k=k_start+dk;
  dsup = fabs(geometry.height[k] - geometry.height[k-dk]) * zmu;
  dsdn = fabs(geometry.height[k+dk] - geometry.height[k]) * zmu;
  dchi_up= (chi[k] - chi[k-dk])/dsup;
  
  /* --- dchi/ds at central point --                   -------------- */

  dchi_c = cent_deriv(dsup,dsdn,chi[k-dk],chi[k],chi[k+dk]);
  
  /* --- upwind path_length (Bezier3 integration) --   -------------- */

  c1 = max(chi[k] - (dsup/3.0) * dchi_c, 0.0);
  c2 = max(chi[k-dk] + (dsup/3.0) * dchi_up,  0.0);
  dtau_uw =  dsup * (chi[k] + chi[k-dk] + c1 + c2) * 0.25;

  /* dS/dtau at upwind point */

  dS_up = (S[k]-S[k-dk]) / dtau_uw;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    
    if (k != k_end) {

      /* --- Downwind path length --                   -------------- */
      
       dsdn = fabs(geometry.height[k+dk] - geometry.height[k]   ) * zmu;
       
      /* --- dchi/ds at downwind point --              -------------- */
       
       if (abs(k - k_end) > 1) {
	 dsdn2=fabs(geometry.height[k+2*dk] -
		    geometry.height[k+dk]) * zmu;
	 dchi_dn = cent_deriv(dsdn,dsdn2,chi[k],chi[k+dk],chi[k+2*dk]);       
       } else {
	 dchi_dn=(chi[k+dk]-chi[k])/dsdn;
       }
       
       /* --- Make sure that c1 and c2 don't go below zero -- ------- */
    
       c1 = max(chi[k]    + (dsdn/3.0) * dchi_c,  0.0);
       c2 = max(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);

       /* --- Downwind optical path length --          -------------- */

       dtau_dw =  dsdn * (chi[k] + chi[k+dk] + c1 + c2) * 0.25;
       dt=dtau_uw, dt03 = dt/3.0;
      
      /* --- Compute interpolation parameters --       -------------- */
       
       Bezier3_coeffs(dt, &alpha, &beta, &gamma, &theta, &eps);
       
       /* --- dS/dt at central point --                -------------- */
       
       dS_c = cent_deriv(dtau_uw, dtau_dw, S[k-dk], S[k], S[k+dk]);

       /* --- Source function control points --        -------------- */
       
       c1 = max(S[k]    - dt03 * dS_c , 0.0);
       c2 = max(S[k-dk] + dt03 * dS_up, 0.0);       
     
       /* --- Solve integral in this interval --       -------------- */
       
       I[k]= I_upw*eps + alpha*S[k] + beta*S[k-dk] +
	 gamma * c1 + theta * c2; 

       /* --- Diagonal operator --                     -------------- */

       if (Psi) Psi[k] = alpha + gamma;
       
    } else { 
      
      /* --- Piecewise linear integration at end of ray -- ---------- */
      
      dtau_uw = 0.5 * zmu * (chi[k] + chi[k-dk]) *
	fabs(geometry.height[k] - geometry.height[k-dk]);
      
      /* --- Defined negative in Han's implementation -- ------------ */
      
      dS_uw = -(S[k] - S[k-dk]) / dtau_uw;
      w3(dtau_uw, w);
      
      I[k] = (1.0 - w[0])*I_upw + w[0]*S[k] + w[1]*dS_uw;

      /* --- Diagonal operator --                      -------------- */
      
      if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
    }
    
    /* --- Re-use downwind quantities for next upwind position -- --- */
    
    I_upw = I[k];
    dsup=dsdn;
    dchi_up=dchi_c;
    dchi_c=dchi_dn;
    dtau_uw=dtau_dw;
    dS_up = dS_c;
  }
}
/* ------- end ---------------------------- Piecewise_Bezier3.c ----- */


/* ------- begin -------------------------- m4inv.c ----------------- */

void m4inv(double MI[4][4])
{

  /* --- In-place Shipley-Coleman matrix inversion
         Fast, but ... how accurate??
         Pivoting is always done in the diagonal.
         Copied here just in case the SIMD matrix inversion 
         gives troubles. --                            -------------- */
  
  register int k, i, j;
  
  for (k = 0;  k < 4;  k++){

    /* --- The pivot element --                        -------------- */
    
    MI[k][k] = -1.0 / MI[k][k];

    /* --- The pivot column --                         -------------- */
    
    for(i = 0;  i < 4;  ++i) if(i != k) MI[i][k]*=MI[k][k];

    /* --- Elements not in a pivot row or column --    -------------- */
    
    for(i = 0;  i < 4;  ++i) {
      if(i != k)
	for(j = 0;  j < 4;  ++j)
	  if(j != k)
	    MI[i][j] += MI[i][k] * MI[k][j];
    }
    /* --- Elements in a pivot row --                  -------------- */
    
    for(i = 0;  i < 4;  ++i) {
      if(i != k)
	MI[k][i] *= MI[k][k];
    }
  }
  
  for(i = 0;  i < 4;  ++i) {
    for(j = 0;  j < 4;  ++j) MI[i][j] = -MI[i][j];
  }
  return;
}
/* ------- end ---------------------------- m4inv.c ----------------- */
