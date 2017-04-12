/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: polyfit.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "call_qpsolver.h"
#include "polyfit.h"
#include "xgeqp3.h"

/* Function Definitions */

/*
 * Arguments    : const double x[9]
 *                const double y[9]
 *                double p[3]
 *                double S_R[9]
 *                double *S_df
 *                double *S_normr
 * Return Type  : void
 */
void polyfit(const double x[9], const double y[9], double p[3], double S_R[9],
             double *S_df, double *S_normr)
{
  double V[27];
  int k;
  double A[27];
  double tau[3];
  int jpvt[3];
  double p1[3];
  double r[9];
  int j;
  double absxk;
  double R[9];
  double scale;
  double t;
  for (k = 0; k < 9; k++) {
    V[18 + k] = 1.0;
    V[9 + k] = x[k];
    V[k] = x[k] * V[9 + k];
  }

  memcpy(&A[0], &V[0], 27U * sizeof(double));
  xgeqp3(A, tau, jpvt);
  for (k = 0; k < 3; k++) {
    p1[k] = 0.0;
  }

  memcpy(&r[0], &y[0], 9U * sizeof(double));
  for (j = 0; j < 3; j++) {
    if (tau[j] != 0.0) {
      absxk = r[j];
      for (k = j + 1; k + 1 < 10; k++) {
        absxk += A[k + 9 * j] * r[k];
      }

      absxk *= tau[j];
      if (absxk != 0.0) {
        r[j] -= absxk;
        for (k = j + 1; k + 1 < 10; k++) {
          r[k] -= A[k + 9 * j] * absxk;
        }
      }
    }
  }

  for (k = 0; k < 3; k++) {
    p1[jpvt[k] - 1] = r[k];
  }

  for (j = 2; j >= 0; j += -1) {
    p1[jpvt[j] - 1] /= A[j + 9 * j];
    for (k = 0; k + 1 <= j; k++) {
      p1[jpvt[k] - 1] -= p1[jpvt[j] - 1] * A[k + 9 * j];
    }
  }

  for (j = 0; j < 3; j++) {
    for (k = 0; k + 1 <= j + 1; k++) {
      R[k + 3 * j] = A[k + 9 * j];
    }

    for (k = j + 1; k + 1 < 4; k++) {
      R[k + 3 * j] = 0.0;
    }
  }

  *S_df = 6.0;
  *S_normr = 0.0;
  scale = 2.2250738585072014E-308;
  for (k = 0; k < 9; k++) {
    absxk = 0.0;
    for (j = 0; j < 3; j++) {
      absxk += V[k + 9 * j] * p1[j];
    }

    r[k] = y[k] - absxk;
    S_R[k] = R[k];
    absxk = fabs(r[k]);
    if (absxk > scale) {
      t = scale / absxk;
      *S_normr = 1.0 + *S_normr * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      *S_normr += t * t;
    }
  }

  *S_normr = scale * sqrt(*S_normr);
  for (j = 0; j < 3; j++) {
    p[j] = p1[j];
  }
}

/*
 * File trailer for polyfit.c
 *
 * [EOF]
 */
