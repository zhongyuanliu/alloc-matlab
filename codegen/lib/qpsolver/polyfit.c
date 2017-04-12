/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: polyfit.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
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
  int i;
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
  for (i = 0; i < 9; i++) {
    V[18 + i] = 1.0;
    V[9 + i] = x[i];
    V[i] = x[i] * V[9 + i];
  }

  memcpy(&A[0], &V[0], 27U * sizeof(double));
  xgeqp3(A, tau, jpvt);
  for (i = 0; i < 3; i++) {
    p1[i] = 0.0;
  }

  memcpy(&r[0], &y[0], 9U * sizeof(double));
  for (j = 0; j < 3; j++) {
    if (tau[j] != 0.0) {
      absxk = r[j];
      for (i = j + 1; i + 1 < 10; i++) {
        absxk += A[i + 9 * j] * r[i];
      }

      absxk *= tau[j];
      if (absxk != 0.0) {
        r[j] -= absxk;
        for (i = j + 1; i + 1 < 10; i++) {
          r[i] -= A[i + 9 * j] * absxk;
        }
      }
    }
  }

  for (i = 0; i < 3; i++) {
    p1[jpvt[i] - 1] = r[i];
  }

  for (j = 2; j >= 0; j += -1) {
    p1[jpvt[j] - 1] /= A[j + 9 * j];
    for (i = 0; i + 1 <= j; i++) {
      p1[jpvt[i] - 1] -= p1[jpvt[j] - 1] * A[i + 9 * j];
    }
  }

  for (j = 0; j < 3; j++) {
    for (i = 0; i + 1 <= j + 1; i++) {
      R[i + 3 * j] = A[i + 9 * j];
    }

    for (i = j + 1; i + 1 < 4; i++) {
      R[i + 3 * j] = 0.0;
    }
  }

  *S_df = 6.0;
  *S_normr = 0.0;
  scale = 2.2250738585072014E-308;
  for (i = 0; i < 9; i++) {
    absxk = 0.0;
    for (j = 0; j < 3; j++) {
      absxk += V[i + 9 * j] * p1[j];
    }

    S_R[i] = R[i];
    absxk = fabs(y[i] - absxk);
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
