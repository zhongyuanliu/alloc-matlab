/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzgetrf.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "xzgetrf.h"
#include "colon.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                int n
 *                emxArray_real_T *A
 *                int lda
 *                emxArray_int32_T *ipiv
 *                int *info
 * Return Type  : void
 */
void xzgetrf(int m, int n, emxArray_real_T *A, int lda, emxArray_int32_T *ipiv,
             int *info)
{
  int b_m;
  int i19;
  int j;
  int mmj;
  int c;
  int iy;
  int ix;
  double smax;
  int jA;
  int i20;
  int jy;
  double s;
  int b_j;
  int ijA;
  if (m <= n) {
    b_m = m;
  } else {
    b_m = n;
  }

  eml_signed_integer_colon(b_m, ipiv);
  *info = 0;
  if ((m < 1) || (n < 1)) {
  } else {
    if (m - 1 <= n) {
      i19 = m - 1;
    } else {
      i19 = n;
    }

    for (j = 0; j + 1 <= i19; j++) {
      mmj = m - j;
      c = j * (lda + 1);
      if (mmj < 1) {
        iy = -1;
      } else {
        iy = 0;
        if (mmj > 1) {
          ix = c;
          smax = fabs(A->data[c]);
          for (jA = 1; jA + 1 <= mmj; jA++) {
            ix++;
            s = fabs(A->data[ix]);
            if (s > smax) {
              iy = jA;
              smax = s;
            }
          }
        }
      }

      if (A->data[c + iy] != 0.0) {
        if (iy != 0) {
          ipiv->data[j] = (j + iy) + 1;
          ix = j;
          iy += j;
          for (jA = 1; jA <= n; jA++) {
            smax = A->data[ix];
            A->data[ix] = A->data[iy];
            A->data[iy] = smax;
            ix += lda;
            iy += lda;
          }
        }

        i20 = c + mmj;
        for (iy = c + 1; iy + 1 <= i20; iy++) {
          A->data[iy] /= A->data[c];
        }
      } else {
        *info = j + 1;
      }

      iy = (n - j) - 1;
      jA = (c + lda) + 1;
      jy = c + lda;
      for (b_j = 1; b_j <= iy; b_j++) {
        smax = A->data[jy];
        if (A->data[jy] != 0.0) {
          ix = c + 1;
          i20 = mmj + jA;
          for (ijA = jA; ijA + 1 < i20; ijA++) {
            A->data[ijA] += A->data[ix] * -smax;
            ix++;
          }
        }

        jy += lda;
        jA += lda;
      }
    }

    if ((*info == 0) && (m <= n) && (!(A->data[(m + A->size[0] * (m - 1)) - 1]
          != 0.0))) {
      *info = m;
    }
  }
}

/*
 * File trailer for xzgetrf.c
 *
 * [EOF]
 */
