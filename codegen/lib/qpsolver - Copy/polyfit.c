/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: polyfit.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "polyfit.h"
#include "qpsolver_emxutil.h"
#include "xgeqp3.h"

/* Function Definitions */

/*
 * Arguments    : const double x[9]
 *                const double y[9]
 *                double n
 *                emxArray_real_T *p
 *                emxArray_real_T *S_R
 *                double *S_df
 *                double *S_normr
 * Return Type  : void
 */
void polyfit(const double x[9], const double y[9], double n, emxArray_real_T *p,
             emxArray_real_T *S_R, double *S_df, double *S_normr)
{
  emxArray_real_T *V;
  int ic;
  int mn;
  int b_mn;
  emxArray_real_T *A;
  int ia;
  double wj;
  int ar;
  emxArray_real_T *p1;
  emxArray_int32_T *jpvt;
  double tau_data[9];
  int tau_size[1];
  double r[9];
  emxArray_real_T *R;
  double scale;
  double absxk;
  double t;
  emxInit_real_T(&V, 2);
  ic = V->size[0] * V->size[1];
  V->size[0] = 9;
  V->size[1] = (int)(n + 1.0);
  emxEnsureCapacity((emxArray__common *)V, ic, (int)sizeof(double));
  if (V->size[1] == 0) {
  } else {
    for (mn = 0; mn < 9; mn++) {
      V->data[mn + V->size[0] * ((int)(n + 1.0) - 1)] = 1.0;
    }

    if (n < 1.0) {
    } else {
      for (mn = 0; mn < 9; mn++) {
        V->data[mn + V->size[0] * ((int)n - 1)] = x[mn];
      }

      for (ia = 0; ia < (int)-(1.0 + (-1.0 - (n - 1.0))); ia++) {
        wj = (n - 1.0) + -(double)ia;
        for (mn = 0; mn < 9; mn++) {
          V->data[mn + V->size[0] * ((int)wj - 1)] = x[mn] * V->data[mn +
            V->size[0] * ((int)(wj + 1.0) - 1)];
        }
      }
    }
  }

  b_mn = V->size[1];
  if (9 <= b_mn) {
    b_mn = 9;
  }

  emxInit_real_T(&A, 2);
  ic = A->size[0] * A->size[1];
  A->size[0] = 9;
  A->size[1] = V->size[1];
  emxEnsureCapacity((emxArray__common *)A, ic, (int)sizeof(double));
  ar = V->size[0] * V->size[1];
  for (ic = 0; ic < ar; ic++) {
    A->data[ic] = V->data[ic];
  }

  emxInit_real_T1(&p1, 1);
  emxInit_int32_T(&jpvt, 2);
  xgeqp3(A, tau_data, tau_size, jpvt);
  mn = A->size[1];
  ic = p1->size[0];
  p1->size[0] = mn;
  emxEnsureCapacity((emxArray__common *)p1, ic, (int)sizeof(double));
  for (ic = 0; ic < mn; ic++) {
    p1->data[ic] = 0.0;
  }

  memcpy(&r[0], &y[0], 9U * sizeof(double));
  mn = A->size[1];
  if (9 <= mn) {
    mn = 9;
  }

  for (ia = 0; ia + 1 <= mn; ia++) {
    if (tau_data[ia] != 0.0) {
      wj = r[ia];
      for (ar = ia + 1; ar + 1 < 10; ar++) {
        wj += A->data[ar + A->size[0] * ia] * r[ar];
      }

      wj *= tau_data[ia];
      if (wj != 0.0) {
        r[ia] -= wj;
        for (ar = ia + 1; ar + 1 < 10; ar++) {
          r[ar] -= A->data[ar + A->size[0] * ia] * wj;
        }
      }
    }
  }

  for (ar = 0; ar + 1 <= b_mn; ar++) {
    p1->data[jpvt->data[ar] - 1] = r[ar];
  }

  for (ia = b_mn - 1; ia + 1 > 0; ia--) {
    p1->data[jpvt->data[ia] - 1] /= A->data[ia + A->size[0] * ia];
    for (ar = 0; ar + 1 <= ia; ar++) {
      p1->data[jpvt->data[ar] - 1] -= p1->data[jpvt->data[ia] - 1] * A->data[ar
        + A->size[0] * ia];
    }
  }

  emxFree_int32_T(&jpvt);
  emxInit_real_T(&R, 2);
  mn = V->size[1];
  ic = R->size[0] * R->size[1];
  R->size[0] = b_mn;
  R->size[1] = mn;
  emxEnsureCapacity((emxArray__common *)R, ic, (int)sizeof(double));
  for (ia = 0; ia + 1 <= V->size[1]; ia++) {
    if (ia + 1 <= b_mn) {
      ic = ia + 1;
    } else {
      ic = b_mn;
    }

    for (ar = 0; ar + 1 <= ic; ar++) {
      R->data[ar + R->size[0] * ia] = A->data[ar + A->size[0] * ia];
    }

    for (ar = ia + 1; ar + 1 <= b_mn; ar++) {
      R->data[ar + R->size[0] * ia] = 0.0;
    }
  }

  emxFree_real_T(&A);
  if ((V->size[1] == 1) || (p1->size[0] == 1)) {
    for (ic = 0; ic < 9; ic++) {
      r[ic] = 0.0;
      ar = V->size[1];
      for (mn = 0; mn < ar; mn++) {
        wj = r[ic] + V->data[ic + V->size[0] * mn] * p1->data[mn];
        r[ic] = wj;
      }
    }
  } else {
    memset(&r[0], 0, 9U * sizeof(double));
    ar = -1;
    for (mn = 0; mn + 1 <= V->size[1]; mn++) {
      if (p1->data[mn] != 0.0) {
        ia = ar;
        for (ic = 0; ic < 9; ic++) {
          ia++;
          wj = r[ic] + p1->data[mn] * V->data[ia];
          r[ic] = wj;
        }
      }

      ar += 9;
    }
  }

  emxFree_real_T(&V);
  for (ic = 0; ic < 9; ic++) {
    r[ic] = y[ic] - r[ic];
  }

  ic = S_R->size[0] * S_R->size[1];
  S_R->size[0] = R->size[0];
  S_R->size[1] = R->size[1];
  emxEnsureCapacity((emxArray__common *)S_R, ic, (int)sizeof(double));
  ar = R->size[0] * R->size[1];
  for (ic = 0; ic < ar; ic++) {
    S_R->data[ic] = R->data[ic];
  }

  emxFree_real_T(&R);
  if ((0.0 >= 9.0 - (n + 1.0)) || rtIsNaN(9.0 - (n + 1.0))) {
    *S_df = 0.0;
  } else {
    *S_df = 9.0 - (n + 1.0);
  }

  wj = 0.0;
  scale = 2.2250738585072014E-308;
  for (mn = 0; mn < 9; mn++) {
    absxk = fabs(r[mn]);
    if (absxk > scale) {
      t = scale / absxk;
      wj = 1.0 + wj * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      wj += t * t;
    }
  }

  wj = scale * sqrt(wj);
  *S_normr = wj;
  ic = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = p1->size[0];
  emxEnsureCapacity((emxArray__common *)p, ic, (int)sizeof(double));
  ar = p1->size[0];
  for (ic = 0; ic < ar; ic++) {
    p->data[p->size[0] * ic] = p1->data[ic];
  }

  emxFree_real_T(&p1);
}

/*
 * File trailer for polyfit.c
 *
 * [EOF]
 */
