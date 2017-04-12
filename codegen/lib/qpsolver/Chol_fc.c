/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Chol_fc.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "Chol_fc.h"
#include "xscal.h"
#include "qpsolver_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *b_H
 *                emxArray_real_T *L
 *                double *p
 * Return Type  : void
 */
void Chol_fc(const emxArray_real_T *b_H, emxArray_real_T *L, double *p)
{
  int n;
  int i0;
  int loop_ub;
  int b_n;
  int info;
  int j;
  boolean_T exitg1;
  int jj;
  int jmax;
  double ajj;
  int ix;
  int iy;
  int b_loop_ub;
  emxArray_real_T *b_L;
  double c;
  int i1;
  n = b_H->size[1];
  i0 = L->size[0] * L->size[1];
  L->size[0] = b_H->size[0];
  L->size[1] = b_H->size[1];
  emxEnsureCapacity((emxArray__common *)L, i0, (int)sizeof(double));
  loop_ub = b_H->size[0] * b_H->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    L->data[i0] = b_H->data[i0];
  }

  b_n = b_H->size[0];
  info = -1;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j + 1 <= b_n)) {
    jj = j + j * n;
    ajj = 0.0;
    if (!(j < 1)) {
      ix = j;
      iy = j;
      for (jmax = 1; jmax <= j; jmax++) {
        ajj += L->data[ix] * L->data[iy];
        ix += n;
        iy += n;
      }
    }

    ajj = L->data[jj] - ajj;
    if (ajj > 0.0) {
      ajj = sqrt(ajj);
      L->data[jj] = ajj;
      if (j + 1 < b_n) {
        jmax = (b_n - j) - 1;
        if ((jmax == 0) || (j == 0)) {
        } else {
          ix = j;
          i0 = (j + n * (j - 1)) + 2;
          for (b_loop_ub = j + 2; b_loop_ub <= i0; b_loop_ub += n) {
            c = -L->data[ix];
            iy = jj + 1;
            i1 = (b_loop_ub + jmax) - 1;
            for (loop_ub = b_loop_ub; loop_ub <= i1; loop_ub++) {
              L->data[iy] += L->data[loop_ub - 1] * c;
              iy++;
            }

            ix += n;
          }
        }

        xscal(jmax, 1.0 / ajj, L, jj + 2);
      }

      j++;
    } else {
      L->data[jj] = ajj;
      info = j;
      exitg1 = true;
    }
  }

  if (info + 1 == 0) {
    jmax = b_H->size[1];
  } else {
    jmax = info;
  }

  for (j = 1; j + 1 <= jmax; j++) {
    for (b_loop_ub = 1; b_loop_ub <= j; b_loop_ub++) {
      L->data[(b_loop_ub + L->size[0] * j) - 1] = 0.0;
    }
  }

  if (1 > jmax) {
    loop_ub = 0;
    b_loop_ub = 0;
  } else {
    loop_ub = jmax;
    b_loop_ub = jmax;
  }

  emxInit_real_T1(&b_L, 2);
  i0 = b_L->size[0] * b_L->size[1];
  b_L->size[0] = loop_ub;
  b_L->size[1] = b_loop_ub;
  emxEnsureCapacity((emxArray__common *)b_L, i0, (int)sizeof(double));
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_L->data[i1 + b_L->size[0] * i0] = L->data[i1 + L->size[0] * i0];
    }
  }

  i0 = L->size[0] * L->size[1];
  L->size[0] = b_L->size[0];
  L->size[1] = b_L->size[1];
  emxEnsureCapacity((emxArray__common *)L, i0, (int)sizeof(double));
  loop_ub = b_L->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = b_L->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      L->data[i1 + L->size[0] * i0] = b_L->data[i1 + b_L->size[0] * i0];
    }
  }

  emxFree_real_T(&b_L);
  *p = (double)info + 1.0;
}

/*
 * File trailer for Chol_fc.c
 *
 * [EOF]
 */
