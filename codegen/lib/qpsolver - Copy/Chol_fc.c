/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Chol_fc.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "Chol_fc.h"
#include "xscal.h"
#include "qpsolver_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *H
 *                emxArray_real_T *L
 *                double *p
 * Return Type  : void
 */
void Chol_fc(const emxArray_real_T *H, emxArray_real_T *L, double *p)
{
  int n;
  int i3;
  int ia;
  int b_n;
  int info;
  int j;
  boolean_T exitg1;
  int jj;
  int jmax;
  double ajj;
  int ix;
  int iy;
  int nmj;
  emxArray_real_T *b_L;
  double c;
  int i4;
  n = H->size[1];
  i3 = L->size[0] * L->size[1];
  L->size[0] = H->size[0];
  L->size[1] = H->size[1];
  emxEnsureCapacity((emxArray__common *)L, i3, (int)sizeof(double));
  ia = H->size[0] * H->size[1];
  for (i3 = 0; i3 < ia; i3++) {
    L->data[i3] = H->data[i3];
  }

  b_n = H->size[0];
  info = -1;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j + 1 <= b_n)) {
    jj = j + j * n;
    ajj = 0.0;
    if (j < 1) {
    } else {
      ix = j;
      iy = j;
      for (nmj = 1; nmj <= j; nmj++) {
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
        nmj = (b_n - j) - 1;
        if ((nmj == 0) || (j == 0)) {
        } else {
          ix = j;
          i3 = (j + n * (j - 1)) + 2;
          for (jmax = j + 2; jmax <= i3; jmax += n) {
            c = -L->data[ix];
            iy = jj + 1;
            i4 = (jmax + nmj) - 1;
            for (ia = jmax; ia <= i4; ia++) {
              L->data[iy] += L->data[ia - 1] * c;
              iy++;
            }

            ix += n;
          }
        }

        b_xscal(nmj, 1.0 / ajj, L, jj + 2);
      }

      j++;
    } else {
      L->data[jj] = ajj;
      info = j;
      exitg1 = true;
    }
  }

  if (info + 1 == 0) {
    jmax = H->size[1];
  } else {
    jmax = info;
  }

  for (j = 1; j + 1 <= jmax; j++) {
    for (nmj = 1; nmj <= j; nmj++) {
      L->data[(nmj + L->size[0] * j) - 1] = 0.0;
    }
  }

  if (1 > jmax) {
    ia = 0;
  } else {
    ia = jmax;
  }

  if (1 > jmax) {
    nmj = 0;
  } else {
    nmj = jmax;
  }

  emxInit_real_T(&b_L, 2);
  i3 = b_L->size[0] * b_L->size[1];
  b_L->size[0] = ia;
  b_L->size[1] = nmj;
  emxEnsureCapacity((emxArray__common *)b_L, i3, (int)sizeof(double));
  for (i3 = 0; i3 < nmj; i3++) {
    for (i4 = 0; i4 < ia; i4++) {
      b_L->data[i4 + b_L->size[0] * i3] = L->data[i4 + L->size[0] * i3];
    }
  }

  i3 = L->size[0] * L->size[1];
  L->size[0] = b_L->size[0];
  L->size[1] = b_L->size[1];
  emxEnsureCapacity((emxArray__common *)L, i3, (int)sizeof(double));
  ia = b_L->size[1];
  for (i3 = 0; i3 < ia; i3++) {
    nmj = b_L->size[0];
    for (i4 = 0; i4 < nmj; i4++) {
      L->data[i4 + L->size[0] * i3] = b_L->data[i4 + b_L->size[0] * i3];
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
