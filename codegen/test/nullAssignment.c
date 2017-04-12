/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: nullAssignment.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "nullAssignment.h"
#include "qpsolver_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *x
 *                int idx
 * Return Type  : void
 */
void b_nullAssignment(emxArray_real_T *x, int idx)
{
  int nxin;
  int nrowx;
  int ncolx;
  int k;
  emxArray_real_T *b_x;
  emxArray_real_T *c_x;
  nxin = x->size[0] * x->size[1] - 1;
  nrowx = x->size[0];
  ncolx = x->size[1];
  for (k = idx; k <= nxin; k++) {
    x->data[k - 1] = x->data[k];
  }

  emxInit_real_T(&b_x, 2);
  emxInit_real_T1(&c_x, 1);
  if ((nrowx != 1) && (ncolx == 1)) {
    if (1 > nxin) {
      nxin = 0;
    }

    nrowx = c_x->size[0];
    c_x->size[0] = nxin;
    emxEnsureCapacity((emxArray__common *)c_x, nrowx, (int)sizeof(double));
    for (nrowx = 0; nrowx < nxin; nrowx++) {
      c_x->data[nrowx] = x->data[nrowx];
    }

    nrowx = x->size[0] * x->size[1];
    x->size[0] = nxin;
    x->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)x, nrowx, (int)sizeof(double));
    for (nrowx = 0; nrowx < 1; nrowx++) {
      for (k = 0; k < nxin; k++) {
        x->data[k] = c_x->data[k];
      }
    }
  } else {
    if (1 > nxin) {
      nxin = 0;
    }

    nrowx = b_x->size[0] * b_x->size[1];
    b_x->size[0] = 1;
    b_x->size[1] = nxin;
    emxEnsureCapacity((emxArray__common *)b_x, nrowx, (int)sizeof(double));
    for (nrowx = 0; nrowx < nxin; nrowx++) {
      b_x->data[b_x->size[0] * nrowx] = x->data[nrowx];
    }

    nrowx = x->size[0] * x->size[1];
    x->size[0] = b_x->size[0];
    x->size[1] = b_x->size[1];
    emxEnsureCapacity((emxArray__common *)x, nrowx, (int)sizeof(double));
    nxin = b_x->size[1];
    for (nrowx = 0; nrowx < nxin; nrowx++) {
      ncolx = b_x->size[0];
      for (k = 0; k < ncolx; k++) {
        x->data[k + x->size[0] * nrowx] = b_x->data[k + b_x->size[0] * nrowx];
      }
    }
  }

  emxFree_real_T(&c_x);
  emxFree_real_T(&b_x);
}

/*
 * Arguments    : emxArray_real_T *x
 *                int idx
 * Return Type  : void
 */
void c_nullAssignment(emxArray_real_T *x, int idx)
{
  int nxin;
  int k;
  emxArray_real_T *b_x;
  nxin = x->size[0] - 1;
  for (k = idx; k <= nxin; k++) {
    x->data[k - 1] = x->data[k];
  }

  if (1 > nxin) {
    nxin = 0;
  }

  emxInit_real_T1(&b_x, 1);
  k = b_x->size[0];
  b_x->size[0] = nxin;
  emxEnsureCapacity((emxArray__common *)b_x, k, (int)sizeof(double));
  for (k = 0; k < nxin; k++) {
    b_x->data[k] = x->data[k];
  }

  k = x->size[0];
  x->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)x, k, (int)sizeof(double));
  nxin = b_x->size[0];
  for (k = 0; k < nxin; k++) {
    x->data[k] = b_x->data[k];
  }

  emxFree_real_T(&b_x);
}

/*
 * Arguments    : emxArray_real_T *x
 *                int idx
 * Return Type  : void
 */
void nullAssignment(emxArray_real_T *x, int idx)
{
  int nrowx;
  int ncolx;
  int j;
  int i;
  emxArray_real_T *b_x;
  nrowx = x->size[0] - 1;
  ncolx = x->size[1];
  for (j = 0; j + 1 <= ncolx; j++) {
    for (i = idx; i <= nrowx; i++) {
      x->data[(i + x->size[0] * j) - 1] = x->data[i + x->size[0] * j];
    }
  }

  if (1 > nrowx) {
    ncolx = 0;
  } else {
    ncolx = nrowx;
  }

  emxInit_real_T(&b_x, 2);
  nrowx = x->size[1];
  j = b_x->size[0] * b_x->size[1];
  b_x->size[0] = ncolx;
  b_x->size[1] = nrowx;
  emxEnsureCapacity((emxArray__common *)b_x, j, (int)sizeof(double));
  for (j = 0; j < nrowx; j++) {
    for (i = 0; i < ncolx; i++) {
      b_x->data[i + b_x->size[0] * j] = x->data[i + x->size[0] * j];
    }
  }

  j = x->size[0] * x->size[1];
  x->size[0] = b_x->size[0];
  x->size[1] = b_x->size[1];
  emxEnsureCapacity((emxArray__common *)x, j, (int)sizeof(double));
  ncolx = b_x->size[1];
  for (j = 0; j < ncolx; j++) {
    nrowx = b_x->size[0];
    for (i = 0; i < nrowx; i++) {
      x->data[i + x->size[0] * j] = b_x->data[i + b_x->size[0] * j];
    }
  }

  emxFree_real_T(&b_x);
}

/*
 * File trailer for nullAssignment.c
 *
 * [EOF]
 */
