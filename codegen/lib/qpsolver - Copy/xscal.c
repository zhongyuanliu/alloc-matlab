/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xscal.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "xscal.h"

/* Function Definitions */

/*
 * Arguments    : int n
 *                double a
 *                emxArray_real_T *x
 *                int ix0
 * Return Type  : void
 */
void b_xscal(int n, double a, emxArray_real_T *x, int ix0)
{
  int i18;
  int k;
  i18 = (ix0 + n) - 1;
  for (k = ix0; k <= i18; k++) {
    x->data[k - 1] *= a;
  }
}

/*
 * Arguments    : int n
 *                double a
 *                emxArray_real_T *x
 *                int ix0
 * Return Type  : void
 */
void c_xscal(int n, double a, emxArray_real_T *x, int ix0)
{
  int i29;
  int k;
  i29 = (ix0 + n) - 1;
  for (k = ix0; k <= i29; k++) {
    x->data[k - 1] *= a;
  }
}

/*
 * Arguments    : int n
 *                double a
 *                emxArray_real_T *x
 *                int ix0
 * Return Type  : void
 */
void xscal(int n, double a, emxArray_real_T *x, int ix0)
{
  int i16;
  int k;
  i16 = (ix0 + n) - 1;
  for (k = ix0; k <= i16; k++) {
    x->data[k - 1] *= a;
  }
}

/*
 * File trailer for xscal.c
 *
 * [EOF]
 */
