/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xscal.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "call_qpsolver.h"
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
  int i28;
  int k;
  i28 = (ix0 + n) - 1;
  for (k = ix0; k <= i28; k++) {
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
  int i17;
  int k;
  i17 = (ix0 + n) - 1;
  for (k = ix0; k <= i17; k++) {
    x->data[k - 1] *= a;
  }
}

/*
 * File trailer for xscal.c
 *
 * [EOF]
 */
