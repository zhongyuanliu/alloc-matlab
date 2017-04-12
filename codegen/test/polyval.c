/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: polyval.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "polyval.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *p
 *                double x
 * Return Type  : double
 */
double polyval(const emxArray_real_T *p, double x)
{
  double y;
  boolean_T b0;
  int k;
  b0 = (p->size[1] == 0);
  if (!b0) {
    y = p->data[0];
    for (k = 0; k <= p->size[1] - 2; k++) {
      y = x * y + p->data[k + 1];
    }
  }

  return y;
}

/*
 * File trailer for polyval.c
 *
 * [EOF]
 */
