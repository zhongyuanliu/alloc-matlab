/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: polyval.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "call_qpsolver.h"
#include "polyval.h"

/* Function Definitions */

/*
 * Arguments    : const double p[3]
 *                double x
 * Return Type  : double
 */
double polyval(const double p[3], double x)
{
  double y;
  int k;
  y = p[0];
  for (k = 0; k < 2; k++) {
    y = x * y + p[k + 1];
  }

  return y;
}

/*
 * File trailer for polyval.c
 *
 * [EOF]
 */
