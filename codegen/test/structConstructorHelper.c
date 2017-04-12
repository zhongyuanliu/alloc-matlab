/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: structConstructorHelper.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "structConstructorHelper.h"

/* Function Definitions */

/*
 * Arguments    : struct2_T s[8]
 * Return Type  : void
 */
void structConstructorHelper(struct2_T s[8])
{
  int j;
  for (j = 0; j < 8; j++) {
    s[j].label = 0.0;
    s[j].enable = 1.0;
    s[j].TSP = 0.0;
    s[j].ASP = 0.0;
    s[j].Tx = 0.0;
    s[j].Ty = 0.0;
    s[j].T = 0.0;
    s[j].phi = 0.0;
    s[j].Tm = 0.0;
  }
}

/*
 * File trailer for structConstructorHelper.c
 *
 * [EOF]
 */
