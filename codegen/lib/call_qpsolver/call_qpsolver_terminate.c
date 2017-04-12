/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: call_qpsolver_terminate.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "call_qpsolver.h"
#include "call_qpsolver_terminate.h"
#include "qpsolver.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void call_qpsolver_terminate(void)
{
  qpsolver_free();
}

/*
 * File trailer for call_qpsolver_terminate.c
 *
 * [EOF]
 */
