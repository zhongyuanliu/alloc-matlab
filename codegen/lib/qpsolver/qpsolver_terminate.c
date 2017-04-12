/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpsolver_terminate.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "qpsolver_terminate.h"
#include "pre_qp.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_terminate(void)
{
  pre_qp_free();
}

/*
 * File trailer for qpsolver_terminate.c
 *
 * [EOF]
 */
