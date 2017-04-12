/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpsolver_initialize.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "qpsolver_initialize.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_initialize(void)
{
  rt_InitInfAndNaN(8U);
  rudder_dat_init_not_empty_init();
  qpsolver_init();
}

/*
 * File trailer for qpsolver_initialize.c
 *
 * [EOF]
 */
