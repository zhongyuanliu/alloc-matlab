/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpsolver_initialize.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "qpsolver_initialize.h"
#include "pre_qp.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_initialize(void)
{
  rt_InitInfAndNaN(8U);
  N_Aeq_not_empty_init();
  b_fpp_not_empty_init();
  A_fpp_not_empty_init();
  NA_azi_not_empty_init();
  A_azi_not_empty_init();
  rudder_dat_init_not_empty_init();
  qpsolver_init();
  pre_qp_init();
}

/*
 * File trailer for qpsolver_initialize.c
 *
 * [EOF]
 */
