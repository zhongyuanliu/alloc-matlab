/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: call_qpsolver.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

#ifndef CALL_QPSOLVER_H
#define CALL_QPSOLVER_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "call_qpsolver_types.h"

/* Function Declarations */
extern void call_qpsolver(const struct0_T *T_r, double method, double N_max_thr,
  emxArray_real_T *solution, struct1_T *alloc_out, double *status);

#endif

/*
 * File trailer for call_qpsolver.h
 *
 * [EOF]
 */
