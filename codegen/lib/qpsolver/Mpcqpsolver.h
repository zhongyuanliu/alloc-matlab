/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Mpcqpsolver.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

#ifndef MPCQPSOLVER_H
#define MPCQPSOLVER_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "qpsolver_types.h"

/* Function Declarations */
extern void Mpcqpsolver(const emxArray_real_T *Linv, const emxArray_real_T *b_f,
  const emxArray_real_T *A, const emxArray_real_T *b, const emxArray_real_T
  *b_Aeq, const emxArray_real_T *b_beq, const emxArray_boolean_T *iA0,
  emxArray_real_T *solution, double *status, emxArray_boolean_T *iA);

#endif

/*
 * File trailer for Mpcqpsolver.h
 *
 * [EOF]
 */
