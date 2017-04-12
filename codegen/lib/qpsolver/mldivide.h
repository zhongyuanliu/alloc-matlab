/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mldivide.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

#ifndef MLDIVIDE_H
#define MLDIVIDE_H

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
extern void b_mldivide(const emxArray_real_T *A, const emxArray_real_T *B,
  emxArray_real_T *Y);
extern void mldivide(const emxArray_real_T *A, const double B[3],
                     emxArray_real_T *Y);

#endif

/*
 * File trailer for mldivide.h
 *
 * [EOF]
 */
