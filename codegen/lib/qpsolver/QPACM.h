/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: QPACM.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

#ifndef QPACM_H
#define QPACM_H

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
extern double QPACM(const emxArray_real_T *B, const emxArray_real_T *b, const
                    emxArray_real_T *C, const emxArray_real_T *c, const
                    emxArray_real_T *A, const emxArray_real_T *d,
                    emxArray_real_T *x0);

#endif

/*
 * File trailer for QPACM.h
 *
 * [EOF]
 */
