/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgeqp3.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

#ifndef XGEQP3_H
#define XGEQP3_H

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
extern void b_xgeqp3(emxArray_real_T *A, double tau_data[], int tau_size[1],
                     emxArray_int32_T *jpvt);
extern void c_xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, emxArray_int32_T *
                     jpvt);
extern void xgeqp3(emxArray_real_T *A, double tau_data[], int tau_size[1],
                   emxArray_int32_T *jpvt);

#endif

/*
 * File trailer for xgeqp3.h
 *
 * [EOF]
 */
