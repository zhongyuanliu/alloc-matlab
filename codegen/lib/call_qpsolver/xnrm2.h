/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xnrm2.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

#ifndef XNRM2_H
#define XNRM2_H

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
extern double b_xnrm2(int n, const emxArray_real_T *x, int ix0);
extern double c_xnrm2(int n, const emxArray_real_T *x, int ix0);
extern double xnrm2(int n, const double x[27], int ix0);

#endif

/*
 * File trailer for xnrm2.h
 *
 * [EOF]
 */
