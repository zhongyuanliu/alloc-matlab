/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: polyfit.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

#ifndef POLYFIT_H
#define POLYFIT_H

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
extern void polyfit(const double x[9], const double y[9], double p[3], double
                    S_R[9], double *S_df, double *S_normr);

#endif

/*
 * File trailer for polyfit.h
 *
 * [EOF]
 */
