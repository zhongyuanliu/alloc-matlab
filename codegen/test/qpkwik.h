/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpkwik.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

#ifndef QPKWIK_H
#define QPKWIK_H

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
extern void qpkwik(const emxArray_real_T *Linv, const emxArray_real_T *Hinv,
                   const emxArray_real_T *f, const emxArray_real_T *Ac, const
                   emxArray_real_T *b, emxArray_int16_T *iA, short m, short n,
                   short meq, emxArray_real_T *x, emxArray_real_T *lambda,
                   double *status);

#endif

/*
 * File trailer for qpkwik.h
 *
 * [EOF]
 */
