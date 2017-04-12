/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_call_qpsolver_api.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

#ifndef _CODER_CALL_QPSOLVER_API_H
#define _CODER_CALL_QPSOLVER_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_call_qpsolver_api.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  real_T Tx;
  real_T Ty;
  real_T Tm;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  real_T label;
  real_T enable;
  real_T TSP;
  real_T ASP;
  real_T Tx;
  real_T Ty;
  real_T T;
  real_T phi;
  real_T Tm;
} struct1_T;

#endif                                 /*typedef_struct1_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void call_qpsolver(struct0_T *T_r, real_T method, real_T N_max_thr,
  emxArray_real_T *solution, struct1_T *alloc_out, real_T *status);
extern void call_qpsolver_api(const mxArray * const prhs[3], const mxArray *
  plhs[3]);
extern void call_qpsolver_atexit(void);
extern void call_qpsolver_initialize(void);
extern void call_qpsolver_terminate(void);
extern void call_qpsolver_xil_terminate(void);

#endif

/*
 * File trailer for _coder_call_qpsolver_api.h
 *
 * [EOF]
 */
