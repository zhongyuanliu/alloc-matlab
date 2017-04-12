/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_Polyfif_api.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 16-Aug-2016 09:21:48
 */

#ifndef _CODER_POLYFIF_API_H
#define _CODER_POLYFIF_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_Polyfif_api.h"

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

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void Polyfif(real_T x[9], real_T y[9], real_T N_poly, real_T x_new,
                    emxArray_real_T *r, real_T *f);
extern void Polyfif_api(const mxArray *prhs[4], const mxArray *plhs[2]);
extern void Polyfif_atexit(void);
extern void Polyfif_initialize(void);
extern void Polyfif_terminate(void);
extern void Polyfif_xil_terminate(void);

#endif

/*
 * File trailer for _coder_Polyfif_api.h
 *
 * [EOF]
 */
