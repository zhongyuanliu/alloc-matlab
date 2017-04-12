/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_qpsolver_api.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

#ifndef _CODER_QPSOLVER_API_H
#define _CODER_QPSOLVER_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_qpsolver_api.h"

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
  real_T label;
  real_T type;
  real_T enable;
  real_T constr_angle;
  real_T dt;
  real_T x0;
  real_T y0;
  real_T x;
  real_T y;
  real_T T;
  real_T phi;
  real_T Tmax;
  real_T Tmin;
  real_T dTmax;
  real_T dphi_max;
  real_T Tplus;
  real_T T_;
  real_T phiPlus;
  real_T phi_;
  real_T phi_min;
  real_T phi_max;
  real_T weight[2];
  real_T weight_s[3];
} struct0_T;

#endif                                 /*typedef_struct0_T*/

#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  real_T Tx;
  real_T Ty;
  real_T Tm;
} struct1_T;

#endif                                 /*typedef_struct1_T*/

#ifndef typedef_struct2_T
#define typedef_struct2_T

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
} struct2_T;

#endif                                 /*typedef_struct2_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void qpsolver(struct0_T thruster_data[8], struct1_T *T_r, real_T
                     N_enabled_thruster, real_T rudder_table0[45], real_T
                     no_azi_angle_constr, real_T method, emxArray_real_T
                     *solution, struct2_T alloc_out[8], real_T *status);
extern void qpsolver_api(const mxArray *prhs[6], const mxArray *plhs[3]);
extern void qpsolver_atexit(void);
extern void qpsolver_initialize(void);
extern void qpsolver_terminate(void);
extern void qpsolver_xil_terminate(void);

#endif

/*
 * File trailer for _coder_qpsolver_api.h
 *
 * [EOF]
 */
