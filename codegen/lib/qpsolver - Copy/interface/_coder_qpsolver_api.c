/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_qpsolver_api.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_qpsolver_api.h"
#include "_coder_qpsolver_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true, false, 131434U, NULL, "qpsolver", NULL,
  false, { 2045744189U, 2170104910U, 2743257031U, 4284093946U }, NULL };

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T y[8]);
static const mxArray *b_emlrt_marshallOut(const struct2_T u[8]);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static const mxArray *c_emlrt_marshallOut(const real_T u);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[2]);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[3]);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *thruster_data,
  const char_T *identifier, struct0_T y[8]);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static struct1_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *T_r,
  const char_T *identifier);
static struct1_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId);
static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *N_enabled_thruster, const char_T *identifier);
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *rudder_table0, const char_T *identifier))[45];
static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[45];
static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[2]);
static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[3]);
static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[45];

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                struct0_T y[8]
 * Return Type  : void
 */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T y[8])
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[23] = { "label", "type", "enable",
    "constr_angle", "dt", "x0", "y0", "x", "y", "T", "phi", "Tmax", "Tmin",
    "dTmax", "dphi_max", "Tplus", "T_", "phiPlus", "phi_", "phi_min", "phi_max",
    "weight", "weight_s" };

  static const int32_T dims[2] = { 1, 8 };

  int32_T i0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 23, fieldNames, 2U, dims);
  for (i0 = 0; i0 < 8; i0++) {
    thisId.fIdentifier = "label";
    y[i0].label = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "label")), &thisId);
    thisId.fIdentifier = "type";
    y[i0].type = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "type")), &thisId);
    thisId.fIdentifier = "enable";
    y[i0].enable = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "enable")), &thisId);
    thisId.fIdentifier = "constr_angle";
    y[i0].constr_angle = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a
      (sp, u, i0, "constr_angle")), &thisId);
    thisId.fIdentifier = "dt";
    y[i0].dt = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "dt")), &thisId);
    thisId.fIdentifier = "x0";
    y[i0].x0 = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "x0")), &thisId);
    thisId.fIdentifier = "y0";
    y[i0].y0 = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "y0")), &thisId);
    thisId.fIdentifier = "x";
    y[i0].x = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "x")), &thisId);
    thisId.fIdentifier = "y";
    y[i0].y = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "y")), &thisId);
    thisId.fIdentifier = "T";
    y[i0].T = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "T")), &thisId);
    thisId.fIdentifier = "phi";
    y[i0].phi = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "phi")), &thisId);
    thisId.fIdentifier = "Tmax";
    y[i0].Tmax = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "Tmax")), &thisId);
    thisId.fIdentifier = "Tmin";
    y[i0].Tmin = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "Tmin")), &thisId);
    thisId.fIdentifier = "dTmax";
    y[i0].dTmax = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "dTmax")), &thisId);
    thisId.fIdentifier = "dphi_max";
    y[i0].dphi_max = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "dphi_max")), &thisId);
    thisId.fIdentifier = "Tplus";
    y[i0].Tplus = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "Tplus")), &thisId);
    thisId.fIdentifier = "T_";
    y[i0].T_ = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "T_")), &thisId);
    thisId.fIdentifier = "phiPlus";
    y[i0].phiPlus = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "phiPlus")), &thisId);
    thisId.fIdentifier = "phi_";
    y[i0].phi_ = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "phi_")), &thisId);
    thisId.fIdentifier = "phi_min";
    y[i0].phi_min = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "phi_min")), &thisId);
    thisId.fIdentifier = "phi_max";
    y[i0].phi_max = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "phi_max")), &thisId);
    thisId.fIdentifier = "weight";
    d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0, "weight")),
                       &thisId, y[i0].weight);
    thisId.fIdentifier = "weight_s";
    e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0, "weight_s")),
                       &thisId, y[i0].weight_s);
  }

  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const struct2_T u[8]
 * Return Type  : const mxArray *
 */
static const mxArray *b_emlrt_marshallOut(const struct2_T u[8])
{
  const mxArray *y;
  int32_T iv1[2];
  int32_T i1;
  const struct2_T *r0;
  y = NULL;
  for (i1 = 0; i1 < 2; i1++) {
    iv1[i1] = 1 + 7 * i1;
  }

  emlrtAssign(&y, emlrtCreateStructArray(2, iv1, 0, NULL));
  for (i1 = 0; i1 < 8; i1++) {
    r0 = &u[i1];
    emlrtAddField(y, c_emlrt_marshallOut(r0->label), "label", i1);
    emlrtAddField(y, c_emlrt_marshallOut(r0->enable), "enable", i1);
    emlrtAddField(y, c_emlrt_marshallOut(r0->TSP), "TSP", i1);
    emlrtAddField(y, c_emlrt_marshallOut(r0->ASP), "ASP", i1);
    emlrtAddField(y, c_emlrt_marshallOut(r0->Tx), "Tx", i1);
    emlrtAddField(y, c_emlrt_marshallOut(r0->Ty), "Ty", i1);
    emlrtAddField(y, c_emlrt_marshallOut(r0->T), "T", i1);
    emlrtAddField(y, c_emlrt_marshallOut(r0->phi), "phi", i1);
    emlrtAddField(y, c_emlrt_marshallOut(r0->Tm), "Tm", i1);
  }

  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const real_T u
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m1;
  y = NULL;
  m1 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m1);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[2]
 * Return Type  : void
 */
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[2])
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[3]
 * Return Type  : void
 */
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[3])
{
  m_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *thruster_data
 *                const char_T *identifier
 *                struct0_T y[8]
 * Return Type  : void
 */
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *thruster_data,
  const char_T *identifier, struct0_T y[8])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(thruster_data), &thisId, y);
  emlrtDestroyArray(&thruster_data);
}

/*
 * Arguments    : const emxArray_real_T *u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[1] = { 0 };

  y = NULL;
  m0 = emlrtCreateNumericArray(1, iv0, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u->data);
  emlrtSetDimensions((mxArray *)m0, u->size, 1);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (real_T *)NULL) && (*pEmxArray)->canFreeData) {
      emlrtFreeMex((void *)(*pEmxArray)->data);
    }

    emlrtFreeMex((void *)(*pEmxArray)->size);
    emlrtFreeMex((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_real_T **pEmxArray
 *                int32_T numDimensions
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)emlrtMallocMex(sizeof(emxArray_real_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2012b(sp, (void *)pEmxArray, (void (*)(void *))
      emxFree_real_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex((uint32_T)(sizeof(int32_T)
    * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *T_r
 *                const char_T *identifier
 * Return Type  : struct1_T
 */
static struct1_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *T_r,
  const char_T *identifier)
{
  struct1_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = g_emlrt_marshallIn(sp, emlrtAlias(T_r), &thisId);
  emlrtDestroyArray(&T_r);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : struct1_T
 */
static struct1_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId)
{
  struct1_T y;
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[3] = { "Tx", "Ty", "Tm" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 3, fieldNames, 0U, &dims);
  thisId.fIdentifier = "Tx";
  y.Tx = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Tx")),
    &thisId);
  thisId.fIdentifier = "Ty";
  y.Ty = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Ty")),
    &thisId);
  thisId.fIdentifier = "Tm";
  y.Tm = c_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Tm")),
    &thisId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *N_enabled_thruster
 *                const char_T *identifier
 * Return Type  : real_T
 */
static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *N_enabled_thruster, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = c_emlrt_marshallIn(sp, emlrtAlias(N_enabled_thruster), &thisId);
  emlrtDestroyArray(&N_enabled_thruster);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *rudder_table0
 *                const char_T *identifier
 * Return Type  : real_T (*)[45]
 */
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *rudder_table0, const char_T *identifier))[45]
{
  real_T (*y)[45];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(rudder_table0), &thisId);
  emlrtDestroyArray(&rudder_table0);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[45]
 */
  static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[45]
{
  real_T (*y)[45];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[2]
 * Return Type  : void
 */
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[2])
{
  static const int32_T dims[2] = { 1, 2 };

  int32_T i2;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  for (i2 = 0; i2 < 2; i2++) {
    ret[i2] = (*(real_T (*)[2])mxGetData(src))[i2];
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[3]
 * Return Type  : void
 */
static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[3])
{
  static const int32_T dims[1] = { 3 };

  int32_T i3;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  for (i3 = 0; i3 < 3; i3++) {
    ret[i3] = (*(real_T (*)[3])mxGetData(src))[i3];
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[45]
 */
static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[45]
{
  real_T (*ret)[45];
  static const int32_T dims[2] = { 5, 9 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[45])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const mxArray *prhs[6]
 *                const mxArray *plhs[3]
 * Return Type  : void
 */
  void qpsolver_api(const mxArray *prhs[6], const mxArray *plhs[3])
{
  emxArray_real_T *solution;
  struct0_T thruster_data[8];
  struct1_T T_r;
  real_T N_enabled_thruster;
  real_T (*rudder_table0)[45];
  real_T no_azi_angle_constr;
  real_T method;
  struct2_T alloc_out[8];
  real_T status;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &solution, 1, true);
  prhs[3] = emlrtProtectR2012b(prhs[3], 3, false, -1);

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "thruster_data", thruster_data);
  T_r = f_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "T_r");
  N_enabled_thruster = h_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]),
    "N_enabled_thruster");
  rudder_table0 = i_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "rudder_table0");
  no_azi_angle_constr = h_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]),
    "no_azi_angle_constr");
  method = h_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "method");

  /* Invoke the target function */
  qpsolver(thruster_data, &T_r, N_enabled_thruster, *rudder_table0,
           no_azi_angle_constr, method, solution, alloc_out, &status);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(solution);
  plhs[1] = b_emlrt_marshallOut(alloc_out);
  plhs[2] = c_emlrt_marshallOut(status);
  solution->canFreeData = false;
  emxFree_real_T(&solution);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  qpsolver_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_qpsolver_api.c
 *
 * [EOF]
 */
