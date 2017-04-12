/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_qpsolver_api.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_qpsolver_api.h"
#include "_coder_qpsolver_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131435U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "qpsolver",                          /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static boolean_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId);
static const mxArray *b_emlrt_marshallOut(const struct2_T u[8]);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *thruster_data, const char_T *identifier, struct0_T y[8]);
static const mxArray *c_emlrt_marshallOut(const real_T u);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T y[8]);
static const mxArray *d_emlrt_marshallOut(const struct0_T u[8]);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static boolean_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *init,
  const char_T *identifier);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[2]);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[3]);
static struct1_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *T_r,
  const char_T *identifier);
static struct1_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId);
static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *N_enabled_thruster, const char_T *identifier);
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *rudder_table0, const char_T *identifier))[45];
static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[45];
static boolean_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static real_T n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[2]);
static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[3]);
static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[45];

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : boolean_T
 */
static boolean_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId)
{
  boolean_T y;
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
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
 *                const mxArray *thruster_data
 *                const char_T *identifier
 *                struct0_T y[8]
 * Return Type  : void
 */
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *thruster_data, const char_T *identifier, struct0_T y[8])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(thruster_data), &thisId, y);
  emlrtDestroyArray(&thruster_data);
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
 *                struct0_T y[8]
 * Return Type  : void
 */
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T y[8])
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[26] = { "label", "type", "enable", "constr",
    "dt", "x0", "y0", "x", "y", "T", "T_reserve", "phi", "phi_reserve", "Tmax",
    "Tmin", "base", "dTmax", "dphi_max", "Tplus", "T_", "phiPlus", "phi_",
    "phi_min", "phi_max", "weight", "weight_s" };

  static const int32_T dims[2] = { 1, 8 };

  int32_T i0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 26, fieldNames, 2U, dims);
  for (i0 = 0; i0 < 8; i0++) {
    thisId.fIdentifier = "label";
    y[i0].label = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "label")), &thisId);
    thisId.fIdentifier = "type";
    y[i0].type = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "type")), &thisId);
    thisId.fIdentifier = "enable";
    y[i0].enable = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "enable")), &thisId);
    thisId.fIdentifier = "constr";
    y[i0].constr = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "constr")), &thisId);
    thisId.fIdentifier = "dt";
    y[i0].dt = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "dt")), &thisId);
    thisId.fIdentifier = "x0";
    y[i0].x0 = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "x0")), &thisId);
    thisId.fIdentifier = "y0";
    y[i0].y0 = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "y0")), &thisId);
    thisId.fIdentifier = "x";
    y[i0].x = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "x")), &thisId);
    thisId.fIdentifier = "y";
    y[i0].y = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "y")), &thisId);
    thisId.fIdentifier = "T";
    y[i0].T = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "T")), &thisId);
    thisId.fIdentifier = "T_reserve";
    y[i0].T_reserve = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp,
      u, i0, "T_reserve")), &thisId);
    thisId.fIdentifier = "phi";
    y[i0].phi = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "phi")), &thisId);
    thisId.fIdentifier = "phi_reserve";
    y[i0].phi_reserve = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp,
      u, i0, "phi_reserve")), &thisId);
    thisId.fIdentifier = "Tmax";
    y[i0].Tmax = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "Tmax")), &thisId);
    thisId.fIdentifier = "Tmin";
    y[i0].Tmin = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "Tmin")), &thisId);
    thisId.fIdentifier = "base";
    f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0, "base")),
                       &thisId, y[i0].base);
    thisId.fIdentifier = "dTmax";
    y[i0].dTmax = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "dTmax")), &thisId);
    thisId.fIdentifier = "dphi_max";
    y[i0].dphi_max = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "dphi_max")), &thisId);
    thisId.fIdentifier = "Tplus";
    y[i0].Tplus = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "Tplus")), &thisId);
    thisId.fIdentifier = "T_";
    y[i0].T_ = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "T_")), &thisId);
    thisId.fIdentifier = "phiPlus";
    y[i0].phiPlus = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "phiPlus")), &thisId);
    thisId.fIdentifier = "phi_";
    y[i0].phi_ = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0,
      "phi_")), &thisId);
    thisId.fIdentifier = "phi_min";
    y[i0].phi_min = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "phi_min")), &thisId);
    thisId.fIdentifier = "phi_max";
    y[i0].phi_max = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u,
      i0, "phi_max")), &thisId);
    thisId.fIdentifier = "weight";
    f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0, "weight")),
                       &thisId, y[i0].weight);
    thisId.fIdentifier = "weight_s";
    g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, i0, "weight_s")),
                       &thisId, y[i0].weight_s);
  }

  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const struct0_T u[8]
 * Return Type  : const mxArray *
 */
static const mxArray *d_emlrt_marshallOut(const struct0_T u[8])
{
  const mxArray *y;
  int32_T iv2[2];
  int32_T i2;
  const struct0_T *r1;
  real_T b_u[2];
  int32_T i;
  const mxArray *b_y;
  const mxArray *m2;
  static const int32_T iv3[2] = { 1, 2 };

  real_T *pData;
  const mxArray *c_y;
  static const int32_T iv4[2] = { 1, 2 };

  real_T c_u[3];
  const mxArray *d_y;
  static const int32_T iv5[1] = { 3 };

  y = NULL;
  for (i2 = 0; i2 < 2; i2++) {
    iv2[i2] = 1 + 7 * i2;
  }

  emlrtAssign(&y, emlrtCreateStructArray(2, iv2, 0, NULL));
  for (i2 = 0; i2 < 8; i2++) {
    r1 = &u[i2];
    emlrtAddField(y, c_emlrt_marshallOut(r1->label), "label", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->type), "type", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->enable), "enable", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->constr), "constr", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->dt), "dt", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->x0), "x0", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->y0), "y0", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->x), "x", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->y), "y", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->T), "T", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->T_reserve), "T_reserve", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->phi), "phi", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->phi_reserve), "phi_reserve", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->Tmax), "Tmax", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->Tmin), "Tmin", i2);
    for (i = 0; i < 2; i++) {
      b_u[i] = r1->base[i];
    }

    b_y = NULL;
    m2 = emlrtCreateNumericArray(2, iv3, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T *)mxGetPr(m2);
    for (i = 0; i < 2; i++) {
      pData[i] = b_u[i];
    }

    emlrtAssign(&b_y, m2);
    emlrtAddField(y, b_y, "base", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->dTmax), "dTmax", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->dphi_max), "dphi_max", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->Tplus), "Tplus", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->T_), "T_", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->phiPlus), "phiPlus", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->phi_), "phi_", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->phi_min), "phi_min", i2);
    emlrtAddField(y, c_emlrt_marshallOut(r1->phi_max), "phi_max", i2);
    for (i = 0; i < 2; i++) {
      b_u[i] = r1->weight[i];
    }

    c_y = NULL;
    m2 = emlrtCreateNumericArray(2, iv4, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T *)mxGetPr(m2);
    for (i = 0; i < 2; i++) {
      pData[i] = b_u[i];
    }

    emlrtAssign(&c_y, m2);
    emlrtAddField(y, c_y, "weight", i2);
    for (i = 0; i < 3; i++) {
      c_u[i] = r1->weight_s[i];
    }

    d_y = NULL;
    m2 = emlrtCreateNumericArray(1, iv5, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T *)mxGetPr(m2);
    for (i = 0; i < 3; i++) {
      pData[i] = c_u[i];
    }

    emlrtAssign(&d_y, m2);
    emlrtAddField(y, d_y, "weight_s", i2);
  }

  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *init
 *                const char_T *identifier
 * Return Type  : boolean_T
 */
static boolean_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *init,
  const char_T *identifier)
{
  boolean_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(init), &thisId);
  emlrtDestroyArray(&init);
  return y;
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
  mxSetData((mxArray *)m0, (void *)&u->data[0]);
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
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[2]
 * Return Type  : void
 */
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[2])
{
  o_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[3]
 * Return Type  : void
 */
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[3])
{
  p_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *T_r
 *                const char_T *identifier
 * Return Type  : struct1_T
 */
static struct1_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *T_r,
  const char_T *identifier)
{
  struct1_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = i_emlrt_marshallIn(sp, emlrtAlias(T_r), &thisId);
  emlrtDestroyArray(&T_r);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : struct1_T
 */
static struct1_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
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
  y.Tx = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Tx")),
    &thisId);
  thisId.fIdentifier = "Ty";
  y.Ty = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Ty")),
    &thisId);
  thisId.fIdentifier = "Tm";
  y.Tm = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2013a(sp, u, 0, "Tm")),
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
static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *N_enabled_thruster, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = e_emlrt_marshallIn(sp, emlrtAlias(N_enabled_thruster), &thisId);
  emlrtDestroyArray(&N_enabled_thruster);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *rudder_table0
 *                const char_T *identifier
 * Return Type  : real_T (*)[45]
 */
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *rudder_table0, const char_T *identifier))[45]
{
  real_T (*y)[45];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = l_emlrt_marshallIn(sp, emlrtAlias(rudder_table0), &thisId);
  emlrtDestroyArray(&rudder_table0);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[45]
 */
  static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[45]
{
  real_T (*y)[45];
  y = q_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : boolean_T
 */
static boolean_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  boolean_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "logical", false, 0U, &dims);
  ret = *mxGetLogicals(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
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
static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[2])
{
  static const int32_T dims[2] = { 1, 2 };

  int32_T i3;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  for (i3 = 0; i3 < 2; i3++) {
    ret[i3] = (*(real_T (*)[2])mxGetData(src))[i3];
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
static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[3])
{
  static const int32_T dims[1] = { 3 };

  int32_T i4;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  for (i4 = 0; i4 < 3; i4++) {
    ret[i4] = (*(real_T (*)[3])mxGetData(src))[i4];
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[45]
 */
static real_T (*q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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
 * Arguments    : const mxArray *prhs[7]
 *                const mxArray *plhs[4]
 * Return Type  : void
 */
  void qpsolver_api(const mxArray *prhs[7], const mxArray *plhs[4])
{
  emxArray_real_T *solution;
  boolean_T init;
  struct0_T thruster_data[8];
  struct1_T T_r;
  real_T N_enabled_thruster;
  real_T (*rudder_table0)[45];
  real_T no_azi_angle_constr;
  real_T method;
  struct2_T alloc_out[8];
  real_T status;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &solution, 1, true);
  prhs[4] = emlrtProtectR2012b(prhs[4], 4, false, -1);

  /* Marshall function inputs */
  init = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "init");
  c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "thruster_data", thruster_data);
  T_r = h_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "T_r");
  N_enabled_thruster = j_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]),
    "N_enabled_thruster");
  rudder_table0 = k_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "rudder_table0");
  no_azi_angle_constr = j_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]),
    "no_azi_angle_constr");
  method = j_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "method");

  /* Invoke the target function */
  qpsolver(init, thruster_data, T_r, N_enabled_thruster, *rudder_table0,
           no_azi_angle_constr, method, solution, alloc_out, &status);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(solution);
  plhs[1] = b_emlrt_marshallOut(alloc_out);
  plhs[2] = d_emlrt_marshallOut(thruster_data);
  plhs[3] = c_emlrt_marshallOut(status);
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
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

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
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

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
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_qpsolver_api.c
 *
 * [EOF]
 */
