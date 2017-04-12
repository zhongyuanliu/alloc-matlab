/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: structConstructorHelper.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "structConstructorHelper.h"

/* Function Declarations */
static void c_structConstructorHelper(const double varargin_2[8], const
  cell_wrap_1 varargin_6[8], struct_T s[8]);

/* Function Definitions */

/*
 * Arguments    : const double varargin_2[8]
 *                const cell_wrap_1 varargin_6[8]
 *                struct_T s[8]
 * Return Type  : void
 */
static void c_structConstructorHelper(const double varargin_2[8], const
  cell_wrap_1 varargin_6[8], struct_T s[8])
{
  int j;
  for (j = 0; j < 8; j++) {
    s[j].label = varargin_2[j];
  }

  for (j = 0; j < 8; j++) {
    s[j].reserve = 0.0;
  }

  for (j = 0; j < 8; j++) {
    memcpy(&s[j].value[0], &varargin_6[j].f1[0], 74U * sizeof(double));
  }
}

/*
 * Arguments    : const double varargin_2[8]
 *                const cell_wrap_1 varargin_6[8]
 *                struct_T s[8]
 * Return Type  : void
 */
void b_structConstructorHelper(const double varargin_2[8], const cell_wrap_1
  varargin_6[8], struct_T s[8])
{
  c_structConstructorHelper(varargin_2, varargin_6, s);
}

/*
 * Arguments    : struct2_T s[8]
 * Return Type  : void
 */
void structConstructorHelper(struct2_T s[8])
{
  memset(&s[0], 0, sizeof(struct2_T) << 3);
}

/*
 * File trailer for structConstructorHelper.c
 *
 * [EOF]
 */
