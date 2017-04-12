/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pre_qp.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

#ifndef PRE_QP_H
#define PRE_QP_H

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
extern void pre_qp(struct0_T thruster_data[8], double T_r_Tx, double T_r_Ty,
                   double T_r_Tm, double N_enabled_thruster, const double
                   rudder_table[45], const emxArray_real_T *b_Rangle_T, const
                   emxArray_real_T *b_Rangle_Fangle, double no_azi_angle_constr,
                   double *Nvar, emxArray_real_T *H, emxArray_real_T *f,
                   emxArray_real_T *A, emxArray_real_T *b, emxArray_real_T *Aeq,
                   emxArray_real_T *beq, emxArray_real_T *x0);

#endif

/*
 * File trailer for pre_qp.h
 *
 * [EOF]
 */
