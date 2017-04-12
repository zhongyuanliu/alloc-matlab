/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pre_qp.h
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
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
extern void A_azi_not_empty_init(void);
extern void A_fpp_not_empty_init(void);
extern void NA_azi_not_empty_init(void);
extern void N_Aeq_not_empty_init(void);
extern void b_fpp_not_empty_init(void);
extern void pre_qp(boolean_T init, struct0_T thruster_data[8], double T_r_Tx,
                   double T_r_Ty, double T_r_Tm, double N_enabled_thruster,
                   const struct_T Ce[8], const double rudder_table[45], const
                   double b_Rangle_T[3], const double b_Rangle_Fangle[3], double
                   *Nvar, emxArray_real_T *Ho, emxArray_real_T *fo,
                   emxArray_real_T *Ao, emxArray_real_T *bo, emxArray_real_T
                   *Aeqo, emxArray_real_T *beqo, emxArray_real_T *x0);
extern void pre_qp_free(void);
extern void pre_qp_init(void);

#endif

/*
 * File trailer for pre_qp.h
 *
 * [EOF]
 */
