/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpsolver.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

#ifndef QPSOLVER_H
#define QPSOLVER_H

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
extern void qpsolver(struct_T thruster_data_data[], double T_r_Tx, double T_r_Ty,
                     double T_r_Tm, double N_enabled_thruster, double method,
                     emxArray_real_T *solution, double *alloc_out_label, double *
                     alloc_out_enable, double *alloc_out_TSP, double
                     *alloc_out_ASP, double *alloc_out_Tx, double *alloc_out_Ty,
                     double *alloc_out_T, double *alloc_out_phi, double
                     *alloc_out_Tm, double *status);
extern void qpsolver_free(void);
extern void qpsolver_init(void);
extern void rudder_dat_init_not_empty_init(void);

#endif

/*
 * File trailer for qpsolver.h
 *
 * [EOF]
 */
