/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpsolver.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
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
#include "qpsolver_types.h"

/* Function Declarations */
extern void qpsolver(struct0_T thruster_data[8], const struct1_T *T_r, double
                     N_enabled_thruster, const double rudder_table0[45], double
                     no_azi_angle_constr, double method, emxArray_real_T
                     *solution, struct2_T alloc_out[8], double *status);
extern void qpsolver_free(void);
extern void qpsolver_init(void);
extern void rudder_dat_init_not_empty_init(void);

#endif

/*
 * File trailer for qpsolver.h
 *
 * [EOF]
 */
