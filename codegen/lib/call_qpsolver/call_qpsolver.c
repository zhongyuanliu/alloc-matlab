/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: call_qpsolver.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "call_qpsolver.h"
#include "qpsolver.h"

/* Function Definitions */

/*
 * Arguments    : const struct0_T *T_r
 *                double method
 *                double N_max_thr
 *                emxArray_real_T *solution
 *                struct1_T *alloc_out
 *                double *status
 * Return Type  : void
 */
void call_qpsolver(const struct0_T *T_r, double method, double N_max_thr,
                   emxArray_real_T *solution, struct1_T *alloc_out, double
                   *status)
{
  struct_T thruster_data_in_data[1];
  static const struct_T r0 = { 1.0, 4.0, 1.0, 0.0, 0.5, 0.0, 0.0, 10.0, 4.0,
    5435.0, 1.5707963267948966, 71500.0, 0.0, 7150.0, 0.1, 0.0, 0.0, 0.0, 0.0,
    -1.04, -2.09, { 1.0, 1.0 }, { 1.0E+6, 1.0E+6, 1.0E+8 } };

  double N_enabled_thruster;
  int i;
  static const struct_T rv0[8] = { { 1.0, 4.0, 1.0, 0.0, 0.5, 0.0, 0.0, 10.0,
      4.0, 5435.0, 1.5707963267948966, 71500.0, 0.0, 7150.0, 0.1, 0.0, 0.0, 0.0,
      0.0, -1.04, -2.09, { 1.0, 1.0 }, { 1.0E+6, 1.0E+6, 1.0E+8 } }, { 2.0, 3.0,
      1.0, 0.0, 0.5, 0.0, 0.0, -30.0, -3.5, 3510.0, 0.0, 71500.0, 0.0, 7150.0,
      0.1, 0.0, 0.0, 0.0, 0.0, -1.04, -2.09, { 1.0, 1.0 }, { 1.0E+6, 1.0E+6,
        1.0E+8 } }, { 3.0, 2.0, 1.0, 0.0, 0.5, 0.0, 0.0, -30.0, 3.5, 3510.0, 0.0,
      71500.0, 0.0, 7150.0, 0.1, 0.0, 0.0, 0.0, 0.0, 2.09, 1.04, { 1.0, 1.0 }, {
        1.0E+6, 1.0E+6, 1.0E+8 } }, { 4.0, 2.0, 1.0, 0.0, 0.5, 0.0, 0.0, 21.0,
      0.0, 70000.0, 0.0, 71500.0, 0.0, 7150.0, 0.1, 0.0, 0.0, 0.0, 0.0, -1.04,
      -2.09, { 1.0, 1.0 }, { 1.0E+6, 1.0E+6, 1.0E+8 } }, { 5.0, 2.0, 1.0, 0.0,
      0.5, 0.0, 0.0, 23.0, 0.0, 3211.0, 0.0, 71500.0, 0.0, 7150.0, 0.1, 0.0, 0.0,
      0.0, 0.0, -1.0, 2.0, { 1.0, 1.0 }, { 1.0E+6, 1.0E+6, 1.0E+8 } }, { 6.0,
      3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 18.0, 0.0, 231.0, 0.0, 71500.0, 0.0, 7150.0,
      0.1, 0.0, 0.0, 0.0, 0.0, -0.69, 0.69, { 1.0, 1.0 }, { 1.0E+6, 1.0E+6,
        1.0E+8 } }, { 7.0, 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 5.0, 4.0, 42200.0, 0.0,
      71500.0, 0.0, 7150.0, 0.1, 0.0, 0.0, 0.0, 0.0, -2.4, 2.4, { 1.0, 1.0 }, {
        1.0E+6, 1.0E+6, 1.0E+8 } }, { 8.0, 3.0, 1.0, 0.0, 0.5, 0.0, 0.0, 5.0,
      -4.0, 5866.0, 0.0, 71500.0, 0.0, 7150.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.01,
      6.0, { 1.0, 1.0 }, { 1.0E+6, 1.0E+6, 1.0E+8 } } };

  thruster_data_in_data[0] = r0;
  N_enabled_thruster = 0.0;
  for (i = 0; i < (int)N_max_thr; i++) {
    N_enabled_thruster++;
    thruster_data_in_data[1] = rv0[i];
  }

  qpsolver(thruster_data_in_data, T_r->Tx, T_r->Ty, T_r->Tm, N_enabled_thruster,
           method, solution, &alloc_out->label, &alloc_out->enable,
           &alloc_out->TSP, &alloc_out->ASP, &alloc_out->Tx, &alloc_out->Ty,
           &alloc_out->T, &alloc_out->phi, &alloc_out->Tm, status);
}

/*
 * File trailer for call_qpsolver.c
 *
 * [EOF]
 */
