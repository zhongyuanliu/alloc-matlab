/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Polyfif.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 16-Aug-2016 09:21:48
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Polyfif.h"
#include "Polyfif_emxutil.h"
#include "polyfit.h"

/* Function Definitions */

/*
 * Arguments    : const double x[9]
 *                const double y[9]
 *                double N_poly
 *                double x_new
 *                emxArray_real_T *r
 *                double *f
 * Return Type  : void
 */
void Polyfif(const double x[9], const double y[9], double N_poly, double x_new,
             emxArray_real_T *r, double *f)
{
  emxArray_real_T *expl_temp;
  double b_expl_temp;
  double c_expl_temp;
  boolean_T b0;
  int k;
  emxInit_real_T(&expl_temp, 2);
  polyfit(x, y, N_poly, r, expl_temp, &b_expl_temp, &c_expl_temp);
  emxFree_real_T(&expl_temp);
  b0 = (r->size[1] == 0);
  if (!b0) {
    *f = r->data[0];
    for (k = 0; k <= r->size[1] - 2; k++) {
      *f = x_new * *f + r->data[k + 1];
    }
  }
}

/*
 * File trailer for Polyfif.c
 *
 * [EOF]
 */
