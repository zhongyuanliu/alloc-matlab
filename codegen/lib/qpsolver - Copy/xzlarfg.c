/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzlarfg.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "xzlarfg.h"
#include "xscal.h"
#include "xnrm2.h"
#include "qpsolver_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : int n
 *                double *alpha1
 *                emxArray_real_T *x
 *                int ix0
 * Return Type  : double
 */
double xzlarfg(int n, double *alpha1, emxArray_real_T *x, int ix0)
{
  double tau;
  double xnorm;
  int knt;
  int k;
  tau = 0.0;
  if (n <= 0) {
  } else {
    xnorm = b_xnrm2(n - 1, x, ix0);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          b_xscal(n - 1, 9.9792015476736E+291, x, ix0);
          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

        xnorm = b_xnrm2(n - 1, x, ix0);
        xnorm = rt_hypotd_snf(*alpha1, xnorm);
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        tau = (xnorm - *alpha1) / xnorm;
        b_xscal(n - 1, 1.0 / (*alpha1 - xnorm), x, ix0);
        for (k = 1; k <= knt; k++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        tau = (xnorm - *alpha1) / xnorm;
        b_xscal(n - 1, 1.0 / (*alpha1 - xnorm), x, ix0);
        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

/*
 * File trailer for xzlarfg.c
 *
 * [EOF]
 */
