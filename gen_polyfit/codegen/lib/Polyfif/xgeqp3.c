/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgeqp3.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 16-Aug-2016 09:21:48
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Polyfif.h"
#include "xgeqp3.h"
#include "xscal.h"
#include "xnrm2.h"
#include "Polyfif_emxutil.h"

/* Function Declarations */
static double rt_hypotd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

/*
 * Arguments    : emxArray_real_T *A
 *                double tau_data[]
 *                int tau_size[1]
 *                emxArray_int32_T *jpvt
 * Return Type  : void
 */
void xgeqp3(emxArray_real_T *A, double tau_data[], int tau_size[1],
            emxArray_int32_T *jpvt)
{
  int n;
  int mn;
  int b_n;
  int i0;
  int yk;
  emxArray_real_T *work;
  int k;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int nmi;
  int i;
  double smax;
  double s;
  int i_i;
  double absxk;
  double t;
  int ix;
  int iy;
  int i_ip1;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int exitg1;
  n = A->size[1];
  if (9 <= A->size[1]) {
    mn = 9;
  } else {
    mn = A->size[1];
  }

  tau_size[0] = mn;
  if (A->size[1] < 1) {
    b_n = 0;
  } else {
    b_n = A->size[1];
  }

  i0 = jpvt->size[0] * jpvt->size[1];
  jpvt->size[0] = 1;
  jpvt->size[1] = b_n;
  emxEnsureCapacity((emxArray__common *)jpvt, i0, (int)sizeof(int));
  if (b_n > 0) {
    jpvt->data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      jpvt->data[k - 1] = yk;
    }
  }

  if (A->size[1] == 0) {
  } else {
    emxInit_real_T1(&work, 1);
    yk = A->size[1];
    i0 = work->size[0];
    work->size[0] = yk;
    emxEnsureCapacity((emxArray__common *)work, i0, (int)sizeof(double));
    for (i0 = 0; i0 < yk; i0++) {
      work->data[i0] = 0.0;
    }

    emxInit_real_T1(&vn1, 1);
    emxInit_real_T1(&vn2, 1);
    yk = A->size[1];
    i0 = vn1->size[0];
    vn1->size[0] = yk;
    emxEnsureCapacity((emxArray__common *)vn1, i0, (int)sizeof(double));
    i0 = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity((emxArray__common *)vn2, i0, (int)sizeof(double));
    k = 1;
    for (nmi = 0; nmi + 1 <= n; nmi++) {
      smax = 0.0;
      s = 2.2250738585072014E-308;
      for (yk = k; yk <= k + 8; yk++) {
        absxk = fabs(A->data[yk - 1]);
        if (absxk > s) {
          t = s / absxk;
          smax = 1.0 + smax * t * t;
          s = absxk;
        } else {
          t = absxk / s;
          smax += t * t;
        }
      }

      smax = s * sqrt(smax);
      vn1->data[nmi] = smax;
      vn2->data[nmi] = vn1->data[nmi];
      k += 9;
    }

    for (i = 1; i <= mn; i++) {
      i_i = (i + (i - 1) * 9) - 1;
      nmi = (n - i) + 1;
      if (nmi < 1) {
        yk = 0;
      } else {
        yk = 1;
        if (nmi > 1) {
          ix = i - 1;
          smax = vn1->data[i - 1];
          for (k = 2; k <= nmi; k++) {
            ix++;
            s = vn1->data[ix];
            if (s > smax) {
              yk = k;
              smax = s;
            }
          }
        }
      }

      b_n = (i + yk) - 2;
      if (b_n + 1 != i) {
        ix = 9 * b_n;
        iy = 9 * (i - 1);
        for (k = 0; k < 9; k++) {
          smax = A->data[ix];
          A->data[ix] = A->data[iy];
          A->data[iy] = smax;
          ix++;
          iy++;
        }

        yk = jpvt->data[b_n];
        jpvt->data[b_n] = jpvt->data[i - 1];
        jpvt->data[i - 1] = yk;
        vn1->data[b_n] = vn1->data[i - 1];
        vn2->data[b_n] = vn2->data[i - 1];
      }

      if (i < 9) {
        absxk = A->data[i_i];
        s = 0.0;
        smax = xnrm2(9 - i, A, i_i + 2);
        if (smax != 0.0) {
          smax = rt_hypotd_snf(A->data[i_i], smax);
          if (A->data[i_i] >= 0.0) {
            smax = -smax;
          }

          if (fabs(smax) < 1.0020841800044864E-292) {
            yk = 0;
            do {
              yk++;
              xscal(9 - i, 9.9792015476736E+291, A, i_i + 2);
              smax *= 9.9792015476736E+291;
              absxk *= 9.9792015476736E+291;
            } while (!(fabs(smax) >= 1.0020841800044864E-292));

            smax = xnrm2(9 - i, A, i_i + 2);
            smax = rt_hypotd_snf(absxk, smax);
            if (absxk >= 0.0) {
              smax = -smax;
            }

            s = (smax - absxk) / smax;
            xscal(9 - i, 1.0 / (absxk - smax), A, i_i + 2);
            for (k = 1; k <= yk; k++) {
              smax *= 1.0020841800044864E-292;
            }

            absxk = smax;
          } else {
            s = (smax - A->data[i_i]) / smax;
            xscal(9 - i, 1.0 / (A->data[i_i] - smax), A, i_i + 2);
            absxk = smax;
          }
        }

        tau_data[i - 1] = s;
        A->data[i_i] = absxk;
      } else {
        tau_data[8] = 0.0;
      }

      if (i < n) {
        absxk = A->data[i_i];
        A->data[i_i] = 1.0;
        i_ip1 = i + i * 9;
        if (tau_data[i - 1] != 0.0) {
          lastv = 10 - i;
          yk = (i_i - i) + 9;
          while ((lastv > 0) && (A->data[yk] == 0.0)) {
            lastv--;
            yk--;
          }

          lastc = nmi - 1;
          exitg2 = false;
          while ((!exitg2) && (lastc > 0)) {
            yk = i_ip1 + (lastc - 1) * 9;
            k = yk;
            do {
              exitg1 = 0;
              if (k <= (yk + lastv) - 1) {
                if (A->data[k - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  k++;
                }
              } else {
                lastc--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = 0;
          lastc = 0;
        }

        if (lastv > 0) {
          if (lastc == 0) {
          } else {
            for (iy = 1; iy <= lastc; iy++) {
              work->data[iy - 1] = 0.0;
            }

            iy = 0;
            i0 = i_ip1 + 9 * (lastc - 1);
            for (yk = i_ip1; yk <= i0; yk += 9) {
              ix = i_i;
              smax = 0.0;
              b_n = (yk + lastv) - 1;
              for (k = yk; k <= b_n; k++) {
                smax += A->data[k - 1] * A->data[ix];
                ix++;
              }

              work->data[iy] += smax;
              iy++;
            }
          }

          if (-tau_data[i - 1] == 0.0) {
          } else {
            yk = i_ip1 - 1;
            b_n = 0;
            for (nmi = 1; nmi <= lastc; nmi++) {
              if (work->data[b_n] != 0.0) {
                smax = work->data[b_n] * -tau_data[i - 1];
                ix = i_i;
                i0 = lastv + yk;
                for (k = yk; k + 1 <= i0; k++) {
                  A->data[k] += A->data[ix] * smax;
                  ix++;
                }
              }

              b_n++;
              yk += 9;
            }
          }
        }

        A->data[i_i] = absxk;
      }

      for (nmi = i; nmi + 1 <= n; nmi++) {
        yk = i + 9 * nmi;
        if (vn1->data[nmi] != 0.0) {
          smax = fabs(A->data[(i + A->size[0] * nmi) - 1]) / vn1->data[nmi];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          s = vn1->data[nmi] / vn2->data[nmi];
          s = smax * (s * s);
          if (s <= 1.4901161193847656E-8) {
            if (i < 9) {
              smax = 0.0;
              if (9 - i == 1) {
                smax = fabs(A->data[yk]);
              } else {
                s = 2.2250738585072014E-308;
                b_n = (yk - i) + 9;
                while (yk + 1 <= b_n) {
                  absxk = fabs(A->data[yk]);
                  if (absxk > s) {
                    t = s / absxk;
                    smax = 1.0 + smax * t * t;
                    s = absxk;
                  } else {
                    t = absxk / s;
                    smax += t * t;
                  }

                  yk++;
                }

                smax = s * sqrt(smax);
              }

              vn1->data[nmi] = smax;
              vn2->data[nmi] = vn1->data[nmi];
            } else {
              vn1->data[nmi] = 0.0;
              vn2->data[nmi] = 0.0;
            }
          } else {
            vn1->data[nmi] *= sqrt(smax);
          }
        }
      }
    }

    emxFree_real_T(&vn2);
    emxFree_real_T(&vn1);
    emxFree_real_T(&work);
  }
}

/*
 * File trailer for xgeqp3.c
 *
 * [EOF]
 */
