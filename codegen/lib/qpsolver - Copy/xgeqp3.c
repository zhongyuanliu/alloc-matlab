/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgeqp3.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "xgeqp3.h"
#include "xnrm2.h"
#include "xzlarf.h"
#include "xzlarfg.h"
#include "ixamax.h"
#include "qpsolver_emxutil.h"
#include "colon.h"
#include "xscal.h"
#include "qpsolver_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *A
 *                double tau_data[]
 *                int tau_size[1]
 *                emxArray_int32_T *jpvt
 * Return Type  : void
 */
void b_xgeqp3(emxArray_real_T *A, double tau_data[], int tau_size[1],
              emxArray_int32_T *jpvt)
{
  int n;
  int mn;
  emxArray_real_T *work;
  int itemp;
  int i28;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int k;
  int iy;
  int i;
  double xnorm;
  double beta1;
  int i_i;
  int nmi;
  double absxk;
  int pvt;
  double t;
  int ix;
  int i_ip1;
  int lastv;
  boolean_T exitg2;
  int exitg1;
  n = A->size[1];
  if (3 <= A->size[1]) {
    mn = 3;
  } else {
    mn = A->size[1];
  }

  tau_size[0] = mn;
  eml_signed_integer_colon(A->size[1], jpvt);
  if (A->size[1] == 0) {
  } else {
    emxInit_real_T1(&work, 1);
    itemp = A->size[1];
    i28 = work->size[0];
    work->size[0] = itemp;
    emxEnsureCapacity((emxArray__common *)work, i28, (int)sizeof(double));
    for (i28 = 0; i28 < itemp; i28++) {
      work->data[i28] = 0.0;
    }

    emxInit_real_T1(&vn1, 1);
    emxInit_real_T1(&vn2, 1);
    itemp = A->size[1];
    i28 = vn1->size[0];
    vn1->size[0] = itemp;
    emxEnsureCapacity((emxArray__common *)vn1, i28, (int)sizeof(double));
    i28 = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity((emxArray__common *)vn2, i28, (int)sizeof(double));
    k = 1;
    for (iy = 0; iy + 1 <= n; iy++) {
      xnorm = 0.0;
      beta1 = 2.2250738585072014E-308;
      for (itemp = k; itemp <= k + 2; itemp++) {
        absxk = fabs(A->data[itemp - 1]);
        if (absxk > beta1) {
          t = beta1 / absxk;
          xnorm = 1.0 + xnorm * t * t;
          beta1 = absxk;
        } else {
          t = absxk / beta1;
          xnorm += t * t;
        }
      }

      xnorm = beta1 * sqrt(xnorm);
      vn1->data[iy] = xnorm;
      vn2->data[iy] = vn1->data[iy];
      k += 3;
    }

    for (i = 1; i <= mn; i++) {
      i_i = (i + (i - 1) * 3) - 1;
      nmi = n - i;
      itemp = ixamax(1 + nmi, vn1, i);
      pvt = (i + itemp) - 2;
      if (pvt + 1 != i) {
        ix = 3 * pvt;
        iy = 3 * (i - 1);
        for (k = 0; k < 3; k++) {
          xnorm = A->data[ix];
          A->data[ix] = A->data[iy];
          A->data[iy] = xnorm;
          ix++;
          iy++;
        }

        itemp = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i - 1];
        jpvt->data[i - 1] = itemp;
        vn1->data[pvt] = vn1->data[i - 1];
        vn2->data[pvt] = vn2->data[i - 1];
      }

      if (i < 3) {
        t = A->data[i_i];
        absxk = 0.0;
        xnorm = c_xnrm2(3 - i, A, i_i + 2);
        if (xnorm != 0.0) {
          beta1 = rt_hypotd_snf(A->data[i_i], xnorm);
          if (A->data[i_i] >= 0.0) {
            beta1 = -beta1;
          }

          if (fabs(beta1) < 1.0020841800044864E-292) {
            itemp = 0;
            do {
              itemp++;
              c_xscal(3 - i, 9.9792015476736E+291, A, i_i + 2);
              beta1 *= 9.9792015476736E+291;
              t *= 9.9792015476736E+291;
            } while (!(fabs(beta1) >= 1.0020841800044864E-292));

            xnorm = c_xnrm2(3 - i, A, i_i + 2);
            beta1 = rt_hypotd_snf(t, xnorm);
            if (t >= 0.0) {
              beta1 = -beta1;
            }

            absxk = (beta1 - t) / beta1;
            c_xscal(3 - i, 1.0 / (t - beta1), A, i_i + 2);
            for (k = 1; k <= itemp; k++) {
              beta1 *= 1.0020841800044864E-292;
            }

            t = beta1;
          } else {
            absxk = (beta1 - A->data[i_i]) / beta1;
            c_xscal(3 - i, 1.0 / (A->data[i_i] - beta1), A, i_i + 2);
            t = beta1;
          }
        }

        tau_data[i - 1] = absxk;
        A->data[i_i] = t;
      } else {
        tau_data[2] = 0.0;
      }

      if (i < n) {
        t = A->data[i_i];
        A->data[i_i] = 1.0;
        i_ip1 = i + i * 3;
        if (tau_data[i - 1] != 0.0) {
          lastv = 4 - i;
          itemp = (i_i - i) + 3;
          while ((lastv > 0) && (A->data[itemp] == 0.0)) {
            lastv--;
            itemp--;
          }

          exitg2 = false;
          while ((!exitg2) && (nmi > 0)) {
            itemp = i_ip1 + (nmi - 1) * 3;
            k = itemp;
            do {
              exitg1 = 0;
              if (k <= (itemp + lastv) - 1) {
                if (A->data[k - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  k++;
                }
              } else {
                nmi--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = 0;
          nmi = 0;
        }

        if (lastv > 0) {
          if (nmi == 0) {
          } else {
            for (iy = 1; iy <= nmi; iy++) {
              work->data[iy - 1] = 0.0;
            }

            iy = 0;
            i28 = i_ip1 + 3 * (nmi - 1);
            for (itemp = i_ip1; itemp <= i28; itemp += 3) {
              ix = i_i;
              xnorm = 0.0;
              pvt = (itemp + lastv) - 1;
              for (k = itemp; k <= pvt; k++) {
                xnorm += A->data[k - 1] * A->data[ix];
                ix++;
              }

              work->data[iy] += xnorm;
              iy++;
            }
          }

          if (-tau_data[i - 1] == 0.0) {
          } else {
            itemp = i_ip1 - 1;
            pvt = 0;
            for (iy = 1; iy <= nmi; iy++) {
              if (work->data[pvt] != 0.0) {
                xnorm = work->data[pvt] * -tau_data[i - 1];
                ix = i_i;
                i28 = lastv + itemp;
                for (k = itemp; k + 1 <= i28; k++) {
                  A->data[k] += A->data[ix] * xnorm;
                  ix++;
                }
              }

              pvt++;
              itemp += 3;
            }
          }
        }

        A->data[i_i] = t;
      }

      for (iy = i; iy + 1 <= n; iy++) {
        itemp = i + 3 * iy;
        if (vn1->data[iy] != 0.0) {
          xnorm = fabs(A->data[(i + A->size[0] * iy) - 1]) / vn1->data[iy];
          xnorm = 1.0 - xnorm * xnorm;
          if (xnorm < 0.0) {
            xnorm = 0.0;
          }

          beta1 = vn1->data[iy] / vn2->data[iy];
          beta1 = xnorm * (beta1 * beta1);
          if (beta1 <= 1.4901161193847656E-8) {
            if (i < 3) {
              xnorm = 0.0;
              if (3 - i == 1) {
                xnorm = fabs(A->data[itemp]);
              } else {
                beta1 = 2.2250738585072014E-308;
                for (k = itemp; k + 1 <= itemp + 2; k++) {
                  absxk = fabs(A->data[k]);
                  if (absxk > beta1) {
                    t = beta1 / absxk;
                    xnorm = 1.0 + xnorm * t * t;
                    beta1 = absxk;
                  } else {
                    t = absxk / beta1;
                    xnorm += t * t;
                  }
                }

                xnorm = beta1 * sqrt(xnorm);
              }

              vn1->data[iy] = xnorm;
              vn2->data[iy] = vn1->data[iy];
            } else {
              vn1->data[iy] = 0.0;
              vn2->data[iy] = 0.0;
            }
          } else {
            vn1->data[iy] *= sqrt(xnorm);
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
 * Arguments    : emxArray_real_T *A
 *                emxArray_real_T *tau
 *                emxArray_int32_T *jpvt
 * Return Type  : void
 */
void c_xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, emxArray_int32_T *jpvt)
{
  int m;
  int n;
  int mn;
  int ix;
  emxArray_real_T *work;
  int k;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int iy;
  int i;
  int i_i;
  int nmi;
  int mmi;
  int pvt;
  double temp1;
  double temp2;
  double absxk;
  double t;
  m = A->size[0];
  n = A->size[1];
  if (A->size[0] <= A->size[1]) {
    mn = A->size[0];
  } else {
    mn = A->size[1];
  }

  ix = tau->size[0];
  tau->size[0] = mn;
  emxEnsureCapacity((emxArray__common *)tau, ix, (int)sizeof(double));
  eml_signed_integer_colon(A->size[1], jpvt);
  if ((A->size[0] == 0) || (A->size[1] == 0)) {
  } else {
    emxInit_real_T1(&work, 1);
    k = A->size[1];
    ix = work->size[0];
    work->size[0] = k;
    emxEnsureCapacity((emxArray__common *)work, ix, (int)sizeof(double));
    for (ix = 0; ix < k; ix++) {
      work->data[ix] = 0.0;
    }

    emxInit_real_T1(&vn1, 1);
    emxInit_real_T1(&vn2, 1);
    k = A->size[1];
    ix = vn1->size[0];
    vn1->size[0] = k;
    emxEnsureCapacity((emxArray__common *)vn1, ix, (int)sizeof(double));
    ix = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity((emxArray__common *)vn2, ix, (int)sizeof(double));
    k = 1;
    for (iy = 0; iy + 1 <= n; iy++) {
      vn1->data[iy] = b_xnrm2(m, A, k);
      vn2->data[iy] = vn1->data[iy];
      k += m;
    }

    for (i = 1; i <= mn; i++) {
      i_i = (i + (i - 1) * m) - 1;
      nmi = n - i;
      mmi = m - i;
      k = ixamax(1 + nmi, vn1, i);
      pvt = (i + k) - 2;
      if (pvt + 1 != i) {
        ix = m * pvt;
        iy = m * (i - 1);
        for (k = 1; k <= m; k++) {
          temp1 = A->data[ix];
          A->data[ix] = A->data[iy];
          A->data[iy] = temp1;
          ix++;
          iy++;
        }

        k = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i - 1];
        jpvt->data[i - 1] = k;
        vn1->data[pvt] = vn1->data[i - 1];
        vn2->data[pvt] = vn2->data[i - 1];
      }

      if (i < m) {
        temp1 = A->data[i_i];
        tau->data[i - 1] = xzlarfg(1 + mmi, &temp1, A, i_i + 2);
        A->data[i_i] = temp1;
      } else {
        tau->data[i - 1] = 0.0;
      }

      if (i < n) {
        temp1 = A->data[i_i];
        A->data[i_i] = 1.0;
        xzlarf(1 + mmi, nmi, i_i + 1, tau->data[i - 1], A, i + i * m, m, work);
        A->data[i_i] = temp1;
      }

      for (iy = i; iy + 1 <= n; iy++) {
        k = i + m * iy;
        if (vn1->data[iy] != 0.0) {
          temp1 = fabs(A->data[(i + A->size[0] * iy) - 1]) / vn1->data[iy];
          temp1 = 1.0 - temp1 * temp1;
          if (temp1 < 0.0) {
            temp1 = 0.0;
          }

          temp2 = vn1->data[iy] / vn2->data[iy];
          temp2 = temp1 * (temp2 * temp2);
          if (temp2 <= 1.4901161193847656E-8) {
            if (i < m) {
              temp1 = 0.0;
              if (mmi < 1) {
              } else if (mmi == 1) {
                temp1 = fabs(A->data[k]);
              } else {
                temp2 = 2.2250738585072014E-308;
                ix = k + mmi;
                while (k + 1 <= ix) {
                  absxk = fabs(A->data[k]);
                  if (absxk > temp2) {
                    t = temp2 / absxk;
                    temp1 = 1.0 + temp1 * t * t;
                    temp2 = absxk;
                  } else {
                    t = absxk / temp2;
                    temp1 += t * t;
                  }

                  k++;
                }

                temp1 = temp2 * sqrt(temp1);
              }

              vn1->data[iy] = temp1;
              vn2->data[iy] = vn1->data[iy];
            } else {
              vn1->data[iy] = 0.0;
              vn2->data[iy] = 0.0;
            }
          } else {
            vn1->data[iy] *= sqrt(temp1);
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
  emxArray_real_T *work;
  int itemp;
  int i15;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int k;
  int iy;
  int i;
  double xnorm;
  double beta1;
  int i_i;
  int nmi;
  double absxk;
  int pvt;
  double t;
  int ix;
  int i_ip1;
  int lastv;
  boolean_T exitg2;
  int exitg1;
  n = A->size[1];
  if (9 <= A->size[1]) {
    mn = 9;
  } else {
    mn = A->size[1];
  }

  tau_size[0] = mn;
  eml_signed_integer_colon(A->size[1], jpvt);
  if (A->size[1] == 0) {
  } else {
    emxInit_real_T1(&work, 1);
    itemp = A->size[1];
    i15 = work->size[0];
    work->size[0] = itemp;
    emxEnsureCapacity((emxArray__common *)work, i15, (int)sizeof(double));
    for (i15 = 0; i15 < itemp; i15++) {
      work->data[i15] = 0.0;
    }

    emxInit_real_T1(&vn1, 1);
    emxInit_real_T1(&vn2, 1);
    itemp = A->size[1];
    i15 = vn1->size[0];
    vn1->size[0] = itemp;
    emxEnsureCapacity((emxArray__common *)vn1, i15, (int)sizeof(double));
    i15 = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity((emxArray__common *)vn2, i15, (int)sizeof(double));
    k = 1;
    for (iy = 0; iy + 1 <= n; iy++) {
      xnorm = 0.0;
      beta1 = 2.2250738585072014E-308;
      for (itemp = k; itemp <= k + 8; itemp++) {
        absxk = fabs(A->data[itemp - 1]);
        if (absxk > beta1) {
          t = beta1 / absxk;
          xnorm = 1.0 + xnorm * t * t;
          beta1 = absxk;
        } else {
          t = absxk / beta1;
          xnorm += t * t;
        }
      }

      xnorm = beta1 * sqrt(xnorm);
      vn1->data[iy] = xnorm;
      vn2->data[iy] = vn1->data[iy];
      k += 9;
    }

    for (i = 1; i <= mn; i++) {
      i_i = (i + (i - 1) * 9) - 1;
      nmi = n - i;
      itemp = ixamax(1 + nmi, vn1, i);
      pvt = (i + itemp) - 2;
      if (pvt + 1 != i) {
        ix = 9 * pvt;
        iy = 9 * (i - 1);
        for (k = 0; k < 9; k++) {
          xnorm = A->data[ix];
          A->data[ix] = A->data[iy];
          A->data[iy] = xnorm;
          ix++;
          iy++;
        }

        itemp = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i - 1];
        jpvt->data[i - 1] = itemp;
        vn1->data[pvt] = vn1->data[i - 1];
        vn2->data[pvt] = vn2->data[i - 1];
      }

      if (i < 9) {
        t = A->data[i_i];
        absxk = 0.0;
        xnorm = xnrm2(9 - i, A, i_i + 2);
        if (xnorm != 0.0) {
          beta1 = rt_hypotd_snf(A->data[i_i], xnorm);
          if (A->data[i_i] >= 0.0) {
            beta1 = -beta1;
          }

          if (fabs(beta1) < 1.0020841800044864E-292) {
            pvt = 0;
            do {
              pvt++;
              xscal(9 - i, 9.9792015476736E+291, A, i_i + 2);
              beta1 *= 9.9792015476736E+291;
              t *= 9.9792015476736E+291;
            } while (!(fabs(beta1) >= 1.0020841800044864E-292));

            xnorm = xnrm2(9 - i, A, i_i + 2);
            beta1 = rt_hypotd_snf(t, xnorm);
            if (t >= 0.0) {
              beta1 = -beta1;
            }

            absxk = (beta1 - t) / beta1;
            xscal(9 - i, 1.0 / (t - beta1), A, i_i + 2);
            for (k = 1; k <= pvt; k++) {
              beta1 *= 1.0020841800044864E-292;
            }

            t = beta1;
          } else {
            absxk = (beta1 - A->data[i_i]) / beta1;
            xscal(9 - i, 1.0 / (A->data[i_i] - beta1), A, i_i + 2);
            t = beta1;
          }
        }

        tau_data[i - 1] = absxk;
        A->data[i_i] = t;
      } else {
        tau_data[8] = 0.0;
      }

      if (i < n) {
        t = A->data[i_i];
        A->data[i_i] = 1.0;
        i_ip1 = i + i * 9;
        if (tau_data[i - 1] != 0.0) {
          lastv = 10 - i;
          itemp = (i_i - i) + 9;
          while ((lastv > 0) && (A->data[itemp] == 0.0)) {
            lastv--;
            itemp--;
          }

          exitg2 = false;
          while ((!exitg2) && (nmi > 0)) {
            itemp = i_ip1 + (nmi - 1) * 9;
            k = itemp;
            do {
              exitg1 = 0;
              if (k <= (itemp + lastv) - 1) {
                if (A->data[k - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  k++;
                }
              } else {
                nmi--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = 0;
          nmi = 0;
        }

        if (lastv > 0) {
          if (nmi == 0) {
          } else {
            for (iy = 1; iy <= nmi; iy++) {
              work->data[iy - 1] = 0.0;
            }

            iy = 0;
            i15 = i_ip1 + 9 * (nmi - 1);
            for (itemp = i_ip1; itemp <= i15; itemp += 9) {
              ix = i_i;
              xnorm = 0.0;
              pvt = (itemp + lastv) - 1;
              for (k = itemp; k <= pvt; k++) {
                xnorm += A->data[k - 1] * A->data[ix];
                ix++;
              }

              work->data[iy] += xnorm;
              iy++;
            }
          }

          if (-tau_data[i - 1] == 0.0) {
          } else {
            itemp = i_ip1 - 1;
            pvt = 0;
            for (iy = 1; iy <= nmi; iy++) {
              if (work->data[pvt] != 0.0) {
                xnorm = work->data[pvt] * -tau_data[i - 1];
                ix = i_i;
                i15 = lastv + itemp;
                for (k = itemp; k + 1 <= i15; k++) {
                  A->data[k] += A->data[ix] * xnorm;
                  ix++;
                }
              }

              pvt++;
              itemp += 9;
            }
          }
        }

        A->data[i_i] = t;
      }

      for (iy = i; iy + 1 <= n; iy++) {
        itemp = i + 9 * iy;
        if (vn1->data[iy] != 0.0) {
          xnorm = fabs(A->data[(i + A->size[0] * iy) - 1]) / vn1->data[iy];
          xnorm = 1.0 - xnorm * xnorm;
          if (xnorm < 0.0) {
            xnorm = 0.0;
          }

          beta1 = vn1->data[iy] / vn2->data[iy];
          beta1 = xnorm * (beta1 * beta1);
          if (beta1 <= 1.4901161193847656E-8) {
            if (i < 9) {
              xnorm = 0.0;
              if (9 - i == 1) {
                xnorm = fabs(A->data[itemp]);
              } else {
                beta1 = 2.2250738585072014E-308;
                pvt = (itemp - i) + 9;
                while (itemp + 1 <= pvt) {
                  absxk = fabs(A->data[itemp]);
                  if (absxk > beta1) {
                    t = beta1 / absxk;
                    xnorm = 1.0 + xnorm * t * t;
                    beta1 = absxk;
                  } else {
                    t = absxk / beta1;
                    xnorm += t * t;
                  }

                  itemp++;
                }

                xnorm = beta1 * sqrt(xnorm);
              }

              vn1->data[iy] = xnorm;
              vn2->data[iy] = vn1->data[iy];
            } else {
              vn1->data[iy] = 0.0;
              vn2->data[iy] = 0.0;
            }
          } else {
            vn1->data[iy] *= sqrt(xnorm);
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
