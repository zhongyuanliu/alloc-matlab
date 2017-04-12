/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xgeqp3.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
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
  int i27;
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
  if (A->size[1] != 0) {
    emxInit_real_T(&work, 1);
    itemp = A->size[1];
    i27 = work->size[0];
    work->size[0] = itemp;
    emxEnsureCapacity((emxArray__common *)work, i27, (int)sizeof(double));
    for (i27 = 0; i27 < itemp; i27++) {
      work->data[i27] = 0.0;
    }

    emxInit_real_T(&vn1, 1);
    emxInit_real_T(&vn2, 1);
    itemp = A->size[1];
    i27 = vn1->size[0];
    vn1->size[0] = itemp;
    emxEnsureCapacity((emxArray__common *)vn1, i27, (int)sizeof(double));
    i27 = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity((emxArray__common *)vn2, i27, (int)sizeof(double));
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
              b_xscal(3 - i, 9.9792015476736E+291, A, i_i + 2);
              beta1 *= 9.9792015476736E+291;
              t *= 9.9792015476736E+291;
            } while (!(fabs(beta1) >= 1.0020841800044864E-292));

            xnorm = c_xnrm2(3 - i, A, i_i + 2);
            beta1 = rt_hypotd_snf(t, xnorm);
            if (t >= 0.0) {
              beta1 = -beta1;
            }

            absxk = (beta1 - t) / beta1;
            b_xscal(3 - i, 1.0 / (t - beta1), A, i_i + 2);
            for (k = 1; k <= itemp; k++) {
              beta1 *= 1.0020841800044864E-292;
            }

            t = beta1;
          } else {
            absxk = (beta1 - A->data[i_i]) / beta1;
            b_xscal(3 - i, 1.0 / (A->data[i_i] - beta1), A, i_i + 2);
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
          if (nmi != 0) {
            for (iy = 1; iy <= nmi; iy++) {
              work->data[iy - 1] = 0.0;
            }

            iy = 0;
            i27 = i_ip1 + 3 * (nmi - 1);
            for (itemp = i_ip1; itemp <= i27; itemp += 3) {
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

          if (!(-tau_data[i - 1] == 0.0)) {
            itemp = i_ip1 - 1;
            pvt = 0;
            for (iy = 1; iy <= nmi; iy++) {
              if (work->data[pvt] != 0.0) {
                xnorm = work->data[pvt] * -tau_data[i - 1];
                ix = i_i;
                i27 = lastv + itemp;
                for (k = itemp; k + 1 <= i27; k++) {
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
  if (!((A->size[0] == 0) || (A->size[1] == 0))) {
    emxInit_real_T(&work, 1);
    k = A->size[1];
    ix = work->size[0];
    work->size[0] = k;
    emxEnsureCapacity((emxArray__common *)work, ix, (int)sizeof(double));
    for (ix = 0; ix < k; ix++) {
      work->data[ix] = 0.0;
    }

    emxInit_real_T(&vn1, 1);
    emxInit_real_T(&vn2, 1);
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
              if (!(mmi < 1)) {
                if (mmi == 1) {
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
 * Arguments    : double A[27]
 *                double tau[3]
 *                int jpvt[3]
 * Return Type  : void
 */
void xgeqp3(double A[27], double tau[3], int jpvt[3])
{
  int k;
  double vn1[3];
  double vn2[3];
  double work[3];
  int iy;
  int i;
  double smax;
  int i_i;
  double temp2;
  int itemp;
  int ix;
  int pvt;
  double absxk;
  double t;
  int i_ip1;
  int i16;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int exitg1;
  k = 1;
  for (iy = 0; iy < 3; iy++) {
    jpvt[iy] = 1 + iy;
    work[iy] = 0.0;
    smax = 0.0;
    temp2 = 2.2250738585072014E-308;
    for (itemp = k; itemp <= k + 8; itemp++) {
      absxk = fabs(A[itemp - 1]);
      if (absxk > temp2) {
        t = temp2 / absxk;
        smax = 1.0 + smax * t * t;
        temp2 = absxk;
      } else {
        t = absxk / temp2;
        smax += t * t;
      }
    }

    smax = temp2 * sqrt(smax);
    vn1[iy] = smax;
    vn2[iy] = vn1[iy];
    k += 9;
  }

  for (i = 0; i < 3; i++) {
    i_i = i + i * 9;
    itemp = 0;
    if (3 - i > 1) {
      ix = i;
      smax = vn1[i];
      for (k = 1; k + 1 <= 3 - i; k++) {
        ix++;
        if (vn1[ix] > smax) {
          itemp = k;
          smax = vn1[ix];
        }
      }
    }

    pvt = i + itemp;
    if (pvt + 1 != i + 1) {
      ix = 9 * pvt;
      iy = 9 * i;
      for (k = 0; k < 9; k++) {
        smax = A[ix];
        A[ix] = A[iy];
        A[iy] = smax;
        ix++;
        iy++;
      }

      itemp = jpvt[pvt];
      jpvt[pvt] = jpvt[i];
      jpvt[i] = itemp;
      vn1[pvt] = vn1[i];
      vn2[pvt] = vn2[i];
    }

    absxk = A[i_i];
    temp2 = 0.0;
    smax = xnrm2(8 - i, A, i_i + 2);
    if (smax != 0.0) {
      smax = rt_hypotd_snf(A[i_i], smax);
      if (A[i_i] >= 0.0) {
        smax = -smax;
      }

      if (fabs(smax) < 1.0020841800044864E-292) {
        itemp = 0;
        do {
          itemp++;
          i16 = (i_i - i) + 9;
          for (k = i_i + 1; k + 1 <= i16; k++) {
            A[k] *= 9.9792015476736E+291;
          }

          smax *= 9.9792015476736E+291;
          absxk *= 9.9792015476736E+291;
        } while (!(fabs(smax) >= 1.0020841800044864E-292));

        smax = rt_hypotd_snf(absxk, xnrm2(8 - i, A, i_i + 2));
        if (absxk >= 0.0) {
          smax = -smax;
        }

        temp2 = (smax - absxk) / smax;
        absxk = 1.0 / (absxk - smax);
        i16 = (i_i - i) + 9;
        for (k = i_i + 1; k + 1 <= i16; k++) {
          A[k] *= absxk;
        }

        for (k = 1; k <= itemp; k++) {
          smax *= 1.0020841800044864E-292;
        }

        absxk = smax;
      } else {
        temp2 = (smax - A[i_i]) / smax;
        absxk = 1.0 / (A[i_i] - smax);
        i16 = (i_i - i) + 9;
        for (k = i_i + 1; k + 1 <= i16; k++) {
          A[k] *= absxk;
        }

        absxk = smax;
      }
    }

    tau[i] = temp2;
    A[i_i] = absxk;
    if (i + 1 < 3) {
      absxk = A[i_i];
      A[i_i] = 1.0;
      i_ip1 = (i + (i + 1) * 9) + 1;
      if (tau[i] != 0.0) {
        lastv = 9 - i;
        itemp = (i_i - i) + 8;
        while ((lastv > 0) && (A[itemp] == 0.0)) {
          lastv--;
          itemp--;
        }

        lastc = 2 - i;
        exitg2 = false;
        while ((!exitg2) && (lastc > 0)) {
          itemp = i_ip1 + (lastc - 1) * 9;
          k = itemp;
          do {
            exitg1 = 0;
            if (k <= (itemp + lastv) - 1) {
              if (A[k - 1] != 0.0) {
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
        if (lastc != 0) {
          for (iy = 1; iy <= lastc; iy++) {
            work[iy - 1] = 0.0;
          }

          iy = 0;
          i16 = i_ip1 + 9 * (lastc - 1);
          for (itemp = i_ip1; itemp <= i16; itemp += 9) {
            ix = i_i;
            smax = 0.0;
            pvt = (itemp + lastv) - 1;
            for (k = itemp; k <= pvt; k++) {
              smax += A[k - 1] * A[ix];
              ix++;
            }

            work[iy] += smax;
            iy++;
          }
        }

        if (!(-tau[i] == 0.0)) {
          itemp = i_ip1 - 1;
          pvt = 0;
          for (iy = 1; iy <= lastc; iy++) {
            if (work[pvt] != 0.0) {
              smax = work[pvt] * -tau[i];
              ix = i_i;
              i16 = lastv + itemp;
              for (k = itemp; k + 1 <= i16; k++) {
                A[k] += A[ix] * smax;
                ix++;
              }
            }

            pvt++;
            itemp += 9;
          }
        }
      }

      A[i_i] = absxk;
    }

    for (iy = i + 1; iy + 1 < 4; iy++) {
      itemp = (i + 9 * iy) + 2;
      if (vn1[iy] != 0.0) {
        smax = fabs(A[i + 9 * iy]) / vn1[iy];
        smax = 1.0 - smax * smax;
        if (smax < 0.0) {
          smax = 0.0;
        }

        temp2 = vn1[iy] / vn2[iy];
        temp2 = smax * (temp2 * temp2);
        if (temp2 <= 1.4901161193847656E-8) {
          smax = 0.0;
          temp2 = 2.2250738585072014E-308;
          pvt = (itemp - i) + 7;
          while (itemp <= pvt) {
            absxk = fabs(A[itemp - 1]);
            if (absxk > temp2) {
              t = temp2 / absxk;
              smax = 1.0 + smax * t * t;
              temp2 = absxk;
            } else {
              t = absxk / temp2;
              smax += t * t;
            }

            itemp++;
          }

          smax = temp2 * sqrt(smax);
          vn1[iy] = smax;
          vn2[iy] = vn1[iy];
        } else {
          vn1[iy] *= sqrt(smax);
        }
      }
    }
  }
}

/*
 * File trailer for xgeqp3.c
 *
 * [EOF]
 */
