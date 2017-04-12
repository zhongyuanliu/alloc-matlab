/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mldivide.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "mldivide.h"
#include "qpsolver_emxutil.h"
#include "xgetrf.h"
#include "qrsolve.h"
#include "xgeqp3.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *A
 *                const emxArray_real_T *B
 *                emxArray_real_T *Y
 * Return Type  : void
 */
void b_mldivide(const emxArray_real_T *A, const emxArray_real_T *B,
                emxArray_real_T *Y)
{
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  emxArray_real_T *b_B;
  unsigned int unnamed_idx_0;
  int m;
  int mn;
  int j;
  int rankA;
  double wj;
  int i;
  emxInit_real_T(&b_A, 2);
  emxInit_real_T1(&tau, 1);
  emxInit_int32_T(&jpvt, 2);
  emxInit_real_T1(&b_B, 1);
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0)) {
    unnamed_idx_0 = (unsigned int)A->size[1];
    m = Y->size[0];
    Y->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)Y, m, (int)sizeof(double));
    j = (int)unnamed_idx_0;
    for (m = 0; m < j; m++) {
      Y->data[m] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    mn = A->size[1];
    m = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, m, (int)sizeof(double));
    j = A->size[0] * A->size[1];
    for (m = 0; m < j; m++) {
      b_A->data[m] = A->data[m];
    }

    xgetrf(A->size[1], A->size[1], b_A, A->size[1], jpvt, &j);
    m = Y->size[0];
    Y->size[0] = B->size[0];
    emxEnsureCapacity((emxArray__common *)Y, m, (int)sizeof(double));
    j = B->size[0];
    for (m = 0; m < j; m++) {
      Y->data[m] = B->data[m];
    }

    for (j = 0; j + 1 < mn; j++) {
      if (jpvt->data[j] != j + 1) {
        wj = Y->data[j];
        Y->data[j] = Y->data[jpvt->data[j] - 1];
        Y->data[jpvt->data[j] - 1] = wj;
      }
    }

    for (j = 0; j + 1 <= mn; j++) {
      m = mn * j;
      if (Y->data[j] != 0.0) {
        for (i = j + 1; i + 1 <= mn; i++) {
          Y->data[i] -= Y->data[j] * b_A->data[i + m];
        }
      }
    }

    for (j = A->size[1] - 1; j + 1 > 0; j--) {
      m = mn * j;
      if (Y->data[j] != 0.0) {
        wj = b_A->data[j + m];
        Y->data[j] /= wj;
        for (i = 0; i + 1 <= j; i++) {
          Y->data[i] -= Y->data[j] * b_A->data[i + m];
        }
      }
    }
  } else {
    m = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, m, (int)sizeof(double));
    j = A->size[0] * A->size[1];
    for (m = 0; m < j; m++) {
      b_A->data[m] = A->data[m];
    }

    c_xgeqp3(b_A, tau, jpvt);
    rankA = rankFromQR(b_A);
    j = b_A->size[1];
    m = Y->size[0];
    Y->size[0] = j;
    emxEnsureCapacity((emxArray__common *)Y, m, (int)sizeof(double));
    for (m = 0; m < j; m++) {
      Y->data[m] = 0.0;
    }

    m = b_B->size[0];
    b_B->size[0] = B->size[0];
    emxEnsureCapacity((emxArray__common *)b_B, m, (int)sizeof(double));
    j = B->size[0];
    for (m = 0; m < j; m++) {
      b_B->data[m] = B->data[m];
    }

    m = b_A->size[0];
    j = b_A->size[0];
    mn = b_A->size[1];
    if (j <= mn) {
      mn = j;
    }

    for (j = 0; j + 1 <= mn; j++) {
      if (tau->data[j] != 0.0) {
        wj = b_B->data[j];
        for (i = j + 1; i + 1 <= m; i++) {
          wj += b_A->data[i + b_A->size[0] * j] * b_B->data[i];
        }

        wj *= tau->data[j];
        if (wj != 0.0) {
          b_B->data[j] -= wj;
          for (i = j + 1; i + 1 <= m; i++) {
            b_B->data[i] -= b_A->data[i + b_A->size[0] * j] * wj;
          }
        }
      }
    }

    for (i = 0; i + 1 <= rankA; i++) {
      Y->data[jpvt->data[i] - 1] = b_B->data[i];
    }

    for (j = rankA - 1; j + 1 > 0; j--) {
      Y->data[jpvt->data[j] - 1] /= b_A->data[j + b_A->size[0] * j];
      for (i = 0; i + 1 <= j; i++) {
        Y->data[jpvt->data[i] - 1] -= Y->data[jpvt->data[j] - 1] * b_A->data[i +
          b_A->size[0] * j];
      }
    }
  }

  emxFree_real_T(&b_B);
  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
}

/*
 * Arguments    : const emxArray_real_T *A
 *                const double B[3]
 *                emxArray_real_T *Y
 * Return Type  : void
 */
void mldivide(const emxArray_real_T *A, const double B[3], emxArray_real_T *Y)
{
  int i9;
  emxArray_real_T *b_A;
  int minmn;
  signed char ipiv[3];
  double A_data[9];
  int j;
  emxArray_int32_T *jpvt;
  double b_B[3];
  int c;
  double tau_data[3];
  int tau_size[1];
  int maxmn;
  int rankR;
  int ix;
  double tol;
  double s;
  int ijA;
  if (A->size[1] == 0) {
    i9 = Y->size[0];
    Y->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)Y, i9, (int)sizeof(double));
  } else if (3 == A->size[1]) {
    minmn = A->size[0] * A->size[1];
    for (i9 = 0; i9 < minmn; i9++) {
      A_data[i9] = A->data[i9];
    }

    for (i9 = 0; i9 < 3; i9++) {
      ipiv[i9] = (signed char)(1 + i9);
    }

    for (j = 0; j < 2; j++) {
      c = j << 2;
      minmn = 0;
      ix = c;
      tol = fabs(A_data[c]);
      for (rankR = 2; rankR <= 3 - j; rankR++) {
        ix++;
        s = fabs(A_data[ix]);
        if (s > tol) {
          minmn = rankR - 1;
          tol = s;
        }
      }

      if (A_data[c + minmn] != 0.0) {
        if (minmn != 0) {
          ipiv[j] = (signed char)((j + minmn) + 1);
          ix = j;
          minmn += j;
          for (rankR = 0; rankR < 3; rankR++) {
            tol = A_data[ix];
            A_data[ix] = A_data[minmn];
            A_data[minmn] = tol;
            ix += 3;
            minmn += 3;
          }
        }

        i9 = (c - j) + 3;
        for (maxmn = c + 1; maxmn + 1 <= i9; maxmn++) {
          A_data[maxmn] /= A_data[c];
        }
      }

      minmn = c + 4;
      maxmn = c + 3;
      for (rankR = 1; rankR <= 2 - j; rankR++) {
        tol = A_data[maxmn];
        if (A_data[maxmn] != 0.0) {
          ix = c + 1;
          i9 = (minmn - j) + 2;
          for (ijA = minmn; ijA + 1 <= i9; ijA++) {
            A_data[ijA] += A_data[ix] * -tol;
            ix++;
          }
        }

        maxmn += 3;
        minmn += 3;
      }
    }

    for (maxmn = 0; maxmn < 3; maxmn++) {
      b_B[maxmn] = B[maxmn];
    }

    for (minmn = 0; minmn < 2; minmn++) {
      if (ipiv[minmn] != minmn + 1) {
        tol = b_B[minmn];
        b_B[minmn] = b_B[ipiv[minmn] - 1];
        b_B[ipiv[minmn] - 1] = tol;
      }
    }

    for (rankR = 0; rankR < 3; rankR++) {
      minmn = 3 * rankR;
      if (b_B[rankR] != 0.0) {
        for (maxmn = rankR + 1; maxmn + 1 < 4; maxmn++) {
          b_B[maxmn] -= b_B[rankR] * A_data[maxmn + minmn];
        }
      }
    }

    for (rankR = 2; rankR >= 0; rankR += -1) {
      minmn = 3 * rankR;
      if (b_B[rankR] != 0.0) {
        b_B[rankR] /= A_data[rankR + minmn];
        for (maxmn = 0; maxmn + 1 <= rankR; maxmn++) {
          b_B[maxmn] -= b_B[rankR] * A_data[maxmn + minmn];
        }
      }
    }

    i9 = Y->size[0];
    Y->size[0] = 3;
    emxEnsureCapacity((emxArray__common *)Y, i9, (int)sizeof(double));
    for (i9 = 0; i9 < 3; i9++) {
      Y->data[i9] = b_B[i9];
    }
  } else {
    emxInit_real_T(&b_A, 2);
    i9 = b_A->size[0] * b_A->size[1];
    b_A->size[0] = 3;
    b_A->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, i9, (int)sizeof(double));
    minmn = A->size[0] * A->size[1];
    for (i9 = 0; i9 < minmn; i9++) {
      b_A->data[i9] = A->data[i9];
    }

    emxInit_int32_T(&jpvt, 2);
    b_xgeqp3(b_A, tau_data, tau_size, jpvt);
    rankR = 0;
    if (3 < b_A->size[1]) {
      minmn = 3;
      maxmn = b_A->size[1];
    } else {
      minmn = b_A->size[1];
      maxmn = 3;
    }

    if (minmn > 0) {
      tol = (double)maxmn * fabs(b_A->data[0]) * 2.2204460492503131E-16;
      while ((rankR < minmn) && (fabs(b_A->data[rankR + b_A->size[0] * rankR]) >=
              tol)) {
        rankR++;
      }
    }

    minmn = b_A->size[1];
    i9 = Y->size[0];
    Y->size[0] = minmn;
    emxEnsureCapacity((emxArray__common *)Y, i9, (int)sizeof(double));
    for (i9 = 0; i9 < minmn; i9++) {
      Y->data[i9] = 0.0;
    }

    for (maxmn = 0; maxmn < 3; maxmn++) {
      b_B[maxmn] = B[maxmn];
    }

    minmn = b_A->size[1];
    if (3 <= minmn) {
      minmn = 3;
    }

    for (j = 0; j + 1 <= minmn; j++) {
      if (tau_data[j] != 0.0) {
        tol = b_B[j];
        for (maxmn = j + 1; maxmn + 1 < 4; maxmn++) {
          tol += b_A->data[maxmn + b_A->size[0] * j] * b_B[maxmn];
        }

        tol *= tau_data[j];
        if (tol != 0.0) {
          b_B[j] -= tol;
          for (maxmn = j + 1; maxmn + 1 < 4; maxmn++) {
            b_B[maxmn] -= b_A->data[maxmn + b_A->size[0] * j] * tol;
          }
        }
      }
    }

    for (maxmn = 0; maxmn + 1 <= rankR; maxmn++) {
      Y->data[jpvt->data[maxmn] - 1] = b_B[maxmn];
    }

    for (j = rankR - 1; j + 1 > 0; j--) {
      Y->data[jpvt->data[j] - 1] /= b_A->data[j + b_A->size[0] * j];
      for (maxmn = 0; maxmn + 1 <= j; maxmn++) {
        Y->data[jpvt->data[maxmn] - 1] -= Y->data[jpvt->data[j] - 1] * b_A->
          data[maxmn + b_A->size[0] * j];
      }
    }

    emxFree_int32_T(&jpvt);
    emxFree_real_T(&b_A);
  }
}

/*
 * File trailer for mldivide.c
 *
 * [EOF]
 */
