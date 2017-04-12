/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qr.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "qr.h"
#include "qpsolver_emxutil.h"
#include "xscal.h"
#include "xgeqrf.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *A
 *                emxArray_real_T *Q
 *                emxArray_real_T *R
 * Return Type  : void
 */
void qr(const emxArray_real_T *A, emxArray_real_T *Q, emxArray_real_T *R)
{
  int m;
  int n;
  int coltop;
  int i;
  int i4;
  emxArray_real_T *tau;
  emxArray_real_T *work;
  emxArray_real_T *b_A;
  int j;
  int b_i;
  int itau;
  int iaii;
  int c;
  int lastv;
  int lastc;
  boolean_T exitg4;
  boolean_T exitg2;
  int exitg3;
  double b_c;
  int ix;
  int exitg1;
  m = A->size[0];
  n = A->size[1];
  coltop = A->size[0];
  i = A->size[0];
  i4 = Q->size[0] * Q->size[1];
  Q->size[0] = coltop;
  Q->size[1] = i;
  emxEnsureCapacity((emxArray__common *)Q, i4, (int)sizeof(double));
  for (i4 = 0; i4 < 2; i4++) {
    coltop = R->size[0] * R->size[1];
    R->size[i4] = A->size[i4];
    emxEnsureCapacity((emxArray__common *)R, coltop, (int)sizeof(double));
  }

  emxInit_real_T(&tau, 1);
  emxInit_real_T(&work, 1);
  if (A->size[0] > A->size[1]) {
    for (j = 0; j + 1 <= n; j++) {
      for (b_i = 0; b_i + 1 <= m; b_i++) {
        Q->data[b_i + Q->size[0] * j] = A->data[b_i + A->size[0] * j];
      }
    }

    for (j = A->size[1]; j + 1 <= m; j++) {
      for (b_i = 1; b_i <= m; b_i++) {
        Q->data[(b_i + Q->size[0] * j) - 1] = 0.0;
      }
    }

    xgeqrf(Q, tau);
    for (j = 0; j + 1 <= n; j++) {
      for (b_i = 0; b_i + 1 <= j + 1; b_i++) {
        R->data[b_i + R->size[0] * j] = Q->data[b_i + Q->size[0] * j];
      }

      for (b_i = j + 1; b_i + 1 <= m; b_i++) {
        R->data[b_i + R->size[0] * j] = 0.0;
      }
    }

    if (!(A->size[0] < 1)) {
      for (j = A->size[1]; j < m; j++) {
        coltop = j * m;
        for (b_i = 0; b_i < m; b_i++) {
          Q->data[coltop + b_i] = 0.0;
        }

        Q->data[coltop + j] = 1.0;
      }

      itau = A->size[1] - 1;
      coltop = Q->size[1];
      i4 = work->size[0];
      work->size[0] = coltop;
      emxEnsureCapacity((emxArray__common *)work, i4, (int)sizeof(double));
      for (i4 = 0; i4 < coltop; i4++) {
        work->data[i4] = 0.0;
      }

      for (b_i = A->size[1]; b_i >= 1; b_i--) {
        iaii = (b_i + (b_i - 1) * m) - 1;
        if (b_i < m) {
          Q->data[iaii] = 1.0;
          coltop = m - b_i;
          c = (iaii + m) + 1;
          if (tau->data[itau] != 0.0) {
            lastv = coltop + 1;
            i = iaii + coltop;
            while ((lastv > 0) && (Q->data[i] == 0.0)) {
              lastv--;
              i--;
            }

            lastc = m - b_i;
            exitg4 = false;
            while ((!exitg4) && (lastc > 0)) {
              coltop = c + (lastc - 1) * m;
              j = coltop;
              do {
                exitg3 = 0;
                if (j <= (coltop + lastv) - 1) {
                  if (Q->data[j - 1] != 0.0) {
                    exitg3 = 1;
                  } else {
                    j++;
                  }
                } else {
                  lastc--;
                  exitg3 = 2;
                }
              } while (exitg3 == 0);

              if (exitg3 == 1) {
                exitg4 = true;
              }
            }
          } else {
            lastv = 0;
            lastc = 0;
          }

          if (lastv > 0) {
            if (lastc != 0) {
              for (i = 1; i <= lastc; i++) {
                work->data[i - 1] = 0.0;
              }

              i = 0;
              i4 = c + m * (lastc - 1);
              for (n = c; n <= i4; n += m) {
                ix = iaii;
                b_c = 0.0;
                coltop = (n + lastv) - 1;
                for (j = n; j <= coltop; j++) {
                  b_c += Q->data[j - 1] * Q->data[ix];
                  ix++;
                }

                work->data[i] += b_c;
                i++;
              }
            }

            if (!(-tau->data[itau] == 0.0)) {
              i = c - 1;
              n = 0;
              for (j = 1; j <= lastc; j++) {
                if (work->data[n] != 0.0) {
                  b_c = work->data[n] * -tau->data[itau];
                  ix = iaii;
                  i4 = lastv + i;
                  for (coltop = i; coltop + 1 <= i4; coltop++) {
                    Q->data[coltop] += Q->data[ix] * b_c;
                    ix++;
                  }
                }

                n++;
                i += m;
              }
            }
          }

          xscal(m - b_i, -tau->data[itau], Q, iaii + 2);
        }

        Q->data[iaii] = 1.0 - tau->data[itau];
        for (j = 1; j < b_i; j++) {
          Q->data[iaii - j] = 0.0;
        }

        itau--;
      }
    }
  } else {
    emxInit_real_T1(&b_A, 2);
    i4 = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, i4, (int)sizeof(double));
    i = A->size[0] * A->size[1];
    for (i4 = 0; i4 < i; i4++) {
      b_A->data[i4] = A->data[i4];
    }

    xgeqrf(b_A, tau);
    for (j = 0; j + 1 <= m; j++) {
      for (b_i = 0; b_i + 1 <= j + 1; b_i++) {
        R->data[b_i + R->size[0] * j] = b_A->data[b_i + b_A->size[0] * j];
      }

      for (b_i = j + 1; b_i + 1 <= m; b_i++) {
        R->data[b_i + R->size[0] * j] = 0.0;
      }
    }

    for (j = A->size[0]; j + 1 <= n; j++) {
      for (b_i = 0; b_i + 1 <= m; b_i++) {
        R->data[b_i + R->size[0] * j] = b_A->data[b_i + b_A->size[0] * j];
      }
    }

    if (!(A->size[0] < 1)) {
      for (j = A->size[0]; j < m; j++) {
        coltop = j * m;
        for (b_i = 0; b_i < m; b_i++) {
          b_A->data[coltop + b_i] = 0.0;
        }

        b_A->data[coltop + j] = 1.0;
      }

      itau = A->size[0] - 1;
      coltop = b_A->size[1];
      i4 = work->size[0];
      work->size[0] = coltop;
      emxEnsureCapacity((emxArray__common *)work, i4, (int)sizeof(double));
      for (i4 = 0; i4 < coltop; i4++) {
        work->data[i4] = 0.0;
      }

      for (b_i = A->size[0]; b_i >= 1; b_i--) {
        iaii = (b_i + (b_i - 1) * m) - 1;
        if (b_i < m) {
          b_A->data[iaii] = 1.0;
          coltop = m - b_i;
          c = (iaii + m) + 1;
          if (tau->data[itau] != 0.0) {
            lastv = coltop + 1;
            i = iaii + coltop;
            while ((lastv > 0) && (b_A->data[i] == 0.0)) {
              lastv--;
              i--;
            }

            lastc = m - b_i;
            exitg2 = false;
            while ((!exitg2) && (lastc > 0)) {
              coltop = c + (lastc - 1) * m;
              j = coltop;
              do {
                exitg1 = 0;
                if (j <= (coltop + lastv) - 1) {
                  if (b_A->data[j - 1] != 0.0) {
                    exitg1 = 1;
                  } else {
                    j++;
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
              for (i = 1; i <= lastc; i++) {
                work->data[i - 1] = 0.0;
              }

              i = 0;
              i4 = c + m * (lastc - 1);
              for (n = c; n <= i4; n += m) {
                ix = iaii;
                b_c = 0.0;
                coltop = (n + lastv) - 1;
                for (j = n; j <= coltop; j++) {
                  b_c += b_A->data[j - 1] * b_A->data[ix];
                  ix++;
                }

                work->data[i] += b_c;
                i++;
              }
            }

            if (!(-tau->data[itau] == 0.0)) {
              i = c - 1;
              n = 0;
              for (j = 1; j <= lastc; j++) {
                if (work->data[n] != 0.0) {
                  b_c = work->data[n] * -tau->data[itau];
                  ix = iaii;
                  i4 = lastv + i;
                  for (coltop = i; coltop + 1 <= i4; coltop++) {
                    b_A->data[coltop] += b_A->data[ix] * b_c;
                    ix++;
                  }
                }

                n++;
                i += m;
              }
            }
          }

          xscal(m - b_i, -tau->data[itau], b_A, iaii + 2);
        }

        b_A->data[iaii] = 1.0 - tau->data[itau];
        for (j = 1; j < b_i; j++) {
          b_A->data[iaii - j] = 0.0;
        }

        itau--;
      }
    }

    for (j = 0; j + 1 <= m; j++) {
      for (b_i = 0; b_i + 1 <= m; b_i++) {
        Q->data[b_i + Q->size[0] * j] = b_A->data[b_i + b_A->size[0] * j];
      }
    }

    emxFree_real_T(&b_A);
  }

  emxFree_real_T(&work);
  emxFree_real_T(&tau);
}

/*
 * File trailer for qr.c
 *
 * [EOF]
 */
