/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: QPACM.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "QPACM.h"
#include "qpsolver_emxutil.h"
#include "nullAssignment.h"
#include "rdivide.h"
#include "sum.h"
#include "abs.h"
#include "mldivide.h"
#include "xtrsm.h"
#include "xgetrf.h"
#include "qrsolve.h"
#include "xgeqp3.h"
#include "qr.h"

/* Function Declarations */
static void null_space(const emxArray_real_T *G, emxArray_real_T *g,
  emxArray_real_T *A, const emxArray_real_T *b, emxArray_real_T *x,
  emxArray_real_T *u);

/* Function Definitions */

/*
 * %
 *  NULL SPACE sol v e s the eq ual i ty constrained convex QP:
 *  min 1/2x 'Gx+g ' x (G i s required to be postiv e semi
 *  d e f i n i t e )
 *  s . t . A' x = b (A i s required to have f u l l column
 *  rank )
 *  where the number of v a ri a b l es i s n and the number of c o nst ra int s i s m.
 *  The nul l space of the OP i s used to fi nd the so l ut io n .
 *
 *  Call
 *  [ x , u ] = n ul l sp ac e (G,A, g , b)
 *
 *  Input parameters
 *  G : i s the Hessian matrix (nxn ) of the QP.
 *  A : i s the c onstr ai nt matrix (nxm) : every column contains
 *  a from the
 *  e qual i ty : a ' x = b .
 *  g : i s the gradient (nx1 ) of the QP.
 *  b : i s the r i ght hand si de of the c o nst ra i nt s .
 *
 *  Output parameters
 *  x : the sol u ti o n
 *  mu : the lagrangian m u l t i p l i e r s
 *
 *  By : Carsten V\ ¨olcker , s961572 & Esben Lundsager Hansen , s022022 .
 *  Subject : Numerical Methods f or Sequential Quadratic Optimization ,
 *  Master Thesis , IMM, DTU, DK-2800 Lyngby .
 *  Supervisor : John Bagterp Jørgensen , Assistant Pr ofe ssor & Per Grove
 *  Thomsen , Professor .
 *  Date : 08. february 2007.
 * Arguments    : const emxArray_real_T *G
 *                emxArray_real_T *g
 *                emxArray_real_T *A
 *                const emxArray_real_T *b
 *                emxArray_real_T *x
 *                emxArray_real_T *u
 * Return Type  : void
 */
static void null_space(const emxArray_real_T *G, emxArray_real_T *g,
  emxArray_real_T *A, const emxArray_real_T *b, emxArray_real_T *x,
  emxArray_real_T *u)
{
  int i10;
  int kAcol;
  emxArray_real_T *b_A;
  int mn;
  int i11;
  emxArray_real_T *b_G;
  emxArray_real_T *Q1;
  emxArray_real_T *Q2t;
  emxArray_real_T *Q;
  emxArray_real_T *Q2;
  int i12;
  emxArray_real_T *R;
  int i13;
  emxArray_real_T *py;
  emxArray_real_T *pz;
  emxArray_int32_T *jpvt;
  emxArray_real_T *B;
  emxArray_real_T *r2;
  unsigned int unnamed_idx_0;
  int n;
  unsigned int unnamed_idx_1;
  int nb;
  int jBcol;
  int ar;
  int m;
  int j;
  emxArray_real_T *y;
  int k;
  double wj;
  int iy;
  emxArray_real_T *C;
  int ix;
  int b_n;
  emxArray_real_T *gz;
  int ia;
  emxArray_real_T *b_y;
  emxArray_real_T *Gz;
  boolean_T guard1 = false;
  boolean_T exitg1;
  emxArray_real_T *b_gz;
  emxArray_real_T *b_Q2t;
  double c;
  emxArray_real_T *b_pz;
  emxArray_real_T *b_C;
  emxArray_real_T *c_C;
  emxArray_real_T *c_y;
  i10 = g->size[0];
  emxEnsureCapacity((emxArray__common *)g, i10, (int)sizeof(double));
  kAcol = g->size[0];
  for (i10 = 0; i10 < kAcol; i10++) {
    g->data[i10] = -g->data[i10];
  }

  emxInit_real_T(&b_A, 2);
  i10 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[1];
  b_A->size[1] = A->size[0];
  emxEnsureCapacity((emxArray__common *)b_A, i10, (int)sizeof(double));
  kAcol = A->size[0];
  for (i10 = 0; i10 < kAcol; i10++) {
    mn = A->size[1];
    for (i11 = 0; i11 < mn; i11++) {
      b_A->data[i11 + b_A->size[0] * i10] = A->data[i10 + A->size[0] * i11];
    }
  }

  i10 = A->size[0] * A->size[1];
  A->size[0] = b_A->size[0];
  A->size[1] = b_A->size[1];
  emxEnsureCapacity((emxArray__common *)A, i10, (int)sizeof(double));
  kAcol = b_A->size[1];
  for (i10 = 0; i10 < kAcol; i10++) {
    mn = b_A->size[0];
    for (i11 = 0; i11 < mn; i11++) {
      A->data[i11 + A->size[0] * i10] = b_A->data[i11 + b_A->size[0] * i10];
    }
  }

  emxFree_real_T(&b_A);
  if (A->size[1] != 0) {
    emxInit_real_T(&Q1, 2);
    emxInit_real_T(&Q2t, 2);
    emxInit_real_T(&Q, 2);

    /*  for si t u a t i o ns where A i s empty */
    qr(A, Q, Q2t);
    kAcol = Q->size[0];
    mn = A->size[1];
    i10 = Q1->size[0] * Q1->size[1];
    Q1->size[0] = kAcol;
    Q1->size[1] = mn;
    emxEnsureCapacity((emxArray__common *)Q1, i10, (int)sizeof(double));
    for (i10 = 0; i10 < mn; i10++) {
      for (i11 = 0; i11 < kAcol; i11++) {
        Q1->data[i11 + Q1->size[0] * i10] = Q->data[i11 + Q->size[0] * i10];
      }
    }

    if (A->size[1] + 1U > (unsigned int)A->size[0]) {
      i10 = 1;
      i11 = 1;
    } else {
      i10 = (int)(A->size[1] + 1U);
      i11 = A->size[0] + 1;
    }

    emxInit_real_T(&Q2, 2);
    kAcol = Q->size[0];
    i12 = Q2->size[0] * Q2->size[1];
    Q2->size[0] = kAcol;
    Q2->size[1] = i11 - i10;
    emxEnsureCapacity((emxArray__common *)Q2, i12, (int)sizeof(double));
    mn = i11 - i10;
    for (i12 = 0; i12 < mn; i12++) {
      for (i13 = 0; i13 < kAcol; i13++) {
        Q2->data[i13 + Q2->size[0] * i12] = Q->data[i13 + Q->size[0] * ((i10 +
          i12) - 1)];
      }
    }

    emxInit_real_T(&R, 2);
    kAcol = A->size[1];
    mn = Q2t->size[1];
    i12 = R->size[0] * R->size[1];
    R->size[0] = kAcol;
    R->size[1] = mn;
    emxEnsureCapacity((emxArray__common *)R, i12, (int)sizeof(double));
    for (i12 = 0; i12 < mn; i12++) {
      for (i13 = 0; i13 < kAcol; i13++) {
        R->data[i13 + R->size[0] * i12] = Q2t->data[i13 + Q2t->size[0] * i12];
      }
    }

    emxInit_real_T(&b_G, 2);
    i12 = b_G->size[0] * b_G->size[1];
    b_G->size[0] = R->size[1];
    b_G->size[1] = R->size[0];
    emxEnsureCapacity((emxArray__common *)b_G, i12, (int)sizeof(double));
    kAcol = R->size[0];
    for (i12 = 0; i12 < kAcol; i12++) {
      mn = R->size[1];
      for (i13 = 0; i13 < mn; i13++) {
        b_G->data[i13 + b_G->size[0] * i12] = R->data[i12 + R->size[0] * i13];
      }
    }

    emxInit_real_T(&py, 2);
    emxInit_real_T1(&pz, 1);
    emxInit_int32_T(&jpvt, 2);
    emxInit_real_T(&B, 2);
    emxInit_real_T(&r2, 2);
    if ((b_G->size[0] == 0) || ((b->size[0] == 0) || (b->size[1] == 0))) {
      unnamed_idx_0 = (unsigned int)b_G->size[1];
      unnamed_idx_1 = (unsigned int)b->size[1];
      i12 = py->size[0] * py->size[1];
      py->size[0] = (int)unnamed_idx_0;
      py->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)py, i12, (int)sizeof(double));
      kAcol = (int)unnamed_idx_0 * (int)unnamed_idx_1;
      for (i12 = 0; i12 < kAcol; i12++) {
        py->data[i12] = 0.0;
      }
    } else if (b_G->size[0] == b_G->size[1]) {
      n = b_G->size[1];
      i12 = Q2t->size[0] * Q2t->size[1];
      Q2t->size[0] = b_G->size[0];
      Q2t->size[1] = b_G->size[1];
      emxEnsureCapacity((emxArray__common *)Q2t, i12, (int)sizeof(double));
      kAcol = b_G->size[0] * b_G->size[1];
      for (i12 = 0; i12 < kAcol; i12++) {
        Q2t->data[i12] = b_G->data[i12];
      }

      xgetrf(b_G->size[1], b_G->size[1], Q2t, b_G->size[1], jpvt, &ar);
      i12 = py->size[0] * py->size[1];
      py->size[0] = b->size[0];
      py->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)py, i12, (int)sizeof(double));
      kAcol = b->size[0] * b->size[1];
      for (i12 = 0; i12 < kAcol; i12++) {
        py->data[i12] = b->data[i12];
      }

      nb = b->size[1];
      for (jBcol = 0; jBcol + 1 < n; jBcol++) {
        if (jpvt->data[jBcol] != jBcol + 1) {
          mn = jpvt->data[jBcol] - 1;
          for (kAcol = 0; kAcol + 1 <= nb; kAcol++) {
            wj = py->data[jBcol + py->size[0] * kAcol];
            py->data[jBcol + py->size[0] * kAcol] = py->data[mn + py->size[0] *
              kAcol];
            py->data[mn + py->size[0] * kAcol] = wj;
          }
        }
      }

      for (j = 1; j <= nb; j++) {
        jBcol = n * (j - 1);
        for (k = 0; k + 1 <= n; k++) {
          kAcol = n * k;
          if (py->data[k + jBcol] != 0.0) {
            for (ar = k + 1; ar + 1 <= n; ar++) {
              py->data[ar + jBcol] -= py->data[k + jBcol] * Q2t->data[ar + kAcol];
            }
          }
        }
      }

      i12 = r2->size[0] * r2->size[1];
      r2->size[0] = py->size[0];
      r2->size[1] = py->size[1];
      emxEnsureCapacity((emxArray__common *)r2, i12, (int)sizeof(double));
      kAcol = py->size[0] * py->size[1];
      for (i12 = 0; i12 < kAcol; i12++) {
        r2->data[i12] = py->data[i12];
      }

      xtrsm(b_G->size[1], b->size[1], Q2t, b_G->size[1], r2, b_G->size[1]);
      i12 = py->size[0] * py->size[1];
      py->size[0] = r2->size[0];
      py->size[1] = r2->size[1];
      emxEnsureCapacity((emxArray__common *)py, i12, (int)sizeof(double));
      kAcol = r2->size[0] * r2->size[1];
      for (i12 = 0; i12 < kAcol; i12++) {
        py->data[i12] = r2->data[i12];
      }
    } else {
      c_xgeqp3(b_G, pz, jpvt);
      nb = rankFromQR(b_G);
      jBcol = b_G->size[1];
      mn = b->size[1];
      i12 = py->size[0] * py->size[1];
      py->size[0] = jBcol;
      py->size[1] = mn;
      emxEnsureCapacity((emxArray__common *)py, i12, (int)sizeof(double));
      kAcol = jBcol * mn;
      for (i12 = 0; i12 < kAcol; i12++) {
        py->data[i12] = 0.0;
      }

      i12 = B->size[0] * B->size[1];
      B->size[0] = b->size[0];
      B->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)B, i12, (int)sizeof(double));
      kAcol = b->size[0] * b->size[1];
      for (i12 = 0; i12 < kAcol; i12++) {
        B->data[i12] = b->data[i12];
      }

      m = b_G->size[0];
      jBcol = b_G->size[0];
      mn = b_G->size[1];
      if (jBcol <= mn) {
        mn = jBcol;
      }

      for (j = 0; j + 1 <= mn; j++) {
        if (pz->data[j] != 0.0) {
          for (k = 0; k + 1 <= b->size[1]; k++) {
            wj = B->data[j + B->size[0] * k];
            for (ar = j + 1; ar + 1 <= m; ar++) {
              wj += b_G->data[ar + b_G->size[0] * j] * B->data[ar + B->size[0] *
                k];
            }

            wj *= pz->data[j];
            if (wj != 0.0) {
              B->data[j + B->size[0] * k] -= wj;
              for (ar = j + 1; ar + 1 <= m; ar++) {
                B->data[ar + B->size[0] * k] -= b_G->data[ar + b_G->size[0] * j]
                  * wj;
              }
            }
          }
        }
      }

      for (k = 0; k + 1 <= b->size[1]; k++) {
        for (ar = 0; ar + 1 <= nb; ar++) {
          py->data[(jpvt->data[ar] + py->size[0] * k) - 1] = B->data[ar +
            B->size[0] * k];
        }

        for (j = nb - 1; j + 1 > 0; j--) {
          py->data[(jpvt->data[j] + py->size[0] * k) - 1] /= b_G->data[j +
            b_G->size[0] * j];
          for (ar = 0; ar + 1 <= j; ar++) {
            py->data[(jpvt->data[ar] + py->size[0] * k) - 1] -= py->data
              [(jpvt->data[j] + py->size[0] * k) - 1] * b_G->data[ar + b_G->
              size[0] * j];
          }
        }
      }
    }

    emxFree_real_T(&r2);
    emxFree_real_T(&B);
    emxFree_int32_T(&jpvt);
    i12 = Q2t->size[0] * Q2t->size[1];
    Q2t->size[0] = Q2->size[1];
    Q2t->size[1] = Q2->size[0];
    emxEnsureCapacity((emxArray__common *)Q2t, i12, (int)sizeof(double));
    kAcol = Q2->size[0];
    for (i12 = 0; i12 < kAcol; i12++) {
      mn = Q2->size[1];
      for (i13 = 0; i13 < mn; i13++) {
        Q2t->data[i13 + Q2t->size[0] * i12] = Q2->data[i12 + Q2->size[0] * i13];
      }
    }

    i12 = A->size[1];
    emxInit_real_T(&y, 2);
    if ((i12 == 1) || (py->size[0] == 1)) {
      i12 = y->size[0] * y->size[1];
      y->size[0] = Q1->size[0];
      y->size[1] = py->size[1];
      emxEnsureCapacity((emxArray__common *)y, i12, (int)sizeof(double));
      kAcol = Q1->size[0];
      for (i12 = 0; i12 < kAcol; i12++) {
        mn = py->size[1];
        for (i13 = 0; i13 < mn; i13++) {
          y->data[i12 + y->size[0] * i13] = 0.0;
          jBcol = Q1->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            y->data[i12 + y->size[0] * i13] += Q1->data[i12 + Q1->size[0] * iy] *
              py->data[iy + py->size[0] * i13];
          }
        }
      }
    } else {
      i12 = A->size[1];
      i13 = Q->size[0];
      unnamed_idx_1 = (unsigned int)py->size[1];
      iy = y->size[0] * y->size[1];
      y->size[0] = i13;
      y->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(double));
      i13 = Q->size[0];
      iy = y->size[0] * y->size[1];
      emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(double));
      kAcol = y->size[1];
      for (iy = 0; iy < kAcol; iy++) {
        mn = y->size[0];
        for (m = 0; m < mn; m++) {
          y->data[m + y->size[0] * iy] = 0.0;
        }
      }

      iy = Q->size[0];
      if ((iy == 0) || (py->size[1] == 0)) {
      } else {
        iy = Q->size[0];
        jBcol = iy * (py->size[1] - 1);
        kAcol = 0;
        while ((i13 > 0) && (kAcol <= jBcol)) {
          iy = kAcol + i13;
          for (ix = kAcol; ix + 1 <= iy; ix++) {
            y->data[ix] = 0.0;
          }

          kAcol += i13;
        }

        nb = 0;
        kAcol = 0;
        while ((i13 > 0) && (kAcol <= jBcol)) {
          ar = 0;
          iy = nb + i12;
          for (b_n = nb; b_n + 1 <= iy; b_n++) {
            if (py->data[b_n] != 0.0) {
              ia = ar;
              m = kAcol + i13;
              for (ix = kAcol; ix + 1 <= m; ix++) {
                ia++;
                mn = Q->size[0];
                y->data[ix] += py->data[b_n] * Q->data[(ia - 1) % mn + Q->size[0]
                  * ((ia - 1) / mn)];
              }
            }

            ar += i13;
          }

          nb += i12;
          kAcol += i13;
        }
      }
    }

    emxInit_real_T(&C, 2);
    if ((G->size[1] == 1) || (y->size[0] == 1)) {
      i12 = C->size[0] * C->size[1];
      C->size[0] = G->size[0];
      C->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)C, i12, (int)sizeof(double));
      kAcol = G->size[0];
      for (i12 = 0; i12 < kAcol; i12++) {
        mn = y->size[1];
        for (i13 = 0; i13 < mn; i13++) {
          C->data[i12 + C->size[0] * i13] = 0.0;
          jBcol = G->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            C->data[i12 + C->size[0] * i13] += G->data[i12 + G->size[0] * iy] *
              y->data[iy + y->size[0] * i13];
          }
        }
      }
    } else {
      k = G->size[1];
      unnamed_idx_0 = (unsigned int)G->size[0];
      unnamed_idx_1 = (unsigned int)y->size[1];
      i12 = C->size[0] * C->size[1];
      C->size[0] = (int)unnamed_idx_0;
      C->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)C, i12, (int)sizeof(double));
      m = G->size[0];
      i12 = C->size[0] * C->size[1];
      emxEnsureCapacity((emxArray__common *)C, i12, (int)sizeof(double));
      kAcol = C->size[1];
      for (i12 = 0; i12 < kAcol; i12++) {
        mn = C->size[0];
        for (i13 = 0; i13 < mn; i13++) {
          C->data[i13 + C->size[0] * i12] = 0.0;
        }
      }

      if ((G->size[0] == 0) || (y->size[1] == 0)) {
      } else {
        jBcol = G->size[0] * (y->size[1] - 1);
        kAcol = 0;
        while ((m > 0) && (kAcol <= jBcol)) {
          i12 = kAcol + m;
          for (ix = kAcol; ix + 1 <= i12; ix++) {
            C->data[ix] = 0.0;
          }

          kAcol += m;
        }

        nb = 0;
        kAcol = 0;
        while ((m > 0) && (kAcol <= jBcol)) {
          ar = 0;
          i12 = nb + k;
          for (b_n = nb; b_n + 1 <= i12; b_n++) {
            if (y->data[b_n] != 0.0) {
              ia = ar;
              i13 = kAcol + m;
              for (ix = kAcol; ix + 1 <= i13; ix++) {
                ia++;
                C->data[ix] += y->data[b_n] * G->data[ia - 1];
              }
            }

            ar += m;
          }

          nb += k;
          kAcol += m;
        }
      }
    }

    emxFree_real_T(&y);
    i12 = pz->size[0];
    pz->size[0] = C->size[0];
    emxEnsureCapacity((emxArray__common *)pz, i12, (int)sizeof(double));
    kAcol = C->size[0];
    for (i12 = 0; i12 < kAcol; i12++) {
      pz->data[i12] = C->data[i12] + g->data[i12];
    }

    emxFree_real_T(&C);
    emxInit_real_T1(&gz, 1);
    if ((Q2t->size[1] == 1) || (pz->size[0] == 1)) {
      i12 = gz->size[0];
      gz->size[0] = Q2t->size[0];
      emxEnsureCapacity((emxArray__common *)gz, i12, (int)sizeof(double));
      kAcol = Q2t->size[0];
      for (i12 = 0; i12 < kAcol; i12++) {
        gz->data[i12] = 0.0;
        mn = Q2t->size[1];
        for (i13 = 0; i13 < mn; i13++) {
          gz->data[i12] += Q2t->data[i12 + Q2t->size[0] * i13] * pz->data[i13];
        }
      }
    } else {
      k = Q2t->size[1];
      unnamed_idx_0 = (unsigned int)Q2t->size[0];
      i12 = gz->size[0];
      gz->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)gz, i12, (int)sizeof(double));
      m = Q2t->size[0];
      jBcol = gz->size[0];
      i12 = gz->size[0];
      gz->size[0] = jBcol;
      emxEnsureCapacity((emxArray__common *)gz, i12, (int)sizeof(double));
      for (i12 = 0; i12 < jBcol; i12++) {
        gz->data[i12] = 0.0;
      }

      if (Q2t->size[0] == 0) {
      } else {
        kAcol = 0;
        while ((m > 0) && (kAcol <= 0)) {
          for (ix = 1; ix <= m; ix++) {
            gz->data[ix - 1] = 0.0;
          }

          kAcol = m;
        }

        nb = 0;
        kAcol = 0;
        while ((m > 0) && (kAcol <= 0)) {
          ar = 0;
          i12 = nb + k;
          for (b_n = nb; b_n + 1 <= i12; b_n++) {
            if (pz->data[b_n] != 0.0) {
              ia = ar;
              for (ix = 0; ix + 1 <= m; ix++) {
                ia++;
                gz->data[ix] += pz->data[b_n] * Q2t->data[ia - 1];
              }
            }

            ar += m;
          }

          nb += k;
          kAcol = m;
        }
      }
    }

    emxInit_real_T(&b_y, 2);
    if ((Q2t->size[1] == 1) || (G->size[0] == 1)) {
      i12 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = Q2t->size[0];
      b_y->size[1] = G->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, i12, (int)sizeof(double));
      kAcol = Q2t->size[0];
      for (i12 = 0; i12 < kAcol; i12++) {
        mn = G->size[1];
        for (i13 = 0; i13 < mn; i13++) {
          b_y->data[i12 + b_y->size[0] * i13] = 0.0;
          jBcol = Q2t->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            b_y->data[i12 + b_y->size[0] * i13] += Q2t->data[i12 + Q2t->size[0] *
              iy] * G->data[iy + G->size[0] * i13];
          }
        }
      }
    } else {
      k = Q2t->size[1];
      unnamed_idx_0 = (unsigned int)Q2t->size[0];
      unnamed_idx_1 = (unsigned int)G->size[1];
      i12 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = (int)unnamed_idx_0;
      b_y->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)b_y, i12, (int)sizeof(double));
      m = Q2t->size[0];
      i12 = b_y->size[0] * b_y->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, i12, (int)sizeof(double));
      kAcol = b_y->size[1];
      for (i12 = 0; i12 < kAcol; i12++) {
        mn = b_y->size[0];
        for (i13 = 0; i13 < mn; i13++) {
          b_y->data[i13 + b_y->size[0] * i12] = 0.0;
        }
      }

      if ((Q2t->size[0] == 0) || (G->size[1] == 0)) {
      } else {
        jBcol = Q2t->size[0] * (G->size[1] - 1);
        kAcol = 0;
        while ((m > 0) && (kAcol <= jBcol)) {
          i12 = kAcol + m;
          for (ix = kAcol; ix + 1 <= i12; ix++) {
            b_y->data[ix] = 0.0;
          }

          kAcol += m;
        }

        nb = 0;
        kAcol = 0;
        while ((m > 0) && (kAcol <= jBcol)) {
          ar = 0;
          i12 = nb + k;
          for (b_n = nb; b_n + 1 <= i12; b_n++) {
            if (G->data[b_n] != 0.0) {
              ia = ar;
              i13 = kAcol + m;
              for (ix = kAcol; ix + 1 <= i13; ix++) {
                ia++;
                b_y->data[ix] += G->data[b_n] * Q2t->data[ia - 1];
              }
            }

            ar += m;
          }

          nb += k;
          kAcol += m;
        }
      }
    }

    emxInit_real_T(&Gz, 2);
    guard1 = false;
    if (b_y->size[1] == 1) {
      guard1 = true;
    } else {
      i12 = Q->size[0];
      if (i12 == 1) {
        guard1 = true;
      } else {
        k = b_y->size[1];
        unnamed_idx_0 = (unsigned int)b_y->size[0];
        i12 = Gz->size[0] * Gz->size[1];
        Gz->size[0] = (int)unnamed_idx_0;
        Gz->size[1] = i11 - i10;
        emxEnsureCapacity((emxArray__common *)Gz, i12, (int)sizeof(double));
        m = b_y->size[0];
        i12 = Gz->size[0] * Gz->size[1];
        emxEnsureCapacity((emxArray__common *)Gz, i12, (int)sizeof(double));
        kAcol = Gz->size[1];
        for (i12 = 0; i12 < kAcol; i12++) {
          mn = Gz->size[0];
          for (i13 = 0; i13 < mn; i13++) {
            Gz->data[i13 + Gz->size[0] * i12] = 0.0;
          }
        }

        if ((b_y->size[0] == 0) || (i11 - i10 == 0)) {
        } else {
          jBcol = b_y->size[0] * ((i11 - i10) - 1);
          kAcol = 0;
          while ((m > 0) && (kAcol <= jBcol)) {
            i12 = kAcol + m;
            for (ix = kAcol; ix + 1 <= i12; ix++) {
              Gz->data[ix] = 0.0;
            }

            kAcol += m;
          }

          nb = 0;
          kAcol = 0;
          while ((m > 0) && (kAcol <= jBcol)) {
            ar = 0;
            i12 = nb + k;
            for (b_n = nb; b_n + 1 <= i12; b_n++) {
              i13 = Q->size[0];
              if (Q->data[b_n % i13 + Q->size[0] * ((i10 + b_n / i13) - 1)] !=
                  0.0) {
                ia = ar;
                i13 = kAcol + m;
                for (ix = kAcol; ix + 1 <= i13; ix++) {
                  ia++;
                  iy = Q->size[0];
                  Gz->data[ix] += Q->data[b_n % iy + Q->size[0] * ((i10 + b_n /
                    iy) - 1)] * b_y->data[ia - 1];
                }
              }

              ar += m;
            }

            nb += k;
            kAcol += m;
          }
        }
      }
    }

    if (guard1) {
      i12 = Gz->size[0] * Gz->size[1];
      Gz->size[0] = b_y->size[0];
      Gz->size[1] = Q2->size[1];
      emxEnsureCapacity((emxArray__common *)Gz, i12, (int)sizeof(double));
      kAcol = b_y->size[0];
      for (i12 = 0; i12 < kAcol; i12++) {
        mn = Q2->size[1];
        for (i13 = 0; i13 < mn; i13++) {
          Gz->data[i12 + Gz->size[0] * i13] = 0.0;
          jBcol = b_y->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            Gz->data[i12 + Gz->size[0] * i13] += b_y->data[i12 + b_y->size[0] *
              iy] * Q2->data[iy + Q2->size[0] * i13];
          }
        }
      }
    }

    emxFree_real_T(&b_y);
    i12 = b_G->size[0] * b_G->size[1];
    b_G->size[0] = Gz->size[0];
    b_G->size[1] = Gz->size[1];
    emxEnsureCapacity((emxArray__common *)b_G, i12, (int)sizeof(double));
    kAcol = Gz->size[0] * Gz->size[1];
    for (i12 = 0; i12 < kAcol; i12++) {
      b_G->data[i12] = Gz->data[i12];
    }

    n = Gz->size[1];
    if (Gz->size[1] == 0) {
    } else {
      b_n = Gz->size[0];
      ar = -1;
      if (Gz->size[0] == 0) {
      } else {
        m = 0;
        j = 0;
        exitg1 = false;
        while ((!exitg1) && (j + 1 <= b_n)) {
          mn = m + j;
          wj = 0.0;
          if (j < 1) {
          } else {
            ix = m;
            iy = m;
            for (k = 1; k <= j; k++) {
              wj += b_G->data[ix] * b_G->data[iy];
              ix++;
              iy++;
            }
          }

          wj = b_G->data[mn] - wj;
          if (wj > 0.0) {
            wj = sqrt(wj);
            b_G->data[mn] = wj;
            if (j + 1 < b_n) {
              nb = (b_n - j) - 2;
              mn += b_n;
              jBcol = (m + b_n) + 1;
              if ((j == 0) || (nb + 1 == 0)) {
              } else {
                iy = mn;
                i12 = jBcol + n * nb;
                for (kAcol = jBcol; kAcol <= i12; kAcol += n) {
                  ix = m + 1;
                  c = 0.0;
                  i13 = (kAcol + j) - 1;
                  for (ia = kAcol; ia <= i13; ia++) {
                    c += b_G->data[ia - 1] * b_G->data[ix - 1];
                    ix++;
                  }

                  b_G->data[iy] += -c;
                  iy += n;
                }
              }

              wj = 1.0 / wj;
              i12 = (mn + b_n * nb) + 1;
              while (mn + 1 <= i12) {
                b_G->data[mn] *= wj;
                mn += b_n;
              }

              m = jBcol - 1;
            }

            j++;
          } else {
            b_G->data[mn] = wj;
            ar = j;
            exitg1 = true;
          }
        }
      }

      if (ar + 1 == 0) {
        jBcol = Gz->size[1];
      } else {
        jBcol = ar;
      }

      for (j = 0; j + 1 <= jBcol; j++) {
        for (ar = j + 1; ar + 1 <= jBcol; ar++) {
          b_G->data[ar + b_G->size[0] * j] = 0.0;
        }
      }
    }

    emxFree_real_T(&Gz);
    i12 = Q2t->size[0] * Q2t->size[1];
    Q2t->size[0] = b_G->size[1];
    Q2t->size[1] = b_G->size[0];
    emxEnsureCapacity((emxArray__common *)Q2t, i12, (int)sizeof(double));
    kAcol = b_G->size[0];
    for (i12 = 0; i12 < kAcol; i12++) {
      mn = b_G->size[1];
      for (i13 = 0; i13 < mn; i13++) {
        Q2t->data[i13 + Q2t->size[0] * i12] = b_G->data[i12 + b_G->size[0] * i13];
      }
    }

    emxFree_real_T(&b_G);
    emxInit_real_T1(&b_gz, 1);
    i12 = b_gz->size[0];
    b_gz->size[0] = gz->size[0];
    emxEnsureCapacity((emxArray__common *)b_gz, i12, (int)sizeof(double));
    kAcol = gz->size[0];
    for (i12 = 0; i12 < kAcol; i12++) {
      b_gz->data[i12] = -gz->data[i12];
    }

    emxFree_real_T(&gz);
    emxInit_real_T(&b_Q2t, 2);
    b_mldivide(Q2t, b_gz, pz);
    i12 = b_Q2t->size[0] * b_Q2t->size[1];
    b_Q2t->size[0] = Q2t->size[1];
    b_Q2t->size[1] = Q2t->size[0];
    emxEnsureCapacity((emxArray__common *)b_Q2t, i12, (int)sizeof(double));
    kAcol = Q2t->size[0];
    emxFree_real_T(&b_gz);
    for (i12 = 0; i12 < kAcol; i12++) {
      mn = Q2t->size[1];
      for (i13 = 0; i13 < mn; i13++) {
        b_Q2t->data[i13 + b_Q2t->size[0] * i12] = Q2t->data[i12 + Q2t->size[0] *
          i13];
      }
    }

    emxInit_real_T1(&b_pz, 1);
    i12 = b_pz->size[0];
    b_pz->size[0] = pz->size[0];
    emxEnsureCapacity((emxArray__common *)b_pz, i12, (int)sizeof(double));
    kAcol = pz->size[0];
    for (i12 = 0; i12 < kAcol; i12++) {
      b_pz->data[i12] = pz->data[i12];
    }

    b_mldivide(b_Q2t, b_pz, pz);
    emxFree_real_T(&b_pz);
    emxFree_real_T(&b_Q2t);
    if ((i11 - i10 == 1) || (pz->size[0] == 1)) {
      i10 = x->size[0];
      x->size[0] = Q2->size[0];
      emxEnsureCapacity((emxArray__common *)x, i10, (int)sizeof(double));
      kAcol = Q2->size[0];
      for (i10 = 0; i10 < kAcol; i10++) {
        x->data[i10] = 0.0;
        mn = Q2->size[1];
        for (i11 = 0; i11 < mn; i11++) {
          x->data[i10] += Q2->data[i10 + Q2->size[0] * i11] * pz->data[i11];
        }
      }
    } else {
      k = i11 - i10;
      i11 = Q->size[0];
      i12 = x->size[0];
      x->size[0] = i11;
      emxEnsureCapacity((emxArray__common *)x, i12, (int)sizeof(double));
      i11 = Q->size[0];
      jBcol = x->size[0];
      i12 = x->size[0];
      x->size[0] = jBcol;
      emxEnsureCapacity((emxArray__common *)x, i12, (int)sizeof(double));
      for (i12 = 0; i12 < jBcol; i12++) {
        x->data[i12] = 0.0;
      }

      i12 = Q->size[0];
      if (i12 == 0) {
      } else {
        kAcol = 0;
        while ((i11 > 0) && (kAcol <= 0)) {
          for (ix = 1; ix <= i11; ix++) {
            x->data[ix - 1] = 0.0;
          }

          kAcol = i11;
        }

        nb = 0;
        kAcol = 0;
        while ((i11 > 0) && (kAcol <= 0)) {
          ar = 0;
          i12 = nb + k;
          for (b_n = nb; b_n + 1 <= i12; b_n++) {
            if (pz->data[b_n] != 0.0) {
              ia = ar;
              for (ix = 0; ix + 1 <= i11; ix++) {
                ia++;
                i13 = Q->size[0];
                x->data[ix] += pz->data[b_n] * Q->data[(ia - 1) % i13 + Q->size
                  [0] * ((i10 + (ia - 1) / i13) - 1)];
              }
            }

            ar += i11;
          }

          nb += k;
          kAcol = i11;
        }
      }
    }

    emxFree_real_T(&Q2);
    i10 = A->size[1];
    emxInit_real_T(&b_C, 2);
    if ((i10 == 1) || (py->size[0] == 1)) {
      i10 = b_C->size[0] * b_C->size[1];
      b_C->size[0] = Q1->size[0];
      b_C->size[1] = py->size[1];
      emxEnsureCapacity((emxArray__common *)b_C, i10, (int)sizeof(double));
      kAcol = Q1->size[0];
      for (i10 = 0; i10 < kAcol; i10++) {
        mn = py->size[1];
        for (i11 = 0; i11 < mn; i11++) {
          b_C->data[i10 + b_C->size[0] * i11] = 0.0;
          jBcol = Q1->size[1];
          for (i12 = 0; i12 < jBcol; i12++) {
            b_C->data[i10 + b_C->size[0] * i11] += Q1->data[i10 + Q1->size[0] *
              i12] * py->data[i12 + py->size[0] * i11];
          }
        }
      }
    } else {
      i10 = A->size[1];
      i11 = Q->size[0];
      unnamed_idx_1 = (unsigned int)py->size[1];
      i12 = b_C->size[0] * b_C->size[1];
      b_C->size[0] = i11;
      b_C->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)b_C, i12, (int)sizeof(double));
      i11 = Q->size[0];
      i12 = b_C->size[0] * b_C->size[1];
      emxEnsureCapacity((emxArray__common *)b_C, i12, (int)sizeof(double));
      kAcol = b_C->size[1];
      for (i12 = 0; i12 < kAcol; i12++) {
        mn = b_C->size[0];
        for (i13 = 0; i13 < mn; i13++) {
          b_C->data[i13 + b_C->size[0] * i12] = 0.0;
        }
      }

      i12 = Q->size[0];
      if ((i12 == 0) || (py->size[1] == 0)) {
      } else {
        i12 = Q->size[0];
        jBcol = i12 * (py->size[1] - 1);
        kAcol = 0;
        while ((i11 > 0) && (kAcol <= jBcol)) {
          i12 = kAcol + i11;
          for (ix = kAcol; ix + 1 <= i12; ix++) {
            b_C->data[ix] = 0.0;
          }

          kAcol += i11;
        }

        nb = 0;
        kAcol = 0;
        while ((i11 > 0) && (kAcol <= jBcol)) {
          ar = 0;
          i12 = nb + i10;
          for (b_n = nb; b_n + 1 <= i12; b_n++) {
            if (py->data[b_n] != 0.0) {
              ia = ar;
              i13 = kAcol + i11;
              for (ix = kAcol; ix + 1 <= i13; ix++) {
                ia++;
                iy = Q->size[0];
                b_C->data[ix] += py->data[b_n] * Q->data[(ia - 1) % iy + Q->
                  size[0] * ((ia - 1) / iy)];
              }
            }

            ar += i11;
          }

          nb += i10;
          kAcol += i11;
        }
      }
    }

    emxFree_real_T(&Q);
    emxFree_real_T(&py);
    i10 = x->size[0];
    x->size[0] = b_C->size[0];
    emxEnsureCapacity((emxArray__common *)x, i10, (int)sizeof(double));
    kAcol = b_C->size[0];
    for (i10 = 0; i10 < kAcol; i10++) {
      x->data[i10] += b_C->data[i10];
    }

    emxFree_real_T(&b_C);
    emxInit_real_T1(&c_C, 1);
    if ((G->size[1] == 1) || (x->size[0] == 1)) {
      i10 = c_C->size[0];
      c_C->size[0] = G->size[0];
      emxEnsureCapacity((emxArray__common *)c_C, i10, (int)sizeof(double));
      kAcol = G->size[0];
      for (i10 = 0; i10 < kAcol; i10++) {
        c_C->data[i10] = 0.0;
        mn = G->size[1];
        for (i11 = 0; i11 < mn; i11++) {
          c_C->data[i10] += G->data[i10 + G->size[0] * i11] * x->data[i11];
        }
      }
    } else {
      k = G->size[1];
      unnamed_idx_0 = (unsigned int)G->size[0];
      i10 = c_C->size[0];
      c_C->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)c_C, i10, (int)sizeof(double));
      m = G->size[0];
      jBcol = c_C->size[0];
      i10 = c_C->size[0];
      c_C->size[0] = jBcol;
      emxEnsureCapacity((emxArray__common *)c_C, i10, (int)sizeof(double));
      for (i10 = 0; i10 < jBcol; i10++) {
        c_C->data[i10] = 0.0;
      }

      if (G->size[0] == 0) {
      } else {
        kAcol = 0;
        while ((m > 0) && (kAcol <= 0)) {
          for (ix = 1; ix <= m; ix++) {
            c_C->data[ix - 1] = 0.0;
          }

          kAcol = m;
        }

        nb = 0;
        kAcol = 0;
        while ((m > 0) && (kAcol <= 0)) {
          ar = 0;
          i10 = nb + k;
          for (b_n = nb; b_n + 1 <= i10; b_n++) {
            if (x->data[b_n] != 0.0) {
              ia = ar;
              for (ix = 0; ix + 1 <= m; ix++) {
                ia++;
                c_C->data[ix] += x->data[b_n] * G->data[ia - 1];
              }
            }

            ar += m;
          }

          nb += k;
          kAcol = m;
        }
      }
    }

    i10 = Q2t->size[0] * Q2t->size[1];
    Q2t->size[0] = Q1->size[1];
    Q2t->size[1] = Q1->size[0];
    emxEnsureCapacity((emxArray__common *)Q2t, i10, (int)sizeof(double));
    kAcol = Q1->size[0];
    for (i10 = 0; i10 < kAcol; i10++) {
      mn = Q1->size[1];
      for (i11 = 0; i11 < mn; i11++) {
        Q2t->data[i11 + Q2t->size[0] * i10] = Q1->data[i10 + Q1->size[0] * i11];
      }
    }

    emxFree_real_T(&Q1);
    i10 = c_C->size[0];
    emxEnsureCapacity((emxArray__common *)c_C, i10, (int)sizeof(double));
    kAcol = c_C->size[0];
    for (i10 = 0; i10 < kAcol; i10++) {
      c_C->data[i10] += g->data[i10];
    }

    emxInit_real_T1(&c_y, 1);
    if ((Q2t->size[1] == 1) || (c_C->size[0] == 1)) {
      i10 = c_y->size[0];
      c_y->size[0] = Q2t->size[0];
      emxEnsureCapacity((emxArray__common *)c_y, i10, (int)sizeof(double));
      kAcol = Q2t->size[0];
      for (i10 = 0; i10 < kAcol; i10++) {
        c_y->data[i10] = 0.0;
        mn = Q2t->size[1];
        for (i11 = 0; i11 < mn; i11++) {
          c_y->data[i10] += Q2t->data[i10 + Q2t->size[0] * i11] * c_C->data[i11];
        }
      }
    } else {
      k = Q2t->size[1];
      unnamed_idx_0 = (unsigned int)Q2t->size[0];
      i10 = c_y->size[0];
      c_y->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)c_y, i10, (int)sizeof(double));
      m = Q2t->size[0];
      jBcol = c_y->size[0];
      i10 = c_y->size[0];
      c_y->size[0] = jBcol;
      emxEnsureCapacity((emxArray__common *)c_y, i10, (int)sizeof(double));
      for (i10 = 0; i10 < jBcol; i10++) {
        c_y->data[i10] = 0.0;
      }

      kAcol = 0;
      while (kAcol <= 0) {
        for (ix = 1; ix <= m; ix++) {
          c_y->data[ix - 1] = 0.0;
        }

        kAcol = m;
      }

      nb = 0;
      kAcol = 0;
      while (kAcol <= 0) {
        ar = 0;
        i10 = nb + k;
        for (b_n = nb; b_n + 1 <= i10; b_n++) {
          if (c_C->data[b_n] != 0.0) {
            ia = ar;
            for (ix = 0; ix + 1 <= m; ix++) {
              ia++;
              c_y->data[ix] += c_C->data[b_n] * Q2t->data[ia - 1];
            }
          }

          ar += m;
        }

        nb += k;
        kAcol = m;
      }
    }

    emxFree_real_T(&c_C);
    emxFree_real_T(&Q2t);
    b_mldivide(R, c_y, pz);
    i10 = u->size[0] * u->size[1];
    u->size[0] = pz->size[0];
    u->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)u, i10, (int)sizeof(double));
    kAcol = pz->size[0];
    emxFree_real_T(&c_y);
    emxFree_real_T(&R);
    for (i10 = 0; i10 < kAcol; i10++) {
      u->data[i10] = pz->data[i10];
    }

    emxFree_real_T(&pz);
  } else {
    emxInit_real_T(&b_G, 2);
    i10 = b_G->size[0] * b_G->size[1];
    b_G->size[0] = G->size[0];
    b_G->size[1] = G->size[1];
    emxEnsureCapacity((emxArray__common *)b_G, i10, (int)sizeof(double));
    kAcol = G->size[0] * G->size[1];
    for (i10 = 0; i10 < kAcol; i10++) {
      b_G->data[i10] = -G->data[i10];
    }

    b_mldivide(b_G, g, x);
    i10 = u->size[0] * u->size[1];
    u->size[0] = 0;
    u->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)u, i10, (int)sizeof(double));
    emxFree_real_T(&b_G);
  }
}

/*
 * %
 *  min f(x)=1/2*x'*B*x - b'*x  %NB
 *  C*x = c
 *  A * x <= d %NB
 *  assert (M_A <= 1000);
 * Arguments    : const emxArray_real_T *B
 *                const emxArray_real_T *b
 *                const emxArray_real_T *C
 *                const emxArray_real_T *c
 *                const emxArray_real_T *A
 *                const emxArray_real_T *d
 *                emxArray_real_T *x0
 * Return Type  : double
 */
double QPACM(const emxArray_real_T *B, const emxArray_real_T *b, const
             emxArray_real_T *C, const emxArray_real_T *c, const emxArray_real_T
             *A, const emxArray_real_T *d, emxArray_real_T *x0)
{
  double flag;
  emxArray_real_T *W_A;
  emxArray_real_T *W_d;
  emxArray_real_T *CW_A;
  int i30;
  int br;
  emxArray_real_T *CW_d;
  int m_C;
  int i;
  emxArray_real_T *p;
  emxArray_real_T *r;
  emxArray_real_T *I;
  emxArray_real_T *temp;
  emxArray_real_T *b_C;
  emxArray_real_T *c_C;
  emxArray_real_T *y;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  emxArray_real_T *a;
  emxArray_real_T *b_y;
  emxArray_real_T *varargin_2;
  emxArray_real_T *d_C;
  emxArray_real_T *e_C;
  emxArray_real_T *r4;
  emxArray_real_T *f_C;
  emxArray_real_T *b_c;
  emxArray_real_T *b_CW_A;
  emxArray_real_T *b_W_A;
  emxArray_real_T *c_W_A;
  emxArray_real_T *c_CW_A;
  int exitg1;
  int k;
  unsigned int B_idx_0;
  int m;
  int ia;
  int ar;
  int i31;
  boolean_T empty_non_axis_sizes;
  int ic;
  int idx;
  unsigned int uv0[2];
  double mtmp;
  boolean_T exitg3;
  boolean_T exitg2;
  boolean_T exitg5;
  boolean_T guard1 = false;
  double alpha;
  int Imin_alpha;
  boolean_T exitg4;
  emxInit_real_T(&W_A, 2);
  emxInit_real_T(&W_d, 2);
  emxInit_real_T(&CW_A, 2);

  /* !!!!!!!!!very important */
  flag = 0.0;

  /*  W.A = A([1],:); */
  /*  W.d = d([1]); */
  i30 = W_A->size[0] * W_A->size[1];
  W_A->size[0] = 0;
  W_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)W_A, i30, (int)sizeof(double));

  /* empty */
  i30 = W_d->size[0] * W_d->size[1];
  W_d->size[0] = 0;
  W_d->size[1] = 0;
  emxEnsureCapacity((emxArray__common *)W_d, i30, (int)sizeof(double));

  /*  */
  /*  [CW.A,IA] = setdiff(A,W.A,'rows','stable');%set1 - set2 */
  i30 = CW_A->size[0] * CW_A->size[1];
  CW_A->size[0] = A->size[0];
  CW_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)CW_A, i30, (int)sizeof(double));
  br = A->size[0] * A->size[1];
  for (i30 = 0; i30 < br; i30++) {
    CW_A->data[i30] = A->data[i30];
  }

  emxInit_real_T1(&CW_d, 1);
  i30 = CW_d->size[0];
  CW_d->size[0] = d->size[0];
  emxEnsureCapacity((emxArray__common *)CW_d, i30, (int)sizeof(double));
  br = d->size[0];
  for (i30 = 0; i30 < br; i30++) {
    CW_d->data[i30] = d->data[i30];
  }

  /*  CW.d = d(IA); */
  m_C = C->size[0];
  i = 0;
  emxInit_real_T1(&p, 1);
  emxInit_real_T(&r, 2);
  emxInit_real_T1(&I, 1);
  emxInit_real_T1(&temp, 1);
  emxInit_real_T1(&b_C, 1);
  emxInit_real_T1(&c_C, 1);
  emxInit_real_T(&y, 2);
  emxInit_boolean_T(&x, 1);
  emxInit_int32_T1(&ii, 1);
  emxInit_real_T(&a, 2);
  emxInit_real_T1(&b_y, 1);
  emxInit_real_T1(&varargin_2, 1);
  emxInit_real_T1(&d_C, 1);
  emxInit_real_T(&e_C, 2);
  emxInit_real_T(&r4, 2);
  emxInit_real_T1(&f_C, 1);
  emxInit_real_T(&b_c, 2);
  emxInit_real_T(&b_CW_A, 2);
  emxInit_real_T(&b_W_A, 2);
  emxInit_real_T(&c_W_A, 2);
  emxInit_real_T(&c_CW_A, 2);
  do {
    exitg1 = 0;
    if (i < 200) {
      i++;
      flag = i;

      /*  check_rank = (rank([C;W.A])==size([C;W.A],1)); */
      if ((B->size[1] == 1) || (x0->size[0] == 1)) {
        i30 = b_C->size[0];
        b_C->size[0] = B->size[0];
        emxEnsureCapacity((emxArray__common *)b_C, i30, (int)sizeof(double));
        br = B->size[0];
        for (i30 = 0; i30 < br; i30++) {
          b_C->data[i30] = 0.0;
          ar = B->size[1];
          for (i31 = 0; i31 < ar; i31++) {
            b_C->data[i30] += B->data[i30 + B->size[0] * i31] * x0->data[i31];
          }
        }
      } else {
        k = B->size[1];
        B_idx_0 = (unsigned int)B->size[0];
        i30 = b_C->size[0];
        b_C->size[0] = (int)B_idx_0;
        emxEnsureCapacity((emxArray__common *)b_C, i30, (int)sizeof(double));
        m = B->size[0];
        ia = b_C->size[0];
        i30 = b_C->size[0];
        b_C->size[0] = ia;
        emxEnsureCapacity((emxArray__common *)b_C, i30, (int)sizeof(double));
        for (i30 = 0; i30 < ia; i30++) {
          b_C->data[i30] = 0.0;
        }

        if (B->size[0] == 0) {
        } else {
          ia = 0;
          while ((m > 0) && (ia <= 0)) {
            for (ic = 1; ic <= m; ic++) {
              b_C->data[ic - 1] = 0.0;
            }

            ia = m;
          }

          br = 0;
          ia = 0;
          while ((m > 0) && (ia <= 0)) {
            ar = 0;
            i30 = br + k;
            for (idx = br; idx + 1 <= i30; idx++) {
              if (x0->data[idx] != 0.0) {
                ia = ar;
                for (ic = 0; ic + 1 <= m; ic++) {
                  ia++;
                  b_C->data[ic] += x0->data[idx] * B->data[ia - 1];
                }
              }

              ar += m;
            }

            br += k;
            ia = m;
          }
        }
      }

      if (!((C->size[0] == 0) || (C->size[1] == 0))) {
        ia = C->size[1];
      } else if (!((W_A->size[0] == 0) || (W_A->size[1] == 0))) {
        ia = W_A->size[1];
      } else {
        ia = C->size[1];
        if (W_A->size[1] > C->size[1]) {
          ia = W_A->size[1];
        }
      }

      empty_non_axis_sizes = (ia == 0);
      if (empty_non_axis_sizes || (!((C->size[0] == 0) || (C->size[1] == 0)))) {
        ic = C->size[0];
      } else {
        ic = 0;
      }

      if (empty_non_axis_sizes || (!((W_A->size[0] == 0) || (W_A->size[1] == 0))))
      {
        ar = W_A->size[0];
      } else {
        ar = 0;
      }

      if (!((W_d->size[0] == 0) || (W_d->size[1] == 0))) {
        idx = W_d->size[0];
      } else {
        idx = 0;
      }

      i30 = d_C->size[0];
      d_C->size[0] = b_C->size[0];
      emxEnsureCapacity((emxArray__common *)d_C, i30, (int)sizeof(double));
      br = b_C->size[0];
      for (i30 = 0; i30 < br; i30++) {
        d_C->data[i30] = b_C->data[i30] - b->data[i30];
      }

      i30 = e_C->size[0] * e_C->size[1];
      e_C->size[0] = ic + ar;
      e_C->size[1] = ia;
      emxEnsureCapacity((emxArray__common *)e_C, i30, (int)sizeof(double));
      for (i30 = 0; i30 < ia; i30++) {
        for (i31 = 0; i31 < ic; i31++) {
          e_C->data[i31 + e_C->size[0] * i30] = C->data[i31 + ic * i30];
        }
      }

      for (i30 = 0; i30 < ia; i30++) {
        for (i31 = 0; i31 < ar; i31++) {
          e_C->data[(i31 + ic) + e_C->size[0] * i30] = W_A->data[i31 + ar * i30];
        }
      }

      ar = c->size[0];
      i30 = b_c->size[0] * b_c->size[1];
      b_c->size[0] = ar + idx;
      b_c->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_c, i30, (int)sizeof(double));
      for (i30 = 0; i30 < 1; i30++) {
        for (i31 = 0; i31 < ar; i31++) {
          b_c->data[i31] = c->data[i31];
        }
      }

      for (i30 = 0; i30 < 1; i30++) {
        for (i31 = 0; i31 < idx; i31++) {
          b_c->data[i31 + ar] = W_d->data[i31];
        }
      }

      i30 = r4->size[0] * r4->size[1];
      r4->size[0] = b_c->size[0];
      r4->size[1] = b_c->size[1];
      emxEnsureCapacity((emxArray__common *)r4, i30, (int)sizeof(double));
      br = b_c->size[1];
      for (i30 = 0; i30 < br; i30++) {
        ar = b_c->size[0];
        for (i31 = 0; i31 < ar; i31++) {
          r4->data[i31 + r4->size[0] * i30] = 0.0 * b_c->data[i31 + b_c->size[0]
            * i30];
        }
      }

      null_space(B, d_C, e_C, r4, p, r);

      /* null space is better */
      c_abs(p, varargin_2);
      i30 = x->size[0];
      x->size[0] = varargin_2->size[0];
      emxEnsureCapacity((emxArray__common *)x, i30, (int)sizeof(boolean_T));
      br = varargin_2->size[0];
      for (i30 = 0; i30 < br; i30++) {
        x->data[i30] = (varargin_2->data[i30] > 1.0E-8);
      }

      i30 = p->size[0];
      emxEnsureCapacity((emxArray__common *)p, i30, (int)sizeof(double));
      br = p->size[0];
      for (i30 = 0; i30 < br; i30++) {
        p->data[i30] *= (double)x->data[i30];
      }

      for (i30 = 0; i30 < 2; i30++) {
        uv0[i30] = (unsigned int)r->size[i30];
      }

      i30 = y->size[0] * y->size[1];
      y->size[0] = (int)uv0[0];
      y->size[1] = (int)uv0[1];
      emxEnsureCapacity((emxArray__common *)y, i30, (int)sizeof(double));
      ia = r->size[0] * r->size[1];
      for (k = 0; k + 1 <= ia; k++) {
        y->data[k] = fabs(r->data[k]);
      }

      i30 = r->size[0] * r->size[1];
      emxEnsureCapacity((emxArray__common *)r, i30, (int)sizeof(double));
      ia = r->size[0];
      ar = r->size[1];
      br = ia * ar;
      for (i30 = 0; i30 < br; i30++) {
        r->data[i30] *= (double)(y->data[i30] > 1.0E-8);
      }

      /*  lambda = r(1:m_C); */
      i30 = r->size[0] * r->size[1];
      if (m_C + 1U > (unsigned int)i30) {
        i31 = 0;
        i30 = 1;
      } else {
        i31 = m_C;
        i30++;
      }

      c_abs(p, varargin_2);
      if (sum(varargin_2) == 0.0) {
        ar = 1;
        ia = (i30 - i31) - 1;
        mtmp = r->data[i31];
        if ((i30 - i31) - 1 > 1) {
          if (rtIsNaN(mtmp)) {
            idx = 2;
            exitg3 = false;
            while ((!exitg3) && (idx <= ia)) {
              ar = idx;
              if (!rtIsNaN(r->data[(i31 + idx) - 1])) {
                mtmp = r->data[(i31 + idx) - 1];
                exitg3 = true;
              } else {
                idx++;
              }
            }
          }

          if (ar < (i30 - i31) - 1) {
            while (ar + 1 <= ia) {
              if (r->data[i31 + ar] < mtmp) {
                mtmp = r->data[i31 + ar];
              }

              ar++;
            }
          }
        }

        if (mtmp >= 0.0) {
          exitg1 = 1;
        } else {
          ar = 1;
          ia = (i30 - i31) - 1;
          mtmp = r->data[i31];
          k = 0;
          if ((i30 - i31) - 1 > 1) {
            if (rtIsNaN(mtmp)) {
              idx = 2;
              exitg2 = false;
              while ((!exitg2) && (idx <= ia)) {
                ar = idx;
                if (!rtIsNaN(r->data[(i31 + idx) - 1])) {
                  mtmp = r->data[(i31 + idx) - 1];
                  k = idx - 1;
                  exitg2 = true;
                } else {
                  idx++;
                }
              }
            }

            if (ar < (i30 - i31) - 1) {
              while (ar + 1 <= ia) {
                if (r->data[i31 + ar] < mtmp) {
                  mtmp = r->data[i31 + ar];
                  k = ar;
                }

                ar++;
              }
            }
          }

          if (!((CW_A->size[0] == 0) || (CW_A->size[1] == 0))) {
            ia = CW_A->size[1];
          } else {
            ia = W_A->size[1];
          }

          if ((ia == 0) || (!((CW_A->size[0] == 0) || (CW_A->size[1] == 0)))) {
            ic = CW_A->size[0];
          } else {
            ic = 0;
          }

          br = W_A->size[1];
          i30 = c_W_A->size[0] * c_W_A->size[1];
          c_W_A->size[0] = 1;
          c_W_A->size[1] = br;
          emxEnsureCapacity((emxArray__common *)c_W_A, i30, (int)sizeof(double));
          for (i30 = 0; i30 < br; i30++) {
            c_W_A->data[c_W_A->size[0] * i30] = W_A->data[k + W_A->size[0] * i30];
          }

          i30 = c_CW_A->size[0] * c_CW_A->size[1];
          c_CW_A->size[0] = ic + 1;
          c_CW_A->size[1] = ia;
          emxEnsureCapacity((emxArray__common *)c_CW_A, i30, (int)sizeof(double));
          for (i30 = 0; i30 < ia; i30++) {
            for (i31 = 0; i31 < ic; i31++) {
              c_CW_A->data[i31 + c_CW_A->size[0] * i30] = CW_A->data[i31 + ic *
                i30];
            }
          }

          for (i30 = 0; i30 < ia; i30++) {
            for (i31 = 0; i31 < 1; i31++) {
              c_CW_A->data[ic + c_CW_A->size[0] * i30] = c_W_A->data[i30];
            }
          }

          i30 = CW_A->size[0] * CW_A->size[1];
          CW_A->size[0] = c_CW_A->size[0];
          CW_A->size[1] = c_CW_A->size[1];
          emxEnsureCapacity((emxArray__common *)CW_A, i30, (int)sizeof(double));
          br = c_CW_A->size[1];
          for (i30 = 0; i30 < br; i30++) {
            ar = c_CW_A->size[0];
            for (i31 = 0; i31 < ar; i31++) {
              CW_A->data[i31 + CW_A->size[0] * i30] = c_CW_A->data[i31 +
                c_CW_A->size[0] * i30];
            }
          }

          ar = CW_d->size[0];
          i30 = CW_d->size[0];
          CW_d->size[0] = ar + 1;
          emxEnsureCapacity((emxArray__common *)CW_d, i30, (int)sizeof(double));
          CW_d->data[ar] = W_d->data[k];
          nullAssignment(W_A, k + 1);

          /* delete constr */
          i30 = a->size[0] * a->size[1];
          a->size[0] = W_d->size[0];
          a->size[1] = W_d->size[1];
          emxEnsureCapacity((emxArray__common *)a, i30, (int)sizeof(double));
          br = W_d->size[0] * W_d->size[1];
          for (i30 = 0; i30 < br; i30++) {
            a->data[i30] = W_d->data[i30];
          }

          b_nullAssignment(a, k + 1);
          i30 = W_d->size[0] * W_d->size[1];
          W_d->size[0] = a->size[0];
          W_d->size[1] = a->size[1];
          emxEnsureCapacity((emxArray__common *)W_d, i30, (int)sizeof(double));
          br = a->size[0] * a->size[1];
          for (i30 = 0; i30 < br; i30++) {
            W_d->data[i30] = a->data[i30];
          }
        }
      } else {
        if ((CW_A->size[1] == 1) || (p->size[0] == 1)) {
          i30 = temp->size[0];
          temp->size[0] = CW_A->size[0];
          emxEnsureCapacity((emxArray__common *)temp, i30, (int)sizeof(double));
          br = CW_A->size[0];
          for (i30 = 0; i30 < br; i30++) {
            temp->data[i30] = 0.0;
            ar = CW_A->size[1];
            for (i31 = 0; i31 < ar; i31++) {
              temp->data[i30] += CW_A->data[i30 + CW_A->size[0] * i31] * p->
                data[i31];
            }
          }
        } else {
          k = CW_A->size[1];
          B_idx_0 = (unsigned int)CW_A->size[0];
          i30 = temp->size[0];
          temp->size[0] = (int)B_idx_0;
          emxEnsureCapacity((emxArray__common *)temp, i30, (int)sizeof(double));
          m = CW_A->size[0];
          ia = temp->size[0];
          i30 = temp->size[0];
          temp->size[0] = ia;
          emxEnsureCapacity((emxArray__common *)temp, i30, (int)sizeof(double));
          for (i30 = 0; i30 < ia; i30++) {
            temp->data[i30] = 0.0;
          }

          if (CW_A->size[0] == 0) {
          } else {
            ia = 0;
            while ((m > 0) && (ia <= 0)) {
              for (ic = 1; ic <= m; ic++) {
                temp->data[ic - 1] = 0.0;
              }

              ia = m;
            }

            br = 0;
            ia = 0;
            while ((m > 0) && (ia <= 0)) {
              ar = 0;
              i30 = br + k;
              for (idx = br; idx + 1 <= i30; idx++) {
                if (p->data[idx] != 0.0) {
                  ia = ar;
                  for (ic = 0; ic + 1 <= m; ic++) {
                    ia++;
                    temp->data[ic] += p->data[idx] * CW_A->data[ia - 1];
                  }
                }

                ar += m;
              }

              br += k;
              ia = m;
            }
          }
        }

        i30 = x->size[0];
        x->size[0] = temp->size[0];
        emxEnsureCapacity((emxArray__common *)x, i30, (int)sizeof(boolean_T));
        br = temp->size[0];
        for (i30 = 0; i30 < br; i30++) {
          x->data[i30] = (temp->data[i30] < -1.0E-8);
        }

        ar = x->size[0];
        idx = 0;
        i30 = ii->size[0];
        ii->size[0] = x->size[0];
        emxEnsureCapacity((emxArray__common *)ii, i30, (int)sizeof(int));
        ia = 1;
        exitg5 = false;
        while ((!exitg5) && (ia <= ar)) {
          guard1 = false;
          if (x->data[ia - 1]) {
            idx++;
            ii->data[idx - 1] = ia;
            if (idx >= ar) {
              exitg5 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            ia++;
          }
        }

        if (x->size[0] == 1) {
          if (idx == 0) {
            i30 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i30, (int)sizeof(int));
          }
        } else {
          i30 = ii->size[0];
          if (1 > idx) {
            ii->size[0] = 0;
          } else {
            ii->size[0] = idx;
          }

          emxEnsureCapacity((emxArray__common *)ii, i30, (int)sizeof(int));
        }

        i30 = I->size[0];
        I->size[0] = ii->size[0];
        emxEnsureCapacity((emxArray__common *)I, i30, (int)sizeof(double));
        br = ii->size[0];
        for (i30 = 0; i30 < br; i30++) {
          I->data[i30] = ii->data[i30];
        }

        if (I->size[0] == 0) {
          alpha = 1.0;
        } else {
          br = CW_A->size[1];
          i30 = a->size[0] * a->size[1];
          a->size[0] = I->size[0];
          a->size[1] = br;
          emxEnsureCapacity((emxArray__common *)a, i30, (int)sizeof(double));
          for (i30 = 0; i30 < br; i30++) {
            ar = I->size[0];
            for (i31 = 0; i31 < ar; i31++) {
              a->data[i31 + a->size[0] * i30] = CW_A->data[((int)I->data[i31] +
                CW_A->size[0] * i30) - 1];
            }
          }

          i30 = CW_A->size[1];
          if ((i30 == 1) || (x0->size[0] == 1)) {
            i30 = c_C->size[0];
            c_C->size[0] = a->size[0];
            emxEnsureCapacity((emxArray__common *)c_C, i30, (int)sizeof(double));
            br = a->size[0];
            for (i30 = 0; i30 < br; i30++) {
              c_C->data[i30] = 0.0;
              ar = a->size[1];
              for (i31 = 0; i31 < ar; i31++) {
                c_C->data[i30] += a->data[i30 + a->size[0] * i31] * x0->data[i31];
              }
            }
          } else {
            i30 = CW_A->size[1];
            B_idx_0 = (unsigned int)I->size[0];
            i31 = c_C->size[0];
            c_C->size[0] = (int)B_idx_0;
            emxEnsureCapacity((emxArray__common *)c_C, i31, (int)sizeof(double));
            m = I->size[0];
            ia = c_C->size[0];
            i31 = c_C->size[0];
            c_C->size[0] = ia;
            emxEnsureCapacity((emxArray__common *)c_C, i31, (int)sizeof(double));
            for (i31 = 0; i31 < ia; i31++) {
              c_C->data[i31] = 0.0;
            }

            ia = 0;
            while (ia <= 0) {
              for (ic = 1; ic <= m; ic++) {
                c_C->data[ic - 1] = 0.0;
              }

              ia = m;
            }

            br = 0;
            ia = 0;
            while (ia <= 0) {
              ar = 0;
              i31 = br + i30;
              for (idx = br; idx + 1 <= i31; idx++) {
                if (x0->data[idx] != 0.0) {
                  ia = ar;
                  for (ic = 0; ic + 1 <= m; ic++) {
                    ia++;
                    c_C->data[ic] += x0->data[idx] * a->data[ia - 1];
                  }
                }

                ar += m;
              }

              br += i30;
              ia = m;
            }
          }

          br = CW_A->size[1];
          i30 = a->size[0] * a->size[1];
          a->size[0] = I->size[0];
          a->size[1] = br;
          emxEnsureCapacity((emxArray__common *)a, i30, (int)sizeof(double));
          for (i30 = 0; i30 < br; i30++) {
            ar = I->size[0];
            for (i31 = 0; i31 < ar; i31++) {
              a->data[i31 + a->size[0] * i30] = CW_A->data[((int)I->data[i31] +
                CW_A->size[0] * i30) - 1];
            }
          }

          i30 = CW_A->size[1];
          if ((i30 == 1) || (p->size[0] == 1)) {
            i30 = b_y->size[0];
            b_y->size[0] = a->size[0];
            emxEnsureCapacity((emxArray__common *)b_y, i30, (int)sizeof(double));
            br = a->size[0];
            for (i30 = 0; i30 < br; i30++) {
              b_y->data[i30] = 0.0;
              ar = a->size[1];
              for (i31 = 0; i31 < ar; i31++) {
                b_y->data[i30] += a->data[i30 + a->size[0] * i31] * p->data[i31];
              }
            }
          } else {
            i30 = CW_A->size[1];
            B_idx_0 = (unsigned int)I->size[0];
            i31 = b_y->size[0];
            b_y->size[0] = (int)B_idx_0;
            emxEnsureCapacity((emxArray__common *)b_y, i31, (int)sizeof(double));
            m = I->size[0];
            ar = b_y->size[0];
            i31 = b_y->size[0];
            b_y->size[0] = ar;
            emxEnsureCapacity((emxArray__common *)b_y, i31, (int)sizeof(double));
            for (i31 = 0; i31 < ar; i31++) {
              b_y->data[i31] = 0.0;
            }

            ia = 0;
            while (ia <= 0) {
              for (ic = 1; ic <= m; ic++) {
                b_y->data[ic - 1] = 0.0;
              }

              ia = m;
            }

            br = 0;
            ia = 0;
            while (ia <= 0) {
              ar = 0;
              i31 = br + i30;
              for (idx = br; idx + 1 <= i31; idx++) {
                if (p->data[idx] != 0.0) {
                  ia = ar;
                  for (ic = 0; ic + 1 <= m; ic++) {
                    ia++;
                    b_y->data[ic] += p->data[idx] * a->data[ia - 1];
                  }
                }

                ar += m;
              }

              br += i30;
              ia = m;
            }
          }

          i30 = f_C->size[0];
          f_C->size[0] = c_C->size[0];
          emxEnsureCapacity((emxArray__common *)f_C, i30, (int)sizeof(double));
          br = c_C->size[0];
          for (i30 = 0; i30 < br; i30++) {
            f_C->data[i30] = c_C->data[i30] - CW_d->data[(int)I->data[i30] - 1];
          }

          rdivide(f_C, b_y, varargin_2);
          ar = 1;
          ia = varargin_2->size[0];
          mtmp = varargin_2->data[0];
          k = 0;
          if (varargin_2->size[0] > 1) {
            if (rtIsNaN(varargin_2->data[0])) {
              idx = 2;
              exitg4 = false;
              while ((!exitg4) && (idx <= ia)) {
                ar = idx;
                if (!rtIsNaN(varargin_2->data[idx - 1])) {
                  mtmp = varargin_2->data[idx - 1];
                  k = idx - 1;
                  exitg4 = true;
                } else {
                  idx++;
                }
              }
            }

            if (ar < varargin_2->size[0]) {
              while (ar + 1 <= ia) {
                if (varargin_2->data[ar] < mtmp) {
                  mtmp = varargin_2->data[ar];
                  k = ar;
                }

                ar++;
              }
            }
          }

          Imin_alpha = k;
          if ((1.0 <= mtmp) || rtIsNaN(mtmp)) {
            alpha = 1.0;
          } else {
            alpha = mtmp;
          }
        }

        i30 = x0->size[0];
        emxEnsureCapacity((emxArray__common *)x0, i30, (int)sizeof(double));
        br = x0->size[0];
        for (i30 = 0; i30 < br; i30++) {
          x0->data[i30] -= alpha * p->data[i30];
        }

        /* NB */
        if (alpha < 1.0) {
          /* alpha>0 && */
          /*          [~,Imin_alpha] = min((CW.A(I,:) * x - CW.d(I))./(CW.A(I,:) * p)); */
          if (!((W_A->size[0] == 0) || (W_A->size[1] == 0))) {
            ia = W_A->size[1];
          } else {
            ia = CW_A->size[1];
          }

          if ((ia == 0) || (!((W_A->size[0] == 0) || (W_A->size[1] == 0)))) {
            ic = W_A->size[0];
          } else {
            ic = 0;
          }

          br = CW_A->size[1];
          ar = (int)I->data[Imin_alpha];
          i30 = b_CW_A->size[0] * b_CW_A->size[1];
          b_CW_A->size[0] = 1;
          b_CW_A->size[1] = br;
          emxEnsureCapacity((emxArray__common *)b_CW_A, i30, (int)sizeof(double));
          for (i30 = 0; i30 < br; i30++) {
            b_CW_A->data[b_CW_A->size[0] * i30] = CW_A->data[(ar + CW_A->size[0]
              * i30) - 1];
          }

          i30 = b_W_A->size[0] * b_W_A->size[1];
          b_W_A->size[0] = ic + 1;
          b_W_A->size[1] = ia;
          emxEnsureCapacity((emxArray__common *)b_W_A, i30, (int)sizeof(double));
          for (i30 = 0; i30 < ia; i30++) {
            for (i31 = 0; i31 < ic; i31++) {
              b_W_A->data[i31 + b_W_A->size[0] * i30] = W_A->data[i31 + ic * i30];
            }
          }

          for (i30 = 0; i30 < ia; i30++) {
            for (i31 = 0; i31 < 1; i31++) {
              b_W_A->data[ic + b_W_A->size[0] * i30] = b_CW_A->data[i30];
            }
          }

          i30 = W_A->size[0] * W_A->size[1];
          W_A->size[0] = b_W_A->size[0];
          W_A->size[1] = b_W_A->size[1];
          emxEnsureCapacity((emxArray__common *)W_A, i30, (int)sizeof(double));
          br = b_W_A->size[1];
          for (i30 = 0; i30 < br; i30++) {
            ar = b_W_A->size[0];
            for (i31 = 0; i31 < ar; i31++) {
              W_A->data[i31 + W_A->size[0] * i30] = b_W_A->data[i31 +
                b_W_A->size[0] * i30];
            }
          }

          /* add constr */
          if (!((W_d->size[0] == 0) || (W_d->size[1] == 0))) {
            ia = W_d->size[0];
          } else {
            ia = 0;
          }

          i30 = varargin_2->size[0];
          varargin_2->size[0] = ia + 1;
          emxEnsureCapacity((emxArray__common *)varargin_2, i30, (int)sizeof
                            (double));
          for (i30 = 0; i30 < ia; i30++) {
            varargin_2->data[i30] = W_d->data[i30];
          }

          varargin_2->data[ia] = CW_d->data[(int)I->data[Imin_alpha] - 1];
          i30 = W_d->size[0] * W_d->size[1];
          W_d->size[0] = varargin_2->size[0];
          W_d->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)W_d, i30, (int)sizeof(double));
          br = varargin_2->size[0];
          for (i30 = 0; i30 < br; i30++) {
            W_d->data[i30] = varargin_2->data[i30];
          }

          nullAssignment(CW_A, (int)I->data[Imin_alpha]);

          /* delete */
          c_nullAssignment(CW_d, (int)I->data[Imin_alpha]);
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_real_T(&c_CW_A);
  emxFree_real_T(&c_W_A);
  emxFree_real_T(&b_W_A);
  emxFree_real_T(&b_CW_A);
  emxFree_real_T(&b_c);
  emxFree_real_T(&f_C);
  emxFree_real_T(&r4);
  emxFree_real_T(&e_C);
  emxFree_real_T(&d_C);
  emxFree_real_T(&varargin_2);
  emxFree_real_T(&b_y);
  emxFree_real_T(&a);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&x);
  emxFree_real_T(&y);
  emxFree_real_T(&c_C);
  emxFree_real_T(&b_C);
  emxFree_real_T(&temp);
  emxFree_real_T(&I);
  emxFree_real_T(&r);
  emxFree_real_T(&p);
  emxFree_real_T(&CW_d);
  emxFree_real_T(&CW_A);
  emxFree_real_T(&W_d);
  emxFree_real_T(&W_A);

  /*  1/2*x'*B*x - b'*x; */
  return flag;
}

/*
 * File trailer for QPACM.c
 *
 * [EOF]
 */
