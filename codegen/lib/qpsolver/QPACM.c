/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: QPACM.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
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
 * NULL SPACE sol v e s the eq ual i ty constrained convex QP:
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
  int i7;
  int br;
  emxArray_real_T *b_A;
  int mn;
  int i8;
  emxArray_real_T *Q1;
  emxArray_real_T *Q2t;
  emxArray_real_T *Q;
  emxArray_real_T *Q2;
  int i9;
  emxArray_real_T *R;
  int i10;
  emxArray_real_T *c_A;
  emxArray_real_T *py;
  emxArray_real_T *pz;
  emxArray_int32_T *jpvt;
  emxArray_real_T *B;
  emxArray_real_T *r0;
  unsigned int unnamed_idx_0;
  int n;
  unsigned int unnamed_idx_1;
  int kAcol;
  int jBcol;
  int nb;
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
  i7 = g->size[0];
  emxEnsureCapacity((emxArray__common *)g, i7, (int)sizeof(double));
  br = g->size[0];
  for (i7 = 0; i7 < br; i7++) {
    g->data[i7] = -g->data[i7];
  }

  emxInit_real_T1(&b_A, 2);
  i7 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[1];
  b_A->size[1] = A->size[0];
  emxEnsureCapacity((emxArray__common *)b_A, i7, (int)sizeof(double));
  br = A->size[0];
  for (i7 = 0; i7 < br; i7++) {
    mn = A->size[1];
    for (i8 = 0; i8 < mn; i8++) {
      b_A->data[i8 + b_A->size[0] * i7] = A->data[i7 + A->size[0] * i8];
    }
  }

  i7 = A->size[0] * A->size[1];
  A->size[0] = b_A->size[0];
  A->size[1] = b_A->size[1];
  emxEnsureCapacity((emxArray__common *)A, i7, (int)sizeof(double));
  br = b_A->size[1];
  for (i7 = 0; i7 < br; i7++) {
    mn = b_A->size[0];
    for (i8 = 0; i8 < mn; i8++) {
      A->data[i8 + A->size[0] * i7] = b_A->data[i8 + b_A->size[0] * i7];
    }
  }

  emxFree_real_T(&b_A);
  if (A->size[1] != 0) {
    emxInit_real_T1(&Q1, 2);
    emxInit_real_T1(&Q2t, 2);
    emxInit_real_T1(&Q, 2);

    /*  for si t u a t i o ns where A i s empty */
    qr(A, Q, Q2t);
    br = Q->size[0];
    mn = A->size[1];
    i7 = Q1->size[0] * Q1->size[1];
    Q1->size[0] = br;
    Q1->size[1] = mn;
    emxEnsureCapacity((emxArray__common *)Q1, i7, (int)sizeof(double));
    for (i7 = 0; i7 < mn; i7++) {
      for (i8 = 0; i8 < br; i8++) {
        Q1->data[i8 + Q1->size[0] * i7] = Q->data[i8 + Q->size[0] * i7];
      }
    }

    if (A->size[1] + 1U > (unsigned int)A->size[0]) {
      i7 = 1;
      i8 = 1;
    } else {
      i7 = (int)(A->size[1] + 1U);
      i8 = A->size[0] + 1;
    }

    emxInit_real_T1(&Q2, 2);
    br = Q->size[0];
    i9 = Q2->size[0] * Q2->size[1];
    Q2->size[0] = br;
    Q2->size[1] = i8 - i7;
    emxEnsureCapacity((emxArray__common *)Q2, i9, (int)sizeof(double));
    mn = i8 - i7;
    for (i9 = 0; i9 < mn; i9++) {
      for (i10 = 0; i10 < br; i10++) {
        Q2->data[i10 + Q2->size[0] * i9] = Q->data[i10 + Q->size[0] * ((i7 + i9)
          - 1)];
      }
    }

    emxInit_real_T1(&R, 2);
    br = A->size[1];
    mn = Q2t->size[1];
    i9 = R->size[0] * R->size[1];
    R->size[0] = br;
    R->size[1] = mn;
    emxEnsureCapacity((emxArray__common *)R, i9, (int)sizeof(double));
    for (i9 = 0; i9 < mn; i9++) {
      for (i10 = 0; i10 < br; i10++) {
        R->data[i10 + R->size[0] * i9] = Q2t->data[i10 + Q2t->size[0] * i9];
      }
    }

    emxInit_real_T1(&c_A, 2);
    i9 = c_A->size[0] * c_A->size[1];
    c_A->size[0] = R->size[1];
    c_A->size[1] = R->size[0];
    emxEnsureCapacity((emxArray__common *)c_A, i9, (int)sizeof(double));
    br = R->size[0];
    for (i9 = 0; i9 < br; i9++) {
      mn = R->size[1];
      for (i10 = 0; i10 < mn; i10++) {
        c_A->data[i10 + c_A->size[0] * i9] = R->data[i9 + R->size[0] * i10];
      }
    }

    emxInit_real_T1(&py, 2);
    emxInit_real_T(&pz, 1);
    emxInit_int32_T(&jpvt, 2);
    emxInit_real_T1(&B, 2);
    emxInit_real_T1(&r0, 2);
    if ((c_A->size[0] == 0) || ((b->size[0] == 0) || (b->size[1] == 0))) {
      unnamed_idx_0 = (unsigned int)c_A->size[1];
      unnamed_idx_1 = (unsigned int)b->size[1];
      i9 = py->size[0] * py->size[1];
      py->size[0] = (int)unnamed_idx_0;
      py->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)py, i9, (int)sizeof(double));
      br = (int)unnamed_idx_0 * (int)unnamed_idx_1;
      for (i9 = 0; i9 < br; i9++) {
        py->data[i9] = 0.0;
      }
    } else if (c_A->size[0] == c_A->size[1]) {
      n = c_A->size[1];
      i9 = Q2t->size[0] * Q2t->size[1];
      Q2t->size[0] = c_A->size[0];
      Q2t->size[1] = c_A->size[1];
      emxEnsureCapacity((emxArray__common *)Q2t, i9, (int)sizeof(double));
      br = c_A->size[0] * c_A->size[1];
      for (i9 = 0; i9 < br; i9++) {
        Q2t->data[i9] = c_A->data[i9];
      }

      xgetrf(c_A->size[1], c_A->size[1], Q2t, c_A->size[1], jpvt, &br);
      nb = b->size[1];
      i9 = py->size[0] * py->size[1];
      py->size[0] = b->size[0];
      py->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)py, i9, (int)sizeof(double));
      br = b->size[0] * b->size[1];
      for (i9 = 0; i9 < br; i9++) {
        py->data[i9] = b->data[i9];
      }

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
            for (br = k + 1; br + 1 <= n; br++) {
              py->data[br + jBcol] -= py->data[k + jBcol] * Q2t->data[br + kAcol];
            }
          }
        }
      }

      i9 = r0->size[0] * r0->size[1];
      r0->size[0] = py->size[0];
      r0->size[1] = py->size[1];
      emxEnsureCapacity((emxArray__common *)r0, i9, (int)sizeof(double));
      br = py->size[0] * py->size[1];
      for (i9 = 0; i9 < br; i9++) {
        r0->data[i9] = py->data[i9];
      }

      xtrsm(c_A->size[1], b->size[1], Q2t, c_A->size[1], r0, c_A->size[1]);
      i9 = py->size[0] * py->size[1];
      py->size[0] = r0->size[0];
      py->size[1] = r0->size[1];
      emxEnsureCapacity((emxArray__common *)py, i9, (int)sizeof(double));
      br = r0->size[0] * r0->size[1];
      for (i9 = 0; i9 < br; i9++) {
        py->data[i9] = r0->data[i9];
      }
    } else {
      c_xgeqp3(c_A, pz, jpvt);
      kAcol = rankFromQR(c_A);
      jBcol = c_A->size[1];
      mn = b->size[1];
      i9 = py->size[0] * py->size[1];
      py->size[0] = jBcol;
      py->size[1] = mn;
      emxEnsureCapacity((emxArray__common *)py, i9, (int)sizeof(double));
      br = jBcol * mn;
      for (i9 = 0; i9 < br; i9++) {
        py->data[i9] = 0.0;
      }

      i9 = B->size[0] * B->size[1];
      B->size[0] = b->size[0];
      B->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)B, i9, (int)sizeof(double));
      br = b->size[0] * b->size[1];
      for (i9 = 0; i9 < br; i9++) {
        B->data[i9] = b->data[i9];
      }

      m = c_A->size[0];
      jBcol = c_A->size[0];
      mn = c_A->size[1];
      if (jBcol <= mn) {
        mn = jBcol;
      }

      for (j = 0; j + 1 <= mn; j++) {
        if (pz->data[j] != 0.0) {
          for (k = 0; k + 1 <= b->size[1]; k++) {
            wj = B->data[j + B->size[0] * k];
            for (br = j + 1; br + 1 <= m; br++) {
              wj += c_A->data[br + c_A->size[0] * j] * B->data[br + B->size[0] *
                k];
            }

            wj *= pz->data[j];
            if (wj != 0.0) {
              B->data[j + B->size[0] * k] -= wj;
              for (br = j + 1; br + 1 <= m; br++) {
                B->data[br + B->size[0] * k] -= c_A->data[br + c_A->size[0] * j]
                  * wj;
              }
            }
          }
        }
      }

      for (k = 0; k + 1 <= b->size[1]; k++) {
        for (br = 0; br + 1 <= kAcol; br++) {
          py->data[(jpvt->data[br] + py->size[0] * k) - 1] = B->data[br +
            B->size[0] * k];
        }

        for (j = kAcol - 1; j + 1 > 0; j--) {
          py->data[(jpvt->data[j] + py->size[0] * k) - 1] /= c_A->data[j +
            c_A->size[0] * j];
          for (br = 0; br + 1 <= j; br++) {
            py->data[(jpvt->data[br] + py->size[0] * k) - 1] -= py->data
              [(jpvt->data[j] + py->size[0] * k) - 1] * c_A->data[br + c_A->
              size[0] * j];
          }
        }
      }
    }

    emxFree_real_T(&r0);
    emxFree_real_T(&B);
    emxFree_int32_T(&jpvt);
    i9 = Q2t->size[0] * Q2t->size[1];
    Q2t->size[0] = Q2->size[1];
    Q2t->size[1] = Q2->size[0];
    emxEnsureCapacity((emxArray__common *)Q2t, i9, (int)sizeof(double));
    br = Q2->size[0];
    for (i9 = 0; i9 < br; i9++) {
      mn = Q2->size[1];
      for (i10 = 0; i10 < mn; i10++) {
        Q2t->data[i10 + Q2t->size[0] * i9] = Q2->data[i9 + Q2->size[0] * i10];
      }
    }

    i9 = A->size[1];
    emxInit_real_T1(&y, 2);
    if ((i9 == 1) || (py->size[0] == 1)) {
      i9 = y->size[0] * y->size[1];
      y->size[0] = Q1->size[0];
      y->size[1] = py->size[1];
      emxEnsureCapacity((emxArray__common *)y, i9, (int)sizeof(double));
      br = Q1->size[0];
      for (i9 = 0; i9 < br; i9++) {
        mn = py->size[1];
        for (i10 = 0; i10 < mn; i10++) {
          y->data[i9 + y->size[0] * i10] = 0.0;
          jBcol = Q1->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            y->data[i9 + y->size[0] * i10] += Q1->data[i9 + Q1->size[0] * iy] *
              py->data[iy + py->size[0] * i10];
          }
        }
      }
    } else {
      i9 = A->size[1];
      i10 = Q->size[0];
      unnamed_idx_1 = (unsigned int)py->size[1];
      iy = y->size[0] * y->size[1];
      y->size[0] = i10;
      y->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(double));
      i10 = Q->size[0];
      iy = y->size[0] * y->size[1];
      emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(double));
      br = y->size[1];
      for (iy = 0; iy < br; iy++) {
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
        while ((i10 > 0) && (kAcol <= jBcol)) {
          iy = kAcol + i10;
          for (ix = kAcol; ix + 1 <= iy; ix++) {
            y->data[ix] = 0.0;
          }

          kAcol += i10;
        }

        br = 0;
        kAcol = 0;
        while ((i10 > 0) && (kAcol <= jBcol)) {
          nb = 0;
          iy = br + i9;
          for (b_n = br; b_n + 1 <= iy; b_n++) {
            if (py->data[b_n] != 0.0) {
              ia = nb;
              m = kAcol + i10;
              for (ix = kAcol; ix + 1 <= m; ix++) {
                ia++;
                mn = Q->size[0];
                y->data[ix] += py->data[b_n] * Q->data[(ia - 1) % mn + Q->size[0]
                  * ((ia - 1) / mn)];
              }
            }

            nb += i10;
          }

          br += i9;
          kAcol += i10;
        }
      }
    }

    emxInit_real_T1(&C, 2);
    if ((G->size[1] == 1) || (y->size[0] == 1)) {
      i9 = C->size[0] * C->size[1];
      C->size[0] = G->size[0];
      C->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)C, i9, (int)sizeof(double));
      br = G->size[0];
      for (i9 = 0; i9 < br; i9++) {
        mn = y->size[1];
        for (i10 = 0; i10 < mn; i10++) {
          C->data[i9 + C->size[0] * i10] = 0.0;
          jBcol = G->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            C->data[i9 + C->size[0] * i10] += G->data[i9 + G->size[0] * iy] *
              y->data[iy + y->size[0] * i10];
          }
        }
      }
    } else {
      k = G->size[1];
      unnamed_idx_0 = (unsigned int)G->size[0];
      unnamed_idx_1 = (unsigned int)y->size[1];
      i9 = C->size[0] * C->size[1];
      C->size[0] = (int)unnamed_idx_0;
      C->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)C, i9, (int)sizeof(double));
      m = G->size[0];
      i9 = C->size[0] * C->size[1];
      emxEnsureCapacity((emxArray__common *)C, i9, (int)sizeof(double));
      br = C->size[1];
      for (i9 = 0; i9 < br; i9++) {
        mn = C->size[0];
        for (i10 = 0; i10 < mn; i10++) {
          C->data[i10 + C->size[0] * i9] = 0.0;
        }
      }

      if ((G->size[0] == 0) || (y->size[1] == 0)) {
      } else {
        jBcol = G->size[0] * (y->size[1] - 1);
        kAcol = 0;
        while ((m > 0) && (kAcol <= jBcol)) {
          i9 = kAcol + m;
          for (ix = kAcol; ix + 1 <= i9; ix++) {
            C->data[ix] = 0.0;
          }

          kAcol += m;
        }

        br = 0;
        kAcol = 0;
        while ((m > 0) && (kAcol <= jBcol)) {
          nb = 0;
          i9 = br + k;
          for (b_n = br; b_n + 1 <= i9; b_n++) {
            if (y->data[b_n] != 0.0) {
              ia = nb;
              i10 = kAcol + m;
              for (ix = kAcol; ix + 1 <= i10; ix++) {
                ia++;
                C->data[ix] += y->data[b_n] * G->data[ia - 1];
              }
            }

            nb += m;
          }

          br += k;
          kAcol += m;
        }
      }
    }

    emxFree_real_T(&y);
    i9 = pz->size[0];
    pz->size[0] = C->size[0];
    emxEnsureCapacity((emxArray__common *)pz, i9, (int)sizeof(double));
    br = C->size[0];
    for (i9 = 0; i9 < br; i9++) {
      pz->data[i9] = C->data[i9] + g->data[i9];
    }

    emxFree_real_T(&C);
    emxInit_real_T(&gz, 1);
    if ((Q2t->size[1] == 1) || (pz->size[0] == 1)) {
      i9 = gz->size[0];
      gz->size[0] = Q2t->size[0];
      emxEnsureCapacity((emxArray__common *)gz, i9, (int)sizeof(double));
      br = Q2t->size[0];
      for (i9 = 0; i9 < br; i9++) {
        gz->data[i9] = 0.0;
        mn = Q2t->size[1];
        for (i10 = 0; i10 < mn; i10++) {
          gz->data[i9] += Q2t->data[i9 + Q2t->size[0] * i10] * pz->data[i10];
        }
      }
    } else {
      k = Q2t->size[1];
      unnamed_idx_0 = (unsigned int)Q2t->size[0];
      i9 = gz->size[0];
      gz->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)gz, i9, (int)sizeof(double));
      m = Q2t->size[0];
      jBcol = gz->size[0];
      i9 = gz->size[0];
      gz->size[0] = jBcol;
      emxEnsureCapacity((emxArray__common *)gz, i9, (int)sizeof(double));
      for (i9 = 0; i9 < jBcol; i9++) {
        gz->data[i9] = 0.0;
      }

      if (Q2t->size[0] != 0) {
        kAcol = 0;
        while ((m > 0) && (kAcol <= 0)) {
          for (ix = 1; ix <= m; ix++) {
            gz->data[ix - 1] = 0.0;
          }

          kAcol = m;
        }

        br = 0;
        kAcol = 0;
        while ((m > 0) && (kAcol <= 0)) {
          nb = 0;
          i9 = br + k;
          for (b_n = br; b_n + 1 <= i9; b_n++) {
            if (pz->data[b_n] != 0.0) {
              ia = nb;
              for (ix = 0; ix + 1 <= m; ix++) {
                ia++;
                gz->data[ix] += pz->data[b_n] * Q2t->data[ia - 1];
              }
            }

            nb += m;
          }

          br += k;
          kAcol = m;
        }
      }
    }

    emxInit_real_T1(&b_y, 2);
    if ((Q2t->size[1] == 1) || (G->size[0] == 1)) {
      i9 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = Q2t->size[0];
      b_y->size[1] = G->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, i9, (int)sizeof(double));
      br = Q2t->size[0];
      for (i9 = 0; i9 < br; i9++) {
        mn = G->size[1];
        for (i10 = 0; i10 < mn; i10++) {
          b_y->data[i9 + b_y->size[0] * i10] = 0.0;
          jBcol = Q2t->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            b_y->data[i9 + b_y->size[0] * i10] += Q2t->data[i9 + Q2t->size[0] *
              iy] * G->data[iy + G->size[0] * i10];
          }
        }
      }
    } else {
      k = Q2t->size[1];
      unnamed_idx_0 = (unsigned int)Q2t->size[0];
      unnamed_idx_1 = (unsigned int)G->size[1];
      i9 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = (int)unnamed_idx_0;
      b_y->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)b_y, i9, (int)sizeof(double));
      m = Q2t->size[0];
      i9 = b_y->size[0] * b_y->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, i9, (int)sizeof(double));
      br = b_y->size[1];
      for (i9 = 0; i9 < br; i9++) {
        mn = b_y->size[0];
        for (i10 = 0; i10 < mn; i10++) {
          b_y->data[i10 + b_y->size[0] * i9] = 0.0;
        }
      }

      if ((Q2t->size[0] == 0) || (G->size[1] == 0)) {
      } else {
        jBcol = Q2t->size[0] * (G->size[1] - 1);
        kAcol = 0;
        while ((m > 0) && (kAcol <= jBcol)) {
          i9 = kAcol + m;
          for (ix = kAcol; ix + 1 <= i9; ix++) {
            b_y->data[ix] = 0.0;
          }

          kAcol += m;
        }

        br = 0;
        kAcol = 0;
        while ((m > 0) && (kAcol <= jBcol)) {
          nb = 0;
          i9 = br + k;
          for (b_n = br; b_n + 1 <= i9; b_n++) {
            if (G->data[b_n] != 0.0) {
              ia = nb;
              i10 = kAcol + m;
              for (ix = kAcol; ix + 1 <= i10; ix++) {
                ia++;
                b_y->data[ix] += G->data[b_n] * Q2t->data[ia - 1];
              }
            }

            nb += m;
          }

          br += k;
          kAcol += m;
        }
      }
    }

    emxInit_real_T1(&Gz, 2);
    guard1 = false;
    if (b_y->size[1] == 1) {
      guard1 = true;
    } else {
      i9 = Q->size[0];
      if (i9 == 1) {
        guard1 = true;
      } else {
        k = b_y->size[1];
        unnamed_idx_0 = (unsigned int)b_y->size[0];
        i9 = Gz->size[0] * Gz->size[1];
        Gz->size[0] = (int)unnamed_idx_0;
        Gz->size[1] = i8 - i7;
        emxEnsureCapacity((emxArray__common *)Gz, i9, (int)sizeof(double));
        m = b_y->size[0];
        i9 = Gz->size[0] * Gz->size[1];
        emxEnsureCapacity((emxArray__common *)Gz, i9, (int)sizeof(double));
        br = Gz->size[1];
        for (i9 = 0; i9 < br; i9++) {
          mn = Gz->size[0];
          for (i10 = 0; i10 < mn; i10++) {
            Gz->data[i10 + Gz->size[0] * i9] = 0.0;
          }
        }

        if ((b_y->size[0] == 0) || (i8 - i7 == 0)) {
        } else {
          jBcol = b_y->size[0] * ((i8 - i7) - 1);
          kAcol = 0;
          while ((m > 0) && (kAcol <= jBcol)) {
            i9 = kAcol + m;
            for (ix = kAcol; ix + 1 <= i9; ix++) {
              Gz->data[ix] = 0.0;
            }

            kAcol += m;
          }

          br = 0;
          kAcol = 0;
          while ((m > 0) && (kAcol <= jBcol)) {
            nb = 0;
            i9 = br + k;
            for (b_n = br; b_n + 1 <= i9; b_n++) {
              i10 = Q->size[0];
              if (Q->data[b_n % i10 + Q->size[0] * ((i7 + b_n / i10) - 1)] !=
                  0.0) {
                ia = nb;
                i10 = kAcol + m;
                for (ix = kAcol; ix + 1 <= i10; ix++) {
                  ia++;
                  iy = Q->size[0];
                  Gz->data[ix] += Q->data[b_n % iy + Q->size[0] * ((i7 + b_n /
                    iy) - 1)] * b_y->data[ia - 1];
                }
              }

              nb += m;
            }

            br += k;
            kAcol += m;
          }
        }
      }
    }

    if (guard1) {
      i9 = Gz->size[0] * Gz->size[1];
      Gz->size[0] = b_y->size[0];
      Gz->size[1] = Q2->size[1];
      emxEnsureCapacity((emxArray__common *)Gz, i9, (int)sizeof(double));
      br = b_y->size[0];
      for (i9 = 0; i9 < br; i9++) {
        mn = Q2->size[1];
        for (i10 = 0; i10 < mn; i10++) {
          Gz->data[i9 + Gz->size[0] * i10] = 0.0;
          jBcol = b_y->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            Gz->data[i9 + Gz->size[0] * i10] += b_y->data[i9 + b_y->size[0] * iy]
              * Q2->data[iy + Q2->size[0] * i10];
          }
        }
      }
    }

    emxFree_real_T(&b_y);
    i9 = c_A->size[0] * c_A->size[1];
    c_A->size[0] = Gz->size[0];
    c_A->size[1] = Gz->size[1];
    emxEnsureCapacity((emxArray__common *)c_A, i9, (int)sizeof(double));
    br = Gz->size[0] * Gz->size[1];
    for (i9 = 0; i9 < br; i9++) {
      c_A->data[i9] = Gz->data[i9];
    }

    n = Gz->size[1];
    if (Gz->size[1] != 0) {
      b_n = Gz->size[0];
      br = -1;
      if (Gz->size[0] != 0) {
        m = 0;
        j = 0;
        exitg1 = false;
        while ((!exitg1) && (j + 1 <= b_n)) {
          mn = m + j;
          wj = 0.0;
          if (!(j < 1)) {
            ix = m;
            iy = m;
            for (k = 1; k <= j; k++) {
              wj += c_A->data[ix] * c_A->data[iy];
              ix++;
              iy++;
            }
          }

          wj = c_A->data[mn] - wj;
          if (wj > 0.0) {
            wj = sqrt(wj);
            c_A->data[mn] = wj;
            if (j + 1 < b_n) {
              nb = (b_n - j) - 2;
              mn += b_n;
              jBcol = (m + b_n) + 1;
              if ((j == 0) || (nb + 1 == 0)) {
              } else {
                iy = mn;
                i9 = jBcol + n * nb;
                for (kAcol = jBcol; kAcol <= i9; kAcol += n) {
                  ix = m + 1;
                  c = 0.0;
                  i10 = (kAcol + j) - 1;
                  for (ia = kAcol; ia <= i10; ia++) {
                    c += c_A->data[ia - 1] * c_A->data[ix - 1];
                    ix++;
                  }

                  c_A->data[iy] += -c;
                  iy += n;
                }
              }

              wj = 1.0 / wj;
              i9 = (mn + b_n * nb) + 1;
              while (mn + 1 <= i9) {
                c_A->data[mn] *= wj;
                mn += b_n;
              }

              m = jBcol - 1;
            }

            j++;
          } else {
            c_A->data[mn] = wj;
            br = j;
            exitg1 = true;
          }
        }
      }

      if (br + 1 == 0) {
        jBcol = Gz->size[1];
      } else {
        jBcol = br;
      }

      for (j = 0; j + 1 <= jBcol; j++) {
        for (br = j + 1; br + 1 <= jBcol; br++) {
          c_A->data[br + c_A->size[0] * j] = 0.0;
        }
      }
    }

    emxFree_real_T(&Gz);
    i9 = Q2t->size[0] * Q2t->size[1];
    Q2t->size[0] = c_A->size[1];
    Q2t->size[1] = c_A->size[0];
    emxEnsureCapacity((emxArray__common *)Q2t, i9, (int)sizeof(double));
    br = c_A->size[0];
    for (i9 = 0; i9 < br; i9++) {
      mn = c_A->size[1];
      for (i10 = 0; i10 < mn; i10++) {
        Q2t->data[i10 + Q2t->size[0] * i9] = c_A->data[i9 + c_A->size[0] * i10];
      }
    }

    emxFree_real_T(&c_A);
    emxInit_real_T(&b_gz, 1);
    i9 = b_gz->size[0];
    b_gz->size[0] = gz->size[0];
    emxEnsureCapacity((emxArray__common *)b_gz, i9, (int)sizeof(double));
    br = gz->size[0];
    for (i9 = 0; i9 < br; i9++) {
      b_gz->data[i9] = -gz->data[i9];
    }

    emxFree_real_T(&gz);
    emxInit_real_T1(&b_Q2t, 2);
    b_mldivide(Q2t, b_gz, pz);
    i9 = b_Q2t->size[0] * b_Q2t->size[1];
    b_Q2t->size[0] = Q2t->size[1];
    b_Q2t->size[1] = Q2t->size[0];
    emxEnsureCapacity((emxArray__common *)b_Q2t, i9, (int)sizeof(double));
    br = Q2t->size[0];
    emxFree_real_T(&b_gz);
    for (i9 = 0; i9 < br; i9++) {
      mn = Q2t->size[1];
      for (i10 = 0; i10 < mn; i10++) {
        b_Q2t->data[i10 + b_Q2t->size[0] * i9] = Q2t->data[i9 + Q2t->size[0] *
          i10];
      }
    }

    emxInit_real_T(&b_pz, 1);
    i9 = b_pz->size[0];
    b_pz->size[0] = pz->size[0];
    emxEnsureCapacity((emxArray__common *)b_pz, i9, (int)sizeof(double));
    br = pz->size[0];
    for (i9 = 0; i9 < br; i9++) {
      b_pz->data[i9] = pz->data[i9];
    }

    b_mldivide(b_Q2t, b_pz, pz);
    i9 = A->size[1];
    emxFree_real_T(&b_pz);
    emxFree_real_T(&b_Q2t);
    emxInit_real_T1(&b_C, 2);
    if ((i9 == 1) || (py->size[0] == 1)) {
      i9 = b_C->size[0] * b_C->size[1];
      b_C->size[0] = Q1->size[0];
      b_C->size[1] = py->size[1];
      emxEnsureCapacity((emxArray__common *)b_C, i9, (int)sizeof(double));
      br = Q1->size[0];
      for (i9 = 0; i9 < br; i9++) {
        mn = py->size[1];
        for (i10 = 0; i10 < mn; i10++) {
          b_C->data[i9 + b_C->size[0] * i10] = 0.0;
          jBcol = Q1->size[1];
          for (iy = 0; iy < jBcol; iy++) {
            b_C->data[i9 + b_C->size[0] * i10] += Q1->data[i9 + Q1->size[0] * iy]
              * py->data[iy + py->size[0] * i10];
          }
        }
      }
    } else {
      i9 = A->size[1];
      i10 = Q->size[0];
      unnamed_idx_1 = (unsigned int)py->size[1];
      iy = b_C->size[0] * b_C->size[1];
      b_C->size[0] = i10;
      b_C->size[1] = (int)unnamed_idx_1;
      emxEnsureCapacity((emxArray__common *)b_C, iy, (int)sizeof(double));
      i10 = Q->size[0];
      iy = b_C->size[0] * b_C->size[1];
      emxEnsureCapacity((emxArray__common *)b_C, iy, (int)sizeof(double));
      br = b_C->size[1];
      for (iy = 0; iy < br; iy++) {
        mn = b_C->size[0];
        for (m = 0; m < mn; m++) {
          b_C->data[m + b_C->size[0] * iy] = 0.0;
        }
      }

      iy = Q->size[0];
      if ((iy == 0) || (py->size[1] == 0)) {
      } else {
        iy = Q->size[0];
        jBcol = iy * (py->size[1] - 1);
        kAcol = 0;
        while ((i10 > 0) && (kAcol <= jBcol)) {
          iy = kAcol + i10;
          for (ix = kAcol; ix + 1 <= iy; ix++) {
            b_C->data[ix] = 0.0;
          }

          kAcol += i10;
        }

        br = 0;
        kAcol = 0;
        while ((i10 > 0) && (kAcol <= jBcol)) {
          nb = 0;
          iy = br + i9;
          for (b_n = br; b_n + 1 <= iy; b_n++) {
            if (py->data[b_n] != 0.0) {
              ia = nb;
              m = kAcol + i10;
              for (ix = kAcol; ix + 1 <= m; ix++) {
                ia++;
                mn = Q->size[0];
                b_C->data[ix] += py->data[b_n] * Q->data[(ia - 1) % mn + Q->
                  size[0] * ((ia - 1) / mn)];
              }
            }

            nb += i10;
          }

          br += i9;
          kAcol += i10;
        }
      }
    }

    emxFree_real_T(&py);
    if ((i8 - i7 == 1) || (pz->size[0] == 1)) {
      i7 = x->size[0];
      x->size[0] = Q2->size[0];
      emxEnsureCapacity((emxArray__common *)x, i7, (int)sizeof(double));
      br = Q2->size[0];
      for (i7 = 0; i7 < br; i7++) {
        x->data[i7] = 0.0;
        mn = Q2->size[1];
        for (i8 = 0; i8 < mn; i8++) {
          x->data[i7] += Q2->data[i7 + Q2->size[0] * i8] * pz->data[i8];
        }
      }
    } else {
      k = i8 - i7;
      i8 = Q->size[0];
      i9 = x->size[0];
      x->size[0] = i8;
      emxEnsureCapacity((emxArray__common *)x, i9, (int)sizeof(double));
      i8 = Q->size[0];
      jBcol = x->size[0];
      i9 = x->size[0];
      x->size[0] = jBcol;
      emxEnsureCapacity((emxArray__common *)x, i9, (int)sizeof(double));
      for (i9 = 0; i9 < jBcol; i9++) {
        x->data[i9] = 0.0;
      }

      i9 = Q->size[0];
      if (i9 != 0) {
        kAcol = 0;
        while ((i8 > 0) && (kAcol <= 0)) {
          for (ix = 1; ix <= i8; ix++) {
            x->data[ix - 1] = 0.0;
          }

          kAcol = i8;
        }

        br = 0;
        kAcol = 0;
        while ((i8 > 0) && (kAcol <= 0)) {
          nb = 0;
          i9 = br + k;
          for (b_n = br; b_n + 1 <= i9; b_n++) {
            if (pz->data[b_n] != 0.0) {
              ia = nb;
              for (ix = 0; ix + 1 <= i8; ix++) {
                ia++;
                i10 = Q->size[0];
                x->data[ix] += pz->data[b_n] * Q->data[(ia - 1) % i10 + Q->size
                  [0] * ((i7 + (ia - 1) / i10) - 1)];
              }
            }

            nb += i8;
          }

          br += k;
          kAcol = i8;
        }
      }
    }

    emxFree_real_T(&Q);
    emxFree_real_T(&Q2);
    i7 = x->size[0];
    x->size[0] = b_C->size[0];
    emxEnsureCapacity((emxArray__common *)x, i7, (int)sizeof(double));
    br = b_C->size[0];
    for (i7 = 0; i7 < br; i7++) {
      x->data[i7] += b_C->data[i7];
    }

    emxFree_real_T(&b_C);
    emxInit_real_T(&c_C, 1);
    if ((G->size[1] == 1) || (x->size[0] == 1)) {
      i7 = c_C->size[0];
      c_C->size[0] = G->size[0];
      emxEnsureCapacity((emxArray__common *)c_C, i7, (int)sizeof(double));
      br = G->size[0];
      for (i7 = 0; i7 < br; i7++) {
        c_C->data[i7] = 0.0;
        mn = G->size[1];
        for (i8 = 0; i8 < mn; i8++) {
          c_C->data[i7] += G->data[i7 + G->size[0] * i8] * x->data[i8];
        }
      }
    } else {
      k = G->size[1];
      unnamed_idx_0 = (unsigned int)G->size[0];
      i7 = c_C->size[0];
      c_C->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)c_C, i7, (int)sizeof(double));
      m = G->size[0];
      jBcol = c_C->size[0];
      i7 = c_C->size[0];
      c_C->size[0] = jBcol;
      emxEnsureCapacity((emxArray__common *)c_C, i7, (int)sizeof(double));
      for (i7 = 0; i7 < jBcol; i7++) {
        c_C->data[i7] = 0.0;
      }

      if (G->size[0] != 0) {
        kAcol = 0;
        while ((m > 0) && (kAcol <= 0)) {
          for (ix = 1; ix <= m; ix++) {
            c_C->data[ix - 1] = 0.0;
          }

          kAcol = m;
        }

        br = 0;
        kAcol = 0;
        while ((m > 0) && (kAcol <= 0)) {
          nb = 0;
          i7 = br + k;
          for (b_n = br; b_n + 1 <= i7; b_n++) {
            if (x->data[b_n] != 0.0) {
              ia = nb;
              for (ix = 0; ix + 1 <= m; ix++) {
                ia++;
                c_C->data[ix] += x->data[b_n] * G->data[ia - 1];
              }
            }

            nb += m;
          }

          br += k;
          kAcol = m;
        }
      }
    }

    i7 = Q2t->size[0] * Q2t->size[1];
    Q2t->size[0] = Q1->size[1];
    Q2t->size[1] = Q1->size[0];
    emxEnsureCapacity((emxArray__common *)Q2t, i7, (int)sizeof(double));
    br = Q1->size[0];
    for (i7 = 0; i7 < br; i7++) {
      mn = Q1->size[1];
      for (i8 = 0; i8 < mn; i8++) {
        Q2t->data[i8 + Q2t->size[0] * i7] = Q1->data[i7 + Q1->size[0] * i8];
      }
    }

    emxFree_real_T(&Q1);
    i7 = c_C->size[0];
    emxEnsureCapacity((emxArray__common *)c_C, i7, (int)sizeof(double));
    br = c_C->size[0];
    for (i7 = 0; i7 < br; i7++) {
      c_C->data[i7] += g->data[i7];
    }

    emxInit_real_T(&c_y, 1);
    if ((Q2t->size[1] == 1) || (c_C->size[0] == 1)) {
      i7 = c_y->size[0];
      c_y->size[0] = Q2t->size[0];
      emxEnsureCapacity((emxArray__common *)c_y, i7, (int)sizeof(double));
      br = Q2t->size[0];
      for (i7 = 0; i7 < br; i7++) {
        c_y->data[i7] = 0.0;
        mn = Q2t->size[1];
        for (i8 = 0; i8 < mn; i8++) {
          c_y->data[i7] += Q2t->data[i7 + Q2t->size[0] * i8] * c_C->data[i8];
        }
      }
    } else {
      k = Q2t->size[1];
      unnamed_idx_0 = (unsigned int)Q2t->size[0];
      i7 = c_y->size[0];
      c_y->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)c_y, i7, (int)sizeof(double));
      m = Q2t->size[0];
      jBcol = c_y->size[0];
      i7 = c_y->size[0];
      c_y->size[0] = jBcol;
      emxEnsureCapacity((emxArray__common *)c_y, i7, (int)sizeof(double));
      for (i7 = 0; i7 < jBcol; i7++) {
        c_y->data[i7] = 0.0;
      }

      kAcol = 0;
      while (kAcol <= 0) {
        for (ix = 1; ix <= m; ix++) {
          c_y->data[ix - 1] = 0.0;
        }

        kAcol = m;
      }

      br = 0;
      kAcol = 0;
      while (kAcol <= 0) {
        nb = 0;
        i7 = br + k;
        for (b_n = br; b_n + 1 <= i7; b_n++) {
          if (c_C->data[b_n] != 0.0) {
            ia = nb;
            for (ix = 0; ix + 1 <= m; ix++) {
              ia++;
              c_y->data[ix] += c_C->data[b_n] * Q2t->data[ia - 1];
            }
          }

          nb += m;
        }

        br += k;
        kAcol = m;
      }
    }

    emxFree_real_T(&c_C);
    emxFree_real_T(&Q2t);
    b_mldivide(R, c_y, pz);
    i7 = u->size[0] * u->size[1];
    u->size[0] = pz->size[0];
    u->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)u, i7, (int)sizeof(double));
    br = pz->size[0];
    emxFree_real_T(&c_y);
    emxFree_real_T(&R);
    for (i7 = 0; i7 < br; i7++) {
      u->data[i7] = pz->data[i7];
    }

    emxFree_real_T(&pz);
  } else {
    emxInit_real_T1(&Q1, 2);
    i7 = Q1->size[0] * Q1->size[1];
    Q1->size[0] = G->size[0];
    Q1->size[1] = G->size[1];
    emxEnsureCapacity((emxArray__common *)Q1, i7, (int)sizeof(double));
    br = G->size[0] * G->size[1];
    for (i7 = 0; i7 < br; i7++) {
      Q1->data[i7] = -G->data[i7];
    }

    b_mldivide(Q1, g, x);
    i7 = u->size[0] * u->size[1];
    u->size[0] = 0;
    u->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)u, i7, (int)sizeof(double));
    emxFree_real_T(&Q1);
  }
}

/*
 * min f(x)=1/2*x'*B*x - b'*x  %NB
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
  int i29;
  int br;
  emxArray_real_T *CW_d;
  int m_C;
  int i;
  emxArray_real_T *p;
  emxArray_real_T *r;
  emxArray_int32_T *I;
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
  emxArray_real_T *r13;
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
  int i30;
  boolean_T empty_non_axis_sizes;
  int ic;
  int idx;
  unsigned int uv0[2];
  double alpha;
  boolean_T exitg3;
  boolean_T exitg2;
  boolean_T exitg5;
  boolean_T guard1 = false;
  boolean_T exitg4;
  emxInit_real_T1(&W_A, 2);
  emxInit_real_T1(&W_d, 2);
  emxInit_real_T1(&CW_A, 2);

  /* !!!!!!!!!very important */
  flag = 0.0;

  /*  W.A = A([1],:); */
  /*  W.d = d([1]); */
  i29 = W_A->size[0] * W_A->size[1];
  W_A->size[0] = 0;
  W_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)W_A, i29, (int)sizeof(double));

  /* empty */
  i29 = W_d->size[0] * W_d->size[1];
  W_d->size[0] = 0;
  W_d->size[1] = 0;
  emxEnsureCapacity((emxArray__common *)W_d, i29, (int)sizeof(double));

  /*  */
  /*  [CW.A,IA] = setdiff(A,W.A,'rows','stable');%set1 - set2 */
  i29 = CW_A->size[0] * CW_A->size[1];
  CW_A->size[0] = A->size[0];
  CW_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)CW_A, i29, (int)sizeof(double));
  br = A->size[0] * A->size[1];
  for (i29 = 0; i29 < br; i29++) {
    CW_A->data[i29] = A->data[i29];
  }

  emxInit_real_T(&CW_d, 1);
  i29 = CW_d->size[0];
  CW_d->size[0] = d->size[0];
  emxEnsureCapacity((emxArray__common *)CW_d, i29, (int)sizeof(double));
  br = d->size[0];
  for (i29 = 0; i29 < br; i29++) {
    CW_d->data[i29] = d->data[i29];
  }

  /*  CW.d = d(IA); */
  m_C = C->size[0];
  i = 0;
  emxInit_real_T(&p, 1);
  emxInit_real_T1(&r, 2);
  emxInit_int32_T1(&I, 1);
  emxInit_real_T(&temp, 1);
  emxInit_real_T(&b_C, 1);
  emxInit_real_T(&c_C, 1);
  emxInit_real_T1(&y, 2);
  emxInit_boolean_T(&x, 1);
  emxInit_int32_T1(&ii, 1);
  emxInit_real_T1(&a, 2);
  emxInit_real_T(&b_y, 1);
  emxInit_real_T(&varargin_2, 1);
  emxInit_real_T(&d_C, 1);
  emxInit_real_T1(&e_C, 2);
  emxInit_real_T1(&r13, 2);
  emxInit_real_T(&f_C, 1);
  emxInit_real_T1(&b_c, 2);
  emxInit_real_T1(&b_CW_A, 2);
  emxInit_real_T1(&b_W_A, 2);
  emxInit_real_T1(&c_W_A, 2);
  emxInit_real_T1(&c_CW_A, 2);
  do {
    exitg1 = 0;
    if (i < 200) {
      i++;
      flag = i;

      /*  check_rank = (rank([C;W.A])==size([C;W.A],1)); */
      if ((B->size[1] == 1) || (x0->size[0] == 1)) {
        i29 = b_C->size[0];
        b_C->size[0] = B->size[0];
        emxEnsureCapacity((emxArray__common *)b_C, i29, (int)sizeof(double));
        br = B->size[0];
        for (i29 = 0; i29 < br; i29++) {
          b_C->data[i29] = 0.0;
          ar = B->size[1];
          for (i30 = 0; i30 < ar; i30++) {
            b_C->data[i29] += B->data[i29 + B->size[0] * i30] * x0->data[i30];
          }
        }
      } else {
        k = B->size[1];
        B_idx_0 = (unsigned int)B->size[0];
        i29 = b_C->size[0];
        b_C->size[0] = (int)B_idx_0;
        emxEnsureCapacity((emxArray__common *)b_C, i29, (int)sizeof(double));
        m = B->size[0];
        ia = b_C->size[0];
        i29 = b_C->size[0];
        b_C->size[0] = ia;
        emxEnsureCapacity((emxArray__common *)b_C, i29, (int)sizeof(double));
        for (i29 = 0; i29 < ia; i29++) {
          b_C->data[i29] = 0.0;
        }

        if (B->size[0] != 0) {
          ar = 0;
          while ((m > 0) && (ar <= 0)) {
            for (ic = 1; ic <= m; ic++) {
              b_C->data[ic - 1] = 0.0;
            }

            ar = m;
          }

          br = 0;
          ar = 0;
          while ((m > 0) && (ar <= 0)) {
            ar = 0;
            i29 = br + k;
            for (idx = br; idx + 1 <= i29; idx++) {
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
            ar = m;
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

      i29 = d_C->size[0];
      d_C->size[0] = b_C->size[0];
      emxEnsureCapacity((emxArray__common *)d_C, i29, (int)sizeof(double));
      br = b_C->size[0];
      for (i29 = 0; i29 < br; i29++) {
        d_C->data[i29] = b_C->data[i29] - b->data[i29];
      }

      i29 = e_C->size[0] * e_C->size[1];
      e_C->size[0] = ic + ar;
      e_C->size[1] = ia;
      emxEnsureCapacity((emxArray__common *)e_C, i29, (int)sizeof(double));
      for (i29 = 0; i29 < ia; i29++) {
        for (i30 = 0; i30 < ic; i30++) {
          e_C->data[i30 + e_C->size[0] * i29] = C->data[i30 + ic * i29];
        }
      }

      for (i29 = 0; i29 < ia; i29++) {
        for (i30 = 0; i30 < ar; i30++) {
          e_C->data[(i30 + ic) + e_C->size[0] * i29] = W_A->data[i30 + ar * i29];
        }
      }

      ia = c->size[0];
      i29 = b_c->size[0] * b_c->size[1];
      b_c->size[0] = ia + idx;
      b_c->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_c, i29, (int)sizeof(double));
      for (i29 = 0; i29 < 1; i29++) {
        for (i30 = 0; i30 < ia; i30++) {
          b_c->data[i30] = c->data[i30];
        }
      }

      for (i29 = 0; i29 < 1; i29++) {
        for (i30 = 0; i30 < idx; i30++) {
          b_c->data[i30 + ia] = W_d->data[i30];
        }
      }

      i29 = r13->size[0] * r13->size[1];
      r13->size[0] = b_c->size[0];
      r13->size[1] = b_c->size[1];
      emxEnsureCapacity((emxArray__common *)r13, i29, (int)sizeof(double));
      br = b_c->size[1];
      for (i29 = 0; i29 < br; i29++) {
        ar = b_c->size[0];
        for (i30 = 0; i30 < ar; i30++) {
          r13->data[i30 + r13->size[0] * i29] = 0.0 * b_c->data[i30 + b_c->size
            [0] * i29];
        }
      }

      null_space(B, d_C, e_C, r13, p, r);

      /* null space is better */
      c_abs(p, varargin_2);
      i29 = x->size[0];
      x->size[0] = varargin_2->size[0];
      emxEnsureCapacity((emxArray__common *)x, i29, (int)sizeof(boolean_T));
      br = varargin_2->size[0];
      for (i29 = 0; i29 < br; i29++) {
        x->data[i29] = (varargin_2->data[i29] > 1.0E-8);
      }

      i29 = p->size[0];
      emxEnsureCapacity((emxArray__common *)p, i29, (int)sizeof(double));
      br = p->size[0];
      for (i29 = 0; i29 < br; i29++) {
        p->data[i29] *= (double)x->data[i29];
      }

      for (i29 = 0; i29 < 2; i29++) {
        uv0[i29] = (unsigned int)r->size[i29];
      }

      i29 = y->size[0] * y->size[1];
      y->size[0] = (int)uv0[0];
      y->size[1] = (int)uv0[1];
      emxEnsureCapacity((emxArray__common *)y, i29, (int)sizeof(double));
      ia = r->size[0] * r->size[1];
      for (k = 0; k + 1 <= ia; k++) {
        y->data[k] = fabs(r->data[k]);
      }

      i29 = r->size[0] * r->size[1];
      emxEnsureCapacity((emxArray__common *)r, i29, (int)sizeof(double));
      ia = r->size[0];
      ar = r->size[1];
      br = ia * ar;
      for (i29 = 0; i29 < br; i29++) {
        r->data[i29] *= (double)(y->data[i29] > 1.0E-8);
      }

      /*  lambda = r(1:m_C); */
      i29 = r->size[0] * r->size[1];
      if (m_C + 1U > (unsigned int)i29) {
        i30 = 0;
        i29 = 1;
      } else {
        i30 = m_C;
        i29++;
      }

      c_abs(p, varargin_2);
      if (sum(varargin_2) == 0.0) {
        ar = 1;
        ia = (i29 - i30) - 1;
        alpha = r->data[i30];
        if ((i29 - i30) - 1 > 1) {
          if (rtIsNaN(alpha)) {
            idx = 2;
            exitg3 = false;
            while ((!exitg3) && (idx <= ia)) {
              ar = idx;
              if (!rtIsNaN(r->data[(i30 + idx) - 1])) {
                alpha = r->data[(i30 + idx) - 1];
                exitg3 = true;
              } else {
                idx++;
              }
            }
          }

          if (ar < (i29 - i30) - 1) {
            while (ar + 1 <= ia) {
              if (r->data[i30 + ar] < alpha) {
                alpha = r->data[i30 + ar];
              }

              ar++;
            }
          }
        }

        if (alpha >= 0.0) {
          exitg1 = 1;
        } else {
          ar = 1;
          ia = (i29 - i30) - 1;
          alpha = r->data[i30];
          k = 0;
          if ((i29 - i30) - 1 > 1) {
            if (rtIsNaN(alpha)) {
              idx = 2;
              exitg2 = false;
              while ((!exitg2) && (idx <= ia)) {
                ar = idx;
                if (!rtIsNaN(r->data[(i30 + idx) - 1])) {
                  alpha = r->data[(i30 + idx) - 1];
                  k = idx - 1;
                  exitg2 = true;
                } else {
                  idx++;
                }
              }
            }

            if (ar < (i29 - i30) - 1) {
              while (ar + 1 <= ia) {
                if (r->data[i30 + ar] < alpha) {
                  alpha = r->data[i30 + ar];
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
          i29 = c_W_A->size[0] * c_W_A->size[1];
          c_W_A->size[0] = 1;
          c_W_A->size[1] = br;
          emxEnsureCapacity((emxArray__common *)c_W_A, i29, (int)sizeof(double));
          for (i29 = 0; i29 < br; i29++) {
            c_W_A->data[c_W_A->size[0] * i29] = W_A->data[k + W_A->size[0] * i29];
          }

          i29 = c_CW_A->size[0] * c_CW_A->size[1];
          c_CW_A->size[0] = ic + 1;
          c_CW_A->size[1] = ia;
          emxEnsureCapacity((emxArray__common *)c_CW_A, i29, (int)sizeof(double));
          for (i29 = 0; i29 < ia; i29++) {
            for (i30 = 0; i30 < ic; i30++) {
              c_CW_A->data[i30 + c_CW_A->size[0] * i29] = CW_A->data[i30 + ic *
                i29];
            }
          }

          for (i29 = 0; i29 < ia; i29++) {
            for (i30 = 0; i30 < 1; i30++) {
              c_CW_A->data[ic + c_CW_A->size[0] * i29] = c_W_A->data[i29];
            }
          }

          i29 = CW_A->size[0] * CW_A->size[1];
          CW_A->size[0] = c_CW_A->size[0];
          CW_A->size[1] = c_CW_A->size[1];
          emxEnsureCapacity((emxArray__common *)CW_A, i29, (int)sizeof(double));
          br = c_CW_A->size[1];
          for (i29 = 0; i29 < br; i29++) {
            ar = c_CW_A->size[0];
            for (i30 = 0; i30 < ar; i30++) {
              CW_A->data[i30 + CW_A->size[0] * i29] = c_CW_A->data[i30 +
                c_CW_A->size[0] * i29];
            }
          }

          ia = CW_d->size[0];
          i29 = CW_d->size[0];
          CW_d->size[0] = ia + 1;
          emxEnsureCapacity((emxArray__common *)CW_d, i29, (int)sizeof(double));
          CW_d->data[ia] = W_d->data[k];
          nullAssignment(W_A, k + 1);

          /* delete constr */
          i29 = a->size[0] * a->size[1];
          a->size[0] = W_d->size[0];
          a->size[1] = W_d->size[1];
          emxEnsureCapacity((emxArray__common *)a, i29, (int)sizeof(double));
          br = W_d->size[0] * W_d->size[1];
          for (i29 = 0; i29 < br; i29++) {
            a->data[i29] = W_d->data[i29];
          }

          b_nullAssignment(a, k + 1);
          i29 = W_d->size[0] * W_d->size[1];
          W_d->size[0] = a->size[0];
          W_d->size[1] = a->size[1];
          emxEnsureCapacity((emxArray__common *)W_d, i29, (int)sizeof(double));
          br = a->size[0] * a->size[1];
          for (i29 = 0; i29 < br; i29++) {
            W_d->data[i29] = a->data[i29];
          }
        }
      } else {
        if ((CW_A->size[1] == 1) || (p->size[0] == 1)) {
          i29 = temp->size[0];
          temp->size[0] = CW_A->size[0];
          emxEnsureCapacity((emxArray__common *)temp, i29, (int)sizeof(double));
          br = CW_A->size[0];
          for (i29 = 0; i29 < br; i29++) {
            temp->data[i29] = 0.0;
            ar = CW_A->size[1];
            for (i30 = 0; i30 < ar; i30++) {
              temp->data[i29] += CW_A->data[i29 + CW_A->size[0] * i30] * p->
                data[i30];
            }
          }
        } else {
          k = CW_A->size[1];
          B_idx_0 = (unsigned int)CW_A->size[0];
          i29 = temp->size[0];
          temp->size[0] = (int)B_idx_0;
          emxEnsureCapacity((emxArray__common *)temp, i29, (int)sizeof(double));
          m = CW_A->size[0];
          ar = temp->size[0];
          i29 = temp->size[0];
          temp->size[0] = ar;
          emxEnsureCapacity((emxArray__common *)temp, i29, (int)sizeof(double));
          for (i29 = 0; i29 < ar; i29++) {
            temp->data[i29] = 0.0;
          }

          if (CW_A->size[0] != 0) {
            ar = 0;
            while ((m > 0) && (ar <= 0)) {
              for (ic = 1; ic <= m; ic++) {
                temp->data[ic - 1] = 0.0;
              }

              ar = m;
            }

            br = 0;
            ar = 0;
            while ((m > 0) && (ar <= 0)) {
              ar = 0;
              i29 = br + k;
              for (idx = br; idx + 1 <= i29; idx++) {
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
              ar = m;
            }
          }
        }

        i29 = x->size[0];
        x->size[0] = temp->size[0];
        emxEnsureCapacity((emxArray__common *)x, i29, (int)sizeof(boolean_T));
        br = temp->size[0];
        for (i29 = 0; i29 < br; i29++) {
          x->data[i29] = (temp->data[i29] < -1.0E-8);
        }

        ar = x->size[0];
        idx = 0;
        i29 = ii->size[0];
        ii->size[0] = x->size[0];
        emxEnsureCapacity((emxArray__common *)ii, i29, (int)sizeof(int));
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
            i29 = ii->size[0];
            ii->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)ii, i29, (int)sizeof(int));
          }
        } else {
          i29 = ii->size[0];
          if (1 > idx) {
            ii->size[0] = 0;
          } else {
            ii->size[0] = idx;
          }

          emxEnsureCapacity((emxArray__common *)ii, i29, (int)sizeof(int));
        }

        i29 = I->size[0];
        I->size[0] = ii->size[0];
        emxEnsureCapacity((emxArray__common *)I, i29, (int)sizeof(int));
        br = ii->size[0];
        for (i29 = 0; i29 < br; i29++) {
          I->data[i29] = ii->data[i29];
        }

        if (I->size[0] == 0) {
          alpha = 1.0;
        } else {
          br = CW_A->size[1];
          i29 = a->size[0] * a->size[1];
          a->size[0] = I->size[0];
          a->size[1] = br;
          emxEnsureCapacity((emxArray__common *)a, i29, (int)sizeof(double));
          for (i29 = 0; i29 < br; i29++) {
            ar = I->size[0];
            for (i30 = 0; i30 < ar; i30++) {
              a->data[i30 + a->size[0] * i29] = CW_A->data[(I->data[i30] +
                CW_A->size[0] * i29) - 1];
            }
          }

          i29 = CW_A->size[1];
          if ((i29 == 1) || (x0->size[0] == 1)) {
            i29 = c_C->size[0];
            c_C->size[0] = a->size[0];
            emxEnsureCapacity((emxArray__common *)c_C, i29, (int)sizeof(double));
            br = a->size[0];
            for (i29 = 0; i29 < br; i29++) {
              c_C->data[i29] = 0.0;
              ar = a->size[1];
              for (i30 = 0; i30 < ar; i30++) {
                c_C->data[i29] += a->data[i29 + a->size[0] * i30] * x0->data[i30];
              }
            }
          } else {
            i29 = CW_A->size[1];
            B_idx_0 = (unsigned int)I->size[0];
            i30 = c_C->size[0];
            c_C->size[0] = (int)B_idx_0;
            emxEnsureCapacity((emxArray__common *)c_C, i30, (int)sizeof(double));
            m = I->size[0];
            ia = c_C->size[0];
            i30 = c_C->size[0];
            c_C->size[0] = ia;
            emxEnsureCapacity((emxArray__common *)c_C, i30, (int)sizeof(double));
            for (i30 = 0; i30 < ia; i30++) {
              c_C->data[i30] = 0.0;
            }

            ar = 0;
            while (ar <= 0) {
              for (ic = 1; ic <= m; ic++) {
                c_C->data[ic - 1] = 0.0;
              }

              ar = m;
            }

            br = 0;
            ar = 0;
            while (ar <= 0) {
              ar = 0;
              i30 = br + i29;
              for (idx = br; idx + 1 <= i30; idx++) {
                if (x0->data[idx] != 0.0) {
                  ia = ar;
                  for (ic = 0; ic + 1 <= m; ic++) {
                    ia++;
                    c_C->data[ic] += x0->data[idx] * a->data[ia - 1];
                  }
                }

                ar += m;
              }

              br += i29;
              ar = m;
            }
          }

          br = CW_A->size[1];
          i29 = a->size[0] * a->size[1];
          a->size[0] = I->size[0];
          a->size[1] = br;
          emxEnsureCapacity((emxArray__common *)a, i29, (int)sizeof(double));
          for (i29 = 0; i29 < br; i29++) {
            ar = I->size[0];
            for (i30 = 0; i30 < ar; i30++) {
              a->data[i30 + a->size[0] * i29] = CW_A->data[(I->data[i30] +
                CW_A->size[0] * i29) - 1];
            }
          }

          i29 = CW_A->size[1];
          if ((i29 == 1) || (p->size[0] == 1)) {
            i29 = b_y->size[0];
            b_y->size[0] = a->size[0];
            emxEnsureCapacity((emxArray__common *)b_y, i29, (int)sizeof(double));
            br = a->size[0];
            for (i29 = 0; i29 < br; i29++) {
              b_y->data[i29] = 0.0;
              ar = a->size[1];
              for (i30 = 0; i30 < ar; i30++) {
                b_y->data[i29] += a->data[i29 + a->size[0] * i30] * p->data[i30];
              }
            }
          } else {
            i29 = CW_A->size[1];
            B_idx_0 = (unsigned int)I->size[0];
            i30 = b_y->size[0];
            b_y->size[0] = (int)B_idx_0;
            emxEnsureCapacity((emxArray__common *)b_y, i30, (int)sizeof(double));
            m = I->size[0];
            ia = b_y->size[0];
            i30 = b_y->size[0];
            b_y->size[0] = ia;
            emxEnsureCapacity((emxArray__common *)b_y, i30, (int)sizeof(double));
            for (i30 = 0; i30 < ia; i30++) {
              b_y->data[i30] = 0.0;
            }

            ar = 0;
            while (ar <= 0) {
              for (ic = 1; ic <= m; ic++) {
                b_y->data[ic - 1] = 0.0;
              }

              ar = m;
            }

            br = 0;
            ar = 0;
            while (ar <= 0) {
              ar = 0;
              i30 = br + i29;
              for (idx = br; idx + 1 <= i30; idx++) {
                if (p->data[idx] != 0.0) {
                  ia = ar;
                  for (ic = 0; ic + 1 <= m; ic++) {
                    ia++;
                    b_y->data[ic] += p->data[idx] * a->data[ia - 1];
                  }
                }

                ar += m;
              }

              br += i29;
              ar = m;
            }
          }

          i29 = f_C->size[0];
          f_C->size[0] = c_C->size[0];
          emxEnsureCapacity((emxArray__common *)f_C, i29, (int)sizeof(double));
          br = c_C->size[0];
          for (i29 = 0; i29 < br; i29++) {
            f_C->data[i29] = c_C->data[i29] - CW_d->data[I->data[i29] - 1];
          }

          rdivide(f_C, b_y, varargin_2);
          ar = 1;
          ia = varargin_2->size[0];
          alpha = varargin_2->data[0];
          k = 0;
          if (varargin_2->size[0] > 1) {
            if (rtIsNaN(varargin_2->data[0])) {
              idx = 1;
              exitg4 = false;
              while ((!exitg4) && (idx + 1 <= ia)) {
                ar = idx + 1;
                if (!rtIsNaN(varargin_2->data[idx])) {
                  alpha = varargin_2->data[idx];
                  k = idx;
                  exitg4 = true;
                } else {
                  idx++;
                }
              }
            }

            if (ar < varargin_2->size[0]) {
              while (ar + 1 <= ia) {
                if (varargin_2->data[ar] < alpha) {
                  alpha = varargin_2->data[ar];
                  k = ar;
                }

                ar++;
              }
            }
          }

          if ((1.0 <= alpha) || rtIsNaN(alpha)) {
            alpha = 1.0;
          }

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
            ar = I->data[k];
            i29 = b_CW_A->size[0] * b_CW_A->size[1];
            b_CW_A->size[0] = 1;
            b_CW_A->size[1] = br;
            emxEnsureCapacity((emxArray__common *)b_CW_A, i29, (int)sizeof
                              (double));
            for (i29 = 0; i29 < br; i29++) {
              b_CW_A->data[b_CW_A->size[0] * i29] = CW_A->data[(ar + CW_A->size
                [0] * i29) - 1];
            }

            i29 = b_W_A->size[0] * b_W_A->size[1];
            b_W_A->size[0] = ic + 1;
            b_W_A->size[1] = ia;
            emxEnsureCapacity((emxArray__common *)b_W_A, i29, (int)sizeof(double));
            for (i29 = 0; i29 < ia; i29++) {
              for (i30 = 0; i30 < ic; i30++) {
                b_W_A->data[i30 + b_W_A->size[0] * i29] = W_A->data[i30 + ic *
                  i29];
              }
            }

            for (i29 = 0; i29 < ia; i29++) {
              for (i30 = 0; i30 < 1; i30++) {
                b_W_A->data[ic + b_W_A->size[0] * i29] = b_CW_A->data[i29];
              }
            }

            i29 = W_A->size[0] * W_A->size[1];
            W_A->size[0] = b_W_A->size[0];
            W_A->size[1] = b_W_A->size[1];
            emxEnsureCapacity((emxArray__common *)W_A, i29, (int)sizeof(double));
            br = b_W_A->size[1];
            for (i29 = 0; i29 < br; i29++) {
              ar = b_W_A->size[0];
              for (i30 = 0; i30 < ar; i30++) {
                W_A->data[i30 + W_A->size[0] * i29] = b_W_A->data[i30 +
                  b_W_A->size[0] * i29];
              }
            }

            /* add constr */
            if (!((W_d->size[0] == 0) || (W_d->size[1] == 0))) {
              ia = W_d->size[0];
            } else {
              ia = 0;
            }

            i29 = varargin_2->size[0];
            varargin_2->size[0] = ia + 1;
            emxEnsureCapacity((emxArray__common *)varargin_2, i29, (int)sizeof
                              (double));
            for (i29 = 0; i29 < ia; i29++) {
              varargin_2->data[i29] = W_d->data[i29];
            }

            varargin_2->data[ia] = CW_d->data[I->data[k] - 1];
            i29 = W_d->size[0] * W_d->size[1];
            W_d->size[0] = varargin_2->size[0];
            W_d->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)W_d, i29, (int)sizeof(double));
            br = varargin_2->size[0];
            for (i29 = 0; i29 < br; i29++) {
              W_d->data[i29] = varargin_2->data[i29];
            }

            nullAssignment(CW_A, I->data[k]);

            /* delete */
            c_nullAssignment(CW_d, I->data[k]);
          }
        }

        i29 = x0->size[0];
        emxEnsureCapacity((emxArray__common *)x0, i29, (int)sizeof(double));
        br = x0->size[0];
        for (i29 = 0; i29 < br; i29++) {
          x0->data[i29] -= alpha * p->data[i29];
        }

        /* NB */
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
  emxFree_real_T(&r13);
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
  emxFree_int32_T(&I);
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
