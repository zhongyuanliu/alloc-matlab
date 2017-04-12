/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpkwik.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "call_qpsolver.h"
#include "qpkwik.h"
#include "call_qpsolver_emxutil.h"
#include "qr.h"
#include "abs.h"
#include "norm.h"

/* Function Declarations */
static void DropConstraint(short kDrop, emxArray_int16_T *b_iA, short *nA, short
  iC_data[]);
static double KWIKfactor(const emxArray_real_T *Ac, const short iC_data[], short
  nA, const emxArray_real_T *Linv, emxArray_real_T *RLinv, emxArray_real_T *D,
  emxArray_real_T *H, short n);
static void ResetToColdStart(short m, short meq, short *nA, short iA_data[], int
  iA_size[1], short iC_data[], int iC_size[1]);
static void Unconstrained(const emxArray_real_T *Hinv, const emxArray_real_T *f,
  emxArray_real_T *x, short n);

/* Function Definitions */

/*
 * Arguments    : short kDrop
 *                emxArray_int16_T *b_iA
 *                short *nA
 *                short iC_data[]
 * Return Type  : void
 */
static void DropConstraint(short kDrop, emxArray_int16_T *b_iA, short *nA, short
  iC_data[])
{
  int i25;
  short i26;
  short i;
  b_iA->data[iC_data[kDrop - 1] - 1] = 0;
  if (kDrop < *nA) {
    i25 = *nA - 1;
    if (i25 < -32768) {
      i25 = -32768;
    }

    i26 = (short)i25;
    for (i = kDrop; i <= i26; i++) {
      iC_data[i - 1] = iC_data[i];
    }
  }

  iC_data[*nA - 1] = 0;
  i25 = *nA - 1;
  if (i25 < -32768) {
    i25 = -32768;
  }

  *nA = (short)i25;
}

/*
 * Arguments    : const emxArray_real_T *Ac
 *                const short iC_data[]
 *                short nA
 *                const emxArray_real_T *Linv
 *                emxArray_real_T *RLinv
 *                emxArray_real_T *D
 *                emxArray_real_T *H
 *                short n
 * Return Type  : double
 */
static double KWIKfactor(const emxArray_real_T *Ac, const short iC_data[], short
  nA, const emxArray_real_T *Linv, emxArray_real_T *RLinv, emxArray_real_T *D,
  emxArray_real_T *H, short n)
{
  double Status;
  emxArray_real_T *TL;
  int i22;
  int ib;
  int ia;
  int ar;
  short i;
  emxArray_real_T *C;
  emxArray_real_T *b;
  emxArray_real_T *QQ;
  emxArray_real_T *RR;
  emxArray_real_T *a;
  int exitg1;
  int k;
  unsigned int Linv_idx_0;
  short j;
  int m;
  short b_k;
  int br;
  int ic;
  short i23;
  boolean_T guard1 = false;
  double y;
  emxInit_real_T1(&TL, 2);
  for (i22 = 0; i22 < 2; i22++) {
    ib = TL->size[0] * TL->size[1];
    TL->size[i22] = Linv->size[i22];
    emxEnsureCapacity((emxArray__common *)TL, ib, (int)sizeof(double));
  }

  Status = 1.0;
  ia = RLinv->size[0];
  ar = RLinv->size[1];
  for (i22 = 0; i22 < ar; i22++) {
    for (ib = 0; ib < ia; ib++) {
      RLinv->data[ib + RLinv->size[0] * i22] = 0.0;
    }
  }

  i = 1;
  emxInit_real_T(&C, 1);
  emxInit_real_T(&b, 1);
  while (i <= nA) {
    ia = Ac->size[1];
    i22 = b->size[0];
    b->size[0] = ia;
    emxEnsureCapacity((emxArray__common *)b, i22, (int)sizeof(double));
    for (i22 = 0; i22 < ia; i22++) {
      b->data[i22] = Ac->data[(iC_data[i - 1] + Ac->size[0] * i22) - 1];
    }

    if ((Linv->size[1] == 1) || (b->size[0] == 1)) {
      i22 = C->size[0];
      C->size[0] = Linv->size[0];
      emxEnsureCapacity((emxArray__common *)C, i22, (int)sizeof(double));
      ia = Linv->size[0];
      for (i22 = 0; i22 < ia; i22++) {
        C->data[i22] = 0.0;
        ar = Linv->size[1];
        for (ib = 0; ib < ar; ib++) {
          C->data[i22] += Linv->data[i22 + Linv->size[0] * ib] * b->data[ib];
        }
      }
    } else {
      k = Linv->size[1];
      Linv_idx_0 = (unsigned int)Linv->size[0];
      i22 = C->size[0];
      C->size[0] = (int)Linv_idx_0;
      emxEnsureCapacity((emxArray__common *)C, i22, (int)sizeof(double));
      m = Linv->size[0];
      ia = C->size[0];
      i22 = C->size[0];
      C->size[0] = ia;
      emxEnsureCapacity((emxArray__common *)C, i22, (int)sizeof(double));
      for (i22 = 0; i22 < ia; i22++) {
        C->data[i22] = 0.0;
      }

      if (Linv->size[0] == 0) {
      } else {
        ia = 0;
        while ((m > 0) && (ia <= 0)) {
          for (ic = 1; ic <= m; ic++) {
            C->data[ic - 1] = 0.0;
          }

          ia = m;
        }

        br = 0;
        ia = 0;
        while ((m > 0) && (ia <= 0)) {
          ar = -1;
          i22 = br + k;
          for (ib = br; ib + 1 <= i22; ib++) {
            if (b->data[ib] != 0.0) {
              ia = ar;
              for (ic = 0; ic + 1 <= m; ic++) {
                ia++;
                C->data[ic] += b->data[ib] * Linv->data[ia];
              }
            }

            ar += m;
          }

          br += k;
          ia = m;
        }
      }
    }

    ia = C->size[0];
    for (i22 = 0; i22 < ia; i22++) {
      RLinv->data[i22 + RLinv->size[0] * (i - 1)] = C->data[i22];
    }

    i++;
  }

  emxFree_real_T(&C);
  emxInit_real_T1(&QQ, 2);
  emxInit_real_T1(&RR, 2);
  qr(RLinv, QQ, RR);
  i = 1;
  emxInit_real_T1(&a, 2);
  do {
    exitg1 = 0;
    if (i <= nA) {
      if (fabs(RR->data[(i + RR->size[0] * (i - 1)) - 1]) < 1.0E-12) {
        Status = -2.0;
        exitg1 = 1;
      } else {
        i++;
      }
    } else {
      for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
          ia = Linv->size[0];
          i22 = a->size[0] * a->size[1];
          a->size[0] = 1;
          a->size[1] = ia;
          emxEnsureCapacity((emxArray__common *)a, i22, (int)sizeof(double));
          for (i22 = 0; i22 < ia; i22++) {
            a->data[a->size[0] * i22] = Linv->data[i22 + Linv->size[0] * (i - 1)];
          }

          ia = QQ->size[0];
          i22 = b->size[0];
          b->size[0] = ia;
          emxEnsureCapacity((emxArray__common *)b, i22, (int)sizeof(double));
          for (i22 = 0; i22 < ia; i22++) {
            b->data[i22] = QQ->data[i22 + QQ->size[0] * (j - 1)];
          }

          guard1 = false;
          if (a->size[1] == 1) {
            guard1 = true;
          } else {
            i22 = QQ->size[0];
            if (i22 == 1) {
              guard1 = true;
            } else {
              y = 0.0;
              for (i22 = 0; i22 < a->size[1]; i22++) {
                y += a->data[a->size[0] * i22] * b->data[i22];
              }
            }
          }

          if (guard1) {
            y = 0.0;
            for (i22 = 0; i22 < a->size[1]; i22++) {
              y += a->data[a->size[0] * i22] * b->data[i22];
            }
          }

          TL->data[(i + TL->size[0] * (j - 1)) - 1] = y;
        }
      }

      ia = RLinv->size[0];
      ar = RLinv->size[1];
      for (i22 = 0; i22 < ar; i22++) {
        for (ib = 0; ib < ia; ib++) {
          RLinv->data[ib + RLinv->size[0] * i22] = 0.0;
        }
      }

      for (j = nA; j > 0; j--) {
        RLinv->data[(j + RLinv->size[0] * (j - 1)) - 1] = 1.0;
        for (b_k = j; b_k <= nA; b_k++) {
          RLinv->data[(j + RLinv->size[0] * (b_k - 1)) - 1] /= RR->data[(j +
            RR->size[0] * (j - 1)) - 1];
        }

        if (j > 1) {
          i23 = (short)(j - 1);
          for (i = 1; i <= i23; i++) {
            for (b_k = j; b_k <= nA; b_k++) {
              RLinv->data[(i + RLinv->size[0] * (b_k - 1)) - 1] -= RR->data[(i +
                RR->size[0] * (j - 1)) - 1] * RLinv->data[(j + RLinv->size[0] *
                (b_k - 1)) - 1];
            }
          }
        }
      }

      for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
          H->data[(i + H->size[0] * (j - 1)) - 1] = 0.0;
          for (b_k = (short)(nA + 1); b_k <= n; b_k++) {
            H->data[(i + H->size[0] * (j - 1)) - 1] -= TL->data[(i + TL->size[0]
              * (b_k - 1)) - 1] * TL->data[(j + TL->size[0] * (b_k - 1)) - 1];
          }

          H->data[(j + H->size[0] * (i - 1)) - 1] = H->data[(i + H->size[0] * (j
            - 1)) - 1];
        }
      }

      for (j = 1; j <= nA; j++) {
        for (i = 1; i <= n; i++) {
          D->data[(i + D->size[0] * (j - 1)) - 1] = 0.0;
          for (b_k = j; b_k <= nA; b_k++) {
            D->data[(i + D->size[0] * (j - 1)) - 1] += TL->data[(i + TL->size[0]
              * (b_k - 1)) - 1] * RLinv->data[(j + RLinv->size[0] * (b_k - 1)) -
              1];
          }
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&RR);
  emxFree_real_T(&QQ);
  emxFree_real_T(&TL);
  return Status;
}

/*
 * Arguments    : short m
 *                short meq
 *                short *nA
 *                short iA_data[]
 *                int iA_size[1]
 *                short iC_data[]
 *                int iC_size[1]
 * Return Type  : void
 */
static void ResetToColdStart(short m, short meq, short *nA, short iA_data[], int
  iA_size[1], short iC_data[], int iC_size[1])
{
  int loop_ub;
  int i9;
  short i;
  short ix;
  *nA = meq;
  iA_size[0] = m;
  loop_ub = m;
  for (i9 = 0; i9 < loop_ub; i9++) {
    iA_data[i9] = 0;
  }

  iC_size[0] = m;
  loop_ub = m;
  for (i9 = 0; i9 < loop_ub; i9++) {
    iC_data[i9] = 0;
  }

  for (i = 1; i <= meq; i++) {
    i9 = m - meq;
    if (i9 > 32767) {
      i9 = 32767;
    } else {
      if (i9 < -32768) {
        i9 = -32768;
      }
    }

    i9 += i;
    if (i9 > 32767) {
      i9 = 32767;
    } else {
      if (i9 < -32768) {
        i9 = -32768;
      }
    }

    ix = (short)i9;
    iA_data[ix - 1] = 1;
    iC_data[i - 1] = ix;
  }
}

/*
 * Arguments    : const emxArray_real_T *Hinv
 *                const emxArray_real_T *f
 *                emxArray_real_T *x
 *                short n
 * Return Type  : void
 */
static void Unconstrained(const emxArray_real_T *Hinv, const emxArray_real_T *f,
  emxArray_real_T *x, short n)
{
  short i;
  emxArray_real_T *a;
  int loop_ub;
  int i21;
  double y;
  i = 1;
  emxInit_real_T1(&a, 2);
  while (i <= n) {
    loop_ub = Hinv->size[1];
    i21 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)a, i21, (int)sizeof(double));
    for (i21 = 0; i21 < loop_ub; i21++) {
      a->data[a->size[0] * i21] = -Hinv->data[(i + Hinv->size[0] * i21) - 1];
    }

    if ((a->size[1] == 1) || (f->size[0] == 1)) {
      y = 0.0;
      for (i21 = 0; i21 < a->size[1]; i21++) {
        y += a->data[a->size[0] * i21] * f->data[i21];
      }
    } else {
      y = 0.0;
      for (i21 = 0; i21 < a->size[1]; i21++) {
        y += a->data[a->size[0] * i21] * f->data[i21];
      }
    }

    x->data[i - 1] = y;
    i++;
  }

  emxFree_real_T(&a);
}

/*
 * Arguments    : const emxArray_real_T *Linv
 *                const emxArray_real_T *Hinv
 *                const emxArray_real_T *f
 *                const emxArray_real_T *Ac
 *                const emxArray_real_T *b
 *                emxArray_int16_T *b_iA
 *                short m
 *                short n
 *                short meq
 *                emxArray_real_T *x
 *                emxArray_real_T *lambda
 *                double *status
 * Return Type  : void
 */
void qpkwik(const emxArray_real_T *Linv, const emxArray_real_T *Hinv, const
            emxArray_real_T *f, const emxArray_real_T *Ac, const emxArray_real_T
            *b, emxArray_int16_T *b_iA, short m, short n, short meq,
            emxArray_real_T *x, emxArray_real_T *lambda, double *status)
{
  int i20;
  int ar;
  emxArray_real_T *r;
  emxArray_real_T *RLinv;
  double rMin;
  emxArray_real_T *D;
  int ia;
  emxArray_real_T *H;
  emxArray_real_T *U;
  emxArray_real_T *cTol;
  boolean_T cTolComputed;
  short nA;
  short iC_data[32767];
  short i;
  emxArray_real_T *Opt;
  emxArray_real_T *Rhs;
  emxArray_real_T *AcRow;
  emxArray_real_T *z;
  emxArray_real_T *b_b;
  emxArray_real_T *a;
  emxArray_real_T *varargin_2;
  emxArray_real_T *b_Ac;
  boolean_T guard1 = false;
  short iSave;
  double Xnorm0;
  int exitg2;
  double cMin;
  short kNext;
  short k;
  boolean_T DualFeasible;
  int exitg1;
  double Xnorm;
  unsigned short q;
  unsigned short b_x;
  unsigned int Rhs_idx_0;
  unsigned int b_Rhs_idx_0;
  int b_n;
  boolean_T ColdReset;
  boolean_T guard2 = false;
  double y;
  boolean_T guard3 = false;
  int exitg3;
  int b_k;
  boolean_T exitg5;
  short iA_data[32767];
  int iA_size[1];
  int iC_size[1];
  short kDrop;
  double t1;
  int b_m;
  boolean_T exitg4;
  int br;
  int ic;
  boolean_T b_guard1 = false;
  *status = 1.0;
  i20 = lambda->size[0];
  lambda->size[0] = m;
  emxEnsureCapacity((emxArray__common *)lambda, i20, (int)sizeof(double));
  ar = m;
  for (i20 = 0; i20 < ar; i20++) {
    lambda->data[i20] = 0.0;
  }

  i20 = x->size[0];
  x->size[0] = n;
  emxEnsureCapacity((emxArray__common *)x, i20, (int)sizeof(double));
  ar = n;
  for (i20 = 0; i20 < ar; i20++) {
    x->data[i20] = 0.0;
  }

  emxInit_real_T(&r, 1);
  i20 = r->size[0];
  r->size[0] = n;
  emxEnsureCapacity((emxArray__common *)r, i20, (int)sizeof(double));
  ar = n;
  for (i20 = 0; i20 < ar; i20++) {
    r->data[i20] = 0.0;
  }

  emxInit_real_T1(&RLinv, 2);
  rMin = 0.0;
  for (i20 = 0; i20 < 2; i20++) {
    ia = RLinv->size[0] * RLinv->size[1];
    RLinv->size[i20] = Linv->size[i20];
    emxEnsureCapacity((emxArray__common *)RLinv, ia, (int)sizeof(double));
  }

  emxInit_real_T1(&D, 2);
  for (i20 = 0; i20 < 2; i20++) {
    ia = D->size[0] * D->size[1];
    D->size[i20] = Linv->size[i20];
    emxEnsureCapacity((emxArray__common *)D, ia, (int)sizeof(double));
  }

  emxInit_real_T1(&H, 2);
  for (i20 = 0; i20 < 2; i20++) {
    ia = H->size[0] * H->size[1];
    H->size[i20] = Linv->size[i20];
    emxEnsureCapacity((emxArray__common *)H, ia, (int)sizeof(double));
  }

  emxInit_real_T1(&U, 2);
  for (i20 = 0; i20 < 2; i20++) {
    ia = U->size[0] * U->size[1];
    U->size[i20] = Linv->size[i20];
    emxEnsureCapacity((emxArray__common *)U, ia, (int)sizeof(double));
  }

  emxInit_real_T(&cTol, 1);
  i20 = cTol->size[0];
  cTol->size[0] = m;
  emxEnsureCapacity((emxArray__common *)cTol, i20, (int)sizeof(double));
  ar = m;
  for (i20 = 0; i20 < ar; i20++) {
    cTol->data[i20] = 1.0;
  }

  cTolComputed = false;
  ar = m;
  for (i20 = 0; i20 < ar; i20++) {
    iC_data[i20] = 0;
  }

  nA = 0;
  for (i = 1; i <= m; i++) {
    if (b_iA->data[i - 1] == 1) {
      i20 = nA + 1;
      if (i20 > 32767) {
        i20 = 32767;
      }

      nA = (short)i20;
      iC_data[nA - 1] = i;
    }
  }

  emxInit_real_T(&Opt, 1);
  emxInit_real_T(&Rhs, 1);
  emxInit_real_T1(&AcRow, 2);
  emxInit_real_T(&z, 1);
  emxInit_real_T(&b_b, 1);
  emxInit_real_T1(&a, 2);
  emxInit_real_T1(&varargin_2, 2);
  emxInit_real_T1(&b_Ac, 2);
  guard1 = false;
  if (nA > 0) {
    if (n > 16383) {
      iSave = MAX_int16_T;
    } else if (n <= -16384) {
      iSave = MIN_int16_T;
    } else {
      iSave = (short)(n << 1);
    }

    i20 = Opt->size[0];
    Opt->size[0] = iSave;
    emxEnsureCapacity((emxArray__common *)Opt, i20, (int)sizeof(double));
    ar = iSave;
    for (i20 = 0; i20 < ar; i20++) {
      Opt->data[i20] = 0.0;
    }

    i20 = Rhs->size[0];
    Rhs->size[0] = f->size[0] + n;
    emxEnsureCapacity((emxArray__common *)Rhs, i20, (int)sizeof(double));
    ar = f->size[0];
    for (i20 = 0; i20 < ar; i20++) {
      Rhs->data[i20] = f->data[i20];
    }

    ar = n;
    for (i20 = 0; i20 < ar; i20++) {
      Rhs->data[i20 + f->size[0]] = 0.0;
    }

    DualFeasible = false;
    i20 = 3 * nA;
    if (i20 > 32767) {
      i20 = 32767;
    }

    iSave = (short)i20;
    if (iSave >= 50) {
    } else {
      iSave = 50;
    }

    q = (unsigned short)(iSave / 10U);
    b_x = (unsigned short)((unsigned int)iSave - q * 10);
    if ((b_x > 0) && (b_x >= 5)) {
      q++;
    }

    ColdReset = false;
    do {
      exitg3 = 0;
      if ((!DualFeasible) && (nA > 0) && (*status <= 100.0)) {
        Xnorm = KWIKfactor(Ac, iC_data, nA, Linv, RLinv, D, H, n);
        if (Xnorm < 0.0) {
          if (ColdReset) {
            *status = -2.0;
            exitg3 = 2;
          } else {
            ResetToColdStart(m, meq, &nA, iA_data, iA_size, iC_data, iC_size);
            i20 = b_iA->size[0];
            b_iA->size[0] = iA_size[0];
            emxEnsureCapacity((emxArray__common *)b_iA, i20, (int)sizeof(short));
            ar = iA_size[0];
            for (i20 = 0; i20 < ar; i20++) {
              b_iA->data[i20] = iA_data[i20];
            }

            ColdReset = true;
          }
        } else {
          for (iSave = 1; iSave <= nA; iSave++) {
            i20 = n + iSave;
            if (i20 > 32767) {
              i20 = 32767;
            } else {
              if (i20 < -32768) {
                i20 = -32768;
              }
            }

            Rhs->data[i20 - 1] = b->data[iC_data[iSave - 1] - 1];
            for (i = iSave; i <= nA; i++) {
              U->data[(i + U->size[0] * (iSave - 1)) - 1] = 0.0;
              for (k = 1; k <= nA; k++) {
                U->data[(i + U->size[0] * (iSave - 1)) - 1] += RLinv->data[(i +
                  RLinv->size[0] * (k - 1)) - 1] * RLinv->data[(iSave +
                  RLinv->size[0] * (k - 1)) - 1];
              }

              U->data[(iSave + U->size[0] * (i - 1)) - 1] = U->data[(i + U->
                size[0] * (iSave - 1)) - 1];
            }
          }

          for (i = 1; i <= n; i++) {
            if (1 > n) {
              ar = 0;
            } else {
              ar = n;
            }

            b_n = H->size[1];
            i20 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = b_n;
            emxEnsureCapacity((emxArray__common *)a, i20, (int)sizeof(double));
            for (i20 = 0; i20 < b_n; i20++) {
              a->data[a->size[0] * i20] = H->data[(i + H->size[0] * i20) - 1];
            }

            i20 = b_b->size[0];
            b_b->size[0] = ar;
            emxEnsureCapacity((emxArray__common *)b_b, i20, (int)sizeof(double));
            for (i20 = 0; i20 < ar; i20++) {
              b_b->data[i20] = Rhs->data[i20];
            }

            i20 = H->size[1];
            if ((i20 == 1) || (ar == 1)) {
              y = 0.0;
              for (i20 = 0; i20 < a->size[1]; i20++) {
                y += a->data[a->size[0] * i20] * b_b->data[i20];
              }
            } else {
              y = 0.0;
              for (i20 = 0; i20 < a->size[1]; i20++) {
                y += a->data[a->size[0] * i20] * b_b->data[i20];
              }
            }

            Opt->data[i - 1] = y;
            for (k = 1; k <= nA; k++) {
              i20 = n + k;
              if (i20 > 32767) {
                i20 = 32767;
              }

              Opt->data[i - 1] += D->data[(i + D->size[0] * (k - 1)) - 1] *
                Rhs->data[i20 - 1];
            }
          }

          for (i = 1; i <= nA; i++) {
            if (1 > n) {
              ar = 0;
            } else {
              ar = n;
            }

            b_n = D->size[0];
            i20 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = b_n;
            emxEnsureCapacity((emxArray__common *)a, i20, (int)sizeof(double));
            for (i20 = 0; i20 < b_n; i20++) {
              a->data[a->size[0] * i20] = D->data[i20 + D->size[0] * (i - 1)];
            }

            i20 = b_b->size[0];
            b_b->size[0] = ar;
            emxEnsureCapacity((emxArray__common *)b_b, i20, (int)sizeof(double));
            for (i20 = 0; i20 < ar; i20++) {
              b_b->data[i20] = Rhs->data[i20];
            }

            if ((a->size[1] == 1) || (ar == 1)) {
              y = 0.0;
              for (i20 = 0; i20 < a->size[1]; i20++) {
                y += a->data[a->size[0] * i20] * b_b->data[i20];
              }
            } else {
              y = 0.0;
              for (i20 = 0; i20 < a->size[1]; i20++) {
                y += a->data[a->size[0] * i20] * b_b->data[i20];
              }
            }

            i20 = n + i;
            if (i20 > 32767) {
              i20 = 32767;
            }

            Opt->data[i20 - 1] = y;
            for (k = 1; k <= nA; k++) {
              i20 = n + i;
              if (i20 > 32767) {
                i20 = 32767;
              }

              ia = n + i;
              if (ia > 32767) {
                ia = 32767;
              }

              ar = n + k;
              if (ar > 32767) {
                ar = 32767;
              }

              Opt->data[i20 - 1] = Opt->data[ia - 1] + U->data[(i + U->size[0] *
                (k - 1)) - 1] * Rhs->data[ar - 1];
            }
          }

          Xnorm = -1.0E-12;
          kDrop = 0;
          for (i = 1; i <= nA; i++) {
            i20 = n + i;
            if (i20 > 32767) {
              i20 = 32767;
            }

            lambda->data[iC_data[i - 1] - 1] = Opt->data[i20 - 1];
            i20 = n + i;
            if (i20 > 32767) {
              i20 = 32767;
            }

            if (Opt->data[i20 - 1] < Xnorm) {
              i20 = nA - meq;
              if (i20 > 32767) {
                i20 = 32767;
              }

              if (i <= i20) {
                kDrop = i;
                i20 = n + i;
                if (i20 > 32767) {
                  i20 = 32767;
                }

                Xnorm = Opt->data[i20 - 1];
              }
            }
          }

          if (kDrop <= 0) {
            DualFeasible = true;
            if (1 > n) {
              ar = 0;
            } else {
              ar = n;
            }

            i20 = x->size[0];
            x->size[0] = ar;
            emxEnsureCapacity((emxArray__common *)x, i20, (int)sizeof(double));
            for (i20 = 0; i20 < ar; i20++) {
              x->data[i20] = Opt->data[i20];
            }
          } else {
            (*status)++;
            if ((int)*status > q) {
              ResetToColdStart(m, meq, &nA, iA_data, iA_size, iC_data, iC_size);
              i20 = b_iA->size[0];
              b_iA->size[0] = iA_size[0];
              emxEnsureCapacity((emxArray__common *)b_iA, i20, (int)sizeof(short));
              ar = iA_size[0];
              for (i20 = 0; i20 < ar; i20++) {
                b_iA->data[i20] = iA_data[i20];
              }

              ColdReset = true;
            } else {
              lambda->data[iC_data[kDrop - 1] - 1] = 0.0;
              DropConstraint(kDrop, b_iA, &nA, iC_data);
            }
          }
        }
      } else {
        if (nA <= 0) {
          i20 = lambda->size[0];
          lambda->size[0] = m;
          emxEnsureCapacity((emxArray__common *)lambda, i20, (int)sizeof(double));
          ar = m;
          for (i20 = 0; i20 < ar; i20++) {
            lambda->data[i20] = 0.0;
          }

          Unconstrained(Hinv, f, x, n);
        }

        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (exitg3 == 1) {
      guard1 = true;
    }
  } else {
    Unconstrained(Hinv, f, x, n);
    guard1 = true;
  }

  if (guard1) {
    Xnorm0 = norm(x);
    do {
      exitg2 = 0;
      if (*status <= 100.0) {
        cMin = -1.0E-12;
        kNext = 0;
        i20 = m - meq;
        if (i20 > 32767) {
          i20 = 32767;
        } else {
          if (i20 < -32768) {
            i20 = -32768;
          }
        }

        k = (short)i20;
        for (i = 1; i <= k; i++) {
          if (!cTolComputed) {
            ar = Ac->size[1];
            i20 = b_Ac->size[0] * b_Ac->size[1];
            b_Ac->size[0] = 1;
            b_Ac->size[1] = ar;
            emxEnsureCapacity((emxArray__common *)b_Ac, i20, (int)sizeof(double));
            for (i20 = 0; i20 < ar; i20++) {
              b_Ac->data[b_Ac->size[0] * i20] = Ac->data[(i + Ac->size[0] * i20)
                - 1] * x->data[i20];
            }

            b_abs(b_Ac, varargin_2);
            ar = 1;
            b_n = varargin_2->size[1];
            Xnorm = varargin_2->data[0];
            if (rtIsNaN(varargin_2->data[0])) {
              ia = 2;
              exitg5 = false;
              while ((!exitg5) && (ia <= b_n)) {
                ar = ia;
                if (!rtIsNaN(varargin_2->data[ia - 1])) {
                  Xnorm = varargin_2->data[ia - 1];
                  exitg5 = true;
                } else {
                  ia++;
                }
              }
            }

            if (ar < varargin_2->size[1]) {
              while (ar + 1 <= b_n) {
                if (varargin_2->data[ar] > Xnorm) {
                  Xnorm = varargin_2->data[ar];
                }

                ar++;
              }
            }

            if ((cTol->data[i - 1] >= Xnorm) || rtIsNaN(Xnorm)) {
              Xnorm = cTol->data[i - 1];
            }

            cTol->data[i - 1] = Xnorm;
          }

          if (b_iA->data[i - 1] == 0) {
            ar = Ac->size[1];
            i20 = a->size[0] * a->size[1];
            a->size[0] = 1;
            a->size[1] = ar;
            emxEnsureCapacity((emxArray__common *)a, i20, (int)sizeof(double));
            for (i20 = 0; i20 < ar; i20++) {
              a->data[a->size[0] * i20] = Ac->data[(i + Ac->size[0] * i20) - 1];
            }

            i20 = Ac->size[1];
            if ((i20 == 1) || (x->size[0] == 1)) {
              y = 0.0;
              for (i20 = 0; i20 < a->size[1]; i20++) {
                y += a->data[a->size[0] * i20] * x->data[i20];
              }
            } else {
              y = 0.0;
              for (i20 = 0; i20 < a->size[1]; i20++) {
                y += a->data[a->size[0] * i20] * x->data[i20];
              }
            }

            Xnorm = (y - b->data[i - 1]) / cTol->data[i - 1];
            if (Xnorm < cMin) {
              cMin = Xnorm;
              kNext = i;
            }
          }
        }

        cTolComputed = true;
        if (kNext <= 0) {
          exitg2 = 1;
        } else {
          do {
            exitg1 = 0;
            if ((kNext > 0) && (*status <= 100.0)) {
              ar = Ac->size[1];
              i20 = AcRow->size[0] * AcRow->size[1];
              AcRow->size[0] = 1;
              AcRow->size[1] = ar;
              emxEnsureCapacity((emxArray__common *)AcRow, i20, (int)sizeof
                                (double));
              for (i20 = 0; i20 < ar; i20++) {
                AcRow->data[AcRow->size[0] * i20] = Ac->data[(kNext + Ac->size[0]
                  * i20) - 1];
              }

              guard2 = false;
              guard3 = false;
              if (nA == 0) {
                i20 = b_b->size[0];
                b_b->size[0] = AcRow->size[1];
                emxEnsureCapacity((emxArray__common *)b_b, i20, (int)sizeof
                                  (double));
                ar = AcRow->size[1];
                for (i20 = 0; i20 < ar; i20++) {
                  b_b->data[i20] = AcRow->data[AcRow->size[0] * i20];
                }

                if ((Hinv->size[1] == 1) || (b_b->size[0] == 1)) {
                  i20 = z->size[0];
                  z->size[0] = Hinv->size[0];
                  emxEnsureCapacity((emxArray__common *)z, i20, (int)sizeof
                                    (double));
                  ar = Hinv->size[0];
                  for (i20 = 0; i20 < ar; i20++) {
                    z->data[i20] = 0.0;
                    b_n = Hinv->size[1];
                    for (ia = 0; ia < b_n; ia++) {
                      z->data[i20] += Hinv->data[i20 + Hinv->size[0] * ia] *
                        b_b->data[ia];
                    }
                  }
                } else {
                  b_k = Hinv->size[1];
                  Rhs_idx_0 = (unsigned int)Hinv->size[0];
                  i20 = z->size[0];
                  z->size[0] = (int)Rhs_idx_0;
                  emxEnsureCapacity((emxArray__common *)z, i20, (int)sizeof
                                    (double));
                  b_m = Hinv->size[0];
                  ar = z->size[0];
                  i20 = z->size[0];
                  z->size[0] = ar;
                  emxEnsureCapacity((emxArray__common *)z, i20, (int)sizeof
                                    (double));
                  for (i20 = 0; i20 < ar; i20++) {
                    z->data[i20] = 0.0;
                  }

                  if (Hinv->size[0] == 0) {
                  } else {
                    ar = 0;
                    while ((b_m > 0) && (ar <= 0)) {
                      for (ic = 1; ic <= b_m; ic++) {
                        z->data[ic - 1] = 0.0;
                      }

                      ar = b_m;
                    }

                    br = 0;
                    ar = 0;
                    while ((b_m > 0) && (ar <= 0)) {
                      ar = 0;
                      i20 = br + b_k;
                      for (b_n = br; b_n + 1 <= i20; b_n++) {
                        if (b_b->data[b_n] != 0.0) {
                          ia = ar;
                          for (ic = 0; ic + 1 <= b_m; ic++) {
                            ia++;
                            z->data[ic] += b_b->data[b_n] * Hinv->data[ia - 1];
                          }
                        }

                        ar += b_m;
                      }

                      br += b_k;
                      ar = b_m;
                    }
                  }
                }

                guard3 = true;
              } else {
                Xnorm = KWIKfactor(Ac, iC_data, nA, Linv, RLinv, D, H, n);
                if (Xnorm <= 0.0) {
                  *status = -2.0;
                  exitg1 = 1;
                } else {
                  i20 = U->size[0] * U->size[1];
                  U->size[0] = H->size[0];
                  U->size[1] = H->size[1];
                  emxEnsureCapacity((emxArray__common *)U, i20, (int)sizeof
                                    (double));
                  ar = H->size[0] * H->size[1];
                  for (i20 = 0; i20 < ar; i20++) {
                    U->data[i20] = -H->data[i20];
                  }

                  i20 = b_b->size[0];
                  b_b->size[0] = AcRow->size[1];
                  emxEnsureCapacity((emxArray__common *)b_b, i20, (int)sizeof
                                    (double));
                  ar = AcRow->size[1];
                  for (i20 = 0; i20 < ar; i20++) {
                    b_b->data[i20] = AcRow->data[AcRow->size[0] * i20];
                  }

                  if ((U->size[1] == 1) || (b_b->size[0] == 1)) {
                    i20 = z->size[0];
                    z->size[0] = U->size[0];
                    emxEnsureCapacity((emxArray__common *)z, i20, (int)sizeof
                                      (double));
                    ar = U->size[0];
                    for (i20 = 0; i20 < ar; i20++) {
                      z->data[i20] = 0.0;
                      b_n = U->size[1];
                      for (ia = 0; ia < b_n; ia++) {
                        z->data[i20] += U->data[i20 + U->size[0] * ia] *
                          b_b->data[ia];
                      }
                    }
                  } else {
                    b_k = U->size[1];
                    Rhs_idx_0 = (unsigned int)U->size[0];
                    i20 = z->size[0];
                    z->size[0] = (int)Rhs_idx_0;
                    emxEnsureCapacity((emxArray__common *)z, i20, (int)sizeof
                                      (double));
                    b_m = U->size[0];
                    ar = z->size[0];
                    i20 = z->size[0];
                    z->size[0] = ar;
                    emxEnsureCapacity((emxArray__common *)z, i20, (int)sizeof
                                      (double));
                    for (i20 = 0; i20 < ar; i20++) {
                      z->data[i20] = 0.0;
                    }

                    if (U->size[0] == 0) {
                    } else {
                      ar = 0;
                      while ((b_m > 0) && (ar <= 0)) {
                        for (ic = 1; ic <= b_m; ic++) {
                          z->data[ic - 1] = 0.0;
                        }

                        ar = b_m;
                      }

                      br = 0;
                      ar = 0;
                      while ((b_m > 0) && (ar <= 0)) {
                        ar = 0;
                        i20 = br + b_k;
                        for (b_n = br; b_n + 1 <= i20; b_n++) {
                          if (b_b->data[b_n] != 0.0) {
                            ia = ar;
                            for (ic = 0; ic + 1 <= b_m; ic++) {
                              ia++;
                              z->data[ic] += b_b->data[b_n] * U->data[ia - 1];
                            }
                          }

                          ar += b_m;
                        }

                        br += b_k;
                        ar = b_m;
                      }
                    }
                  }

                  for (i = 1; i <= nA; i++) {
                    ar = D->size[0];
                    i20 = b_b->size[0];
                    b_b->size[0] = ar;
                    emxEnsureCapacity((emxArray__common *)b_b, i20, (int)sizeof
                                      (double));
                    for (i20 = 0; i20 < ar; i20++) {
                      b_b->data[i20] = D->data[i20 + D->size[0] * (i - 1)];
                    }

                    i20 = Ac->size[1];
                    b_guard1 = false;
                    if (i20 == 1) {
                      b_guard1 = true;
                    } else {
                      i20 = D->size[0];
                      if (i20 == 1) {
                        b_guard1 = true;
                      } else {
                        y = 0.0;
                        for (i20 = 0; i20 < AcRow->size[1]; i20++) {
                          y += AcRow->data[AcRow->size[0] * i20] * b_b->data[i20];
                        }
                      }
                    }

                    if (b_guard1) {
                      y = 0.0;
                      for (i20 = 0; i20 < AcRow->size[1]; i20++) {
                        y += AcRow->data[AcRow->size[0] * i20] * b_b->data[i20];
                      }
                    }

                    r->data[i - 1] = y;
                  }

                  guard3 = true;
                }
              }

              if (guard3) {
                kDrop = 0;
                t1 = 0.0;
                ColdReset = true;
                DualFeasible = true;
                if (nA > meq) {
                  i20 = nA - meq;
                  if (i20 > 32767) {
                    i20 = 32767;
                  } else {
                    if (i20 < -32768) {
                      i20 = -32768;
                    }
                  }

                  k = (short)i20;
                  iSave = 1;
                  exitg4 = false;
                  while ((!exitg4) && (iSave <= k)) {
                    if (r->data[iSave - 1] >= 1.0E-12) {
                      DualFeasible = false;
                      exitg4 = true;
                    } else {
                      iSave++;
                    }
                  }
                }

                if ((nA == meq) || DualFeasible) {
                  DualFeasible = true;
                } else {
                  DualFeasible = false;
                }

                if (!DualFeasible) {
                  i20 = nA - meq;
                  if (i20 > 32767) {
                    i20 = 32767;
                  } else {
                    if (i20 < -32768) {
                      i20 = -32768;
                    }
                  }

                  k = (short)i20;
                  for (i = 1; i <= k; i++) {
                    if (r->data[i - 1] > 1.0E-12) {
                      Xnorm = lambda->data[iC_data[i - 1] - 1] / r->data[i - 1];
                      if ((kDrop == 0) || (Xnorm < rMin)) {
                        rMin = Xnorm;
                        kDrop = i;
                      }
                    }
                  }

                  if (kDrop > 0) {
                    t1 = rMin;
                    ColdReset = false;
                  }
                }

                i20 = a->size[0] * a->size[1];
                a->size[0] = 1;
                a->size[1] = z->size[0];
                emxEnsureCapacity((emxArray__common *)a, i20, (int)sizeof(double));
                ar = z->size[0];
                for (i20 = 0; i20 < ar; i20++) {
                  a->data[a->size[0] * i20] = z->data[i20];
                }

                i20 = b_b->size[0];
                b_b->size[0] = AcRow->size[1];
                emxEnsureCapacity((emxArray__common *)b_b, i20, (int)sizeof
                                  (double));
                ar = AcRow->size[1];
                for (i20 = 0; i20 < ar; i20++) {
                  b_b->data[i20] = AcRow->data[AcRow->size[0] * i20];
                }

                if ((a->size[1] == 1) || (b_b->size[0] == 1)) {
                  Xnorm = 0.0;
                  for (i20 = 0; i20 < a->size[1]; i20++) {
                    Xnorm += a->data[a->size[0] * i20] * b_b->data[i20];
                  }
                } else {
                  Xnorm = 0.0;
                  for (i20 = 0; i20 < a->size[1]; i20++) {
                    Xnorm += a->data[a->size[0] * i20] * b_b->data[i20];
                  }
                }

                if (Xnorm <= 0.0) {
                  cMin = 0.0;
                  DualFeasible = true;
                } else {
                  i20 = Ac->size[1];
                  if ((i20 == 1) || (x->size[0] == 1)) {
                    y = 0.0;
                    for (i20 = 0; i20 < AcRow->size[1]; i20++) {
                      y += AcRow->data[AcRow->size[0] * i20] * x->data[i20];
                    }
                  } else {
                    y = 0.0;
                    for (i20 = 0; i20 < AcRow->size[1]; i20++) {
                      y += AcRow->data[AcRow->size[0] * i20] * x->data[i20];
                    }
                  }

                  cMin = (b->data[kNext - 1] - y) / Xnorm;
                  DualFeasible = false;
                }

                if (ColdReset && DualFeasible) {
                  *status = -1.0;
                  exitg1 = 1;
                } else {
                  if ((t1 <= cMin) || rtIsNaN(cMin)) {
                    y = t1;
                  } else {
                    y = cMin;
                  }

                  if (DualFeasible) {
                    Xnorm = t1;
                  } else if (ColdReset) {
                    Xnorm = cMin;
                  } else {
                    Xnorm = y;
                  }

                  for (i = 1; i <= nA; i++) {
                    lambda->data[iC_data[i - 1] - 1] -= Xnorm * r->data[i - 1];
                    i20 = m - meq;
                    if (i20 > 32767) {
                      i20 = 32767;
                    } else {
                      if (i20 < -32768) {
                        i20 = -32768;
                      }
                    }

                    if ((iC_data[i - 1] <= i20) && (lambda->data[iC_data[i - 1]
                         - 1] < 0.0)) {
                      lambda->data[iC_data[i - 1] - 1] = 0.0;
                    }
                  }

                  lambda->data[kNext - 1] += Xnorm;
                  if (Xnorm == t1) {
                    DropConstraint(kDrop, b_iA, &nA, iC_data);
                  }

                  if (!DualFeasible) {
                    i20 = x->size[0];
                    emxEnsureCapacity((emxArray__common *)x, i20, (int)sizeof
                                      (double));
                    ar = x->size[0];
                    for (i20 = 0; i20 < ar; i20++) {
                      x->data[i20] += Xnorm * z->data[i20];
                    }

                    if (Xnorm == cMin) {
                      if (nA == n) {
                        *status = -1.0;
                        exitg1 = 1;
                      } else {
                        i20 = nA + 1;
                        if (i20 > 32767) {
                          i20 = 32767;
                        }

                        nA = (short)i20;
                        iC_data[nA - 1] = kNext;
                        i = nA;
                        while ((i > 1) && (!(iC_data[i - 1] > iC_data[i - 2])))
                        {
                          iSave = iC_data[i - 1];
                          iC_data[i - 1] = iC_data[i - 2];
                          iC_data[i - 2] = iSave;
                          i--;
                        }

                        b_iA->data[kNext - 1] = 1;
                        kNext = 0;
                        guard2 = true;
                      }
                    } else {
                      guard2 = true;
                    }
                  } else {
                    guard2 = true;
                  }
                }
              }

              if (guard2) {
                (*status)++;
              }
            } else {
              Xnorm = norm(x);
              if (fabs(Xnorm - Xnorm0) > 0.001) {
                Xnorm0 = Xnorm;
                c_abs(b, Rhs);
                Rhs_idx_0 = (unsigned int)Rhs->size[0];
                b_Rhs_idx_0 = (unsigned int)Rhs->size[0];
                i20 = cTol->size[0];
                cTol->size[0] = (int)b_Rhs_idx_0;
                emxEnsureCapacity((emxArray__common *)cTol, i20, (int)sizeof
                                  (double));
                for (b_k = 0; b_k + 1 <= (int)Rhs_idx_0; b_k++) {
                  if (Rhs->data[b_k] >= 1.0) {
                    Xnorm = Rhs->data[b_k];
                  } else {
                    Xnorm = 1.0;
                  }

                  cTol->data[b_k] = Xnorm;
                }

                cTolComputed = false;
              }

              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = 1;
          }
        }
      } else {
        *status = 0.0;
        exitg2 = 1;
      }
    } while (exitg2 == 0);
  }

  emxFree_real_T(&b_Ac);
  emxFree_real_T(&varargin_2);
  emxFree_real_T(&a);
  emxFree_real_T(&b_b);
  emxFree_real_T(&z);
  emxFree_real_T(&AcRow);
  emxFree_real_T(&Rhs);
  emxFree_real_T(&Opt);
  emxFree_real_T(&cTol);
  emxFree_real_T(&U);
  emxFree_real_T(&H);
  emxFree_real_T(&D);
  emxFree_real_T(&RLinv);
  emxFree_real_T(&r);
}

/*
 * File trailer for qpkwik.c
 *
 * [EOF]
 */
