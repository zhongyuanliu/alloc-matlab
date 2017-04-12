/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Mpcqpsolver.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "Mpcqpsolver.h"
#include "qpsolver_emxutil.h"
#include "qpkwik.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *Linv
 *                const emxArray_real_T *f
 *                const emxArray_real_T *A
 *                const emxArray_real_T *b
 *                const emxArray_real_T *Aeq
 *                const emxArray_real_T *beq
 *                const emxArray_boolean_T *iA0
 *                emxArray_real_T *solution
 *                double *status
 *                emxArray_boolean_T *iA
 * Return Type  : void
 */
void Mpcqpsolver(const emxArray_real_T *Linv, const emxArray_real_T *f, const
                 emxArray_real_T *A, const emxArray_real_T *b, const
                 emxArray_real_T *Aeq, const emxArray_real_T *beq, const
                 emxArray_boolean_T *iA0, emxArray_real_T *solution, double
                 *status, emxArray_boolean_T *iA)
{
  emxArray_real_T *a;
  int i5;
  int ar;
  emxArray_real_T *Hinv;
  int br;
  int i6;
  int k;
  unsigned int unnamed_idx_0;
  unsigned int unnamed_idx_1;
  int m;
  int c;
  boolean_T empty_non_axis_sizes;
  int cr;
  int ic;
  emxArray_int16_T *iA1;
  int ib;
  int ia;
  emxArray_real_T *b_A;
  emxArray_real_T *b_b;
  emxArray_real_T *lam;
  double b_status;
  emxInit_real_T(&a, 2);
  i5 = a->size[0] * a->size[1];
  a->size[0] = Linv->size[1];
  a->size[1] = Linv->size[0];
  emxEnsureCapacity((emxArray__common *)a, i5, (int)sizeof(double));
  ar = Linv->size[0];
  for (i5 = 0; i5 < ar; i5++) {
    br = Linv->size[1];
    for (i6 = 0; i6 < br; i6++) {
      a->data[i6 + a->size[0] * i5] = Linv->data[i5 + Linv->size[0] * i6];
    }
  }

  emxInit_real_T(&Hinv, 2);
  if ((a->size[1] == 1) || (Linv->size[0] == 1)) {
    i5 = Hinv->size[0] * Hinv->size[1];
    Hinv->size[0] = a->size[0];
    Hinv->size[1] = Linv->size[1];
    emxEnsureCapacity((emxArray__common *)Hinv, i5, (int)sizeof(double));
    ar = a->size[0];
    for (i5 = 0; i5 < ar; i5++) {
      br = Linv->size[1];
      for (i6 = 0; i6 < br; i6++) {
        Hinv->data[i5 + Hinv->size[0] * i6] = 0.0;
        c = a->size[1];
        for (cr = 0; cr < c; cr++) {
          Hinv->data[i5 + Hinv->size[0] * i6] += a->data[i5 + a->size[0] * cr] *
            Linv->data[cr + Linv->size[0] * i6];
        }
      }
    }
  } else {
    k = a->size[1];
    unnamed_idx_0 = (unsigned int)a->size[0];
    unnamed_idx_1 = (unsigned int)Linv->size[1];
    i5 = Hinv->size[0] * Hinv->size[1];
    Hinv->size[0] = (int)unnamed_idx_0;
    Hinv->size[1] = (int)unnamed_idx_1;
    emxEnsureCapacity((emxArray__common *)Hinv, i5, (int)sizeof(double));
    m = a->size[0];
    i5 = Hinv->size[0] * Hinv->size[1];
    emxEnsureCapacity((emxArray__common *)Hinv, i5, (int)sizeof(double));
    ar = Hinv->size[1];
    for (i5 = 0; i5 < ar; i5++) {
      br = Hinv->size[0];
      for (i6 = 0; i6 < br; i6++) {
        Hinv->data[i6 + Hinv->size[0] * i5] = 0.0;
      }
    }

    if ((a->size[0] == 0) || (Linv->size[1] == 0)) {
    } else {
      c = a->size[0] * (Linv->size[1] - 1);
      cr = 0;
      while ((m > 0) && (cr <= c)) {
        i5 = cr + m;
        for (ic = cr; ic + 1 <= i5; ic++) {
          Hinv->data[ic] = 0.0;
        }

        cr += m;
      }

      br = 0;
      cr = 0;
      while ((m > 0) && (cr <= c)) {
        ar = -1;
        i5 = br + k;
        for (ib = br; ib + 1 <= i5; ib++) {
          if (Linv->data[ib] != 0.0) {
            ia = ar;
            i6 = cr + m;
            for (ic = cr; ic + 1 <= i6; ic++) {
              ia++;
              Hinv->data[ic] += Linv->data[ib] * a->data[ia];
            }
          }

          ar += m;
        }

        br += k;
        cr += m;
      }
    }
  }

  emxFree_real_T(&a);
  if (!((A->size[0] == 0) || (A->size[1] == 0))) {
    c = A->size[1];
  } else if (!((Aeq->size[0] == 0) || (Aeq->size[1] == 0))) {
    c = Aeq->size[1];
  } else {
    c = A->size[1];
    if (Aeq->size[1] > A->size[1]) {
      c = Aeq->size[1];
    }
  }

  empty_non_axis_sizes = (c == 0);
  if (empty_non_axis_sizes || (!((A->size[0] == 0) || (A->size[1] == 0)))) {
    cr = A->size[0];
  } else {
    cr = 0;
  }

  if (empty_non_axis_sizes || (!((Aeq->size[0] == 0) || (Aeq->size[1] == 0)))) {
    br = Aeq->size[0];
  } else {
    br = 0;
  }

  emxInit_int16_T(&iA1, 1);
  i5 = iA1->size[0];
  iA1->size[0] = iA0->size[0] + Aeq->size[0];
  emxEnsureCapacity((emxArray__common *)iA1, i5, (int)sizeof(short));
  ar = iA0->size[0];
  for (i5 = 0; i5 < ar; i5++) {
    iA1->data[i5] = iA0->data[i5];
  }

  ar = Aeq->size[0];
  for (i5 = 0; i5 < ar; i5++) {
    iA1->data[i5 + iA0->size[0]] = 1;
  }

  emxInit_real_T(&b_A, 2);
  i5 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = cr + br;
  b_A->size[1] = c;
  emxEnsureCapacity((emxArray__common *)b_A, i5, (int)sizeof(double));
  for (i5 = 0; i5 < c; i5++) {
    for (i6 = 0; i6 < cr; i6++) {
      b_A->data[i6 + b_A->size[0] * i5] = A->data[i6 + cr * i5];
    }
  }

  for (i5 = 0; i5 < c; i5++) {
    for (i6 = 0; i6 < br; i6++) {
      b_A->data[(i6 + cr) + b_A->size[0] * i5] = Aeq->data[i6 + br * i5];
    }
  }

  emxInit_real_T1(&b_b, 1);
  i5 = b_b->size[0];
  b_b->size[0] = b->size[0] + beq->size[0];
  emxEnsureCapacity((emxArray__common *)b_b, i5, (int)sizeof(double));
  ar = b->size[0];
  for (i5 = 0; i5 < ar; i5++) {
    b_b->data[i5] = b->data[i5];
  }

  ar = beq->size[0];
  for (i5 = 0; i5 < ar; i5++) {
    b_b->data[i5 + b->size[0]] = beq->data[i5];
  }

  emxInit_real_T1(&lam, 1);
  unnamed_idx_0 = (unsigned int)A->size[0] + Aeq->size[0];
  if (unnamed_idx_0 > 32767U) {
    unnamed_idx_0 = 32767U;
  }

  i5 = Linv->size[0];
  if (i5 > 32767) {
    i5 = 32767;
  } else {
    if (i5 < -32768) {
      i5 = -32768;
    }
  }

  i6 = Aeq->size[0];
  if (i6 > 32767) {
    i6 = 32767;
  } else {
    if (i6 < -32768) {
      i6 = -32768;
    }
  }

  qpkwik(Linv, Hinv, f, b_A, b_b, iA1, (short)unnamed_idx_0, (short)i5, (short)
         i6, solution, lam, &b_status);
  emxFree_real_T(&b_b);
  emxFree_real_T(&b_A);
  emxFree_real_T(&lam);
  emxFree_real_T(&Hinv);
  if (A->size[0] == 0) {
    i5 = iA->size[0];
    iA->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)iA, i5, (int)sizeof(boolean_T));
  } else {
    ar = A->size[0];
    i5 = iA->size[0];
    iA->size[0] = ar;
    emxEnsureCapacity((emxArray__common *)iA, i5, (int)sizeof(boolean_T));
    for (i5 = 0; i5 < ar; i5++) {
      iA->data[i5] = (iA1->data[i5] != 0);
    }
  }

  emxFree_int16_T(&iA1);
  *status = b_status;
}

/*
 * File trailer for Mpcqpsolver.c
 *
 * [EOF]
 */
