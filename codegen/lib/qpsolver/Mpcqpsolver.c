/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Mpcqpsolver.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
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
 *                const emxArray_real_T *b_f
 *                const emxArray_real_T *A
 *                const emxArray_real_T *b
 *                const emxArray_real_T *b_Aeq
 *                const emxArray_real_T *b_beq
 *                const emxArray_boolean_T *iA0
 *                emxArray_real_T *solution
 *                double *status
 *                emxArray_boolean_T *iA
 * Return Type  : void
 */
void Mpcqpsolver(const emxArray_real_T *Linv, const emxArray_real_T *b_f, const
                 emxArray_real_T *A, const emxArray_real_T *b, const
                 emxArray_real_T *b_Aeq, const emxArray_real_T *b_beq, const
                 emxArray_boolean_T *iA0, emxArray_real_T *solution, double
                 *status, emxArray_boolean_T *iA)
{
  emxArray_real_T *a;
  int i2;
  int ar;
  emxArray_real_T *Hinv;
  int br;
  int i3;
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
  emxInit_real_T1(&a, 2);
  i2 = a->size[0] * a->size[1];
  a->size[0] = Linv->size[1];
  a->size[1] = Linv->size[0];
  emxEnsureCapacity((emxArray__common *)a, i2, (int)sizeof(double));
  ar = Linv->size[0];
  for (i2 = 0; i2 < ar; i2++) {
    br = Linv->size[1];
    for (i3 = 0; i3 < br; i3++) {
      a->data[i3 + a->size[0] * i2] = Linv->data[i2 + Linv->size[0] * i3];
    }
  }

  emxInit_real_T1(&Hinv, 2);
  if ((a->size[1] == 1) || (Linv->size[0] == 1)) {
    i2 = Hinv->size[0] * Hinv->size[1];
    Hinv->size[0] = a->size[0];
    Hinv->size[1] = Linv->size[1];
    emxEnsureCapacity((emxArray__common *)Hinv, i2, (int)sizeof(double));
    ar = a->size[0];
    for (i2 = 0; i2 < ar; i2++) {
      br = Linv->size[1];
      for (i3 = 0; i3 < br; i3++) {
        Hinv->data[i2 + Hinv->size[0] * i3] = 0.0;
        c = a->size[1];
        for (cr = 0; cr < c; cr++) {
          Hinv->data[i2 + Hinv->size[0] * i3] += a->data[i2 + a->size[0] * cr] *
            Linv->data[cr + Linv->size[0] * i3];
        }
      }
    }
  } else {
    k = a->size[1];
    unnamed_idx_0 = (unsigned int)a->size[0];
    unnamed_idx_1 = (unsigned int)Linv->size[1];
    i2 = Hinv->size[0] * Hinv->size[1];
    Hinv->size[0] = (int)unnamed_idx_0;
    Hinv->size[1] = (int)unnamed_idx_1;
    emxEnsureCapacity((emxArray__common *)Hinv, i2, (int)sizeof(double));
    m = a->size[0];
    i2 = Hinv->size[0] * Hinv->size[1];
    emxEnsureCapacity((emxArray__common *)Hinv, i2, (int)sizeof(double));
    ar = Hinv->size[1];
    for (i2 = 0; i2 < ar; i2++) {
      br = Hinv->size[0];
      for (i3 = 0; i3 < br; i3++) {
        Hinv->data[i3 + Hinv->size[0] * i2] = 0.0;
      }
    }

    if ((a->size[0] == 0) || (Linv->size[1] == 0)) {
    } else {
      c = a->size[0] * (Linv->size[1] - 1);
      cr = 0;
      while ((m > 0) && (cr <= c)) {
        i2 = cr + m;
        for (ic = cr; ic + 1 <= i2; ic++) {
          Hinv->data[ic] = 0.0;
        }

        cr += m;
      }

      br = 0;
      cr = 0;
      while ((m > 0) && (cr <= c)) {
        ar = -1;
        i2 = br + k;
        for (ib = br; ib + 1 <= i2; ib++) {
          if (Linv->data[ib] != 0.0) {
            ia = ar;
            i3 = cr + m;
            for (ic = cr; ic + 1 <= i3; ic++) {
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
  } else if (!((b_Aeq->size[0] == 0) || (b_Aeq->size[1] == 0))) {
    c = b_Aeq->size[1];
  } else {
    c = A->size[1];
    if (b_Aeq->size[1] > A->size[1]) {
      c = b_Aeq->size[1];
    }
  }

  empty_non_axis_sizes = (c == 0);
  if (empty_non_axis_sizes || (!((A->size[0] == 0) || (A->size[1] == 0)))) {
    cr = A->size[0];
  } else {
    cr = 0;
  }

  if (empty_non_axis_sizes || (!((b_Aeq->size[0] == 0) || (b_Aeq->size[1] == 0))))
  {
    br = b_Aeq->size[0];
  } else {
    br = 0;
  }

  emxInit_int16_T(&iA1, 1);
  i2 = iA1->size[0];
  iA1->size[0] = iA0->size[0] + b_Aeq->size[0];
  emxEnsureCapacity((emxArray__common *)iA1, i2, (int)sizeof(short));
  ar = iA0->size[0];
  for (i2 = 0; i2 < ar; i2++) {
    iA1->data[i2] = iA0->data[i2];
  }

  ar = b_Aeq->size[0];
  for (i2 = 0; i2 < ar; i2++) {
    iA1->data[i2 + iA0->size[0]] = 1;
  }

  emxInit_real_T1(&b_A, 2);
  i2 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = cr + br;
  b_A->size[1] = c;
  emxEnsureCapacity((emxArray__common *)b_A, i2, (int)sizeof(double));
  for (i2 = 0; i2 < c; i2++) {
    for (i3 = 0; i3 < cr; i3++) {
      b_A->data[i3 + b_A->size[0] * i2] = A->data[i3 + cr * i2];
    }
  }

  for (i2 = 0; i2 < c; i2++) {
    for (i3 = 0; i3 < br; i3++) {
      b_A->data[(i3 + cr) + b_A->size[0] * i2] = b_Aeq->data[i3 + br * i2];
    }
  }

  emxInit_real_T(&b_b, 1);
  i2 = b_b->size[0];
  b_b->size[0] = b->size[0] + b_beq->size[0];
  emxEnsureCapacity((emxArray__common *)b_b, i2, (int)sizeof(double));
  ar = b->size[0];
  for (i2 = 0; i2 < ar; i2++) {
    b_b->data[i2] = b->data[i2];
  }

  ar = b_beq->size[0];
  for (i2 = 0; i2 < ar; i2++) {
    b_b->data[i2 + b->size[0]] = b_beq->data[i2];
  }

  emxInit_real_T(&lam, 1);
  unnamed_idx_0 = (unsigned int)A->size[0] + b_Aeq->size[0];
  if (unnamed_idx_0 > 32767U) {
    unnamed_idx_0 = 32767U;
  }

  i2 = Linv->size[0];
  if (i2 > 32767) {
    i2 = 32767;
  } else {
    if (i2 < -32768) {
      i2 = -32768;
    }
  }

  i3 = b_Aeq->size[0];
  if (i3 > 32767) {
    i3 = 32767;
  } else {
    if (i3 < -32768) {
      i3 = -32768;
    }
  }

  qpkwik(Linv, Hinv, b_f, b_A, b_b, iA1, (short)unnamed_idx_0, (short)i2, (short)
         i3, solution, lam, &b_status);
  emxFree_real_T(&b_b);
  emxFree_real_T(&b_A);
  emxFree_real_T(&lam);
  emxFree_real_T(&Hinv);
  if (A->size[0] == 0) {
    i2 = iA->size[0];
    iA->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)iA, i2, (int)sizeof(boolean_T));
  } else {
    ar = A->size[0];
    i2 = iA->size[0];
    iA->size[0] = ar;
    emxEnsureCapacity((emxArray__common *)iA, i2, (int)sizeof(boolean_T));
    for (i2 = 0; i2 < ar; i2++) {
      iA->data[i2] = (iA1->data[i2] != 0);
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
