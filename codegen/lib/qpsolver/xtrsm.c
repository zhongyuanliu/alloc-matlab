/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xtrsm.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "xtrsm.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                int n
 *                const emxArray_real_T *A
 *                int lda
 *                emxArray_real_T *B
 *                int ldb
 * Return Type  : void
 */
void xtrsm(int m, int n, const emxArray_real_T *A, int lda, emxArray_real_T *B,
           int ldb)
{
  int j;
  int jBcol;
  int k;
  int kAcol;
  double x;
  double y;
  int i;
  if ((n == 0) || ((B->size[0] == 0) || (B->size[1] == 0))) {
  } else {
    for (j = 1; j <= n; j++) {
      jBcol = ldb * (j - 1) - 1;
      for (k = m; k > 0; k--) {
        kAcol = lda * (k - 1) - 1;
        if (B->data[k + jBcol] != 0.0) {
          x = B->data[k + jBcol];
          y = A->data[k + kAcol];
          B->data[k + jBcol] = x / y;
          for (i = 1; i < k; i++) {
            B->data[i + jBcol] -= B->data[k + jBcol] * A->data[i + kAcol];
          }
        }
      }
    }
  }
}

/*
 * File trailer for xtrsm.c
 *
 * [EOF]
 */
