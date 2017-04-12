/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzlarf.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "call_qpsolver.h"
#include "xzlarf.h"

/* Function Definitions */

/*
 * Arguments    : int m
 *                int n
 *                int iv0
 *                double tau
 *                emxArray_real_T *C
 *                int ic0
 *                int ldc
 *                emxArray_real_T *work
 * Return Type  : void
 */
void xzlarf(int m, int n, int iv0, double tau, emxArray_real_T *C, int ic0, int
            ldc, emxArray_real_T *work)
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  int jy;
  int i24;
  int j;
  int ia;
  int exitg1;
  int ix;
  double c;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C->data[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      i = ic0 + (lastc - 1) * ldc;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C->data[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
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
    if (lastc == 0) {
    } else {
      for (i = 1; i <= lastc; i++) {
        work->data[i - 1] = 0.0;
      }

      i = 0;
      i24 = ic0 + ldc * (lastc - 1);
      jy = ic0;
      while ((ldc > 0) && (jy <= i24)) {
        ix = iv0;
        c = 0.0;
        j = (jy + lastv) - 1;
        for (ia = jy; ia <= j; ia++) {
          c += C->data[ia - 1] * C->data[ix - 1];
          ix++;
        }

        work->data[i] += c;
        i++;
        jy += ldc;
      }
    }

    if (-tau == 0.0) {
    } else {
      i = ic0 - 1;
      jy = 0;
      for (j = 1; j <= lastc; j++) {
        if (work->data[jy] != 0.0) {
          c = work->data[jy] * -tau;
          ix = iv0;
          i24 = lastv + i;
          for (ia = i; ia + 1 <= i24; ia++) {
            C->data[ia] += C->data[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += ldc;
      }
    }
  }
}

/*
 * File trailer for xzlarf.c
 *
 * [EOF]
 */
