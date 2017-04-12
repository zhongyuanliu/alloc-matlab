/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: azi_eff.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "azi_eff.h"
#include "angle_dist.h"
#include "qpsolver_rtwutil.h"

/* Function Definitions */

/*
 * y = 0;
 *  N = size(Ce,1);
 * Arguments    : const double Ce[74]
 *                double phi
 *                double y_data[]
 *                int y_size[1]
 * Return Type  : void
 */
void azi_eff(const double Ce[74], double phi, double y_data[], int y_size[1])
{
  signed char tmp_data[36];
  int trueCount;
  boolean_T bv0[36];
  boolean_T bv1[36];
  int i;
  int partialTrueCount;
  boolean_T b0;
  boolean_T b1;

  /*  change angle to the range -pi to pi */
  phi /= 6.2831853071795862;
  if (fabs(phi - rt_roundd_snf(phi)) <= 2.2204460492503131E-16 * fabs(phi)) {
    phi = 0.0;
  } else {
    phi = (phi - floor(phi)) * 6.2831853071795862;
  }

  phi *= (double)(fabs(phi) > 1.0E-8);

  /*  if phi > pi */
  /*      y = -2 * pi + phi; */
  /*  else */
  /*      y = phi; */
  /*  end */
  /*  if phi>=Ce(end-1,1) */
  /*      y1 = Ce(end-1,2); */
  /*      y2 = Ce(end,2); */
  /*   */
  /*      y = (phi - Ce(end-1,1)) / (Ce(end,1) - Ce(end-1,1)) * (y2 - y1) + y1; */
  /*  else */
  /*      t= ((phi >= Ce(1:end-1,1) & phi<Ce(2:end,1))); */
  /*      I = find(t==1); */
  /*      y1 = Ce(I,2); */
  /*      y2 = Ce(I + 1,2); */
  /*      y = (phi - Ce(I,1)) / (Ce(I + 1,1) - Ce(I,1)) * (y2 - y1) + y1; */
  /*  end */
  /*  % % % % % % % % % % % % % % % % % % % % % % % %  */
  if (phi >= Ce[35]) {
    y_size[0] = 1;
    y_data[0] = Ce[72];
  } else {
    trueCount = 0;
    for (i = 0; i < 36; i++) {
      b0 = (phi >= Ce[i]);
      b1 = (phi < Ce[1 + i]);
      if (b0 && b1) {
        trueCount++;
      }

      bv0[i] = b0;
      bv1[i] = b1;
    }

    partialTrueCount = 0;
    for (i = 0; i < 36; i++) {
      if (bv0[i] && bv1[i]) {
        tmp_data[partialTrueCount] = (signed char)(i + 1);
        partialTrueCount++;
      }
    }

    y_size[0] = trueCount;
    for (i = 0; i < trueCount; i++) {
      y_data[i] = Ce[tmp_data[i] + 36];
    }
  }

  /*  for i = 1 : N - 1 */
  /*      if (phi>= Ce(i,1) && phi<Ce(i + 1,1)) ... */
  /*          || phi>=Ce(end-1,1) */
  /*          y = Ce(i,2); */
  /*      end */
  /*  end */
}

/*
 * File trailer for azi_eff.c
 *
 * [EOF]
 */
