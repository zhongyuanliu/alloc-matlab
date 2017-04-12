/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: angle_dist.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "angle_dist.h"
#include "qpsolver_rtwutil.h"

/* Function Definitions */

/*
 * angle distance
 *  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise
 * Arguments    : double angle_start
 *                double angle_end
 * Return Type  : double
 */
double angle_dist(double angle_start, double angle_end)
{
  double k;
  double phi;
  double b_k;
  double b_phi;

  /*  change angle to the range -pi to pi */
  k = angle_end / 6.2831853071795862;
  if (fabs(k - rt_roundd_snf(k)) <= 2.2204460492503131E-16 * fabs(k)) {
    k = 0.0;
  } else {
    k = (k - floor(k)) * 6.2831853071795862;
  }

  /*  change angle to the range -pi to pi */
  phi = angle_start / 6.2831853071795862;
  if (fabs(phi - rt_roundd_snf(phi)) <= 2.2204460492503131E-16 * fabs(phi)) {
    phi = 0.0;
  } else {
    phi = (phi - floor(phi)) * 6.2831853071795862;
  }

  if (k > 3.1415926535897931) {
    b_k = -6.2831853071795862 + k;
  } else {
    b_k = k;
  }

  if (phi > 3.1415926535897931) {
    b_phi = -6.2831853071795862 + phi;
  } else {
    b_phi = phi;
  }

  k = b_k - b_phi;
  k *= (double)(fabs(k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  return k * (double)(k >= 0.0) + (6.2831853071795862 + k) * (double)(k < 0.0);
}

/*
 * File trailer for angle_dist.c
 *
 * [EOF]
 */
