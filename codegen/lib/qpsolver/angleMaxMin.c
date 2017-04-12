/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: angleMaxMin.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "angleMaxMin.h"
#include "angle_dist.h"
#include "qpsolver_rtwutil.h"

/* Function Definitions */

/*
 * up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise
 *  phi: current angle
 *  phi1 to phi2
 * !!!!!!!!!!!!!!!!!!!!!watch out for FeasibilityTol
 * Arguments    : double angle_start
 *                double angle_end
 *                double phi
 *                double phi1
 *                double phi2
 *                double *phi_
 *                double *phiPlus
 * Return Type  : void
 */
void angleMaxMin(double angle_start, double angle_end, double phi, double phi1,
                 double phi2, double *phi_, double *phiPlus)
{
  double phi1_end;
  double phi1_phi2;

  /* !!!!!!!!!!!!!check if there are bugs */
  /*  change angle to the range -pi to pi */
  phi1_end = angle_start / 6.2831853071795862;
  if (fabs(phi1_end - rt_roundd_snf(phi1_end)) <= 2.2204460492503131E-16 * fabs
      (phi1_end)) {
    phi1_end = 0.0;
  } else {
    phi1_end = (phi1_end - floor(phi1_end)) * 6.2831853071795862;
  }

  if (phi1_end > 3.1415926535897931) {
    angle_start = -6.2831853071795862 + phi1_end;
  } else {
    angle_start = phi1_end;
  }

  /*  change angle to the range -pi to pi */
  phi1_end = angle_end / 6.2831853071795862;
  if (fabs(phi1_end - rt_roundd_snf(phi1_end)) <= 2.2204460492503131E-16 * fabs
      (phi1_end)) {
    phi1_end = 0.0;
  } else {
    phi1_end = (phi1_end - floor(phi1_end)) * 6.2831853071795862;
  }

  if (phi1_end > 3.1415926535897931) {
    angle_end = -6.2831853071795862 + phi1_end;
  } else {
    angle_end = phi1_end;
  }

  /*  change angle to the range -pi to pi */
  phi1_end = phi / 6.2831853071795862;
  if (fabs(phi1_end - rt_roundd_snf(phi1_end)) <= 2.2204460492503131E-16 * fabs
      (phi1_end)) {
    phi1_end = 0.0;
  } else {
    phi1_end = (phi1_end - floor(phi1_end)) * 6.2831853071795862;
  }

  if (phi1_end > 3.1415926535897931) {
    phi = -6.2831853071795862 + phi1_end;
  } else {
    phi = phi1_end;
  }

  /*  change angle to the range -pi to pi */
  phi1_end = phi1 / 6.2831853071795862;
  if (fabs(phi1_end - rt_roundd_snf(phi1_end)) <= 2.2204460492503131E-16 * fabs
      (phi1_end)) {
    phi1_end = 0.0;
  } else {
    phi1_end = (phi1_end - floor(phi1_end)) * 6.2831853071795862;
  }

  if (phi1_end > 3.1415926535897931) {
    phi1 = -6.2831853071795862 + phi1_end;
  } else {
    phi1 = phi1_end;
  }

  /*  change angle to the range -pi to pi */
  phi1_end = phi2 / 6.2831853071795862;
  if (fabs(phi1_end - rt_roundd_snf(phi1_end)) <= 2.2204460492503131E-16 * fabs
      (phi1_end)) {
    phi1_end = 0.0;
  } else {
    phi1_end = (phi1_end - floor(phi1_end)) * 6.2831853071795862;
  }

  if (phi1_end > 3.1415926535897931) {
    phi2 = -6.2831853071795862 + phi1_end;
  } else {
    phi2 = phi1_end;
  }

  phi1_end = angle_dist(phi1, angle_end);
  phi1_phi2 = angle_dist(phi1, phi2);
  *phi_ = 0.0;
  *phiPlus = 0.0;
  if (fabs(angle_start - angle_end) < 1.0E-8) {
    /* angle_start == angle_end means no forbidden zone */
    *phi_ = phi1;
    *phiPlus = phi2;
  } else if (angle_dist(angle_start, angle_end) >= angle_dist(angle_start, phi1))
  {
    /* phi1 is at the range from start to end */
    if (phi1_end >= phi1_phi2) {
      /* phi2 is at the range from phi1 to end */
      *phi_ = phi1;
      *phiPlus = phi2;
    } else if ((phi1_phi2 > phi1_end) && (phi1_phi2 <= angle_dist(phi1,
                 angle_start))) {
      *phi_ = phi1;
      *phiPlus = angle_end;
    } else {
      /* phi2 is at the range from start to phi1 */
      if (angle_dist(phi1, phi) <= phi1_end) {
        /* phi is at the range from phi1 to end */
        *phi_ = phi1;
        *phiPlus = angle_end;
      } else {
        if (angle_dist(angle_start, phi) <= angle_dist(angle_start, phi2)) {
          /* phi is at the range from start to phi2 */
          *phi_ = angle_start;
          *phiPlus = phi2;
        }
      }
    }
  } else {
    /* phi1 is at the range from end to start */
    *phi_ = angle_start;
    *phiPlus = phi2;
  }
}

/*
 * File trailer for angleMaxMin.c
 *
 * [EOF]
 */
