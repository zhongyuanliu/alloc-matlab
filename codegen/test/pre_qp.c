/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pre_qp.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 12:18:18
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "pre_qp.h"
#include "polyval.h"
#include "sign.h"
#include "qpsolver_emxutil.h"
#include "eye.h"

/* Function Declarations */
static void N_sided_linear(double A[46]);
static void N_sided_rudder(const double points[18], double A[36], double b[18]);
static void angleMaxMin(double angle_start, double angle_end, double phi, double
  phi1, double phi2, double *phi_, double *phiPlus);
static double rt_roundd_snf(double u);

/* Function Definitions */

/*
 * %
 * Arguments    : double A[46]
 * Return Type  : void
 */
static void N_sided_linear(double A[46])
{
  int i;
  double phi;

  /*      r = - R * cos(pi / N);%change <= to >= */
  for (i = 0; i < 23; i++) {
    phi = (2.0 * (double)i + 1.0) * 3.1415926535897931 / 23.0;
    A[i] = -cos(phi);

    /* change <= to >= */
    A[23 + i] = -sin(phi);

    /* change <= to >= */
  }
}

/*
 * x,y
 * %
 *  points =[0.980,0;0.965,0.0800;0.956,0.163;0.920,0.200;0.895,0.287;0.840,0.300
 *      0.811,0.341;0.770,0.345;0.743,0.347;0.720,0.300];
 * Arguments    : const double points[18]
 *                double A[36]
 *                double b[18]
 * Return Type  : void
 */
static void N_sided_rudder(const double points[18], double A[36], double b[18])
{
  int j;
  double temp[16];
  int i;
  double xtmp;
  double b_points[38];
  for (j = 0; j < 2; j++) {
    memcpy(&temp[j << 3], &points[j * 9 + 1], sizeof(double) << 3);
    for (i = 0; i < 4; i++) {
      xtmp = temp[i + (j << 3)];
      temp[i + (j << 3)] = temp[((j << 3) - i) + 7];
      temp[((j << 3) - i) + 7] = xtmp;
    }
  }

  for (j = 0; j < 8; j++) {
    temp[8 + j] = -temp[8 + j];
  }

  for (j = 0; j < 2; j++) {
    b_points[19 * j] = 0.0;
    memcpy(&b_points[j * 19 + 1], &temp[j << 3], sizeof(double) << 3);
    memcpy(&b_points[j * 19 + 9], &points[j * 9], 9U * sizeof(double));
    b_points[18 + 19 * j] = 0.0;
  }

  /*      points(end + 1,:) = points(1,:); */
  for (i = 0; i < 18; i++) {
    A[i] = -b_points[i + 20] + b_points[19 + i];
    A[18 + i] = -b_points[i] + b_points[1 + i];
    b[i] = -b_points[i] * b_points[i + 20] + b_points[1 + i] * b_points[19 + i];
  }

  /*  check_rudder_constr(points,A,b); */
}

/*
 * % up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise
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
static void angleMaxMin(double angle_start, double angle_end, double phi, double
  phi1, double phi2, double *phi_, double *phiPlus)
{
  double b_phi;
  double k;
  double b_k;
  double c_k;
  double phi1_end;
  double phi1_phi2;
  double d_k;
  double e_k;
  double f_k;

  /* !!!!!!!!!!!!!check if there are bugs */
  /* % change angle to the range -pi to pi */
  b_phi = angle_start / 6.2831853071795862;
  if (fabs(b_phi - rt_roundd_snf(b_phi)) <= 2.2204460492503131E-16 * fabs(b_phi))
  {
    b_phi = 0.0;
  } else {
    b_phi = (b_phi - floor(b_phi)) * 6.2831853071795862;
  }

  if (b_phi > 3.1415926535897931) {
    angle_start = -6.2831853071795862 + b_phi;
  } else {
    angle_start = b_phi;
  }

  /* % change angle to the range -pi to pi */
  b_phi = angle_end / 6.2831853071795862;
  if (fabs(b_phi - rt_roundd_snf(b_phi)) <= 2.2204460492503131E-16 * fabs(b_phi))
  {
    b_phi = 0.0;
  } else {
    b_phi = (b_phi - floor(b_phi)) * 6.2831853071795862;
  }

  if (b_phi > 3.1415926535897931) {
    angle_end = -6.2831853071795862 + b_phi;
  } else {
    angle_end = b_phi;
  }

  /* % change angle to the range -pi to pi */
  b_phi = phi / 6.2831853071795862;
  if (fabs(b_phi - rt_roundd_snf(b_phi)) <= 2.2204460492503131E-16 * fabs(b_phi))
  {
    b_phi = 0.0;
  } else {
    b_phi = (b_phi - floor(b_phi)) * 6.2831853071795862;
  }

  if (b_phi > 3.1415926535897931) {
    phi = -6.2831853071795862 + b_phi;
  } else {
    phi = b_phi;
  }

  /* % change angle to the range -pi to pi */
  b_phi = phi1 / 6.2831853071795862;
  if (fabs(b_phi - rt_roundd_snf(b_phi)) <= 2.2204460492503131E-16 * fabs(b_phi))
  {
    b_phi = 0.0;
  } else {
    b_phi = (b_phi - floor(b_phi)) * 6.2831853071795862;
  }

  if (b_phi > 3.1415926535897931) {
    phi1 = -6.2831853071795862 + b_phi;
  } else {
    phi1 = b_phi;
  }

  /* % change angle to the range -pi to pi */
  b_phi = phi2 / 6.2831853071795862;
  if (fabs(b_phi - rt_roundd_snf(b_phi)) <= 2.2204460492503131E-16 * fabs(b_phi))
  {
    b_phi = 0.0;
  } else {
    b_phi = (b_phi - floor(b_phi)) * 6.2831853071795862;
  }

  if (b_phi > 3.1415926535897931) {
    phi2 = -6.2831853071795862 + b_phi;
  } else {
    phi2 = b_phi;
  }

  /* % angle distance */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  k = angle_end - angle_start;
  k *= (double)(fabs(k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  /* % angle distance */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  b_k = phi1 - angle_start;
  b_k *= (double)(fabs(b_k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  /* % angle distance */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  c_k = angle_end - phi1;
  c_k *= (double)(fabs(c_k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  phi1_end = c_k * (double)(c_k >= 0.0) + (6.2831853071795862 + c_k) * (double)
    (c_k < 0.0);

  /* % angle distance */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  c_k = phi2 - phi1;
  c_k *= (double)(fabs(c_k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  phi1_phi2 = c_k * (double)(c_k >= 0.0) + (6.2831853071795862 + c_k) * (double)
    (c_k < 0.0);

  /* % angle distance */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  c_k = angle_start - phi1;
  c_k *= (double)(fabs(c_k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  /* % angle distance */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  d_k = phi - phi1;
  d_k *= (double)(fabs(d_k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  /* % angle distance */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  e_k = phi - angle_start;
  e_k *= (double)(fabs(e_k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  /* % angle distance */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  f_k = phi2 - angle_start;
  f_k *= (double)(fabs(f_k) >= 1.0E-8);

  /* set small value to zero, avoid numerical problem */
  *phi_ = 0.0;
  *phiPlus = 0.0;
  if (k * (double)(k >= 0.0) + (6.2831853071795862 + k) * (double)(k < 0.0) >=
      b_k * (double)(b_k >= 0.0) + (6.2831853071795862 + b_k) * (double)(b_k <
       0.0)) {
    /* phi1 is at the range from start to end */
    if (phi1_end >= phi1_phi2) {
      /* phi2 is at the range from phi1 to end */
      *phi_ = phi1;
      *phiPlus = phi2;
    } else if ((phi1_phi2 > phi1_end) && (phi1_phi2 <= c_k * (double)(c_k >= 0.0)
                + (6.2831853071795862 + c_k) * (double)(c_k < 0.0))) {
      *phi_ = phi1;
      *phiPlus = angle_end;
    } else {
      /* phi2 is at the range from start to phi1 */
      if (d_k * (double)(d_k >= 0.0) + (6.2831853071795862 + d_k) * (double)(d_k
           < 0.0) <= phi1_end) {
        /* phi is at the range from phi1 to end */
        *phi_ = phi1;
        *phiPlus = angle_end;
      } else {
        if (e_k * (double)(e_k >= 0.0) + (6.2831853071795862 + e_k) * (double)
            (e_k < 0.0) <= f_k * (double)(f_k >= 0.0) + (6.2831853071795862 +
             f_k) * (double)(f_k < 0.0)) {
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
 * Arguments    : double u
 * Return Type  : double
 */
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/*
 * % ---------------precalculation------------------
 * Arguments    : struct0_T thruster_data[8]
 *                double T_r_Tx
 *                double T_r_Ty
 *                double T_r_Tm
 *                double N_enabled_thruster
 *                const double rudder_table[45]
 *                const emxArray_real_T *b_Rangle_T
 *                const emxArray_real_T *b_Rangle_Fangle
 *                double no_azi_angle_constr
 *                double *Nvar
 *                emxArray_real_T *H
 *                emxArray_real_T *f
 *                emxArray_real_T *A
 *                emxArray_real_T *b
 *                emxArray_real_T *Aeq
 *                emxArray_real_T *beq
 *                emxArray_real_T *x0
 * Return Type  : void
 */
void pre_qp(struct0_T thruster_data[8], double T_r_Tx, double T_r_Ty, double
            T_r_Tm, double N_enabled_thruster, const double rudder_table[45],
            const emxArray_real_T *b_Rangle_T, const emxArray_real_T
            *b_Rangle_Fangle, double no_azi_angle_constr, double *Nvar,
            emxArray_real_T *H, emxArray_real_T *f, emxArray_real_T *A,
            emxArray_real_T *b, emxArray_real_T *Aeq, emxArray_real_T *beq,
            emxArray_real_T *x0)
{
  double Ntunnel;
  double Nazi;
  double Nfpp;
  double dt;
  double A_azi[46];
  int i17;
  int loop_ub;
  int j;
  double x;
  double N_A;
  double b_x;
  double c_x;
  double temp1;
  double temp2;
  unsigned int i_H;
  double numOfTunnel;
  double numOfAzi;
  double numOfFpp;
  double pos_A;
  unsigned int pos_Aeq;
  double b_thruster_data[18];
  int i;
  int ixy;
  double A_fpp[36];
  double b_fpp[18];
  int i_Azi;
  int i_Fpp;
  Ntunnel = 0.0;
  Nazi = 0.0;
  Nfpp = 0.0;
  dt = thruster_data[0].dt;

  /* !!!!!!!!!!!!!!!!!!!!!!!!!! */
  N_sided_linear(A_azi);

  /*  NA_azi = 0;%not used at the moment!!!!!!!!!!!!!!!! */
  i17 = x0->size[0];
  x0->size[0] = (int)(N_enabled_thruster * 2.0 + 3.0);
  emxEnsureCapacity((emxArray__common *)x0, i17, (int)sizeof(double));
  loop_ub = (int)(N_enabled_thruster * 2.0 + 3.0);
  for (i17 = 0; i17 < loop_ub; i17++) {
    x0->data[i17] = 0.0;
  }

  for (j = 0; j < (int)N_enabled_thruster; j++) {
    x0->data[(int)(2.0 * (1.0 + (double)j) - 1.0) - 1] = thruster_data[j].T *
      cos(thruster_data[j].phi);
    x0->data[(int)(2.0 * (1.0 + (double)j)) - 1] = thruster_data[j].T * sin
      (thruster_data[j].phi);
    switch ((int)thruster_data[j].type) {
     case 4:
      /* tunnel */
      x = thruster_data[j].T + dt * thruster_data[j].dTmax;
      if ((x <= thruster_data[j].Tmax) || rtIsNaN(thruster_data[j].Tmax)) {
        thruster_data[j].Tplus = x;
      } else {
        thruster_data[j].Tplus = thruster_data[j].Tmax;
      }

      x = thruster_data[j].T - dt * thruster_data[j].dTmax;
      if ((x >= thruster_data[j].Tmin) || rtIsNaN(thruster_data[j].Tmin)) {
        thruster_data[j].T_ = x;
      } else {
        thruster_data[j].T_ = thruster_data[j].Tmin;
      }

      Ntunnel++;
      break;

     case 3:
      /* azi */
      x = thruster_data[j].T + dt * thruster_data[j].dTmax;
      if ((x <= thruster_data[j].Tmax) || rtIsNaN(thruster_data[j].Tmax)) {
        thruster_data[j].Tplus = x;
      } else {
        thruster_data[j].Tplus = thruster_data[j].Tmax;
      }

      x = thruster_data[j].T - dt * thruster_data[j].dTmax;
      if ((x >= thruster_data[j].Tmin) || rtIsNaN(thruster_data[j].Tmin)) {
        thruster_data[j].T_ = x;
      } else {
        thruster_data[j].T_ = thruster_data[j].Tmin;
      }

      if ((thruster_data[j].constr_angle == 1.0) && (fabs(thruster_data[j].
            phi_max - thruster_data[j].phi_min) > 1.0E-8)) {
        /* phi_min == phi_max means no angle constraint */
        angleMaxMin(thruster_data[j].phi_min, thruster_data[j].phi_max,
                    thruster_data[j].phi, thruster_data[j].phi - dt *
                    thruster_data[j].dphi_max, thruster_data[j].phi + dt *
                    thruster_data[j].dphi_max, &thruster_data[j].phi_,
                    &thruster_data[j].phiPlus);
      } else {
        thruster_data[j].phi_ = thruster_data[j].phi - dt * thruster_data[j].
          dphi_max;
        thruster_data[j].phiPlus = thruster_data[j].phi + dt * thruster_data[j].
          dphi_max;

        /*                  thruster_data(j).constr_angle = 1; */
      }

      /*              thruster_data(j).phiPlus = min(thruster_data(j).phi + ... */
      /*                  dt * thruster_data(j).dphi_max, thruster_data(j).phi_max); */
      /*  */
      /*              thruster_data(j).phi_ = max(thruster_data(j).phi - ... */
      /*                  dt * thruster_data(j).dphi_max, thruster_data(j).phi_min); */
      Nazi++;
      break;

     case 2:
      /* fpp */
      x = thruster_data[j].T + dt * thruster_data[j].dTmax;
      if ((x <= thruster_data[j].Tmax) || rtIsNaN(thruster_data[j].Tmax)) {
        b_x = x;
      } else {
        b_x = thruster_data[j].Tmax;
      }

      thruster_data[j].Tplus = b_x * polyval(b_Rangle_T, fabs(thruster_data[j].
        phi));
      x = thruster_data[j].T - dt * thruster_data[j].dTmax;
      if ((x >= thruster_data[j].Tmin) || rtIsNaN(thruster_data[j].Tmin)) {
        c_x = x;
      } else {
        c_x = thruster_data[j].Tmin;
      }

      thruster_data[j].T_ = c_x * polyval(b_Rangle_T, fabs(thruster_data[j].phi));
      angleMaxMin(thruster_data[j].phi_min, thruster_data[j].phi_max,
                  thruster_data[j].phi, thruster_data[j].phi - dt *
                  thruster_data[j].dphi_max, thruster_data[j].phi + dt *
                  thruster_data[j].dphi_max, &temp1, &temp2);
      x = temp2;
      b_sign(&x);
      thruster_data[j].phiPlus = x * polyval(b_Rangle_Fangle, fabs(temp2));
      x = temp1;
      b_sign(&x);
      thruster_data[j].phi_ = x * polyval(b_Rangle_Fangle, fabs(temp1));

      /*              temp = (min(thruster_data(j).phi + ... %force angle */
      /*                  dt * thruster_data(j).dphi_max, thruster_data(j).phi_max)); */
      /*              thruster_data(j).phiPlus = sign(temp) * polyval(Rangle_Fangle,abs(temp)); */
      /*              temp = (max(thruster_data(j).phi - ... %force angle */
      /*                  dt * thruster_data(j).dphi_max, thruster_data(j).phi_min)); */
      /*              thruster_data(j).phi_ = sign(temp) * polyval(Rangle_Fangle,abs(temp)); */
      Nfpp++;
      break;
    }
  }

  /* !!!!!to be improved */
  /*  chck_sol = [[thruster_data(:).phi_min]' [thruster_data(:).phi_]' [thruster_data.phi]' [thruster_data.phiPlus]' [thruster_data(:).phi_max]']; */
  /*  check_angle(chck_sol(:,1),chck_sol(:,2),chck_sol(:,3),chck_sol(:,4),chck_sol(:,5)) */
  *Nvar = ((2.0 * Ntunnel + 2.0 * Nazi) + 2.0 * Nfpp) + 3.0;

  /* number of variables in quad */
  if (no_azi_angle_constr == 0.0) {
    /*  there is angle constraint */
    N_A = (Ntunnel * 2.0 + Nazi * 28.0) + Nfpp * 23.0;
  } else {
    N_A = (Ntunnel * 2.0 + Nazi * 23.0) + Nfpp * 23.0;
  }

  /* % ------------------------------------------------------ */
  eye(*Nvar, H);
  i17 = f->size[0];
  f->size[0] = (int)*Nvar;
  emxEnsureCapacity((emxArray__common *)f, i17, (int)sizeof(double));
  loop_ub = (int)*Nvar;
  for (i17 = 0; i17 < loop_ub; i17++) {
    f->data[i17] = 0.0;
  }

  i17 = Aeq->size[0] * Aeq->size[1];
  Aeq->size[0] = (int)(3.0 + Ntunnel);
  Aeq->size[1] = (int)*Nvar;
  emxEnsureCapacity((emxArray__common *)Aeq, i17, (int)sizeof(double));
  loop_ub = (int)(3.0 + Ntunnel) * (int)*Nvar;
  for (i17 = 0; i17 < loop_ub; i17++) {
    Aeq->data[i17] = 0.0;
  }

  i17 = beq->size[0];
  beq->size[0] = (int)(3.0 + Ntunnel);
  emxEnsureCapacity((emxArray__common *)beq, i17, (int)sizeof(double));
  loop_ub = (int)(3.0 + Ntunnel);
  for (i17 = 0; i17 < loop_ub; i17++) {
    beq->data[i17] = 0.0;
  }

  beq->data[0] = T_r_Tx;
  beq->data[1] = T_r_Ty;
  beq->data[2] = T_r_Tm;
  i17 = A->size[0] * A->size[1];
  A->size[0] = (int)N_A;
  A->size[1] = (int)*Nvar;
  emxEnsureCapacity((emxArray__common *)A, i17, (int)sizeof(double));
  loop_ub = (int)N_A * (int)*Nvar;
  for (i17 = 0; i17 < loop_ub; i17++) {
    A->data[i17] = 0.0;
  }

  i17 = b->size[0];
  b->size[0] = (int)N_A;
  emxEnsureCapacity((emxArray__common *)b, i17, (int)sizeof(double));
  loop_ub = (int)N_A;
  for (i17 = 0; i17 < loop_ub; i17++) {
    b->data[i17] = 0.0;
  }

  /*  ------------------------------------------------------ */
  i_H = 2U;
  numOfTunnel = 0.0;
  numOfAzi = 0.0;
  numOfFpp = 0.0;
  pos_A = 1.0;
  pos_Aeq = 4U;
  for (j = 0; j < (int)N_enabled_thruster; j++) {
    switch ((int)thruster_data[j].type) {
     case 4:
      /* tunnel */
      H->data[((int)i_H + H->size[0] * ((int)i_H - 2)) - 2] = thruster_data[j].
        weight[0];
      H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = thruster_data[j].
        weight[1];
      i_H += 2U;
      for (i = 0; i < 3; i++) {
        /* R1 R2 R3 */
        switch (i + 1) {
         case 1:
          Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
            2.0) + numOfFpp * 2.0) + 1.0) - 1)] = 1.0;
          Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
            2.0) + numOfFpp * 2.0) + 2.0) - 1)] = 0.0;
          break;

         case 2:
          Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
            2.0) + numOfFpp * 2.0) + 1.0) - 1)] = 0.0;
          Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
            2.0) + numOfFpp * 2.0) + 2.0) - 1)] = 1.0;
          break;

         case 3:
          Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
            2.0) + numOfFpp * 2.0) + 1.0) - 1)] = -thruster_data[j].y +
            thruster_data[j].y0;
          Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
            2.0) + numOfFpp * 2.0) + 2.0) - 1)] = thruster_data[j].x -
            thruster_data[j].x0;
          break;
        }
      }

      Aeq->data[((int)pos_Aeq + Aeq->size[0] * ((int)(((numOfAzi * 2.0 +
        numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) - 1] = sin
        (thruster_data[j].phi);
      Aeq->data[((int)pos_Aeq + Aeq->size[0] * ((int)(((numOfAzi * 2.0 +
        numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) - 1] = -cos
        (thruster_data[j].phi);
      beq->data[(int)pos_Aeq - 1] = 0.0;
      A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
        2.0) + numOfFpp * 2.0) + 1.0) - 1)) - 1] = -cos(thruster_data[j].phi);
      A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
        2.0) + numOfFpp * 2.0) + 2.0) - 1)) - 1] = -sin(thruster_data[j].phi);
      b->data[(int)pos_A - 1] = -thruster_data[j].Tplus;
      A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
        numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = cos(thruster_data[j]
        .phi);
      A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
        numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = sin(thruster_data[j]
        .phi);
      b->data[(int)(unsigned int)pos_A] = thruster_data[j].T_;

      /*              check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 1, ... */
      /*                  numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*                  b(pos_A:pos_A + 1)); */
      pos_A += 2.0;
      pos_Aeq++;
      numOfTunnel++;
      break;

     case 3:
      /* azi */
      H->data[((int)i_H + H->size[0] * ((int)i_H - 2)) - 2] = thruster_data[j].
        weight[0];
      H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = thruster_data[j].
        weight[1];
      i_H += 2U;
      for (ixy = 0; ixy < 2; ixy++) {
        /*  xthurst and ythrust */
        for (i = 0; i < 3; i++) {
          /* R1 R2 R3 */
          if (1 + ixy == 1) {
            switch (i + 1) {
             case 1:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = 1.0;
              break;

             case 2:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = 0.0;
              break;

             case 3:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = -thruster_data[j].y +
                thruster_data[j].y0;
              break;
            }
          } else {
            switch (i + 1) {
             case 1:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = 0.0;
              break;

             case 2:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = 1.0;
              break;

             case 3:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = thruster_data[j].x -
                thruster_data[j].x0;
              break;
            }
          }
        }

        if (1 + ixy == 1) {
          if (no_azi_angle_constr == 0.0) {
            /*  there is angle constraint */
            A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
              numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) - 1] = -sin
              (thruster_data[j].phi_);

            /* (3.15a) */
            b->data[(int)pos_A - 1] = -0.0 * fabs(cos(thruster_data[j].phi_)) *
              (thruster_data[j].Tplus - thruster_data[j].T_);
            A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
              2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = sin
              (thruster_data[j].phiPlus);

            /* (3.15b) */
            b->data[(int)(unsigned int)pos_A] = -0.0 * fabs(cos(thruster_data[j]
              .phiPlus)) * (thruster_data[j].Tplus - thruster_data[j].T_);
            A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
              2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 1] = cos
              (thruster_data[j].phi);

            /* (3.16a) */
            b->data[(int)(unsigned int)pos_A + 1] = thruster_data[j].T_;
            A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
              2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 2] =
              -cos(thruster_data[j].phi_);

            /* (3.16b) */
            b->data[(int)(unsigned int)pos_A + 2] = -cos(thruster_data[j].phi_ -
              thruster_data[j].phi) * thruster_data[j].Tplus;
            A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
              2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 3] =
              -cos(thruster_data[j].phiPlus);

            /* (3.16c) */
            b->data[(int)(unsigned int)pos_A + 3] = -cos(thruster_data[j].
              phiPlus - thruster_data[j].phi) * thruster_data[j].Tplus;
            for (i_Azi = 0; i_Azi < 23; i_Azi++) {
              /* linearized thrust zone */
              A->data[((int)((unsigned int)pos_A + (1 + i_Azi)) + A->size[0] *
                       ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) + numOfFpp *
                               2.0) + 1.0) - 1)) + 3] = A_azi[i_Azi];
              b->data[(int)((unsigned int)pos_A + (1 + i_Azi)) + 3] =
                -thruster_data[j].Tmax * 0.99068594603633076;
            }
          } else {
            for (i_Azi = 0; i_Azi < 23; i_Azi++) {
              /* linearized thrust zone */
              A->data[((int)((pos_A - 1.0) + (1.0 + (double)i_Azi)) + A->size[0]
                       * ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) + numOfFpp
                                 * 2.0) + 1.0) - 1)) - 1] = A_azi[i_Azi];
              b->data[(int)((pos_A - 1.0) + (1.0 + (double)i_Azi)) - 1] =
                -thruster_data[j].Tmax * 0.99068594603633076;
            }
          }
        } else if (no_azi_angle_constr == 0.0) {
          /*  there is angle constraint */
          A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
            numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) - 1] = cos
            (thruster_data[j].phi_);

          /* (3.15a) */
          A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = -cos
            (thruster_data[j].phiPlus);

          /* (3.15b) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 1] = sin
            (thruster_data[j].phi);

          /* (3.16a) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 2] = -sin
            (thruster_data[j].phi_);

          /* (3.16b) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 3] = -sin
            (thruster_data[j].phiPlus);

          /* (3.16c) */
          for (i_Azi = 0; i_Azi < 23; i_Azi++) {
            /* linearized thrust zone */
            A->data[((int)((unsigned int)pos_A + (1 + i_Azi)) + A->size[0] *
                     ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) + numOfFpp *
                             2.0) + 2.0) - 1)) + 3] = A_azi[23 + i_Azi];
          }
        } else {
          for (i_Azi = 0; i_Azi < 23; i_Azi++) {
            /* linearized thrust zone */
            A->data[((int)((pos_A - 1.0) + (1.0 + (double)i_Azi)) + A->size[0] *
                     ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) + numOfFpp *
                             2.0) + 2.0) - 1)) - 1] = A_azi[23 + i_Azi];
          }
        }
      }

      /*              check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 4, ... */
      /*                  numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*                  b(pos_A:pos_A + 4));%azi */
      /*  */
      /*                      check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A + 4 + 1:pos_A + 4 + i_Azi, ... */
      /*                          numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*                          b(pos_A + 4 + 1:pos_A + 4 + i_Azi)); */
      /*          check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 4 + i_Azi, ... */
      /*              numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*              b(pos_A:pos_A + 4 + i_Azi)); */
      if (no_azi_angle_constr == 0.0) {
        /*  there is angle constraint */
        pos_A = (pos_A + 5.0) + 23.0;
      } else {
        pos_A += 23.0;
      }

      numOfAzi++;
      break;

     case 2:
      /* fpp */
      for (i17 = 0; i17 < 2; i17++) {
        for (loop_ub = 0; loop_ub < 9; loop_ub++) {
          b_thruster_data[loop_ub + 9 * i17] = thruster_data[j].Tmax *
            rudder_table[loop_ub + 9 * (1 + i17)];
        }
      }

      N_sided_rudder(b_thruster_data, A_fpp, b_fpp);

      /* x,y */
      H->data[((int)i_H + H->size[0] * ((int)i_H - 2)) - 2] = thruster_data[j].
        weight[0];
      H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = thruster_data[j].
        weight[1];
      i_H += 2U;
      for (ixy = 0; ixy < 2; ixy++) {
        /*  xthurst and ythrust */
        for (i = 0; i < 3; i++) {
          /* R1 R2 R3 */
          if (1 + ixy == 1) {
            switch (i + 1) {
             case 1:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = 1.0;
              break;

             case 2:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = 0.0;
              break;

             case 3:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = -thruster_data[j].y +
                thruster_data[j].y0;
              break;
            }
          } else {
            switch (i + 1) {
             case 1:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = 0.0;
              break;

             case 2:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = 1.0;
              break;

             case 3:
              Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel
                * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = thruster_data[j].x -
                thruster_data[j].x0;
              break;
            }
          }
        }

        if (1 + ixy == 1) {
          A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
            numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) - 1] = -sin
            (thruster_data[j].phi_);

          /* (3.15a) */
          b->data[(int)pos_A - 1] = -0.0 * fabs(cos(thruster_data[j].phi_)) *
            (thruster_data[j].Tplus - thruster_data[j].T_);
          A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = sin
            (thruster_data[j].phiPlus);

          /* (3.15b) */
          b->data[(int)(unsigned int)pos_A] = -0.0 * fabs(cos(thruster_data[j].
            phiPlus)) * (thruster_data[j].Tplus - thruster_data[j].T_);
          x = polyval(b_Rangle_Fangle, fabs(thruster_data[j].phi));
          x = cos(x);
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 1] = x;

          /* (3.16a) */
          b->data[(int)(unsigned int)pos_A + 1] = thruster_data[j].T_;
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 2] = -cos
            (thruster_data[j].phi_);

          /* (3.16b) */
          x = thruster_data[j].phi;
          b_sign(&x);
          x = thruster_data[j].phi_ - x * polyval(b_Rangle_Fangle, fabs
            (thruster_data[j].phi));
          b->data[(int)(unsigned int)pos_A + 2] = -cos(x) * thruster_data[j].
            Tplus;
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 3] = -cos
            (thruster_data[j].phiPlus);

          /* (3.16c) */
          x = thruster_data[j].phi;
          b_sign(&x);
          x = thruster_data[j].phiPlus - x * polyval(b_Rangle_Fangle, fabs
            (thruster_data[j].phi));
          b->data[(int)(unsigned int)pos_A + 3] = -cos(x) * thruster_data[j].
            Tplus;
          for (i_Fpp = 0; i_Fpp < 18; i_Fpp++) {
            /* linearized thrust zone */
            A->data[((int)((unsigned int)pos_A + (1 + i_Fpp)) + A->size[0] *
                     ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) + numOfFpp *
                             2.0) + 1.0) - 1)) + 3] = A_fpp[i_Fpp];
            b->data[(int)((unsigned int)pos_A + (1 + i_Fpp)) + 3] = b_fpp[i_Fpp];
          }
        } else {
          A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
            numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) - 1] = cos
            (thruster_data[j].phi_);

          /* (3.15a) */
          A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = -cos
            (thruster_data[j].phiPlus);

          /* (3.15b) */
          x = thruster_data[j].phi;
          b_sign(&x);
          x *= polyval(b_Rangle_Fangle, fabs(thruster_data[j].phi));
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 1] = sin(x);

          /* (3.16a) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 2] = -sin
            (thruster_data[j].phi_);

          /* (3.16b) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 3] = -sin
            (thruster_data[j].phiPlus);

          /* (3.16c) */
          for (i_Fpp = 0; i_Fpp < 18; i_Fpp++) {
            A->data[((int)((unsigned int)pos_A + (1 + i_Fpp)) + A->size[0] *
                     ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) + numOfFpp *
                             2.0) + 2.0) - 1)) + 3] = A_fpp[18 + i_Fpp];
          }
        }
      }

      /*          check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 4, ... */
      /*              numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*              b(pos_A:pos_A + 4));%fpp */
      /*  */
      /*          check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A + 4 + 1:pos_A + 4 + i_Fpp, ... */
      /*              numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*              b(pos_A + 4 + 1:pos_A + 4 + i_Fpp)); */
      /*              check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 4 + i_Fpp, ... */
      /*                  numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*                  b(pos_A:pos_A + 4 + i_Fpp)); */
      pos_A = (pos_A + 5.0) + 18.0;
      numOfFpp++;
      break;
    }
  }

  /* !!-----------------add bound for s1 s2 s3---does not help------------ */
  /*  A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = 1; */
  /*  A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = 1; */
  /*  A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = 1; */
  /*  A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = -1; */
  /*  A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = -1; */
  /*  A(end + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = -1; */
  /*  b(end + 1) = -1000000; */
  /*  b(end + 1) = -1000000; */
  /*  b(end + 1) = -1000000; */
  /*  b(end + 1) = 1000000; */
  /*  b(end + 1) = 1000000; */
  /*  b(end + 1) = 1000000; */
  /* ------------------------------------------------- */
  Aeq->data[Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) +
    numOfFpp * 2.0) + 1.0) - 1)] = -1.0;

  /* s1 */
  Aeq->data[1 + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) +
    numOfFpp * 2.0) + 2.0) - 1)] = -1.0;

  /* s2 */
  Aeq->data[2 + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) +
    numOfFpp * 2.0) + 3.0) - 1)] = -1.0;

  /* s3 */
  for (i = 0; i < 3; i++) {
    H->data[((int)(i_H + (1 + i)) + H->size[0] * ((int)(i_H + (1 + i)) - 3)) - 3]
      = thruster_data[0].weight_s[i];

    /* !!!thruster_data(1) must be available */
  }
}

/*
 * File trailer for pre_qp.c
 *
 * [EOF]
 */
