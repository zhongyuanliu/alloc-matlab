/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: pre_qp.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "pre_qp.h"
#include "qpsolver_emxutil.h"
#include "polyval.h"
#include "sign.h"
#include "angleMaxMin.h"
#include "azi_eff.h"
#include "eye.h"

/* Variable Definitions */
static double A_azi[20];
static boolean_T A_azi_not_empty;
static boolean_T NA_azi_not_empty;
static boolean_T A_fpp_not_empty;
static boolean_T b_fpp_not_empty;
static emxArray_real_T *H;
static boolean_T H_not_empty;
static emxArray_real_T *f;
static boolean_T f_not_empty;
static boolean_T N_Aeq_not_empty;
static emxArray_real_T *Aeq;
static boolean_T Aeq_not_empty;
static emxArray_real_T *beq;
static boolean_T beq_not_empty;

/* Function Declarations */
static void N_sided_linear(double A[20]);
static void N_sided_rudder(const double points[18], double A[36], double b[18]);
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double A[20]
 * Return Type  : void
 */
static void N_sided_linear(double A[20])
{
  int i;
  double phi;

  /*      r = - R * cos(pi / N);%change <= to >= */
  for (i = 0; i < 10; i++) {
    phi = (2.0 * (double)i + 1.0) * 3.1415926535897931 / 10.0;
    A[i] = -cos(phi);

    /* change <= to >= */
    A[10 + i] = -sin(phi);

    /* change <= to >= */
  }
}

/*
 * x,y
 *
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
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = rtNaN;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void A_azi_not_empty_init(void)
{
  A_azi_not_empty = false;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void A_fpp_not_empty_init(void)
{
  A_fpp_not_empty = false;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void NA_azi_not_empty_init(void)
{
  NA_azi_not_empty = false;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void N_Aeq_not_empty_init(void)
{
  N_Aeq_not_empty = false;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void b_fpp_not_empty_init(void)
{
  b_fpp_not_empty = false;
}

/*
 * 1.04 thruster_data(i).phi is replaced by thruster_data(i).phi_min for tunnel thruster
 *   1.06 guarantee initial x feasible
 *   4.06 u is scaled. c_scale * Tmax^2 == 1 |13-03-2017 by Liu
 *  ---------------precalculation------------------
 * Arguments    : boolean_T init
 *                struct0_T thruster_data[8]
 *                double T_r_Tx
 *                double T_r_Ty
 *                double T_r_Tm
 *                double N_enabled_thruster
 *                const struct_T Ce[8]
 *                const double rudder_table[45]
 *                const double b_Rangle_T[3]
 *                const double b_Rangle_Fangle[3]
 *                double *Nvar
 *                emxArray_real_T *Ho
 *                emxArray_real_T *fo
 *                emxArray_real_T *Ao
 *                emxArray_real_T *bo
 *                emxArray_real_T *Aeqo
 *                emxArray_real_T *beqo
 *                emxArray_real_T *x0
 * Return Type  : void
 */
void pre_qp(boolean_T init, struct0_T thruster_data[8], double T_r_Tx, double
            T_r_Ty, double T_r_Tm, double N_enabled_thruster, const struct_T Ce
            [8], const double rudder_table[45], const double b_Rangle_T[3],
            const double b_Rangle_Fangle[3], double *Nvar, emxArray_real_T *Ho,
            emxArray_real_T *fo, emxArray_real_T *Ao, emxArray_real_T *bo,
            emxArray_real_T *Aeqo, emxArray_real_T *beqo, emxArray_real_T *x0)
{
  double Ntunnel;
  double Nazi;
  double Nfpp;
  double dt;
  int ix;
  int ixstart;
  int n;
  double N_A;
  double b_N_A;
  double c_N_A;
  double c;
  double i_H;
  emxArray_real_T *eff_all;
  double pos_Aeq;
  double A_fpp[36];
  int tmp_size[1];
  int itmp;
  boolean_T exitg1;
  double b_thruster_data[18];
  double b_fpp[18];
  double b_eff_all[6];
  double dv3[6];
  double dv4[2];
  double dv5[10];
  double dv6[2];
  int tmp_data[10];
  double dv7[10];
  double dv8[2];
  int iv0[2];
  int b_tmp_data[10];

  /* !!!avoid numerical problem */
  Ntunnel = 0.0;
  Nazi = 0.0;
  Nfpp = 0.0;
  dt = thruster_data[0].dt;

  /* !!!!!!!!!!!!!!!!!!!!!!!!!! */
  if ((!A_azi_not_empty) || (!NA_azi_not_empty) || init) {
    N_sided_linear(A_azi);
    A_azi_not_empty = true;
    NA_azi_not_empty = true;
  }

  if ((!A_fpp_not_empty) || (!b_fpp_not_empty) || init) {
    A_fpp_not_empty = true;
    b_fpp_not_empty = true;

    /* !!what if there is no fpp!!!!!!!!!! */
  }

  /*  NA_azi = 0;%not used at the moment!!!!!!!!!!!!!!!! */
  /*  init = true; */
  ix = x0->size[0];
  x0->size[0] = (int)(N_enabled_thruster * 2.0 + 3.0);
  emxEnsureCapacity((emxArray__common *)x0, ix, (int)sizeof(double));
  ixstart = (int)(N_enabled_thruster * 2.0 + 3.0);
  for (ix = 0; ix < ixstart; ix++) {
    x0->data[ix] = 0.0;
  }

  for (n = 0; n < (int)N_enabled_thruster; n++) {
    /*      check if phi meets constraint, if check_phi == 0, phi is not between */
    /*      max and min, then it should be changed. VERY IMPORTANT */
    /*      check_phi = check_angle(thruster_data(j).phi_min,thruster_data(j).phi_min, ... */
    /*          thruster_data(j).phi,thruster_data(j).phi_max,thruster_data(j).phi_max); */
    /*       */
    /*      if check_phi == 0 % guarantee initial phi feasible */
    /*          thruster_data(j).phi = ... */
    /*              thruster_data(j).phi_min + ... */
    /*              angle_dist(thruster_data(j).phi_min, thruster_data(j).phi_max) / 2; */
    /*   */
    /*          thruster_data(j).phi = Pi_toPi(thruster_data(j).phi); */
    /*      end */
    if ((thruster_data[n].type == 4.0) && (fabs(thruster_data[n].phi_min -
          thruster_data[n].phi) > 1.0E-6)) {
      /* added, avoid tunnel infeasible angle input, v1.07 */
      thruster_data[n].phi = thruster_data[n].phi_min;
    }

    if (fabs(thruster_data[n].Tmax - thruster_data[n].Tmin) < 1.0E-6) {
      /* avoid Tmax == Tmin */
      thruster_data[n].Tmax = 1.0E+7;

      /* big enough */
      if (thruster_data[n].type == 4.0) {
        thruster_data[n].Tmin = -1.0E+7;
      } else {
        thruster_data[n].Tmin = 0.0;
      }
    }

    if ((thruster_data[n].T > thruster_data[n].Tmax) || (thruster_data[n].T <
         thruster_data[n].Tmin)) {
      /*  guarantee initial T feasible  */
      if (thruster_data[n].type == 4.0) {
        thruster_data[n].T = 0.0;
      } else {
        thruster_data[n].T = thruster_data[n].Tmin;
      }
    }

    if ((thruster_data[n].T_reserve > thruster_data[n].Tmax) || (thruster_data[n]
         .T_reserve < thruster_data[n].Tmin)) {
      /*  guarantee initial T feasible  */
      if (thruster_data[n].type == 4.0) {
        thruster_data[n].T_reserve = 0.0;
      } else {
        thruster_data[n].T_reserve = thruster_data[n].Tmin;
      }
    }

    /*  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %      */
    x0->data[(int)(2.0 * (1.0 + (double)n) - 1.0) - 1] = thruster_data[n].T *
      cos(thruster_data[n].phi);
    x0->data[(int)(2.0 * (1.0 + (double)n)) - 1] = thruster_data[n].T * sin
      (thruster_data[n].phi);
    switch ((int)thruster_data[n].type) {
     case 4:
      /* tunnel */
      N_A = thruster_data[n].T + dt * thruster_data[n].dTmax;
      if ((N_A <= thruster_data[n].Tmax) || rtIsNaN(thruster_data[n].Tmax)) {
        thruster_data[n].Tplus = N_A;
      } else {
        thruster_data[n].Tplus = thruster_data[n].Tmax;
      }

      N_A = thruster_data[n].T - dt * thruster_data[n].dTmax;
      if ((N_A >= thruster_data[n].Tmin) || rtIsNaN(thruster_data[n].Tmin)) {
        thruster_data[n].T_ = N_A;
      } else {
        thruster_data[n].T_ = thruster_data[n].Tmin;
      }

      Ntunnel++;
      break;

     case 3:
      /* azi */
      /*              thruster_data(j).Tplus = min(thruster_data(j).T + ... */
      /*                  dt * thruster_data(j).dTmax, thruster_data(j).Tmax); */
      /*               */
      /*              thruster_data(j).T_ = max(thruster_data(j).T - ... */
      /*                  dt * thruster_data(j).dTmax, thruster_data(j).Tmin); */
      /*              ----------------------------- */
      /*              if thruster_data(j).constr == 2  */
      /*                  thruster_data(j).phi_ = thruster_data(j).phi;%!!!!!!add -tol or not */
      /*                  thruster_data(j).phiPlus = thruster_data(j).phi;%!!!!!add +tol or not */
      /*              elseif thruster_data(j).constr == 3  */
      /*                   */
      /*              end */
      /*              thruster_data(j).phi_ = thruster_data(j).phi - dt * thruster_data(j).dphi_max; */
      /*                  thruster_data(j).phiPlus = thruster_data(j).phi + dt * thruster_data(j).dphi_max; */
      /*              ----------------------------- */
      /*              thruster_data(j).phiPlus = min(thruster_data(j).phi + ... */
      /*                  dt * thruster_data(j).dphi_max, thruster_data(j).phi_max); */
      /*  */
      /*              thruster_data(j).phi_ = max(thruster_data(j).phi - ... */
      /*                  dt * thruster_data(j).dphi_max, thruster_data(j).phi_min); */
      Nazi++;
      break;

     case 2:
      /* fpp */
      N_A = thruster_data[n].T + dt * thruster_data[n].dTmax;
      if ((N_A <= thruster_data[n].Tmax) || rtIsNaN(thruster_data[n].Tmax)) {
        b_N_A = N_A;
      } else {
        b_N_A = thruster_data[n].Tmax;
      }

      thruster_data[n].Tplus = b_N_A * polyval(b_Rangle_T, fabs(thruster_data[n]
        .phi));
      N_A = thruster_data[n].T - dt * thruster_data[n].dTmax;
      if ((N_A >= thruster_data[n].Tmin) || rtIsNaN(thruster_data[n].Tmin)) {
        c_N_A = N_A;
      } else {
        c_N_A = thruster_data[n].Tmin;
      }

      thruster_data[n].T_ = c_N_A * polyval(b_Rangle_T, fabs(thruster_data[n].
        phi));
      angleMaxMin(thruster_data[n].phi_min, thruster_data[n].phi_max,
                  thruster_data[n].phi, thruster_data[n].phi - dt *
                  thruster_data[n].dphi_max, thruster_data[n].phi + dt *
                  thruster_data[n].dphi_max, &N_A, &c);
      i_H = c;
      b_sign(&i_H);
      thruster_data[n].phiPlus = i_H * polyval(b_Rangle_Fangle, fabs(c));
      i_H = N_A;
      b_sign(&i_H);
      thruster_data[n].phi_ = i_H * polyval(b_Rangle_Fangle, fabs(N_A));

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
  /*      for j=1:N_enabled_thruster%!!!!!the above Tplus is overwitten */
  /*          thruster_data(j).Tplus = thruster_data(j).Tmax; */
  /*          thruster_data(j).T_ = thruster_data(j).Tmin; */
  /*          thruster_data(j).phi_ = thruster_data(j).phi_min;%not correct. if phi_max - phi_min >pi, constraint falls */
  /*          thruster_data(j).phiPlus = thruster_data(j).phi_max; */
  /*      end */
  /*  chck_sol = [[thruster_data(:).phi_min]' [thruster_data(:).phi_]' [thruster_data.phi]' [thruster_data.phiPlus]' [thruster_data(:).phi_max]']; */
  /*  check_angle(chck_sol(:,1),chck_sol(:,2),chck_sol(:,3),chck_sol(:,4),chck_sol(:,5)) */
  *Nvar = ((2.0 * Ntunnel + 2.0 * Nazi) + 2.0 * Nfpp) + 3.0;

  /* number of variables in quad */
  /*  if no_azi_angle_constr == 0 % there is angle constraint */
  /*      N_A = Ntunnel * 2 + Nazi * (5 + NA_azi) + Nfpp * (5 + NA_fpp); */
  /*  else */
  /*      N_A = Ntunnel * 2 + Nazi *  NA_azi + Nfpp * (5 + NA_fpp); */
  /*  end */
  N_A = 0.0;

  /* added v1.08 */
  for (n = 0; n < (int)N_enabled_thruster; n++) {
    /* added v1.08 */
    switch ((int)thruster_data[n].type) {
     case 4:
      N_A += 2.0;
      break;

     case 3:
      if ((thruster_data[n].constr == 2.0) || (thruster_data[n].constr == 3.0))
      {
        N_A += 15.0;
      } else {
        N_A += 10.0;
      }
      break;

     case 2:
      N_A += 23.0;
      break;
    }
  }

  /*  ------------------------------------------------------ */
  if ((!H_not_empty) || (!f_not_empty) || (!N_Aeq_not_empty) || (!Aeq_not_empty)
      || (!beq_not_empty) || init) {
    eye(*Nvar, H);
    H_not_empty = true;
    ix = f->size[0];
    f->size[0] = (int)*Nvar;
    emxEnsureCapacity((emxArray__common *)f, ix, (int)sizeof(double));
    ixstart = (int)*Nvar;
    for (ix = 0; ix < ixstart; ix++) {
      f->data[ix] = 0.0;
    }

    f_not_empty = true;
    N_Aeq_not_empty = true;
    ix = Aeq->size[0] * Aeq->size[1];
    Aeq->size[0] = (int)(3.0 + Ntunnel);
    Aeq->size[1] = (int)*Nvar;
    emxEnsureCapacity((emxArray__common *)Aeq, ix, (int)sizeof(double));
    ixstart = (int)(3.0 + Ntunnel) * (int)*Nvar;
    for (ix = 0; ix < ixstart; ix++) {
      Aeq->data[ix] = 0.0;
    }

    Aeq_not_empty = true;
    ix = beq->size[0];
    beq->size[0] = (int)(3.0 + Ntunnel);
    emxEnsureCapacity((emxArray__common *)beq, ix, (int)sizeof(double));
    ixstart = (int)(3.0 + Ntunnel);
    for (ix = 0; ix < ixstart; ix++) {
      beq->data[ix] = 0.0;
    }

    beq_not_empty = true;
  }

  ix = Ao->size[0] * Ao->size[1];
  Ao->size[0] = (int)N_A;
  Ao->size[1] = (int)*Nvar;
  emxEnsureCapacity((emxArray__common *)Ao, ix, (int)sizeof(double));
  ixstart = (int)N_A * (int)*Nvar;
  for (ix = 0; ix < ixstart; ix++) {
    Ao->data[ix] = 0.0;
  }

  ix = bo->size[0];
  bo->size[0] = (int)N_A;
  emxEnsureCapacity((emxArray__common *)bo, ix, (int)sizeof(double));
  ixstart = (int)N_A;
  for (ix = 0; ix < ixstart; ix++) {
    bo->data[ix] = 0.0;
  }

  emxInit_real_T1(&eff_all, 2);
  beq->data[0] = T_r_Tx;
  beq->data[1] = T_r_Ty;
  beq->data[2] = T_r_Tm;

  /*  ------------------------------------------------------ */
  i_H = 1.0;
  Nazi = 0.0;
  Nfpp = 0.0;
  dt = 0.0;
  Ntunnel = 1.0;
  pos_Aeq = 4.0;
  ix = eff_all->size[0] * eff_all->size[1];
  eff_all->size[0] = 1;
  eff_all->size[1] = (int)N_enabled_thruster;
  emxEnsureCapacity((emxArray__common *)eff_all, ix, (int)sizeof(double));
  ixstart = (int)N_enabled_thruster;
  for (ix = 0; ix < ixstart; ix++) {
    eff_all->data[ix] = 1.0;
  }

  for (n = 0; n < (int)N_enabled_thruster; n++) {
    if (thruster_data[n].type == 3.0) {
      azi_eff(Ce[n].value, thruster_data[n].phi, A_fpp, tmp_size);
      eff_all->data[n] = A_fpp[0];
    }
  }

  ixstart = 1;
  n = eff_all->size[1];
  N_A = eff_all->data[0];
  itmp = 0;
  if (eff_all->size[1] > 1) {
    if (rtIsNaN(eff_all->data[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix <= n)) {
        ixstart = ix;
        if (!rtIsNaN(eff_all->data[ix - 1])) {
          N_A = eff_all->data[ix - 1];
          itmp = ix - 1;
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < eff_all->size[1]) {
      while (ixstart + 1 <= n) {
        if (eff_all->data[ixstart] < N_A) {
          N_A = eff_all->data[ixstart];
          itmp = ixstart;
        }

        ixstart++;
      }
    }
  }

  ix = eff_all->size[0] * eff_all->size[1];
  eff_all->size[0] = 1;
  eff_all->size[1] = (int)N_enabled_thruster;
  emxEnsureCapacity((emxArray__common *)eff_all, ix, (int)sizeof(double));
  ixstart = (int)N_enabled_thruster;
  for (ix = 0; ix < ixstart; ix++) {
    eff_all->data[ix] = 1.0;
  }

  eff_all->data[itmp] = N_A;
  for (n = 0; n < (int)N_enabled_thruster; n++) {
    switch ((int)thruster_data[n].type) {
     case 4:
      /* tunnel */
      N_A = (Nfpp * 2.0 + Nazi * 2.0) + dt * 2.0;
      if (init) {
        H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = thruster_data[n]
          .weight[0] * rt_powd_snf(thruster_data[n].Tmax, 4.0);
        H->data[(int)(unsigned int)i_H + H->size[0] * (int)(unsigned int)i_H] =
          thruster_data[n].weight[1] * rt_powd_snf(thruster_data[n].Tmax, 4.0);
        c = thruster_data[n].Tmax * thruster_data[n].Tmax;
        for (ix = 0; ix < 2; ix++) {
          dv3[3 * ix] = 1.0 - (double)ix;
          dv3[1 + 3 * ix] = ix;
        }

        dv3[2] = -thruster_data[n].y + thruster_data[n].y0;
        dv3[5] = thruster_data[n].x - thruster_data[n].x0;
        for (ix = 0; ix < 2; ix++) {
          for (ixstart = 0; ixstart < 3; ixstart++) {
            Aeq->data[ixstart + Aeq->size[0] * ((int)(N_A + (1.0 + (double)ix))
              - 1)] = dv3[ixstart + 3 * ix] * c;
          }
        }

        /* R1 R2 R3 */
        /*              Aeq(pos_Aeq,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = sin(thruster_data(j).phi_min); */
        /*              Aeq(pos_Aeq,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2) = -cos(thruster_data(j).phi_min); */
        c = thruster_data[n].Tmax * thruster_data[n].Tmax;
        dv8[0] = sin(thruster_data[n].phi_min);
        dv8[1] = -cos(thruster_data[n].phi_min);
        for (ix = 0; ix < 2; ix++) {
          Aeq->data[((int)pos_Aeq + Aeq->size[0] * ((int)(N_A + (1.0 + (double)
            ix)) - 1)) - 1] = dv8[ix] * c;
        }

        beq->data[(int)pos_Aeq - 1] = 0.0;
      }

      /*              A(pos_A,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = -cos(thruster_data(j).phi_min); */
      /*              A(pos_A,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2) = -sin(thruster_data(j).phi_min); */
      c = thruster_data[n].Tmax * thruster_data[n].Tmax;
      dv4[0] = -cos(thruster_data[n].phi_min);
      dv4[1] = -sin(thruster_data[n].phi_min);
      for (ix = 0; ix < 2; ix++) {
        Ao->data[((int)Ntunnel + Ao->size[0] * ((int)(N_A + (1.0 + (double)ix))
                   - 1)) - 1] = dv4[ix] * c;
      }

      /*              b(pos_A) = -thruster_data(j).Tplus; */
      /*              A(pos_A + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1) = cos(thruster_data(j).phi_min); */
      /*              A(pos_A + 1,numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2) = sin(thruster_data(j).phi_min); */
      c = thruster_data[n].Tmax * thruster_data[n].Tmax;
      dv6[0] = cos(thruster_data[n].phi_min);
      dv6[1] = sin(thruster_data[n].phi_min);
      ixstart = (int)((unsigned int)Ntunnel + 1U);
      for (ix = 0; ix < 2; ix++) {
        Ao->data[(ixstart + Ao->size[0] * ((int)(N_A + (1.0 + (double)ix)) - 1))
          - 1] = dv6[ix] * c;
      }

      /*              b(pos_A + 1) = thruster_data(j).T_; */
      bo->data[(int)(unsigned int)Ntunnel - 1] = -thruster_data[n].Tplus;
      bo->data[(int)(unsigned int)Ntunnel] = thruster_data[n].T_;

      /*              check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 1, ... */
      /*                  numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*                  b(pos_A:pos_A + 1)); */
      i_H += 2.0;
      Ntunnel += 2.0;
      pos_Aeq++;
      Nazi++;
      break;

     case 3:
      /* azi */
      if (init) {
        H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = thruster_data[n]
          .weight[0] * rt_powd_snf(thruster_data[n].Tmax, 4.0);
        H->data[(int)(unsigned int)i_H + H->size[0] * (int)(unsigned int)i_H] =
          thruster_data[n].weight[1] * rt_powd_snf(thruster_data[n].Tmax, 4.0);
      }

      N_A = (Nfpp * 2.0 + Nazi * 2.0) + dt * 2.0;
      c = thruster_data[n].Tmax * thruster_data[n].Tmax;
      b_eff_all[0] = eff_all->data[n];
      b_eff_all[3] = 0.0;
      b_eff_all[1] = 0.0;
      b_eff_all[4] = eff_all->data[n];
      b_eff_all[2] = -thruster_data[n].y + thruster_data[n].y0;
      b_eff_all[5] = thruster_data[n].x - thruster_data[n].x0;
      for (ix = 0; ix < 2; ix++) {
        for (ixstart = 0; ixstart < 3; ixstart++) {
          Aeq->data[ixstart + Aeq->size[0] * ((int)(N_A + (1.0 + (double)ix)) -
            1)] = b_eff_all[ixstart + 3 * ix] * c;
        }
      }

      if ((thruster_data[n].constr == 2.0) || (thruster_data[n].constr == 3.0))
      {
        /*  there is constraint */
        c = thruster_data[n].Tmax * thruster_data[n].Tmax;
        dv7[0] = -sin(thruster_data[n].phi_);
        dv7[5] = cos(thruster_data[n].phi_);
        dv7[1] = sin(thruster_data[n].phiPlus);
        dv7[6] = -cos(thruster_data[n].phiPlus);
        dv7[2] = cos(thruster_data[n].phi);
        dv7[7] = sin(thruster_data[n].phi);
        dv7[3] = -cos(thruster_data[n].phi_);
        dv7[8] = -sin(thruster_data[n].phi_);
        dv7[4] = -cos(thruster_data[n].phiPlus);
        dv7[9] = -sin(thruster_data[n].phiPlus);
        for (ix = 0; ix < 2; ix++) {
          for (ixstart = 0; ixstart < 5; ixstart++) {
            Ao->data[((int)(Ntunnel + (double)ixstart) + Ao->size[0] * ((int)
                       (N_A + (1.0 + (double)ix)) - 1)) - 1] = dv7[ixstart + 5 *
              ix] * c;
          }
        }

        bo->data[(int)Ntunnel - 1] = -0.0 * fabs(cos(thruster_data[n].phi_)) *
          (thruster_data[n].Tplus - thruster_data[n].T_);
        bo->data[(int)(Ntunnel + 1.0) - 1] = -0.0 * fabs(cos(thruster_data[n].
          phiPlus)) * (thruster_data[n].Tplus - thruster_data[n].T_);
        bo->data[(int)(Ntunnel + 2.0) - 1] = thruster_data[n].T_;
        bo->data[(int)(Ntunnel + 3.0) - 1] = -cos(thruster_data[n].phi_ -
          thruster_data[n].phi) * thruster_data[n].Tplus;
        bo->data[(int)(Ntunnel + 4.0) - 1] = -cos(thruster_data[n].phiPlus -
          thruster_data[n].phi) * thruster_data[n].Tplus;

        /*                  if init == true */
        for (ix = 0; ix < 10; ix++) {
          tmp_data[ix] = (int)((Ntunnel + 4.0) + (1.0 + (double)ix)) - 1;
        }

        for (ix = 0; ix < 2; ix++) {
          iv0[ix] = (int)(N_A + (1.0 + (double)ix)) - 1;
        }

        c = thruster_data[n].Tmax * thruster_data[n].Tmax;
        for (ix = 0; ix < 2; ix++) {
          for (ixstart = 0; ixstart < 10; ixstart++) {
            Ao->data[tmp_data[ixstart] + Ao->size[0] * iv0[ix]] = A_azi[ixstart
              + 10 * ix] * c;
          }
        }

        for (ix = 0; ix < 10; ix++) {
          b_tmp_data[ix] = (int)((Ntunnel + 4.0) + (1.0 + (double)ix));
        }

        for (ix = 0; ix < 10; ix++) {
          bo->data[b_tmp_data[ix] - 1] = -thruster_data[n].Tmax *
            0.95105651629515353;
        }

        /*                  end */
      } else {
        /*                  if init == true */
        for (ix = 0; ix < 10; ix++) {
          tmp_data[ix] = (int)((Ntunnel - 1.0) + (1.0 + (double)ix)) - 1;
        }

        for (ix = 0; ix < 2; ix++) {
          iv0[ix] = (int)(N_A + (1.0 + (double)ix)) - 1;
        }

        c = thruster_data[n].Tmax * thruster_data[n].Tmax;
        for (ix = 0; ix < 2; ix++) {
          for (ixstart = 0; ixstart < 10; ixstart++) {
            Ao->data[tmp_data[ixstart] + Ao->size[0] * iv0[ix]] = A_azi[ixstart
              + 10 * ix] * c;
          }
        }

        for (ix = 0; ix < 10; ix++) {
          b_tmp_data[ix] = (int)((Ntunnel - 1.0) + (1.0 + (double)ix));
        }

        for (ix = 0; ix < 10; ix++) {
          bo->data[b_tmp_data[ix] - 1] = -thruster_data[n].Tmax *
            0.95105651629515353;
        }

        /*                  end */
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
      if ((thruster_data[n].constr == 2.0) || (thruster_data[n].constr == 3.0))
      {
        /*  there is constraint  */
        Ntunnel = (Ntunnel + 5.0) + 10.0;
      } else {
        Ntunnel += 10.0;
      }

      i_H += 2.0;
      Nfpp++;
      break;

     case 2:
      /* fpp */
      for (ix = 0; ix < 2; ix++) {
        for (ixstart = 0; ixstart < 9; ixstart++) {
          b_thruster_data[ixstart + 9 * ix] = thruster_data[n].Tmax *
            rudder_table[ixstart + 9 * (1 + ix)];
        }
      }

      N_sided_rudder(b_thruster_data, A_fpp, b_fpp);

      /* x,y */
      if (init) {
        H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = thruster_data[n]
          .weight[0] * rt_powd_snf(thruster_data[n].Tmax, 4.0);
        H->data[(int)(unsigned int)i_H + H->size[0] * (int)(unsigned int)i_H] =
          thruster_data[n].weight[1] * rt_powd_snf(thruster_data[n].Tmax, 4.0);
      }

      N_A = (Nfpp * 2.0 + Nazi * 2.0) + dt * 2.0;
      if (init) {
        c = thruster_data[n].Tmax * thruster_data[n].Tmax;
        for (ix = 0; ix < 2; ix++) {
          dv3[3 * ix] = 1.0 - (double)ix;
          dv3[1 + 3 * ix] = ix;
        }

        dv3[2] = -thruster_data[n].y + thruster_data[n].y0;
        dv3[5] = thruster_data[n].x - thruster_data[n].x0;
        for (ix = 0; ix < 2; ix++) {
          for (ixstart = 0; ixstart < 3; ixstart++) {
            Aeq->data[ixstart + Aeq->size[0] * ((int)(N_A + (1.0 + (double)ix))
              - 1)] = dv3[ixstart + 3 * ix] * c;
          }
        }
      }

      c = thruster_data[n].Tmax * thruster_data[n].Tmax;
      dv5[0] = -sin(thruster_data[n].phi_);
      dv5[5] = cos(thruster_data[n].phi_);
      dv5[1] = sin(thruster_data[n].phiPlus);
      dv5[6] = -cos(thruster_data[n].phiPlus);
      dv5[2] = cos(thruster_data[n].phi);
      dv5[7] = sin(thruster_data[n].phi);
      dv5[3] = -cos(thruster_data[n].phi_);
      dv5[8] = -sin(thruster_data[n].phi_);
      dv5[4] = -cos(thruster_data[n].phiPlus);
      dv5[9] = -sin(thruster_data[n].phiPlus);
      for (ix = 0; ix < 2; ix++) {
        for (ixstart = 0; ixstart < 5; ixstart++) {
          Ao->data[((int)(Ntunnel + (double)ixstart) + Ao->size[0] * ((int)(N_A
                      + (1.0 + (double)ix)) - 1)) - 1] = dv5[ixstart + 5 * ix] *
            c;
        }
      }

      bo->data[(int)Ntunnel - 1] = -0.0 * fabs(cos(thruster_data[n].phi_)) *
        (thruster_data[n].Tplus - thruster_data[n].T_);
      bo->data[(int)(Ntunnel + 1.0) - 1] = -0.0 * fabs(cos(thruster_data[n].
        phiPlus)) * (thruster_data[n].Tplus - thruster_data[n].T_);
      bo->data[(int)(Ntunnel + 2.0) - 1] = thruster_data[n].T_;
      bo->data[(int)(Ntunnel + 3.0) - 1] = -cos(thruster_data[n].phi_ -
        thruster_data[n].phi) * thruster_data[n].Tplus;
      bo->data[(int)(Ntunnel + 4.0) - 1] = -cos(thruster_data[n].phiPlus -
        thruster_data[n].phi) * thruster_data[n].Tplus;
      if (init) {
        c = thruster_data[n].Tmax * thruster_data[n].Tmax;
        for (ix = 0; ix < 2; ix++) {
          for (ixstart = 0; ixstart < 18; ixstart++) {
            Ao->data[((int)((Ntunnel + 4.0) + (1.0 + (double)ixstart)) +
                      Ao->size[0] * ((int)(N_A + (1.0 + (double)ix)) - 1)) - 1] =
              A_fpp[ixstart + 18 * ix] * c;
          }
        }

        for (ix = 0; ix < 18; ix++) {
          bo->data[(int)((Ntunnel + 4.0) + (1.0 + (double)ix)) - 1] = b_fpp[ix];
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
      i_H += 2.0;
      Ntunnel = (Ntunnel + 5.0) + 18.0;
      dt++;
      break;
    }
  }

  emxFree_real_T(&eff_all);

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
  if (init) {
    N_A = (Nfpp * 2.0 + Nazi * 2.0) + dt * 2.0;
    Aeq->data[Aeq->size[0] * ((int)(N_A + 1.0) - 1)] = -1.0;

    /* s1 */
    Aeq->data[1 + Aeq->size[0] * ((int)(N_A + 2.0) - 1)] = -1.0;

    /* s2 */
    Aeq->data[2 + Aeq->size[0] * ((int)(N_A + 3.0) - 1)] = -1.0;

    /* s3 */
    for (ixstart = 0; ixstart < 3; ixstart++) {
      H->data[((int)((i_H + (1.0 + (double)ixstart)) - 1.0) + H->size[0] * ((int)
                ((i_H + (1.0 + (double)ixstart)) - 1.0) - 1)) - 1] =
        thruster_data[0].weight_s[ixstart];

      /* !!!thruster_data(1) must be available */
    }
  }

  ix = Ho->size[0] * Ho->size[1];
  Ho->size[0] = H->size[0];
  Ho->size[1] = H->size[1];
  emxEnsureCapacity((emxArray__common *)Ho, ix, (int)sizeof(double));
  ixstart = H->size[0] * H->size[1];
  for (ix = 0; ix < ixstart; ix++) {
    Ho->data[ix] = H->data[ix];
  }

  ix = fo->size[0];
  fo->size[0] = f->size[0];
  emxEnsureCapacity((emxArray__common *)fo, ix, (int)sizeof(double));
  ixstart = f->size[0];
  for (ix = 0; ix < ixstart; ix++) {
    fo->data[ix] = f->data[ix];
  }

  ix = Aeqo->size[0] * Aeqo->size[1];
  Aeqo->size[0] = Aeq->size[0];
  Aeqo->size[1] = Aeq->size[1];
  emxEnsureCapacity((emxArray__common *)Aeqo, ix, (int)sizeof(double));
  ixstart = Aeq->size[0] * Aeq->size[1];
  for (ix = 0; ix < ixstart; ix++) {
    Aeqo->data[ix] = Aeq->data[ix];
  }

  ix = beqo->size[0];
  beqo->size[0] = beq->size[0];
  emxEnsureCapacity((emxArray__common *)beqo, ix, (int)sizeof(double));
  ixstart = beq->size[0];
  for (ix = 0; ix < ixstart; ix++) {
    beqo->data[ix] = beq->data[ix];
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void pre_qp_free(void)
{
  emxFree_real_T(&beq);
  emxFree_real_T(&Aeq);
  emxFree_real_T(&f);
  emxFree_real_T(&H);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void pre_qp_init(void)
{
  emxInit_real_T(&beq, 1);
  emxInit_real_T1(&Aeq, 2);
  emxInit_real_T(&f, 1);
  emxInit_real_T1(&H, 2);
  beq_not_empty = false;
  Aeq_not_empty = false;
  f_not_empty = false;
  H_not_empty = false;
}

/*
 * File trailer for pre_qp.c
 *
 * [EOF]
 */
