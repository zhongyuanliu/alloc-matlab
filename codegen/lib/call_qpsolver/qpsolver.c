/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpsolver.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "call_qpsolver.h"
#include "qpsolver.h"
#include "call_qpsolver_emxutil.h"
#include "polyval.h"
#include "sign.h"
#include "QPACM.h"
#include "mldivide.h"
#include "Mpcqpsolver.h"
#include "inv.h"
#include "Chol_fc.h"
#include "eye.h"
#include "polyfit.h"

/* Variable Definitions */
static emxArray_boolean_T *iA;
static boolean_T iA_not_empty;
static double Rangle_T[3];
static double Rangle_Fangle[3];
static double Fangle_Rangle[3];
static boolean_T rudder_dat_init_not_empty;

/* Function Declarations */
static void N_sided_linear(double A[46]);
static void N_sided_rudder(const double points[18], double A[36], double b[18]);
static void angleMaxMin(double angle_start, double angle_end, double phi, double
  phi1, double phi2, double *phi_, double *phiPlus);
static double rt_atan2d_snf(double u0, double u1);
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
  double effective_dist;
  double start_phi1;
  double phi1_end;
  double phi1_phi2;
  double phi1_start;
  double phi1_phi;
  double start_phi;
  double start_phi2;

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

  /* % angle distance  */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  effective_dist = angle_end - angle_start;
  if (effective_dist < 0.0) {
    effective_dist += 6.2831853071795862;
  }

  /* % angle distance  */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  start_phi1 = phi1 - angle_start;
  if (start_phi1 < 0.0) {
    start_phi1 += 6.2831853071795862;
  }

  /* % angle distance  */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  phi1_end = angle_end - phi1;
  if (phi1_end < 0.0) {
    phi1_end += 6.2831853071795862;
  }

  /* % angle distance  */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  phi1_phi2 = phi2 - phi1;
  if (phi1_phi2 < 0.0) {
    phi1_phi2 += 6.2831853071795862;
  }

  /* % angle distance  */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  phi1_start = angle_start - phi1;
  if (phi1_start < 0.0) {
    phi1_start += 6.2831853071795862;
  }

  /* % angle distance  */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  phi1_phi = phi - phi1;
  if (phi1_phi < 0.0) {
    phi1_phi += 6.2831853071795862;
  }

  /* % angle distance  */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  start_phi = phi - angle_start;
  if (start_phi < 0.0) {
    start_phi += 6.2831853071795862;
  }

  /* % angle distance  */
  /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
  start_phi2 = phi2 - angle_start;
  if (start_phi2 < 0.0) {
    start_phi2 += 6.2831853071795862;
  }

  *phi_ = 0.0;
  *phiPlus = 0.0;
  if (effective_dist >= start_phi1) {
    /* phi1 is at the range from start to end */
    if (phi1_end >= phi1_phi2) {
      /* phi2 is at the range from phi1 to end */
      *phi_ = phi1;
      *phiPlus = phi2;
    } else if ((phi1_phi2 > phi1_end) && (phi1_phi2 <= phi1_start)) {
      *phi_ = phi1;
      *phiPlus = angle_end;
    } else {
      /* phi2 is at the range from start to phi1  */
      if (phi1_phi <= phi1_end) {
        /* phi is at the range from phi1 to end  */
        *phi_ = phi1;
        *phiPlus = angle_end;
      } else {
        if (start_phi <= start_phi2) {
          /* phi is at the range from start to phi2  */
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
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2(b_u0, b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
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
 * method == 1 qpkwik
 *  method == 2 primal active set
 * Arguments    : struct_T thruster_data_data[]
 *                double T_r_Tx
 *                double T_r_Ty
 *                double T_r_Tm
 *                double N_enabled_thruster
 *                double method
 *                emxArray_real_T *solution
 *                double *alloc_out_label
 *                double *alloc_out_enable
 *                double *alloc_out_TSP
 *                double *alloc_out_ASP
 *                double *alloc_out_Tx
 *                double *alloc_out_Ty
 *                double *alloc_out_T
 *                double *alloc_out_phi
 *                double *alloc_out_Tm
 *                double *status
 * Return Type  : void
 */
void qpsolver(struct_T thruster_data_data[], double T_r_Tx, double T_r_Ty,
              double T_r_Tm, double N_enabled_thruster, double method,
              emxArray_real_T *solution, double *alloc_out_label, double
              *alloc_out_enable, double *alloc_out_TSP, double *alloc_out_ASP,
              double *alloc_out_Tx, double *alloc_out_Ty, double *alloc_out_T,
              double *alloc_out_phi, double *alloc_out_Tm, double *status)
{
  emxArray_real_T *x0;
  static const double dv0[9] = { 0.0, 0.0872664625997165, 0.174532925199433,
    0.261799387799149, 0.349065850398866, 0.436332312998582, 0.523598775598299,
    0.610865238198015, 0.698131700797732 };

  static const double dv1[9] = { 0.98, 0.965, 0.956, 0.92, 0.895, 0.84, 0.811,
    0.77, 0.743 };

  double Rangle_x[3];
  double expl_temp[9];
  double b_expl_temp;
  double c_expl_temp;
  static const double dv2[9] = { 0.0, 0.08, 0.163, 0.2, 0.287, 0.3, 0.341, 0.345,
    0.347 };

  static const double dv3[9] = { 0.98, 0.968310384122777, 0.969796370378854,
    0.941488183675186, 0.9398904191447, 0.89196412483911, 0.879773834573409,
    0.843756481456587, 0.820035365091043 };

  static const double dv4[9] = { 0.0, 0.082712415448173, 0.168878105792867,
    0.214060683563822, 0.310310945764607, 0.343023940420703, 0.398026222516986,
    0.421232744089294, 0.436921840812821 };

  double Ntunnel;
  double Nazi;
  double Nfpp;
  double A_azi[46];
  int i0;
  int loop_ub;
  int j;
  emxArray_real_T *H;
  emxArray_real_T *Aeq;
  double Nvar;
  double N_A;
  double b_thruster_data_data;
  double c_thruster_data_data;
  double temp1;
  double temp2;
  emxArray_real_T *beq;
  emxArray_real_T *A;
  emxArray_real_T *b;
  unsigned int i_H;
  double numOfTunnel;
  double numOfAzi;
  double numOfFpp;
  double pos_A;
  unsigned int pos_Aeq;
  static const double y[18] = { 70070.0, 68997.5, 68354.0, 65780.0, 63992.5,
    60060.0, 57986.500000000007, 55055.0, 53124.5, 0.0, 5720.0, 11654.5, 14300.0,
    20520.5, 21450.0, 24381.5, 24667.499999999996, 24810.5 };

  double A_fpp[36];
  double b_fpp[18];
  int i;
  int ixy;
  int i_Azi;
  emxArray_real_T *L;
  double p;
  int i_Fpp;
  unsigned int unnamed_idx_0;
  int i1;
  emxArray_real_T *r1;
  int ic;
  emxArray_boolean_T *r2;
  int i2;
  int i3;
  emxArray_real_T *b_solution;
  emxArray_real_T *b_x0;
  double b_y[3];
  int ar;
  emxArray_real_T *d_expl_temp;
  emxArray_real_T *e_expl_temp;
  int ib;
  int ia;
  emxArray_real_T *b_Aeq;
  double b_beq[3];
  emxArray_int32_T *r3;
  emxArray_real_T *r4;
  emxArray_real_T *b_A;
  emxArray_real_T *b_b;

  /* % rudder table */
  if (!rudder_dat_init_not_empty) {
    /* assume this is only one type of rudder */
    rudder_dat_init_not_empty = true;
    polyfit(dv0, dv1, Rangle_x, expl_temp, &b_expl_temp, &c_expl_temp);

    /* rudder angle vs Fx */
    polyfit(dv0, dv2, Rangle_x, expl_temp, &b_expl_temp, &c_expl_temp);

    /* rudder angle vs Fy */
    polyfit(dv0, dv3, Rangle_T, expl_temp, &b_expl_temp, &c_expl_temp);

    /* rudder angle vs T */
    polyfit(dv3, dv0, Rangle_x, expl_temp, &b_expl_temp, &c_expl_temp);

    /* T vs rudder angle */
    polyfit(dv0, dv4, Rangle_Fangle, expl_temp, &b_expl_temp, &c_expl_temp);

    /* rudder angle vs F angle */
    polyfit(dv0, dv4, Fangle_Rangle, expl_temp, &b_expl_temp, &c_expl_temp);

    /* F angle vs rudder angle */
  }

  emxInit_real_T(&x0, 1);
  *alloc_out_label = 1.0;
  *alloc_out_enable = 1.0;
  *alloc_out_TSP = 0.0;
  *alloc_out_ASP = 0.0;
  *alloc_out_Tx = 0.0;
  *alloc_out_Ty = 0.0;
  *alloc_out_T = 0.0;
  *alloc_out_phi = 0.0;
  *alloc_out_Tm = 0.0;

  /*  phi_f=atan2(rudder_table(:,3),rudder_table(:,2));%angle of Fx Fy */
  /*  rudder_table(:,end+1)=phi_f; */
  /* % ---------------precalculation------------------ */
  Ntunnel = 0.0;
  Nazi = 0.0;
  Nfpp = 0.0;
  N_sided_linear(A_azi);

  /*  NA_azi = 0;%not used at the moment!!!!!!!!!!!!!!!! */
  i0 = x0->size[0];
  x0->size[0] = (int)(N_enabled_thruster * 2.0 + 3.0);
  emxEnsureCapacity((emxArray__common *)x0, i0, (int)sizeof(double));
  loop_ub = (int)(N_enabled_thruster * 2.0 + 3.0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    x0->data[i0] = 0.0;
  }

  i0 = solution->size[0];
  solution->size[0] = (int)(N_enabled_thruster * 2.0 + 3.0);
  emxEnsureCapacity((emxArray__common *)solution, i0, (int)sizeof(double));
  loop_ub = (int)(N_enabled_thruster * 2.0 + 3.0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    solution->data[i0] = 0.0;
  }

  *status = 0.0;
  for (j = 0; j < (int)N_enabled_thruster; j++) {
    x0->data[(int)(2.0 * (1.0 + (double)j)) - 1] = thruster_data_data[0].T * cos
      (thruster_data_data[0].phi);
    x0->data[2] = thruster_data_data[0].T * sin(thruster_data_data[0].phi);
    switch ((int)thruster_data_data[0].type) {
     case 4:
      /* tunnel */
      if (thruster_data_data[0].T + 3575.0 <= 71500.0) {
        thruster_data_data[0].Tplus = thruster_data_data[0].T + 3575.0;
      } else {
        thruster_data_data[0].Tplus = 71500.0;
      }

      if (thruster_data_data[0].T - 3575.0 >= 0.0) {
        thruster_data_data[0].T_ = thruster_data_data[0].T - 3575.0;
      } else {
        thruster_data_data[0].T_ = 0.0;
      }

      Ntunnel++;
      break;

     case 3:
      /* azi */
      if (thruster_data_data[0].T + 3575.0 <= 71500.0) {
        thruster_data_data[0].Tplus = thruster_data_data[0].T + 3575.0;
      } else {
        thruster_data_data[0].Tplus = 71500.0;
      }

      if (thruster_data_data[0].T - 3575.0 >= 0.0) {
        thruster_data_data[0].T_ = thruster_data_data[0].T - 3575.0;
      } else {
        thruster_data_data[0].T_ = 0.0;
      }

      angleMaxMin(thruster_data_data[0].phi_min, thruster_data_data[0].phi_max,
                  thruster_data_data[0].phi, thruster_data_data[0].phi - 0.05,
                  thruster_data_data[0].phi + 0.05, &c_expl_temp, &b_expl_temp);
      thruster_data_data[0].phi_ = c_expl_temp;
      thruster_data_data[0].phiPlus = b_expl_temp;

      /*              thruster_data(j).phiPlus = min(thruster_data(j).phi + ... */
      /*                  dt * thruster_data(j).dphi_max, thruster_data(j).phi_max); */
      /*               */
      /*              thruster_data(j).phi_ = max(thruster_data(j).phi - ... */
      /*                  dt * thruster_data(j).dphi_max, thruster_data(j).phi_min); */
      Nazi++;
      break;

     case 2:
      /* fpp */
      if (thruster_data_data[0].T + 3575.0 <= 71500.0) {
        b_thruster_data_data = thruster_data_data[0].T + 3575.0;
      } else {
        b_thruster_data_data = 71500.0;
      }

      thruster_data_data[0].Tplus = b_thruster_data_data * polyval(Rangle_T,
        fabs(thruster_data_data[0].phi));
      if (thruster_data_data[0].T - 3575.0 >= 0.0) {
        c_thruster_data_data = thruster_data_data[0].T - 3575.0;
      } else {
        c_thruster_data_data = 0.0;
      }

      thruster_data_data[0].T_ = c_thruster_data_data * polyval(Rangle_T, fabs
        (thruster_data_data[0].phi));
      angleMaxMin(thruster_data_data[0].phi_min, thruster_data_data[0].phi_max,
                  thruster_data_data[0].phi, thruster_data_data[0].phi - 0.05,
                  thruster_data_data[0].phi + 0.05, &temp1, &temp2);
      c_expl_temp = temp2;
      b_sign(&c_expl_temp);
      thruster_data_data[0].phiPlus = c_expl_temp * polyval(Rangle_Fangle, fabs
        (temp2));
      c_expl_temp = temp1;
      b_sign(&c_expl_temp);
      thruster_data_data[0].phi_ = c_expl_temp * polyval(Rangle_Fangle, fabs
        (temp1));

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

  emxInit_real_T1(&H, 2);
  emxInit_real_T1(&Aeq, 2);
  Nvar = ((2.0 * Ntunnel + 2.0 * Nazi) + 2.0 * Nfpp) + 3.0;

  /* number of variables in quad */
  N_A = (Ntunnel * 2.0 + Nazi * 28.0) + Nfpp * 23.0;

  /* % ------------------------------------------------------ */
  eye(Nvar, H);
  i0 = Aeq->size[0] * Aeq->size[1];
  Aeq->size[0] = (int)(3.0 + Ntunnel);
  Aeq->size[1] = (int)Nvar;
  emxEnsureCapacity((emxArray__common *)Aeq, i0, (int)sizeof(double));
  loop_ub = (int)(3.0 + Ntunnel) * (int)Nvar;
  for (i0 = 0; i0 < loop_ub; i0++) {
    Aeq->data[i0] = 0.0;
  }

  emxInit_real_T(&beq, 1);
  i0 = beq->size[0];
  beq->size[0] = (int)(3.0 + Ntunnel);
  emxEnsureCapacity((emxArray__common *)beq, i0, (int)sizeof(double));
  loop_ub = (int)(3.0 + Ntunnel);
  for (i0 = 0; i0 < loop_ub; i0++) {
    beq->data[i0] = 0.0;
  }

  emxInit_real_T1(&A, 2);
  beq->data[0] = T_r_Tx;
  beq->data[1] = T_r_Ty;
  beq->data[2] = T_r_Tm;
  i0 = A->size[0] * A->size[1];
  A->size[0] = (int)N_A;
  A->size[1] = (int)Nvar;
  emxEnsureCapacity((emxArray__common *)A, i0, (int)sizeof(double));
  loop_ub = (int)N_A * (int)Nvar;
  for (i0 = 0; i0 < loop_ub; i0++) {
    A->data[i0] = 0.0;
  }

  emxInit_real_T(&b, 1);
  i0 = b->size[0];
  b->size[0] = (int)N_A;
  emxEnsureCapacity((emxArray__common *)b, i0, (int)sizeof(double));
  loop_ub = (int)N_A;
  for (i0 = 0; i0 < loop_ub; i0++) {
    b->data[i0] = 0.0;
  }

  /*  ------------------------------------------------------ */
  i_H = 2U;
  numOfTunnel = 0.0;
  numOfAzi = 0.0;
  numOfFpp = 0.0;
  pos_A = 1.0;
  pos_Aeq = 4U;
  for (j = 0; j < (int)N_enabled_thruster; j++) {
    switch ((int)thruster_data_data[0].type) {
     case 4:
      /* tunnel */
      H->data[((int)i_H + H->size[0] * ((int)i_H - 2)) - 2] = 1.0;
      H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = 1.0;
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
            2.0) + numOfFpp * 2.0) + 1.0) - 1)] = -thruster_data_data[0].y;
          Aeq->data[i + Aeq->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
            2.0) + numOfFpp * 2.0) + 2.0) - 1)] = thruster_data_data[0].x;
          break;
        }
      }

      Aeq->data[((int)pos_Aeq + Aeq->size[0] * ((int)(((numOfAzi * 2.0 +
        numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) - 1] = sin
        (thruster_data_data[0].phi);
      Aeq->data[((int)pos_Aeq + Aeq->size[0] * ((int)(((numOfAzi * 2.0 +
        numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) - 1] = -cos
        (thruster_data_data[0].phi);
      beq->data[(int)pos_Aeq - 1] = 0.0;
      A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
        2.0) + numOfFpp * 2.0) + 1.0) - 1)) - 1] = -cos(thruster_data_data[0].
        phi);
      A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 + numOfTunnel *
        2.0) + numOfFpp * 2.0) + 2.0) - 1)) - 1] = -sin(thruster_data_data[0].
        phi);
      b->data[(int)pos_A - 1] = -thruster_data_data[0].Tplus;
      A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
        numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = cos
        (thruster_data_data[0].phi);
      A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
        numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = sin
        (thruster_data_data[0].phi);
      b->data[(int)(unsigned int)pos_A] = thruster_data_data[0].T_;

      /*              check_rudder_constr(thruster_data(j).Tmax*[-1 -1;-1 1;1 1;1 -1],A(pos_A:pos_A + 1, ... */
      /*                  numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 1:numOfAzi * 2 + numOfTunnel * 2 + numOfFpp * 2 + 2), ... */
      /*                  b(pos_A:pos_A + 1)); */
      pos_A += 2.0;
      pos_Aeq++;
      numOfTunnel++;
      break;

     case 3:
      /* azi */
      H->data[((int)i_H + H->size[0] * ((int)i_H - 2)) - 2] = 1.0;
      H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = 1.0;
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
                * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = -thruster_data_data[0].y;
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
                * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = thruster_data_data[0].x;
              break;
            }
          }
        }

        if (1 + ixy == 1) {
          A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
            numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) - 1] = -sin
            (thruster_data_data[0].phi_);

          /* (3.15a) */
          b->data[(int)pos_A - 1] = -0.001 * fabs(cos(thruster_data_data[0].phi_))
            * (thruster_data_data[0].Tplus - thruster_data_data[0].T_);
          A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = sin
            (thruster_data_data[0].phiPlus);

          /* (3.15b) */
          b->data[(int)(unsigned int)pos_A] = -0.001 * fabs(cos
            (thruster_data_data[0].phiPlus)) * (thruster_data_data[0].Tplus -
            thruster_data_data[0].T_);
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 1] = cos
            (thruster_data_data[0].phi);

          /* (3.16a) */
          b->data[(int)(unsigned int)pos_A + 1] = thruster_data_data[0].T_;
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 2] = -cos
            (thruster_data_data[0].phi_);

          /* (3.16b) */
          b->data[(int)(unsigned int)pos_A + 2] = -cos(thruster_data_data[0].
            phi_ - thruster_data_data[0].phi) * thruster_data_data[0].Tplus;
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 3] = -cos
            (thruster_data_data[0].phiPlus);

          /* (3.16c) */
          b->data[(int)(unsigned int)pos_A + 3] = -cos(thruster_data_data[0].
            phiPlus - thruster_data_data[0].phi) * thruster_data_data[0].Tplus;
          for (i_Azi = 0; i_Azi < 23; i_Azi++) {
            /* linearized thrust zone */
            A->data[((int)((unsigned int)pos_A + (1 + i_Azi)) + A->size[0] *
                     ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) + numOfFpp *
                             2.0) + 1.0) - 1)) + 3] = A_azi[i_Azi];
            b->data[(int)((unsigned int)pos_A + (1 + i_Azi)) + 3] =
              -70834.045141597642;
          }
        } else {
          A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
            numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) - 1] = cos
            (thruster_data_data[0].phi_);

          /* (3.15a) */
          A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = -cos
            (thruster_data_data[0].phiPlus);

          /* (3.15b) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 1] = sin
            (thruster_data_data[0].phi);

          /* (3.16a) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 2] = -sin
            (thruster_data_data[0].phi_);

          /* (3.16b) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 3] = -sin
            (thruster_data_data[0].phiPlus);

          /* (3.16c) */
          for (i_Azi = 0; i_Azi < 23; i_Azi++) {
            /* linearized thrust zone */
            A->data[((int)((unsigned int)pos_A + (1 + i_Azi)) + A->size[0] *
                     ((int)(((numOfAzi * 2.0 + numOfTunnel * 2.0) + numOfFpp *
                             2.0) + 2.0) - 1)) + 3] = A_azi[23 + i_Azi];
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
      pos_A = (pos_A + 5.0) + 23.0;
      numOfAzi++;
      break;

     case 2:
      /* fpp */
      N_sided_rudder(y, A_fpp, b_fpp);

      /* x,y */
      H->data[((int)i_H + H->size[0] * ((int)i_H - 2)) - 2] = 1.0;
      H->data[((int)i_H + H->size[0] * ((int)i_H - 1)) - 1] = 1.0;
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
                * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = -thruster_data_data[0].y;
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
                * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = thruster_data_data[0].x;
              break;
            }
          }
        }

        if (1 + ixy == 1) {
          A->data[((int)pos_A + A->size[0] * ((int)(((numOfAzi * 2.0 +
            numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) - 1] = -sin
            (thruster_data_data[0].phi_);

          /* (3.15a) */
          b->data[(int)pos_A - 1] = -0.001 * fabs(cos(thruster_data_data[0].phi_))
            * (thruster_data_data[0].Tplus - thruster_data_data[0].T_);
          A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)] = sin
            (thruster_data_data[0].phiPlus);

          /* (3.15b) */
          b->data[(int)(unsigned int)pos_A] = -0.001 * fabs(cos
            (thruster_data_data[0].phiPlus)) * (thruster_data_data[0].Tplus -
            thruster_data_data[0].T_);
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 1] = cos
            (polyval(Rangle_Fangle, fabs(thruster_data_data[0].phi)));

          /* (3.16a) */
          b->data[(int)(unsigned int)pos_A + 1] = thruster_data_data[0].T_;
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 2] = -cos
            (thruster_data_data[0].phi_);

          /* (3.16b) */
          c_expl_temp = thruster_data_data[0].phi;
          b_sign(&c_expl_temp);
          b->data[(int)(unsigned int)pos_A + 2] = -cos(thruster_data_data[0].
            phi_ - c_expl_temp * polyval(Rangle_Fangle, fabs(thruster_data_data
            [0].phi))) * thruster_data_data[0].Tplus;
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 1.0) - 1)) + 3] = -cos
            (thruster_data_data[0].phiPlus);

          /* (3.16c) */
          c_expl_temp = thruster_data_data[0].phi;
          b_sign(&c_expl_temp);
          b->data[(int)(unsigned int)pos_A + 3] = -cos(thruster_data_data[0].
            phiPlus - c_expl_temp * polyval(Rangle_Fangle, fabs
            (thruster_data_data[0].phi))) * thruster_data_data[0].Tplus;
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
            (thruster_data_data[0].phi_);

          /* (3.15a) */
          A->data[(int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)] = -cos
            (thruster_data_data[0].phiPlus);

          /* (3.15b) */
          c_expl_temp = thruster_data_data[0].phi;
          b_sign(&c_expl_temp);
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 1] = sin
            (c_expl_temp * polyval(Rangle_Fangle, fabs(thruster_data_data[0].phi)));

          /* (3.16a) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 2] = -sin
            (thruster_data_data[0].phi_);

          /* (3.16b) */
          A->data[((int)(unsigned int)pos_A + A->size[0] * ((int)(((numOfAzi *
            2.0 + numOfTunnel * 2.0) + numOfFpp * 2.0) + 2.0) - 1)) + 3] = -sin
            (thruster_data_data[0].phiPlus);

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
      = thruster_data_data[0].weight_s[i];

    /* !!!thruster_data(1) must be available */
  }

  emxInit_real_T1(&L, 2);
  Chol_fc(H, L, &p);

  /*  before chol */
  if (p == 0.0) {
    switch ((int)method) {
     case 1:
      inv(L, H);
      unnamed_idx_0 = (unsigned int)b->size[0];
      if (!iA_not_empty) {
        i0 = iA->size[0];
        iA->size[0] = (int)unnamed_idx_0;
        emxEnsureCapacity((emxArray__common *)iA, i0, (int)sizeof(boolean_T));
        loop_ub = (int)unnamed_idx_0;
        for (i0 = 0; i0 < loop_ub; i0++) {
          iA->data[i0] = false;
        }

        iA_not_empty = !(iA->size[0] == 0);
      }

      emxInit_real_T(&r1, 1);
      i0 = r1->size[0];
      r1->size[0] = (int)Nvar;
      emxEnsureCapacity((emxArray__common *)r1, i0, (int)sizeof(double));
      loop_ub = (int)Nvar;
      for (i0 = 0; i0 < loop_ub; i0++) {
        r1->data[i0] = 0.0;
      }

      emxInit_boolean_T(&r2, 1);
      i0 = r2->size[0];
      r2->size[0] = (int)unnamed_idx_0;
      emxEnsureCapacity((emxArray__common *)r2, i0, (int)sizeof(boolean_T));
      loop_ub = (int)unnamed_idx_0;
      for (i0 = 0; i0 < loop_ub; i0++) {
        r2->data[i0] = false;
      }

      emxInit_real_T(&b_solution, 1);
      emxInit_real_T(&d_expl_temp, 1);
      emxInit_real_T(&e_expl_temp, 1);
      Mpcqpsolver(H, r1, A, b, Aeq, beq, r2, b_solution, status, iA, d_expl_temp,
                  e_expl_temp);
      i0 = solution->size[0];
      solution->size[0] = b_solution->size[0];
      emxEnsureCapacity((emxArray__common *)solution, i0, (int)sizeof(double));
      loop_ub = b_solution->size[0];
      emxFree_boolean_T(&r2);
      emxFree_real_T(&r1);
      emxFree_real_T(&e_expl_temp);
      emxFree_real_T(&d_expl_temp);
      for (i0 = 0; i0 < loop_ub; i0++) {
        solution->data[i0] = b_solution->data[i0];
      }

      emxFree_real_T(&b_solution);
      iA_not_empty = !(iA->size[0] == 0);

      /* faster to use iA instead of iA0 */
      /*      [solution1,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq) */
      break;

     case 2:
      /*      x0=[x0;0;0;0]; */
      if (Aeq->size[1] - 2 > Aeq->size[1]) {
        i0 = 1;
        i1 = 0;
      } else {
        i0 = Aeq->size[1] - 2;
        i1 = Aeq->size[1];
      }

      if (1 > Aeq->size[1] - 3) {
        loop_ub = 0;
      } else {
        loop_ub = Aeq->size[1] - 3;
      }

      if (1 > Aeq->size[1] - 3) {
        ic = 0;
      } else {
        ic = Aeq->size[1] - 3;
      }

      if (Aeq->size[1] - 2 > x0->size[0]) {
        i2 = 1;
        i3 = 0;
      } else {
        i2 = Aeq->size[1] - 2;
        i3 = x0->size[0];
      }

      emxInit_real_T1(&b_solution, 2);
      emxInit_real_T(&b_x0, 1);
      if ((loop_ub == 1) || (ic == 1)) {
        ar = b_solution->size[0] * b_solution->size[1];
        b_solution->size[0] = 3;
        b_solution->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)b_solution, ar, (int)sizeof(double));
        for (ar = 0; ar < loop_ub; ar++) {
          for (ib = 0; ib < 3; ib++) {
            b_solution->data[ib + b_solution->size[0] * ar] = Aeq->data[ib +
              Aeq->size[0] * ar];
          }
        }

        ar = b_x0->size[0];
        b_x0->size[0] = ic;
        emxEnsureCapacity((emxArray__common *)b_x0, ar, (int)sizeof(double));
        for (ar = 0; ar < ic; ar++) {
          b_x0->data[ar] = x0->data[ar];
        }

        for (ar = 0; ar < 3; ar++) {
          b_y[ar] = 0.0;
          loop_ub = b_solution->size[1];
          for (ib = 0; ib < loop_ub; ib++) {
            b_expl_temp = b_y[ar] + b_solution->data[ar + b_solution->size[0] *
              ib] * b_x0->data[ib];
            b_y[ar] = b_expl_temp;
          }
        }
      } else {
        for (ic = 0; ic < 3; ic++) {
          b_y[ic] = 0.0;
        }

        ar = -1;
        for (ib = 0; ib + 1 <= loop_ub; ib++) {
          if (x0->data[ib] != 0.0) {
            ia = ar;
            for (ic = 0; ic < 3; ic++) {
              ia++;
              b_expl_temp = b_y[ic] + x0->data[ib] * Aeq->data[ia % 3 +
                Aeq->size[0] * (ia / 3)];
              b_y[ic] = b_expl_temp;
            }
          }

          ar += 3;
        }
      }

      emxFree_real_T(&b_x0);
      emxFree_real_T(&b_solution);
      emxInit_real_T1(&b_Aeq, 2);
      ar = b_Aeq->size[0] * b_Aeq->size[1];
      b_Aeq->size[0] = 3;
      b_Aeq->size[1] = (i1 - i0) + 1;
      emxEnsureCapacity((emxArray__common *)b_Aeq, ar, (int)sizeof(double));
      loop_ub = (i1 - i0) + 1;
      for (i1 = 0; i1 < loop_ub; i1++) {
        for (ar = 0; ar < 3; ar++) {
          b_Aeq->data[ar + b_Aeq->size[0] * i1] = Aeq->data[ar + Aeq->size[0] *
            ((i0 + i1) - 1)];
        }
      }

      for (i0 = 0; i0 < 3; i0++) {
        b_beq[i0] = beq->data[i0] - b_y[i0];
      }

      emxInit_real_T(&r1, 1);
      emxInit_int32_T(&r3, 2);
      mldivide(b_Aeq, b_beq, r1);
      i0 = r3->size[0] * r3->size[1];
      r3->size[0] = 1;
      r3->size[1] = (i3 - i2) + 1;
      emxEnsureCapacity((emxArray__common *)r3, i0, (int)sizeof(int));
      loop_ub = (i3 - i2) + 1;
      emxFree_real_T(&b_Aeq);
      for (i0 = 0; i0 < loop_ub; i0++) {
        r3->data[r3->size[0] * i0] = (i2 + i0) - 1;
      }

      loop_ub = r3->size[0] * r3->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        x0->data[r3->data[i0]] = r1->data[i0];
      }

      emxFree_int32_T(&r3);
      emxFree_real_T(&r1);

      /* does not work if there is tunnel!!!!!!!!! */
      i0 = solution->size[0];
      solution->size[0] = x0->size[0];
      emxEnsureCapacity((emxArray__common *)solution, i0, (int)sizeof(double));
      loop_ub = x0->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        solution->data[i0] = x0->data[i0];
      }

      emxInit_real_T(&r4, 1);
      i0 = r4->size[0];
      r4->size[0] = (int)Nvar;
      emxEnsureCapacity((emxArray__common *)r4, i0, (int)sizeof(double));
      loop_ub = (int)Nvar;
      for (i0 = 0; i0 < loop_ub; i0++) {
        r4->data[i0] = -0.0;
      }

      emxInit_real_T1(&b_A, 2);
      i0 = b_A->size[0] * b_A->size[1];
      b_A->size[0] = A->size[0];
      b_A->size[1] = A->size[1];
      emxEnsureCapacity((emxArray__common *)b_A, i0, (int)sizeof(double));
      loop_ub = A->size[0] * A->size[1];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_A->data[i0] = -A->data[i0];
      }

      emxInit_real_T(&b_b, 1);
      i0 = b_b->size[0];
      b_b->size[0] = b->size[0];
      emxEnsureCapacity((emxArray__common *)b_b, i0, (int)sizeof(double));
      loop_ub = b->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        b_b->data[i0] = -b->data[i0];
      }

      *status = QPACM(H, r4, Aeq, beq, b_A, b_b, solution);

      /*  x = quadprog(H,f,A,b,Aeq,beq); */
      emxFree_real_T(&b_b);
      emxFree_real_T(&b_A);
      emxFree_real_T(&r4);
      break;
    }
  } else {
    i0 = solution->size[0];
    solution->size[0] = (int)Nvar;
    emxEnsureCapacity((emxArray__common *)solution, i0, (int)sizeof(double));
    loop_ub = (int)Nvar;
    for (i0 = 0; i0 < loop_ub; i0++) {
      solution->data[i0] = 0.0;
    }

    *status = -2.0;
  }

  emxFree_real_T(&L);
  emxFree_real_T(&b);
  emxFree_real_T(&A);
  emxFree_real_T(&beq);
  emxFree_real_T(&Aeq);
  emxFree_real_T(&H);
  emxFree_real_T(&x0);

  /*  output */
  numOfTunnel = 0.0;
  numOfAzi = 0.0;
  numOfFpp = 0.0;
  for (i = 0; i < (int)N_enabled_thruster; i++) {
    *alloc_out_label = thruster_data_data[0].label;
    switch ((int)thruster_data_data[0].type) {
     case 4:
      *alloc_out_Tx = 0.0;
      *alloc_out_Ty = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 1.0) - 1];
      *alloc_out_T = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 1.0) - 1];
      *alloc_out_Tm = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 1.0) - 1] * thruster_data_data[0].x + 0.0 * (0.0 -
        thruster_data_data[0].y);
      numOfTunnel++;
      break;

     case 3:
      *alloc_out_Tx = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 1.0) - 1];
      *alloc_out_Ty = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 2.0) - 1];
      *alloc_out_T = sqrt(solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi *
        2.0) + numOfFpp * 2.0) + 1.0) - 1] * solution->data[(int)(((numOfTunnel *
        2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0) - 1] + solution->data
                          [(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0) +
        numOfFpp * 2.0) + 2.0) - 1] * solution->data[(int)(((numOfTunnel * 2.0 +
        numOfAzi * 2.0) + numOfFpp * 2.0) + 2.0) - 1]);
      *alloc_out_phi = rt_atan2d_snf(solution->data[(int)(((numOfTunnel * 2.0 +
        numOfAzi * 2.0) + numOfFpp * 2.0) + 2.0) - 1], solution->data[(int)
        (((numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0) - 1]);
      *alloc_out_Tm = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 2.0) - 1] * thruster_data_data[0].x + solution->
        data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0)
        - 1] * (0.0 - thruster_data_data[0].y);
      numOfAzi++;
      break;

     case 2:
      *alloc_out_Tx = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 1.0) - 1];
      *alloc_out_Ty = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 2.0) - 1];
      *alloc_out_T = sqrt(solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi *
        2.0) + numOfFpp * 2.0) + 1.0) - 1] * solution->data[(int)(((numOfTunnel *
        2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0) - 1] + solution->data
                          [(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0) +
        numOfFpp * 2.0) + 2.0) - 1] * solution->data[(int)(((numOfTunnel * 2.0 +
        numOfAzi * 2.0) + numOfFpp * 2.0) + 2.0) - 1]);
      *alloc_out_phi = rt_atan2d_snf(solution->data[(int)(((numOfTunnel * 2.0 +
        numOfAzi * 2.0) + numOfFpp * 2.0) + 2.0) - 1], solution->data[(int)
        (((numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0) - 1]);
      c_expl_temp = *alloc_out_phi;
      b_sign(&c_expl_temp);
      *alloc_out_phi = c_expl_temp * polyval(Fangle_Rangle, fabs(*alloc_out_phi));
      *alloc_out_Tm = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0)
        + numOfFpp * 2.0) + 2.0) - 1] * thruster_data_data[0].x + solution->
        data[(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0)
        - 1] * (0.0 - thruster_data_data[0].y);
      numOfFpp++;
      break;
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_free(void)
{
  emxFree_boolean_T(&iA);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_init(void)
{
  emxInit_boolean_T(&iA, 1);
  iA_not_empty = false;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void rudder_dat_init_not_empty_init(void)
{
  rudder_dat_init_not_empty = false;
}

/*
 * File trailer for qpsolver.c
 *
 * [EOF]
 */
