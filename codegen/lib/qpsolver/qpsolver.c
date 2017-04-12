/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpsolver.c
 *
 * MATLAB Coder version            : 3.2
 * C/C++ source code generated on  : 07-Apr-2017 15:39:37
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "polyval.h"
#include "sign.h"
#include "angle_dist.h"
#include "mod.h"
#include "angleMaxMin.h"
#include "QPACM.h"
#include "qpsolver_emxutil.h"
#include "mldivide.h"
#include "Mpcqpsolver.h"
#include "inv.h"
#include "Chol_fc.h"
#include "pre_qp.h"
#include "structConstructorHelper.h"
#include "polyfit.h"

/* Variable Definitions */
static double Rangle_T[3];
static double Rangle_Fangle[3];
static double Fangle_Rangle[3];
static boolean_T rudder_dat_init_not_empty;
static double on_border[8];
static double phi_guider_old[8];

/* Function Declarations */
static void angle_point_to(double s, double phi, double dphi_max, double dt,
  double phi_min, double phi_max, double *phi_0, double *phiPlus0, double
  *inside_fz);
static void baseThrust(double *T, double phi, double T_base, double phi_base,
  double phi_min, double phi_max, double *new, double *in_fb);
static double rt_atan2d_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double s
 *                double phi
 *                double dphi_max
 *                double dt
 *                double phi_min
 *                double phi_max
 *                double *phi_0
 *                double *phiPlus0
 *                double *inside_fz
 * Return Type  : void
 */
static void angle_point_to(double s, double phi, double dphi_max, double dt,
  double phi_min, double phi_max, double *phi_0, double *phiPlus0, double
  *inside_fz)
{
  double move;
  double varargin_3;
  double b_move;

  /* s: setpoint from global solution */
  *inside_fz = 0.0;

  /* phi, phi_0 and phiPlus0 is outside forbidden zone */
  move = angle_dist(phi, s);
  varargin_3 = angle_dist(s, phi);
  if ((move <= varargin_3) || rtIsNaN(varargin_3)) {
    b_move = move;
  } else {
    b_move = varargin_3;
  }

  if (b_move > dphi_max * dt) {
    move = dphi_max * dt / 2.0;
  } else {
    move = 0.0;
  }

  /*  find out if phi, phi_0 or phiPlus0 is inside fb zone */
  if ((fabs(phi_max - phi_min) > 1.0E-6) && (angle_dist(phi_min, phi) >=
       angle_dist(phi_min, phi_max))) {
    /* means there is forbidden zone */
    /*  ... */
    /*          || angle_dist(phi_min,phi_0) >= angle_dist(phi_min,phi_max) ... */
    /*          || angle_dist(phi_min,phiPlus0) >= angle_dist(phi_min,phi_max) */
    /*        commented at v4.10 */
    *inside_fz = 1.0;

    /* phi, phi_0 or phiPlus0 is inside forbidden zone */
  }

  /*  % % % % % % % % */
  if (fabs(phi_max - phi_min) > 1.0E-6) {
    /* means there is forbidden zone */
    if (*inside_fz == 0.0) {
      /* phi is outside forbidden zone */
      if (angle_dist(phi, s) < 3.1415926535897931) {
        *phi_0 = phi + move;
        *phiPlus0 = phi + dphi_max * dt;
      } else {
        *phiPlus0 = phi - move;
        *phi_0 = phi - dphi_max * dt;
      }
    } else {
      /* phi is inside forbidden zone */
      if (angle_dist(phi_min, s) <= angle_dist(phi_min, phi_max)) {
        /* s is outside forbidden zone */
        if (angle_dist(phi, s) < 3.1415926535897931) {
          *phi_0 = phi + move;
          *phiPlus0 = phi + dphi_max * dt;
        } else {
          *phiPlus0 = phi - move;
          *phi_0 = phi - dphi_max * dt;
        }
      } else {
        /*  s is inside forbidden zone */
        if (angle_dist(phi_max, phi) < angle_dist(phi, phi_min)) {
          /* phi is closer to phi_max */
          *phiPlus0 = phi - dphi_max * dt / 2.0;
          *phi_0 = phi - dphi_max * dt;
        } else {
          /* phi is closer to phi_min */
          *phi_0 = phi + dphi_max * dt / 2.0;
          *phiPlus0 = phi + dphi_max * dt;
        }
      }
    }
  } else {
    /* there is no forbidden zone */
    if (angle_dist(phi, s) < 3.1415926535897931) {
      *phi_0 = phi + move;
      *phiPlus0 = phi + dphi_max * dt;
    } else {
      *phiPlus0 = phi - move;
      *phi_0 = phi - dphi_max * dt;
    }
  }
}

/*
 * Arguments    : double *T
 *                double phi
 *                double T_base
 *                double phi_base
 *                double phi_min
 *                double phi_max
 *                double *new
 *                double *in_fb
 * Return Type  : void
 */
static void baseThrust(double *T, double phi, double T_base, double phi_base,
  double phi_min, double phi_max, double *new, double *in_fb)
{
  double x;
  double y;
  double b_x[2];
  double scale;
  int k;
  double absxk;
  double t;
  x = *T * cos(phi) + T_base * cos(phi_base);
  y = *T * sin(phi) + T_base * sin(phi_base);
  b_x[0] = x;
  b_x[1] = y;
  *T = 0.0;
  scale = 2.2250738585072014E-308;
  for (k = 0; k < 2; k++) {
    absxk = fabs(b_x[k]);
    if (absxk > scale) {
      t = scale / absxk;
      *T = 1.0 + *T * t * t;
      scale = absxk;
    } else {
      t = absxk / scale;
      *T += t * t;
    }
  }

  *T = scale * sqrt(*T);
  *new = rt_atan2d_snf(y, x);
  *in_fb = 0.0;
  if ((fabs(phi_max - phi_min) > 1.0E-6) && (angle_dist(phi_min, *new) >
       angle_dist(phi_min, phi_max))) {
    /* means there is forbidden zone */
    /* inside fbzone */
    *in_fb = 1.0;
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
 * method == 1 qpkwik
 *  method == 2 sef made primal active set
 * Arguments    : boolean_T init
 *                struct0_T thruster_data[8]
 *                const struct1_T *T_r
 *                double N_enabled_thruster
 *                const double rudder_table0[45]
 *                double no_azi_angle_constr
 *                double method
 *                emxArray_real_T *solution
 *                struct2_T alloc_out[8]
 *                double *status
 * Return Type  : void
 */
void qpsolver(boolean_T init, struct0_T thruster_data[8], const struct1_T *T_r,
              double N_enabled_thruster, const double rudder_table0[45], double
              no_azi_angle_constr, double method, emxArray_real_T *solution,
              struct2_T alloc_out[8], double *status)
{
  double rudder_table[45];
  int i12;
  int i13;
  double expl_temp[9];
  double diff_phi;
  double numOfTunnel;
  double Ce3[74];
  static const double dv1[37] = { 0.0, 0.174532925199433, 0.349065850398866,
    0.523598775598299, 0.698131700797732, 0.872664625997165, 1.0471975511966,
    1.22173047639603, 1.39626340159546, 1.5707963267949, 1.74532925199433,
    1.91986217719376, 2.0943951023932, 2.26892802759263, 2.44346095279206,
    2.61799387799149, 2.79252680319093, 2.96705972839036, 3.14159265358979,
    3.31612557878923, 3.49065850398866, 3.66519142918809, 3.83972435438753,
    4.01425727958696, 4.18879020478639, 4.36332312998582, 4.53785605518526,
    4.71238898038469, 4.88692190558412, 5.06145483078356, 5.23598775598299,
    5.41052068118242, 5.58505360638185, 5.75958653158129, 5.93411945678072,
    6.10865238198015, 0.0 };

  cell_wrap_1 r2;
  cell_wrap_1 r3;
  cell_wrap_1 r4;
  cell_wrap_1 r5;
  cell_wrap_1 r6;
  cell_wrap_1 r7;
  cell_wrap_1 r8;
  cell_wrap_1 r9;
  double dv2[8];
  cell_wrap_1 rv0[8];
  struct_T Ce[8];
  double dt;
  int loop_ub;
  double phi_0[8];
  double phiPlus0[8];
  double inside_fz[8];
  double alloc_out_T[8];
  double alloc_out_phi[8];
  emxArray_real_T *x0;
  emxArray_real_T *Linv;
  emxArray_real_T *b_f;
  emxArray_real_T *A;
  emxArray_real_T *b;
  emxArray_real_T *b_Aeq;
  emxArray_real_T *b_beq;
  emxArray_real_T *L;
  emxArray_real_T *b_solution;
  emxArray_boolean_T *unusedU3;
  emxArray_real_T *r10;
  emxArray_int32_T *r11;
  emxArray_boolean_T *r12;
  emxArray_real_T *c_f;
  emxArray_real_T *b_A;
  emxArray_real_T *b_b;
  emxArray_real_T *c_Aeq;
  emxArray_real_T *d_Aeq;
  emxArray_real_T *b_x0;
  int step;
  int i;
  double counter;
  double phi0[8];
  boolean_T run;
  int flag_run;
  int count_run;
  boolean_T exitg1;
  boolean_T use_base;
  unsigned int unnamed_idx_0;
  double numOfAzi;
  double numOfFpp;
  int i14;
  int i15;
  double y[3];
  int ar;
  int ib;
  int ia;
  double c_beq[3];
  (void)no_azi_angle_constr;

  /*  solves a quadratic program to determine an optimal solution, x. */
  /*  It minimizes the quadratic objective function, J = 0.5*x'*H*x + f'*x, */
  /*  subject to linear inequality constraints, A*x >= b, and linear equality */
  /*  constraints, Aeq*x = beq, where x is a column vector of length n. */
  /*  */
  /*  H is the n-by-n hessian matrix, which must be symmetric and positive */
  /*  definite.  L is its lower-triangular Cholesky decomposition and Linv is */
  /*  the inverse of L.  Consequently, Linv can be computed from H as follows */
  /*  in MATLAB: */
  /*    L = chol(H, 'lower') */
  /*    Linv = L\eye(size(H,1)) */
  /*  Note that H = L*L' */
  /*  */
  /*  Input arguments (all are mandatory): */
  /*     Linv:  n-by-n (n>0) matrix representing the inverse of L. */
  /*        f:  column vector of length n (n>0). */
  /*        A:  m-by-n matrix of linear inequality constraint coefficients.  If */
  /*            your problem has no inequality constraints, use []. */
  /*        b:  column vector of length m, the right-hand side of A*x >= b.  If */
  /*            your problem has no inequality constraints, use "zeros(0,1)". */
  /*      Aeq:  q-by-n matrix of linear equality constraint coefficients, */
  /*            q <= n.  If your problem has no equality constraints, use []. */
  /*            NOTE:  equality constraints must be linearly independent with */
  /*            rank(Aeq) = q. */
  /*      beq:  column vector of length q, the right-hand side of Aeq*x = beq. */
  /*            If your problem has no equality constraints, use "zeros(0,1)". */
  /*      iA0:  logical vector of length m.  If your problem has no inequality */
  /*            constraints, use "false(0,1)".  For a "cold start": use */
  /*            false(m,1).  For a "warm start": if iA0(i) == true, the */
  /*            algorithm begins with A(i,:)*x = b(i), i.e., with the ith */
  /*            inequality active. */
  /*            Tips: */
  /*            (1) Normally, you should use the optional output argument iA */
  /*            from a previous solution as the input iA0 in the next */
  /*            calculation. */
  /*            (2) if iA0(i) and iA0(j) are true, rows i and j of A should be */
  /*            linearly independent.  Otherwise, the solution may fail (status */
  /*            = -2). This should not happen if you use iA from a previous */
  /*            solution as recommended above. */
  /*  options:  Structure with following fields used by the QP solver. */
  /*  */
  /*                   DataType: string, either 'double' or 'single'.  All the */
  /*                             real input arguments to the "mpcqpsolver" */
  /*                             command must match this data type and it is */
  /*                             used in both simulation and code generation. */
  /*                             Default value is 'double'. */
  /*                    MaxIter: scalar, maximum number of iterations allowed */
  /*                             when computing QP solution.  Default value is */
  /*                             200. */
  /*             FeasibilityTol: scalar, tolerance used to verify that */
  /*                             inequality constraints are satisfied at the */
  /*                             optimal solution.  Increasing this causes */
  /*                             MPCQPSOLVER to allow larger constraint */
  /*                             violations. Default value is 1.0e-6. */
  /*  */
  /*            Use the MPCQPSOLVEROPTIONS command to create such a structure */
  /*            with default values. */
  /*  */
  /*  Output arguments: */
  /*        x: column vector of length n, representing the optimal solution. */
  /*           MPCQPSOLVER will return x in all cases, but it is likely to be */
  /*           sub-optimal and/or infeasible unless status > 0. */
  /*   status: scalar, indicating the validity of the returned x as follows: */
  /*             >  0:  x is optimal, representing the number of iterations */
  /*                    used in optimization. */
  /*            ==  0:  x was obtained when maximum number of iterations is */
  /*                    reached.  It may be sub-optimal and/or violate A*x >=b. */
  /*            == -1:  The problem appears to be infeasible.  In other words, */
  /*                    A*x >= b cannot be satisfied. */
  /*            == -2:  An unrecoverable numerical error occurred (e.g., see */
  /*                    description of iA0). */
  /*       iA: logical vector of length m, indicating the inequalities that */
  /*           are active (at equality) at the solution.  In a series of */
  /*           problems in which A is constant, you can use iA from one */
  /*           solution as the input iA0 to the next ("warm start"). */
  /*   lambda: structure of Lagrange multipliers with two fields: */
  /*            ineqlin: column vector of length m, multipliers of the */
  /*                     inequality constraints.  Must be non-negative at an */
  /*                     optimal solution. */
  /*              eqlin: column vector of length q, multipliers of the equality */
  /*                     constraints.  No sign restriction. */
  /*   command for code generation: */
  /*   cfg = coder.config('lib'); */
  /*   codegen -config cfg -c  -args {init,thruster_data,T_r,N_enabled_thruster,rudder_table,0,1} qpsolver */
  /*   update: */
  /*   1.01 add check_phi in pre_qp.m to make sure initial phi meets */
  /*   constraint. |29-11-2016 by Liu */
  /*   1.02 correct mark8 and add if status < 0 then break to avoid risk of */
  /*   infinite loop. |29-11-2016 by Liu */
  /*   1.04 change alloc_out(i).phi = thruster_data(i).phi;%output is equal to input for tunnel thruster */
  /*   1.05 improve function check_angle(phi_min,phi_,phi,phiPlus,phi_max) and angleMaxMin(angle_start,angle_end,phi,phi1,phi2) */
  /*   1.06 guarantee initial x feasible, add mark8(i) = 2 to avoid infinite */
  /*   loop, add count_run > 50 then break |03-01-2017 by Liu */
  /*   1.07 guarantee tunnel angle input feasible |23-01-2017 by Liu */
  /*   4.01 base thrust added, Tmin of azi must be zero in thrust_config !!!!! */
  /*   and bastT must be half of real Tmin. dt*dphi_max is replaced by a bigger constant |01-02-2017 by Liu */
  /*   4.04 azi use efficiency one by one |03-02-2017 by Liu */
  /*   4.05 azi efficiency is smoothy |06-02-2017 by Liu */
  /*   4.06 u is scaled. c_scale * Tmax^2 == 1 |13-03-2017 by Liu */
  for (i12 = 0; i12 < 5; i12++) {
    for (i13 = 0; i13 < 9; i13++) {
      rudder_table[i13 + 9 * i12] = rudder_table0[i12 + 5 * i13];
    }
  }

  /* !!!!!!!!!!!!!!! */
  /*  no_azi_angle_constr = 1; */
  /*  rudder table processing */
  /*  persistent iA; */
  if (!rudder_dat_init_not_empty) {
    /* assume this is only one type of rudder */
    rudder_dat_init_not_empty = true;
    polyfit(*(double (*)[9])&rudder_table[0], *(double (*)[9])&rudder_table[36],
            Rangle_T, expl_temp, &diff_phi, &numOfTunnel);

    /* rudder angle vs T */
    polyfit(*(double (*)[9])&rudder_table[0], *(double (*)[9])&rudder_table[27],
            Rangle_Fangle, expl_temp, &diff_phi, &numOfTunnel);

    /* rudder angle vs F angle */
    polyfit(*(double (*)[9])&rudder_table[27], *(double (*)[9])&rudder_table[0],
            Fangle_Rangle, expl_temp, &diff_phi, &numOfTunnel);

    /* F angle vs rudder angle */
  }

  /* max half range of phi; */
  /*  coder.varsize('alloc_out');%!!!!!!!!!very important??? */
  structConstructorHelper(alloc_out);
  for (i12 = 0; i12 < 74; i12++) {
    Ce3[i12] = 1.0;
  }

  memcpy(&Ce3[0], &dv1[0], 37U * sizeof(double));

  /*  Ce2 = Ce1; */
  /*  Ce3 = Ce1; */
  memcpy(&r2.f1[0], &Ce3[0], 74U * sizeof(double));
  memcpy(&r3.f1[0], &Ce3[0], 74U * sizeof(double));
  memcpy(&r4.f1[0], &Ce3[0], 74U * sizeof(double));
  memcpy(&r5.f1[0], &Ce3[0], 74U * sizeof(double));
  memcpy(&r6.f1[0], &Ce3[0], 74U * sizeof(double));
  memcpy(&r7.f1[0], &Ce3[0], 74U * sizeof(double));
  memcpy(&r8.f1[0], &Ce3[0], 74U * sizeof(double));
  memcpy(&r9.f1[0], &Ce3[0], 74U * sizeof(double));
  for (i12 = 0; i12 < 8; i12++) {
    dv2[i12] = 1.0 + (double)i12;
  }

  rv0[0] = r2;
  rv0[1] = r3;
  rv0[2] = r4;
  rv0[3] = r5;
  rv0[4] = r6;
  rv0[5] = r7;
  rv0[6] = r8;
  rv0[7] = r9;
  b_structConstructorHelper(dv2, rv0, Ce);
  dt = thruster_data[0].dt;
  i12 = solution->size[0];
  solution->size[0] = (int)(N_enabled_thruster * 2.0 + 3.0);
  emxEnsureCapacity((emxArray__common *)solution, i12, (int)sizeof(double));
  loop_ub = (int)(N_enabled_thruster * 2.0 + 3.0);
  for (i12 = 0; i12 < loop_ub; i12++) {
    solution->data[i12] = 0.0;
  }

  /*  phi_reserve = zeros(8,1); */
  /*  T_reserve = zeros(8,1); */
  memset(&phi_0[0], 0, sizeof(double) << 3);
  memset(&phiPlus0[0], 0, sizeof(double) << 3);
  memset(&inside_fz[0], 0, sizeof(double) << 3);
  memset(&alloc_out_T[0], 0, sizeof(double) << 3);
  memset(&alloc_out_phi[0], 0, sizeof(double) << 3);
  emxInit_real_T(&x0, 1);
  emxInit_real_T1(&Linv, 2);
  emxInit_real_T(&b_f, 1);
  emxInit_real_T1(&A, 2);
  emxInit_real_T(&b, 1);
  emxInit_real_T1(&b_Aeq, 2);
  emxInit_real_T(&b_beq, 1);
  emxInit_real_T1(&L, 2);
  emxInit_real_T(&b_solution, 1);
  emxInit_boolean_T(&unusedU3, 1);
  emxInit_real_T(&r10, 1);
  emxInit_int32_T(&r11, 2);
  emxInit_boolean_T(&r12, 1);
  emxInit_real_T(&c_f, 1);
  emxInit_real_T1(&b_A, 2);
  emxInit_real_T(&b_b, 1);
  emxInit_real_T1(&c_Aeq, 2);
  emxInit_real_T1(&d_Aeq, 2);
  emxInit_real_T(&b_x0, 1);
  for (step = 0; step < 2; step++) {
    /* step 1: global solution, step 2: local solution */
    for (i = 0; i < (int)N_enabled_thruster; i++) {
      /* added v1.09 */
      /*  0 no angle constraint no thrust constraint */
      /*  2 only angle constraint */
      /*  1 only thrust constraint */
      /*  3 means angle must be fixed because solution is inside fb zone */
      if (1 + step == 1) {
        thruster_data[i].constr = 0.0;

        /* no constraint at all */
      } else {
        thruster_data[i].constr = 2.0;

        /* no constraint at all */
      }
    }

    *status = 0.0;
    counter = 0.0;
    memset(&phi0[0], 0, sizeof(double) << 3);
    run = true;
    while (run) {
      counter++;
      flag_run = 1;
      count_run = 0;
      exitg1 = false;
      while ((!exitg1) && (flag_run != 0)) {
        flag_run = 0;
        count_run++;
        pre_qp(init, thruster_data, T_r->Tx, T_r->Ty, T_r->Tm,
               N_enabled_thruster, Ce, rudder_table, Rangle_T, Rangle_Fangle,
               &numOfTunnel, Linv, b_f, A, b, b_Aeq, b_beq, x0);
        init = false;
        Chol_fc(Linv, L, &diff_phi);

        /*  before chol */
        if (diff_phi == 0.0) {
          switch ((int)method) {
           case 1:
            inv(L, Linv);
            unnamed_idx_0 = (unsigned int)b->size[0];

            /*                          if isempty(iA) */
            /*                              iA = iA0; */
            /*                          end */
            i12 = r12->size[0];
            r12->size[0] = (int)unnamed_idx_0;
            emxEnsureCapacity((emxArray__common *)r12, i12, (int)sizeof
                              (boolean_T));
            loop_ub = (int)unnamed_idx_0;
            for (i12 = 0; i12 < loop_ub; i12++) {
              r12->data[i12] = false;
            }

            Mpcqpsolver(Linv, b_f, A, b, b_Aeq, b_beq, r12, b_solution, status,
                        unusedU3);
            i12 = solution->size[0];
            solution->size[0] = b_solution->size[0];
            emxEnsureCapacity((emxArray__common *)solution, i12, (int)sizeof
                              (double));
            loop_ub = b_solution->size[0];
            for (i12 = 0; i12 < loop_ub; i12++) {
              solution->data[i12] = b_solution->data[i12];
            }

            /* faster to use iA instead of iA0 */
            /*      [solution1,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq) */
            break;

           case 2:
            /*      x0=[x0;0;0;0]; */
            /*      ????????????????ux,uy,?????????s1?s2?s3??? */
            if (b_Aeq->size[1] - 2 > b_Aeq->size[1]) {
              i12 = 1;
              i13 = 0;
            } else {
              i12 = b_Aeq->size[1] - 2;
              i13 = b_Aeq->size[1];
            }

            if (1 > b_Aeq->size[1] - 3) {
              loop_ub = 0;
            } else {
              loop_ub = b_Aeq->size[1] - 3;
            }

            if (1 > b_Aeq->size[1] - 3) {
              i = 0;
            } else {
              i = b_Aeq->size[1] - 3;
            }

            if (b_Aeq->size[1] - 2 > x0->size[0]) {
              i14 = 1;
              i15 = 0;
            } else {
              i14 = b_Aeq->size[1] - 2;
              i15 = x0->size[0];
            }

            if ((loop_ub == 1) || (i == 1)) {
              ar = d_Aeq->size[0] * d_Aeq->size[1];
              d_Aeq->size[0] = 3;
              d_Aeq->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)d_Aeq, ar, (int)sizeof
                                (double));
              for (ar = 0; ar < loop_ub; ar++) {
                for (ib = 0; ib < 3; ib++) {
                  d_Aeq->data[ib + d_Aeq->size[0] * ar] = b_Aeq->data[ib +
                    b_Aeq->size[0] * ar];
                }
              }

              ar = b_x0->size[0];
              b_x0->size[0] = i;
              emxEnsureCapacity((emxArray__common *)b_x0, ar, (int)sizeof(double));
              for (ar = 0; ar < i; ar++) {
                b_x0->data[ar] = x0->data[ar];
              }

              for (ar = 0; ar < 3; ar++) {
                y[ar] = 0.0;
                loop_ub = d_Aeq->size[1];
                for (ib = 0; ib < loop_ub; ib++) {
                  diff_phi = y[ar] + d_Aeq->data[ar + d_Aeq->size[0] * ib] *
                    b_x0->data[ib];
                  y[ar] = diff_phi;
                }
              }
            } else {
              for (i = 0; i < 3; i++) {
                y[i] = 0.0;
              }

              ar = -1;
              for (ib = 0; ib + 1 <= loop_ub; ib++) {
                if (x0->data[ib] != 0.0) {
                  ia = ar;
                  for (i = 0; i < 3; i++) {
                    ia++;
                    diff_phi = y[i] + x0->data[ib] * b_Aeq->data[ia % 3 +
                      b_Aeq->size[0] * (ia / 3)];
                    y[i] = diff_phi;
                  }
                }

                ar += 3;
              }
            }

            ar = c_Aeq->size[0] * c_Aeq->size[1];
            c_Aeq->size[0] = 3;
            c_Aeq->size[1] = (i13 - i12) + 1;
            emxEnsureCapacity((emxArray__common *)c_Aeq, ar, (int)sizeof(double));
            loop_ub = (i13 - i12) + 1;
            for (i13 = 0; i13 < loop_ub; i13++) {
              for (ar = 0; ar < 3; ar++) {
                c_Aeq->data[ar + c_Aeq->size[0] * i13] = b_Aeq->data[ar +
                  b_Aeq->size[0] * ((i12 + i13) - 1)];
              }
            }

            for (i12 = 0; i12 < 3; i12++) {
              c_beq[i12] = b_beq->data[i12] - y[i12];
            }

            mldivide(c_Aeq, c_beq, r10);
            i12 = r11->size[0] * r11->size[1];
            r11->size[0] = 1;
            r11->size[1] = (i15 - i14) + 1;
            emxEnsureCapacity((emxArray__common *)r11, i12, (int)sizeof(int));
            loop_ub = (i15 - i14) + 1;
            for (i12 = 0; i12 < loop_ub; i12++) {
              r11->data[r11->size[0] * i12] = (i14 + i12) - 1;
            }

            loop_ub = r11->size[0] * r11->size[1];
            for (i12 = 0; i12 < loop_ub; i12++) {
              x0->data[r11->data[i12]] = r10->data[i12];
            }

            /* does not work if there is tunnel!!!!!!!!! */
            i12 = solution->size[0];
            solution->size[0] = x0->size[0];
            emxEnsureCapacity((emxArray__common *)solution, i12, (int)sizeof
                              (double));
            loop_ub = x0->size[0];
            for (i12 = 0; i12 < loop_ub; i12++) {
              solution->data[i12] = x0->data[i12];
            }

            i12 = c_f->size[0];
            c_f->size[0] = b_f->size[0];
            emxEnsureCapacity((emxArray__common *)c_f, i12, (int)sizeof(double));
            loop_ub = b_f->size[0];
            for (i12 = 0; i12 < loop_ub; i12++) {
              c_f->data[i12] = -b_f->data[i12];
            }

            i12 = b_A->size[0] * b_A->size[1];
            b_A->size[0] = A->size[0];
            b_A->size[1] = A->size[1];
            emxEnsureCapacity((emxArray__common *)b_A, i12, (int)sizeof(double));
            loop_ub = A->size[0] * A->size[1];
            for (i12 = 0; i12 < loop_ub; i12++) {
              b_A->data[i12] = -A->data[i12];
            }

            i12 = b_b->size[0];
            b_b->size[0] = b->size[0];
            emxEnsureCapacity((emxArray__common *)b_b, i12, (int)sizeof(double));
            loop_ub = b->size[0];
            for (i12 = 0; i12 < loop_ub; i12++) {
              b_b->data[i12] = -b->data[i12];
            }

            *status = QPACM(Linv, c_f, b_Aeq, b_beq, b_A, b_b, solution);

            /*  x = quadprog(H,f,A,b,Aeq,beq); */
            break;
          }
        } else {
          i12 = solution->size[0];
          solution->size[0] = (int)numOfTunnel;
          emxEnsureCapacity((emxArray__common *)solution, i12, (int)sizeof
                            (double));
          loop_ub = (int)numOfTunnel;
          for (i12 = 0; i12 < loop_ub; i12++) {
            solution->data[i12] = 0.0;
          }

          *status = -2.0;
        }

        /*              solution = solution.*(abs(solution)>=FeasibilityTol);%!!!!!avoid nuemrical problem */
        /*  output */
        numOfTunnel = 0.0;
        numOfAzi = 0.0;
        numOfFpp = 0.0;
        for (i = 0; i < (int)N_enabled_thruster; i++) {
          alloc_out[i].label = thruster_data[i].label;
          alloc_out[i].enable = thruster_data[i].enable;
          switch ((int)thruster_data[i].type) {
           case 4:
            diff_phi = (numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0;
            alloc_out[i].Tx = solution->data[(int)(diff_phi + 1.0) - 1] *
              (thruster_data[i].Tmax * thruster_data[i].Tmax);
            alloc_out[i].Ty = solution->data[(int)(diff_phi + 2.0) - 1] *
              (thruster_data[i].Tmax * thruster_data[i].Tmax);
            alloc_out[i].T = alloc_out[i].Ty;

            /*                  alloc_out(i).phi = atan2(solution(numOfTunnel * 2 + numOfAzi * 2 + numOfFpp * 2 + 2), ... */
            /*                      solution(numOfTunnel * 2 + numOfAzi * 2 + numOfFpp * 2 + 1)); */
            alloc_out[i].phi = thruster_data[i].phi;
            alloc_out[i].phi = thruster_data[i].phi * (double)(fabs(alloc_out[i]
              .T) <= 1.0E-8) + alloc_out[i].phi * (double)(fabs(alloc_out[i].T) >
              1.0E-8);
            alloc_out[i].Tm = alloc_out[i].Ty * (thruster_data[i].x -
              thruster_data[i].x0) + alloc_out[i].Tx * (thruster_data[i].y0 -
              thruster_data[i].y);
            numOfTunnel++;
            break;

           case 3:
            diff_phi = (numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0;
            alloc_out[i].Tx = solution->data[(int)(diff_phi + 1.0) - 1] *
              (thruster_data[i].Tmax * thruster_data[i].Tmax);
            alloc_out[i].Ty = solution->data[(int)(diff_phi + 2.0) - 1] *
              (thruster_data[i].Tmax * thruster_data[i].Tmax);
            alloc_out[i].T = sqrt(alloc_out[i].Tx * alloc_out[i].Tx +
                                  alloc_out[i].Ty * alloc_out[i].Ty);
            alloc_out[i].phi = rt_atan2d_snf(alloc_out[i].Ty, alloc_out[i].Tx);
            alloc_out[i].phi = thruster_data[i].phi * (double)(fabs(alloc_out[i]
              .T) <= 1.0E-8) + alloc_out[i].phi * (double)(fabs(alloc_out[i].T) >
              1.0E-8);
            alloc_out[i].Tm = alloc_out[i].Ty * (thruster_data[i].x -
              thruster_data[i].x0) + alloc_out[i].Tx * (thruster_data[i].y0 -
              thruster_data[i].y);
            if ((1 + step == 1) && (fabs(thruster_data[i].phi_max -
                  thruster_data[i].phi_min) > 1.0E-8)) {
              /* means there is forbidden zone */
              if (angle_dist(thruster_data[i].phi_min, alloc_out[i].phi) >
                  angle_dist(thruster_data[i].phi_min, thruster_data[i].phi_max)
                  + 1.0E-8) {
                /*  if solution is inside fb zone */
                /*                              thruster_data(i).constr = 3;%means angle is fixed and  solution is inside fb zone */
                thruster_data[i].constr = 2.0;

                /* added v4.00!!!!!!!!!!!!!!! */
                if (on_border[i] == 0.0) {
                  /* esle use previous phi, phi_ and phiPlus */
                  /* %phi closer to min */
                  if (angle_dist(thruster_data[i].phi_max, alloc_out[i].phi) >
                      angle_dist(alloc_out[i].phi, thruster_data[i].phi_min)) {
                    thruster_data[i].phi_ = thruster_data[i].phi_min;
                    thruster_data[i].phiPlus = thruster_data[i].phi_min + 2.0 *
                      dt * thruster_data[i].dphi_max;
                    thruster_data[i].phi = thruster_data[i].phi_min + dt *
                      thruster_data[i].dphi_max;
                  } else {
                    /* phi closer to max */
                    thruster_data[i].phi_ = thruster_data[i].phi_max - 2.0 * dt *
                      thruster_data[i].dphi_max;
                    thruster_data[i].phiPlus = thruster_data[i].phi_max;
                    thruster_data[i].phi = thruster_data[i].phi_max - dt *
                      thruster_data[i].dphi_max;
                  }
                }

                thruster_data[i].T_ = thruster_data[i].Tmin;
                thruster_data[i].Tplus = thruster_data[i].Tmax;
                flag_run = 1;
              }

              if (on_border[i] == 1.0) {
                /*                              if i == 1 */
                /*                                  temp = phi_guider_old(i); */
                /*                                  save('phi_guider_old.txt','temp','-ASCII','-append'); */
                /*                              end */
                thruster_data[i].constr = 2.0;

                /* added!!!!! */
                if (fabs(thruster_data[i].phi - phi_guider_old[i]) > 1.0E-8) {
                  thruster_data[i].phi = phi_guider_old[i];
                  thruster_data[i].phi_ = thruster_data[i].phi;
                  thruster_data[i].phiPlus = thruster_data[i].phi;
                  flag_run = 1;
                } else {
                  thruster_data[i].phi_ = thruster_data[i].phi;
                  thruster_data[i].phiPlus = thruster_data[i].phi;
                }

                thruster_data[i].T_ = thruster_data[i].Tmin;
                thruster_data[i].Tplus = thruster_data[i].Tmax;
              }
            }

            numOfAzi++;
            break;

           case 2:
            diff_phi = (numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0;
            alloc_out[i].Tx = solution->data[(int)(diff_phi + 1.0) - 1] *
              (thruster_data[i].Tmax * thruster_data[i].Tmax);
            alloc_out[i].Ty = solution->data[(int)(diff_phi + 2.0) - 1] *
              (thruster_data[i].Tmax * thruster_data[i].Tmax);
            alloc_out[i].T = sqrt(alloc_out[i].Tx * alloc_out[i].Tx +
                                  alloc_out[i].Ty * alloc_out[i].Ty);
            alloc_out[i].phi = rt_atan2d_snf(alloc_out[i].Ty, alloc_out[i].Tx);
            alloc_out[i].phi = thruster_data[i].phi * (double)(fabs(alloc_out[i]
              .T) <= 1.0E-8) + alloc_out[i].phi * (double)(fabs(alloc_out[i].T) >
              1.0E-8);
            diff_phi = alloc_out[i].phi;
            b_sign(&diff_phi);
            alloc_out[i].phi = diff_phi * polyval(Fangle_Rangle, fabs
              (alloc_out[i].phi));
            alloc_out[i].Tm = alloc_out[i].Ty * (thruster_data[i].x -
              thruster_data[i].x0) + alloc_out[i].Tx * (thruster_data[i].y0 -
              thruster_data[i].y);
            numOfFpp++;
            break;
          }
        }

        if ((*status < 0.0) || (count_run > 50)) {
          /* avoid too many loops */
          exitg1 = true;
        }
      }

      if (1 + step == 1) {
        diff_phi = 0.0;
        for (i = 0; i < (int)N_enabled_thruster; i++) {
          /*          thruster_data(i).T = alloc_out(i).T; */
          if (thruster_data[i].type == 3.0) {
            thruster_data[i].phi = alloc_out[i].phi;
            diff_phi += b_mod(fabs(phi0[i] - thruster_data[i].phi),
                              3.1415926535897931);

            /* check!!!!!!!!! */
            phi0[i] = thruster_data[i].phi;
          }
        }

        if ((diff_phi < 0.01) || (counter > 10.0)) {
          /* this is the last loop in step one */
          run = false;

          /* break the while loop */
          use_base = false;
          for (i = 0; i < (int)N_enabled_thruster; i++) {
            if ((thruster_data[i].type == 3.0) && (alloc_out[i].T <
                 thruster_data[i].base[0] / 2.0)) {
              /* smaller than half base T */
              use_base = true;
            }
          }

          /*                  use_base = true;%use base will increase error when base leads phi in fbzone */
          if (use_base) {
            numOfTunnel = 0.0;
            for (i = 0; i < (int)N_enabled_thruster; i++) {
              if (thruster_data[i].type == 3.0) {
                alloc_out_T[i] = alloc_out[i].T;
                baseThrust(&alloc_out_T[i], alloc_out[i].phi, thruster_data[i].
                           base[0], thruster_data[i].base[1], thruster_data[i].
                           phi_min, thruster_data[i].phi_max, &alloc_out_phi[i],
                           &diff_phi);
                numOfTunnel += diff_phi;
              }
            }

            /*  the drawback is that when in step */
            /*  1 check_base_infb == 1, the result will */
            /*  not be accurate because phi */
            /*  in fbzone from step 1 will guide phi in step 2 to the */
            /*  fbzone border. If phi in step 2 is guided to fbzone, */
            /*  the result is accurate. However it is not allowed */
            if (numOfTunnel == 0.0) {
              /* only use base when it does not bring phi to fbzone */
              for (i = 0; i < (int)N_enabled_thruster; i++) {
                if (thruster_data[i].type == 3.0) {
                  alloc_out[i].T = alloc_out_T[i];
                  alloc_out[i].phi = alloc_out_phi[i];
                }
              }
            }
          }

          for (i = 0; i < (int)N_enabled_thruster; i++) {
            if (thruster_data[i].type == 3.0) {
              angle_point_to(alloc_out[i].phi, thruster_data[i].phi_reserve,
                             thruster_data[i].dphi_max, dt, thruster_data[i].
                             phi_min, thruster_data[i].phi_max, &phi_0[i],
                             &phiPlus0[i], &inside_fz[i]);
              on_border[i] = inside_fz[i];
              if (inside_fz[i] == 0.0) {
                phi_guider_old[i] = alloc_out[i].phi;
              }

              if (inside_fz[i] == 0.0) {
                /* phi is outside forbidden zone */
                angleMaxMin(thruster_data[i].phi_min, thruster_data[i].phi_max,
                            thruster_data[i].phi_reserve, phi_0[i], phiPlus0[i],
                            &thruster_data[i].phi_, &thruster_data[i].phiPlus);
                numOfTunnel = thruster_data[i].T_reserve - dt * thruster_data[i]
                  .dTmax;
                if ((numOfTunnel >= thruster_data[i].Tmin) || rtIsNaN
                    (thruster_data[i].Tmin)) {
                  thruster_data[i].T_ = numOfTunnel;
                } else {
                  thruster_data[i].T_ = thruster_data[i].Tmin;
                }

                numOfTunnel = thruster_data[i].T_reserve + dt * thruster_data[i]
                  .dTmax;
                if ((numOfTunnel <= thruster_data[i].Tmax) || rtIsNaN
                    (thruster_data[i].Tmax)) {
                  thruster_data[i].Tplus = numOfTunnel;
                } else {
                  thruster_data[i].Tplus = thruster_data[i].Tmax;
                }
              } else {
                thruster_data[i].phi_ = phi_0[i];
                thruster_data[i].phiPlus = phiPlus0[i];
                thruster_data[i].T_ = thruster_data[i].Tmin;
                thruster_data[i].Tplus = thruster_data[i].Tmin + 10.0;
              }

              thruster_data[i].phi = thruster_data[i].phi_ + angle_dist
                (thruster_data[i].phi_, thruster_data[i].phiPlus) / 2.0;
              thruster_data[i].T = thruster_data[i].T_reserve;

              /*                          if abs(thruster_data(i).phi_max - thruster_data(i).phi_min)> FeasibilityTol%means there is forbidden zone */
              /*                           */
              /*                             if min([angle_dist(alloc_out(i).phi,thruster_data(i).phi_min) */
              /*                                     angle_dist(thruster_data(i).phi_min,alloc_out(i).phi) */
              /*                                     angle_dist(alloc_out(i).phi,thruster_data(i).phi_max) */
              /*                                     angle_dist(thruster_data(i).phi_max,alloc_out(i).phi)]) ... */
              /*                                     <= FeasibilityTol */
              /*                                 on_border(i) = 1; */
              /*                             else */
              /*                                 on_border(i) = 0; */
              /*                             end */
              /*                              */
              /*                           */
              /*                          end */
            } else {
              thruster_data[i].phi = thruster_data[i].phi_reserve;
              thruster_data[i].T = thruster_data[i].T_reserve;
            }
          }
        }
      } else {
        /* if step == 2 */
        for (i = 0; i < (int)N_enabled_thruster; i++) {
          thruster_data[i].phi_reserve = alloc_out[i].phi;

          /* previous phi */
          thruster_data[i].T_reserve = alloc_out[i].T;

          /* previous T */
          thruster_data[i].phi = alloc_out[i].phi;
          thruster_data[i].T = alloc_out[i].T;
        }

        run = false;

        /* break the while loop */
      }
    }
  }

  emxFree_real_T(&b_x0);
  emxFree_real_T(&d_Aeq);
  emxFree_real_T(&c_Aeq);
  emxFree_real_T(&b_b);
  emxFree_real_T(&b_A);
  emxFree_real_T(&c_f);
  emxFree_boolean_T(&r12);
  emxFree_int32_T(&r11);
  emxFree_real_T(&r10);
  emxFree_boolean_T(&unusedU3);
  emxFree_real_T(&b_solution);
  emxFree_real_T(&L);
  emxFree_real_T(&b_beq);
  emxFree_real_T(&b_Aeq);
  emxFree_real_T(&b);
  emxFree_real_T(&A);
  emxFree_real_T(&b_f);
  emxFree_real_T(&Linv);
  emxFree_real_T(&x0);

  /*  status */
  /*  chck_sol = [[thruster_data(:).phi_min]' [thruster_data(:).phi_]' [alloc_out.phi]' [thruster_data.phiPlus]' [thruster_data(:).phi_max]']; */
  /*  check_angle(chck_sol(:,1),chck_sol(:,2),chck_sol(:,3),chck_sol(:,4),chck_sol(:,5)) */
  /*  solution */
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_init(void)
{
  memset(&on_border[0], 0, sizeof(double) << 3);
  memset(&phi_guider_old[0], 0, sizeof(double) << 3);
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
