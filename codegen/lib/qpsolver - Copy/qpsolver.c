/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: qpsolver.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "qpsolver_emxutil.h"
#include "polyval.h"
#include "abs.h"
#include "QPACM.h"
#include "mldivide.h"
#include "Mpcqpsolver.h"
#include "xtrsm.h"
#include "colon.h"
#include "xzgetrf.h"
#include "Chol_fc.h"
#include "pre_qp.h"
#include "structConstructorHelper.h"
#include "polyfit.h"

/* Variable Definitions */
static emxArray_real_T *Rangle_T;
static emxArray_real_T *Rangle_Fangle;
static emxArray_real_T *Fangle_Rangle;
static boolean_T rudder_dat_init_not_empty;

/* Function Declarations */
static double rt_atan2d_snf(double u0, double u1);

/* Function Definitions */

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
 *  method == 2 primal active set
 * Arguments    : struct0_T thruster_data[8]
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
void qpsolver(struct0_T thruster_data[8], const struct1_T *T_r, double
              N_enabled_thruster, const double rudder_table0[45], double
              no_azi_angle_constr, double method, emxArray_real_T *solution,
              struct2_T alloc_out[8], double *status)
{
  double rudder_table[45];
  int i0;
  int n;
  emxArray_real_T *expl_temp;
  emxArray_real_T *b_expl_temp;
  emxArray_real_T *c_expl_temp;
  int flag_run;
  signed char mark8[8];
  int pipk;
  double d_expl_temp;
  double e_expl_temp;
  int loop_ub;
  emxArray_real_T *x0;
  emxArray_real_T *Linv;
  emxArray_real_T *f;
  emxArray_real_T *A;
  emxArray_real_T *b;
  emxArray_real_T *Aeq;
  emxArray_real_T *beq;
  emxArray_real_T *L;
  emxArray_real_T *b_solution;
  emxArray_boolean_T *lambda;
  emxArray_real_T *r0;
  emxArray_int32_T *p;
  emxArray_real_T *x;
  emxArray_int32_T *ipiv;
  emxArray_boolean_T *r1;
  emxArray_real_T *b_f;
  emxArray_real_T *b_A;
  emxArray_real_T *b_b;
  emxArray_real_T *b_Aeq;
  emxArray_real_T *c_Aeq;
  emxArray_real_T *b_x0;
  double Nvar;
  double b_p;
  unsigned int unnamed_idx_0;
  int i1;
  int i2;
  double numOfTunnel;
  double numOfAzi;
  double y[3];
  double numOfFpp;
  int ia;
  int c;
  int j;
  double b_beq[3];
  double k;
  double b_k;
  for (i0 = 0; i0 < 5; i0++) {
    for (n = 0; n < 9; n++) {
      rudder_table[n + 9 * i0] = rudder_table0[i0 + 5 * n];
    }
  }

  /* !!!!!!!!!!!!!!! */
  /*  no_azi_angle_constr = 1; */
  /* % rudder table */
  if (!rudder_dat_init_not_empty) {
    emxInit_real_T(&expl_temp, 2);
    emxInit_real_T(&b_expl_temp, 2);
    emxInit_real_T(&c_expl_temp, 2);

    /* assume this is only one type of rudder */
    rudder_dat_init_not_empty = true;
    polyfit(*(double (*)[9])&rudder_table[0], *(double (*)[9])&rudder_table[36],
            no_azi_angle_constr + 2.0, Rangle_T, expl_temp, &d_expl_temp,
            &e_expl_temp);

    /* rudder angle vs T */
    polyfit(*(double (*)[9])&rudder_table[0], *(double (*)[9])&rudder_table[27],
            no_azi_angle_constr + 2.0, Rangle_Fangle, b_expl_temp, &d_expl_temp,
            &e_expl_temp);

    /* rudder angle vs F angle */
    polyfit(*(double (*)[9])&rudder_table[27], *(double (*)[9])&rudder_table[0],
            no_azi_angle_constr + 2.0, Fangle_Rangle, c_expl_temp, &d_expl_temp,
            &e_expl_temp);

    /* F angle vs rudder angle */
    emxFree_real_T(&c_expl_temp);
    emxFree_real_T(&b_expl_temp);
    emxFree_real_T(&expl_temp);
  }

  /*  coder.varsize('alloc_out');%!!!!!!!!!very important??? */
  structConstructorHelper(alloc_out);
  flag_run = 1;
  for (pipk = 0; pipk < 8; pipk++) {
    mark8[pipk] = 0;
  }

  i0 = solution->size[0];
  solution->size[0] = (int)(N_enabled_thruster * 2.0 + 3.0);
  emxEnsureCapacity((emxArray__common *)solution, i0, (int)sizeof(double));
  loop_ub = (int)(N_enabled_thruster * 2.0 + 3.0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    solution->data[i0] = 0.0;
  }

  *status = 0.0;
  emxInit_real_T1(&x0, 1);
  emxInit_real_T(&Linv, 2);
  emxInit_real_T1(&f, 1);
  emxInit_real_T(&A, 2);
  emxInit_real_T1(&b, 1);
  emxInit_real_T(&Aeq, 2);
  emxInit_real_T1(&beq, 1);
  emxInit_real_T(&L, 2);
  emxInit_real_T1(&b_solution, 1);
  emxInit_boolean_T(&lambda, 1);
  emxInit_real_T1(&r0, 1);
  emxInit_int32_T(&p, 2);
  emxInit_real_T(&x, 2);
  emxInit_int32_T(&ipiv, 2);
  emxInit_boolean_T(&r1, 1);
  emxInit_real_T1(&b_f, 1);
  emxInit_real_T(&b_A, 2);
  emxInit_real_T1(&b_b, 1);
  emxInit_real_T(&b_Aeq, 2);
  emxInit_real_T(&c_Aeq, 2);
  emxInit_real_T1(&b_x0, 1);
  while (flag_run != 0) {
    flag_run = 0;
    pre_qp(thruster_data, T_r->Tx, T_r->Ty, T_r->Tm, N_enabled_thruster,
           rudder_table, Rangle_T, Rangle_Fangle, no_azi_angle_constr, &Nvar,
           Linv, f, A, b, Aeq, beq, x0);
    Chol_fc(Linv, L, &b_p);

    /*  before chol */
    if (b_p == 0.0) {
      switch ((int)method) {
       case 1:
        if ((L->size[0] == 0) || (L->size[1] == 0)) {
          i0 = Linv->size[0] * Linv->size[1];
          Linv->size[0] = L->size[0];
          Linv->size[1] = L->size[1];
          emxEnsureCapacity((emxArray__common *)Linv, i0, (int)sizeof(double));
          loop_ub = L->size[0] * L->size[1];
          for (i0 = 0; i0 < loop_ub; i0++) {
            Linv->data[i0] = L->data[i0];
          }
        } else {
          n = L->size[0];
          i0 = Linv->size[0] * Linv->size[1];
          Linv->size[0] = L->size[0];
          Linv->size[1] = L->size[1];
          emxEnsureCapacity((emxArray__common *)Linv, i0, (int)sizeof(double));
          loop_ub = L->size[0] * L->size[1];
          for (i0 = 0; i0 < loop_ub; i0++) {
            Linv->data[i0] = 0.0;
          }

          i0 = x->size[0] * x->size[1];
          x->size[0] = L->size[0];
          x->size[1] = L->size[1];
          emxEnsureCapacity((emxArray__common *)x, i0, (int)sizeof(double));
          loop_ub = L->size[0] * L->size[1];
          for (i0 = 0; i0 < loop_ub; i0++) {
            x->data[i0] = L->data[i0];
          }

          xzgetrf(L->size[0], L->size[0], x, L->size[0], ipiv, &pipk);
          eml_signed_integer_colon(L->size[0], p);
          for (ia = 0; ia < ipiv->size[1]; ia++) {
            if (ipiv->data[ia] > 1 + ia) {
              pipk = p->data[ipiv->data[ia] - 1];
              p->data[ipiv->data[ia] - 1] = p->data[ia];
              p->data[ia] = pipk;
            }
          }

          for (ia = 0; ia + 1 <= n; ia++) {
            c = p->data[ia] - 1;
            Linv->data[ia + Linv->size[0] * (p->data[ia] - 1)] = 1.0;
            for (j = ia; j + 1 <= n; j++) {
              if (Linv->data[j + Linv->size[0] * c] != 0.0) {
                for (pipk = j + 1; pipk + 1 <= n; pipk++) {
                  Linv->data[pipk + Linv->size[0] * c] -= Linv->data[j +
                    Linv->size[0] * c] * x->data[pipk + x->size[0] * j];
                }
              }
            }
          }

          xtrsm(L->size[0], L->size[0], x, L->size[0], Linv, L->size[0]);
        }

        unnamed_idx_0 = (unsigned int)b->size[0];

        /*      if isempty(iA) */
        /*          iA = iA0; */
        /*      end   */
        i0 = r1->size[0];
        r1->size[0] = (int)unnamed_idx_0;
        emxEnsureCapacity((emxArray__common *)r1, i0, (int)sizeof(boolean_T));
        loop_ub = (int)unnamed_idx_0;
        for (i0 = 0; i0 < loop_ub; i0++) {
          r1->data[i0] = false;
        }

        Mpcqpsolver(Linv, f, A, b, Aeq, beq, r1, b_solution, status, lambda);
        i0 = solution->size[0];
        solution->size[0] = b_solution->size[0];
        emxEnsureCapacity((emxArray__common *)solution, i0, (int)sizeof(double));
        loop_ub = b_solution->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          solution->data[i0] = b_solution->data[i0];
        }

        /* faster to use iA instead of iA0 */
        /*      [solution1,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq) */
        break;

       case 2:
        /*      x0=[x0;0;0;0]; */
        /*      ????????????????ux,uy,?????????s1?s2?s3??? */
        if (Aeq->size[1] - 2 > Aeq->size[1]) {
          i0 = 1;
          n = 0;
        } else {
          i0 = Aeq->size[1] - 2;
          n = Aeq->size[1];
        }

        if (1 > Aeq->size[1] - 3) {
          loop_ub = 0;
        } else {
          loop_ub = Aeq->size[1] - 3;
        }

        if (1 > Aeq->size[1] - 3) {
          pipk = 0;
        } else {
          pipk = Aeq->size[1] - 3;
        }

        if (Aeq->size[1] - 2 > x0->size[0]) {
          i1 = 1;
          i2 = 0;
        } else {
          i1 = Aeq->size[1] - 2;
          i2 = x0->size[0];
        }

        if ((loop_ub == 1) || (pipk == 1)) {
          c = c_Aeq->size[0] * c_Aeq->size[1];
          c_Aeq->size[0] = 3;
          c_Aeq->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)c_Aeq, c, (int)sizeof(double));
          for (c = 0; c < loop_ub; c++) {
            for (j = 0; j < 3; j++) {
              c_Aeq->data[j + c_Aeq->size[0] * c] = Aeq->data[j + Aeq->size[0] *
                c];
            }
          }

          c = b_x0->size[0];
          b_x0->size[0] = pipk;
          emxEnsureCapacity((emxArray__common *)b_x0, c, (int)sizeof(double));
          for (c = 0; c < pipk; c++) {
            b_x0->data[c] = x0->data[c];
          }

          for (c = 0; c < 3; c++) {
            y[c] = 0.0;
            loop_ub = c_Aeq->size[1];
            for (j = 0; j < loop_ub; j++) {
              d_expl_temp = y[c] + c_Aeq->data[c + c_Aeq->size[0] * j] *
                b_x0->data[j];
              y[c] = d_expl_temp;
            }
          }
        } else {
          for (pipk = 0; pipk < 3; pipk++) {
            y[pipk] = 0.0;
          }

          c = -1;
          for (j = 0; j + 1 <= loop_ub; j++) {
            if (x0->data[j] != 0.0) {
              ia = c;
              for (pipk = 0; pipk < 3; pipk++) {
                ia++;
                d_expl_temp = y[pipk] + x0->data[j] * Aeq->data[ia % 3 +
                  Aeq->size[0] * (ia / 3)];
                y[pipk] = d_expl_temp;
              }
            }

            c += 3;
          }
        }

        c = b_Aeq->size[0] * b_Aeq->size[1];
        b_Aeq->size[0] = 3;
        b_Aeq->size[1] = (n - i0) + 1;
        emxEnsureCapacity((emxArray__common *)b_Aeq, c, (int)sizeof(double));
        loop_ub = (n - i0) + 1;
        for (n = 0; n < loop_ub; n++) {
          for (c = 0; c < 3; c++) {
            b_Aeq->data[c + b_Aeq->size[0] * n] = Aeq->data[c + Aeq->size[0] *
              ((i0 + n) - 1)];
          }
        }

        for (i0 = 0; i0 < 3; i0++) {
          b_beq[i0] = beq->data[i0] - y[i0];
        }

        mldivide(b_Aeq, b_beq, r0);
        i0 = ipiv->size[0] * ipiv->size[1];
        ipiv->size[0] = 1;
        ipiv->size[1] = (i2 - i1) + 1;
        emxEnsureCapacity((emxArray__common *)ipiv, i0, (int)sizeof(int));
        loop_ub = (i2 - i1) + 1;
        for (i0 = 0; i0 < loop_ub; i0++) {
          ipiv->data[ipiv->size[0] * i0] = (i1 + i0) - 1;
        }

        loop_ub = ipiv->size[0] * ipiv->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          x0->data[ipiv->data[i0]] = r0->data[i0];
        }

        /* does not work if there is tunnel!!!!!!!!! */
        i0 = solution->size[0];
        solution->size[0] = x0->size[0];
        emxEnsureCapacity((emxArray__common *)solution, i0, (int)sizeof(double));
        loop_ub = x0->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          solution->data[i0] = x0->data[i0];
        }

        i0 = b_f->size[0];
        b_f->size[0] = f->size[0];
        emxEnsureCapacity((emxArray__common *)b_f, i0, (int)sizeof(double));
        loop_ub = f->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_f->data[i0] = -f->data[i0];
        }

        i0 = b_A->size[0] * b_A->size[1];
        b_A->size[0] = A->size[0];
        b_A->size[1] = A->size[1];
        emxEnsureCapacity((emxArray__common *)b_A, i0, (int)sizeof(double));
        loop_ub = A->size[0] * A->size[1];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_A->data[i0] = -A->data[i0];
        }

        i0 = b_b->size[0];
        b_b->size[0] = b->size[0];
        emxEnsureCapacity((emxArray__common *)b_b, i0, (int)sizeof(double));
        loop_ub = b->size[0];
        for (i0 = 0; i0 < loop_ub; i0++) {
          b_b->data[i0] = -b->data[i0];
        }

        *status = QPACM(Linv, b_f, Aeq, beq, b_A, b_b, solution);

        /*  x = quadprog(H,f,A,b,Aeq,beq); */
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

    c_abs(solution, r0);
    i0 = lambda->size[0];
    lambda->size[0] = r0->size[0];
    emxEnsureCapacity((emxArray__common *)lambda, i0, (int)sizeof(boolean_T));
    loop_ub = r0->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      lambda->data[i0] = (r0->data[i0] >= 1.0E-8);
    }

    i0 = solution->size[0];
    emxEnsureCapacity((emxArray__common *)solution, i0, (int)sizeof(double));
    loop_ub = solution->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      solution->data[i0] *= (double)lambda->data[i0];
    }

    /* !!!!!avoid nuemrical problem */
    /*  output */
    numOfTunnel = 0.0;
    numOfAzi = 0.0;
    numOfFpp = 0.0;
    for (pipk = 0; pipk < (int)N_enabled_thruster; pipk++) {
      alloc_out[pipk].label = thruster_data[pipk].label;
      alloc_out[pipk].enable = thruster_data[pipk].enable;
      switch ((int)thruster_data[pipk].type) {
       case 4:
        alloc_out[pipk].Tx = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi
          * 2.0) + numOfFpp * 2.0) + 1.0) - 1];
        alloc_out[pipk].Ty = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi
          * 2.0) + numOfFpp * 2.0) + 2.0) - 1];
        alloc_out[pipk].T = alloc_out[pipk].Ty;
        alloc_out[pipk].phi = rt_atan2d_snf(solution->data[(int)(((numOfTunnel *
          2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 2.0) - 1], solution->data
          [(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0)
          - 1]);
        alloc_out[pipk].phi = thruster_data[pipk].phi * (double)(fabs
          (alloc_out[pipk].T) <= 1.0E-8) + alloc_out[pipk].phi * (double)(fabs
          (alloc_out[pipk].T) > 1.0E-8);
        alloc_out[pipk].Tm = alloc_out[pipk].Ty * (thruster_data[pipk].x -
          thruster_data[pipk].x0) + alloc_out[pipk].Tx * (thruster_data[pipk].y0
          - thruster_data[pipk].y);
        numOfTunnel++;
        break;

       case 3:
        alloc_out[pipk].Tx = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi
          * 2.0) + numOfFpp * 2.0) + 1.0) - 1];
        alloc_out[pipk].Ty = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi
          * 2.0) + numOfFpp * 2.0) + 2.0) - 1];
        alloc_out[pipk].T = sqrt(alloc_out[pipk].Tx * alloc_out[pipk].Tx +
          alloc_out[pipk].Ty * alloc_out[pipk].Ty);
        alloc_out[pipk].phi = rt_atan2d_snf(solution->data[(int)(((numOfTunnel *
          2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 2.0) - 1], solution->data
          [(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0)
          - 1]);
        alloc_out[pipk].phi = thruster_data[pipk].phi * (double)(fabs
          (alloc_out[pipk].T) <= 1.0E-8) + alloc_out[pipk].phi * (double)(fabs
          (alloc_out[pipk].T) > 1.0E-8);
        alloc_out[pipk].Tm = alloc_out[pipk].Ty * (thruster_data[pipk].x -
          thruster_data[pipk].x0) + alloc_out[pipk].Tx * (thruster_data[pipk].y0
          - thruster_data[pipk].y);
        if (fabs(thruster_data[pipk].phi_max - thruster_data[pipk].phi_min) >
            1.0E-8) {
          /* means the thruster has physical constraint */
          /* % angle distance */
          /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
          k = alloc_out[pipk].phi - thruster_data[pipk].phi_min;
          k *= (double)(fabs(k) >= 1.0E-8);

          /* set small value to zero, avoid numerical problem */
          /* % angle distance */
          /*  up x axis, right y axis, angle starts from x axis, -pi to pi, clockwise */
          b_k = thruster_data[pipk].phi_max - thruster_data[pipk].phi_min;
          b_k *= (double)(fabs(b_k) >= 1.0E-8);

          /* set small value to zero, avoid numerical problem */
          if (k * (double)(k >= 0.0) + (6.2831853071795862 + k) * (double)(k <
               0.0) > (b_k * (double)(b_k >= 0.0) + (6.2831853071795862 + b_k) *
                       (double)(b_k < 0.0)) + 1.0E-8) {
            /*  if solution is inside fb zone */
            thruster_data[pipk].constr_angle = 1.0;
            mark8[pipk] = 1;

            /* when mark8(i)==1, the thruster will be constrained from now on in while loop */
            flag_run = 1;
          }

          if (((fabs(alloc_out[pipk].phi - thruster_data[pipk].phi_min) <=
                1.0E-8) || (fabs(alloc_out[pipk].phi - thruster_data[pipk].
                                 phi_max) <= 1.0E-8)) && (fabs(alloc_out[pipk].T
                - thruster_data[pipk].Tmin) <= 1.0E-8) && (mark8[pipk] == 0)) {
            thruster_data[pipk].constr_angle = 0.0;
            flag_run = 1;

            /*  remove constr and run again */
          }
        }

        numOfAzi++;
        break;

       case 2:
        alloc_out[pipk].Tx = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi
          * 2.0) + numOfFpp * 2.0) + 1.0) - 1];
        alloc_out[pipk].Ty = solution->data[(int)(((numOfTunnel * 2.0 + numOfAzi
          * 2.0) + numOfFpp * 2.0) + 2.0) - 1];
        alloc_out[pipk].T = sqrt(alloc_out[pipk].Tx * alloc_out[pipk].Tx +
          alloc_out[pipk].Ty * alloc_out[pipk].Ty);
        alloc_out[pipk].phi = rt_atan2d_snf(solution->data[(int)(((numOfTunnel *
          2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 2.0) - 1], solution->data
          [(int)(((numOfTunnel * 2.0 + numOfAzi * 2.0) + numOfFpp * 2.0) + 1.0)
          - 1]);
        alloc_out[pipk].phi = thruster_data[pipk].phi * (double)(fabs
          (alloc_out[pipk].T) <= 1.0E-8) + alloc_out[pipk].phi * (double)(fabs
          (alloc_out[pipk].T) > 1.0E-8);
        if (alloc_out[pipk].phi < 0.0) {
          d_expl_temp = -1.0;
        } else if (alloc_out[pipk].phi > 0.0) {
          d_expl_temp = 1.0;
        } else if (alloc_out[pipk].phi == 0.0) {
          d_expl_temp = 0.0;
        } else {
          d_expl_temp = alloc_out[pipk].phi;
        }

        alloc_out[pipk].phi = d_expl_temp * polyval(Fangle_Rangle, fabs
          (alloc_out[pipk].phi));
        alloc_out[pipk].Tm = alloc_out[pipk].Ty * (thruster_data[pipk].x -
          thruster_data[pipk].x0) + alloc_out[pipk].Tx * (thruster_data[pipk].y0
          - thruster_data[pipk].y);
        numOfFpp++;
        break;
      }
    }
  }

  emxFree_real_T(&b_x0);
  emxFree_real_T(&c_Aeq);
  emxFree_real_T(&b_Aeq);
  emxFree_real_T(&b_b);
  emxFree_real_T(&b_A);
  emxFree_real_T(&b_f);
  emxFree_boolean_T(&r1);
  emxFree_int32_T(&ipiv);
  emxFree_real_T(&x);
  emxFree_int32_T(&p);
  emxFree_real_T(&r0);
  emxFree_boolean_T(&lambda);
  emxFree_real_T(&b_solution);
  emxFree_real_T(&L);
  emxFree_real_T(&beq);
  emxFree_real_T(&Aeq);
  emxFree_real_T(&b);
  emxFree_real_T(&A);
  emxFree_real_T(&f);
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
void qpsolver_free(void)
{
  emxFree_real_T(&Fangle_Rangle);
  emxFree_real_T(&Rangle_Fangle);
  emxFree_real_T(&Rangle_T);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void qpsolver_init(void)
{
  emxInit_real_T(&Fangle_Rangle, 2);
  emxInit_real_T(&Rangle_Fangle, 2);
  emxInit_real_T(&Rangle_T, 2);
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
