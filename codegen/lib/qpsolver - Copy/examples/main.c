/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 19-Aug-2016 16:37:49
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "qpsolver.h"
#include "main.h"
#include "qpsolver_terminate.h"
#include "qpsolver_emxAPI.h"
#include "qpsolver_initialize.h"

/* Function Declarations */
static void argInit_1x2_real_T(double result[2]);
static void argInit_1x8_struct0_T(struct0_T result[8]);
static void argInit_3x1_real_T(double result[3]);
static void argInit_5x9_real_T(double result[45]);
static double argInit_real_T(void);
static void argInit_struct0_T(struct0_T *result);
static struct1_T argInit_struct1_T(void);
static void main_qpsolver(void);

/* Function Definitions */

/*
 * Arguments    : double result[2]
 * Return Type  : void
 */
static void argInit_1x2_real_T(double result[2])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 2; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

/*
 * Arguments    : struct0_T result[8]
 * Return Type  : void
 */
static void argInit_1x8_struct0_T(struct0_T result[8])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 8; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    argInit_struct0_T(&result[idx1]);
  }
}

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[45]
 * Return Type  : void
 */
static void argInit_5x9_real_T(double result[45])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 5; idx0++) {
    for (idx1 = 0; idx1 < 9; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 5 * idx1] = argInit_real_T();
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : struct0_T *result
 * Return Type  : void
 */
static void argInit_struct0_T(struct0_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result->label = argInit_real_T();
  result->type = argInit_real_T();
  result->enable = argInit_real_T();
  result->constr_angle = argInit_real_T();
  result->dt = argInit_real_T();
  result->x0 = argInit_real_T();
  result->y0 = argInit_real_T();
  result->x = argInit_real_T();
  result->y = argInit_real_T();
  result->T = argInit_real_T();
  result->phi = argInit_real_T();
  result->Tmax = argInit_real_T();
  result->Tmin = argInit_real_T();
  result->dTmax = argInit_real_T();
  result->dphi_max = argInit_real_T();
  result->Tplus = argInit_real_T();
  result->T_ = argInit_real_T();
  result->phiPlus = argInit_real_T();
  result->phi_ = argInit_real_T();
  result->phi_min = argInit_real_T();
  result->phi_max = argInit_real_T();
  argInit_1x2_real_T(result->weight);
  argInit_3x1_real_T(result->weight_s);
}

/*
 * Arguments    : void
 * Return Type  : struct1_T
 */
static struct1_T argInit_struct1_T(void)
{
  struct1_T result;

  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result.Tx = argInit_real_T();
  result.Ty = argInit_real_T();
  result.Tm = argInit_real_T();
  return result;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_qpsolver(void)
{
  emxArray_real_T *solution;
  struct0_T rv0[8];
  double dv0[45];
  struct1_T r3;
  struct2_T alloc_out[8];
  double status;
  emxInitArray_real_T(&solution, 1);

  /* Initialize function 'qpsolver' input arguments. */
  /* Initialize function input argument 'thruster_data'. */
  /* Initialize function input argument 'T_r'. */
  /* Initialize function input argument 'rudder_table0'. */
  /* Call the entry-point 'qpsolver'. */
  argInit_1x8_struct0_T(rv0);
  argInit_5x9_real_T(dv0);
  r3 = argInit_struct1_T();
  qpsolver(rv0, &r3, argInit_real_T(), dv0, argInit_real_T(), argInit_real_T(),
           solution, alloc_out, &status);
  emxDestroyArray_real_T(solution);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  qpsolver_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_qpsolver();

  /* Terminate the application.
     You do not need to do this more than one time. */
  qpsolver_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
