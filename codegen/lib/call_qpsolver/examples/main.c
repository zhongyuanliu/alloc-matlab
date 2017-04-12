/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 01-Aug-2016 09:19:40
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
#include "call_qpsolver.h"
#include "main.h"
#include "call_qpsolver_terminate.h"
#include "call_qpsolver_emxAPI.h"
#include "call_qpsolver_initialize.h"

/* Function Declarations */
static double argInit_real_T(void);
static struct0_T argInit_struct0_T(void);
static void main_call_qpsolver(void);

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : struct0_T
 */
static struct0_T argInit_struct0_T(void)
{
  struct0_T result;

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
static void main_call_qpsolver(void)
{
  emxArray_real_T *solution;
  struct0_T r7;
  struct1_T alloc_out;
  double status;
  emxInitArray_real_T(&solution, 1);

  /* Initialize function 'call_qpsolver' input arguments. */
  /* Initialize function input argument 'T_r'. */
  /* Call the entry-point 'call_qpsolver'. */
  r7 = argInit_struct0_T();
  call_qpsolver(&r7, argInit_real_T(), argInit_real_T(), solution, &alloc_out,
                &status);
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
  call_qpsolver_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_call_qpsolver();

  /* Terminate the application.
     You do not need to do this more than one time. */
  call_qpsolver_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
