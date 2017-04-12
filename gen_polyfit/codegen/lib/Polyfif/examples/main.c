/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 16-Aug-2016 09:21:48
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
#include "Polyfif.h"
#include "main.h"
#include "Polyfif_terminate.h"
#include "Polyfif_emxAPI.h"
#include "Polyfif_initialize.h"

/* Function Declarations */
static void argInit_1x9_real_T(double result[9]);
static double argInit_real_T(void);
static void main_Polyfif(void);

/* Function Definitions */

/*
 * Arguments    : double result[9]
 * Return Type  : void
 */
static void argInit_1x9_real_T(double result[9])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 9; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
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
 * Arguments    : void
 * Return Type  : void
 */
static void main_Polyfif(void)
{
  emxArray_real_T *r;
  double dv0[9];
  double dv1[9];
  double f;
  emxInitArray_real_T(&r, 2);

  /* Initialize function 'Polyfif' input arguments. */
  /* Initialize function input argument 'x'. */
  /* Initialize function input argument 'y'. */
  /* Call the entry-point 'Polyfif'. */
  argInit_1x9_real_T(dv0);
  argInit_1x9_real_T(dv1);
  Polyfif(dv0, dv1, argInit_real_T(), argInit_real_T(), r, &f);
  emxDestroyArray_real_T(r);
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
  Polyfif_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_Polyfif();

  /* Terminate the application.
     You do not need to do this more than one time. */
  Polyfif_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
