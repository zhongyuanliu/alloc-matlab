/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: Polyfif_initialize.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 16-Aug-2016 09:21:48
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Polyfif.h"
#include "Polyfif_initialize.h"

/* Function Definitions */

/*
 * Arguments    : void
 * Return Type  : void
 */
void Polyfif_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

/*
 * File trailer for Polyfif_initialize.c
 *
 * [EOF]
 */
