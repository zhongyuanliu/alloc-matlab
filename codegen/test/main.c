/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 02-Aug-2016 16:16:16
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
#include "Qpsolver.h"
//#include "main.h"

#include "qpsolver_terminate.h"
#include "qpsolver_emxAPI.h"
#include "qpsolver_initialize.h"
#include <stdio.h>
#include <windows.h>
#include "qpmn.h"
typedef alloation_out (*Dllfun)(thruster_data_in*);
extern void reald_file(struct0_T *thruster_in, double *rudder_table);
extern void FCN_ALLOCATION(thruster_data_in* pInput,alloation_out (*pOutput) );
/* Function Declarations */
static void argInit_1x2_real_T(double result[2]);
static void argInit_1x8_struct0_T(struct0_T result[8]);
static void argInit_3x1_real_T(double result[3]);
static void argInit_9x5_real_T(double result[45]);
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
    result[idx1] =1.0;
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
    result[idx0] = 100000.0;
  }
}

/*
 * Arguments    : double result[45]
 * Return Type  : void
 */
static void argInit_9x5_real_T(double result[45])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 9; idx0++) {
    for (idx1 = 0; idx1 < 5; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 9 * idx1] = argInit_real_T();
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
	int i = 0;
(result+i)->label = i;
  (result+i)->type = 3;
  (result+i)->enable = 1;
  (result+i)->constr_angle = 0;
  (result+i)->dt = 0.5;
  (result+i)->x0 = 0.0;
  (result+i)->y0 = 0.0;
  (result+i)->x = 3;
  (result+i)->y = 5;
  (result+i)->T = 4365;
  (result+i)->phi = 0.0;
  (result+i)->Tmax = 71500;
  (result+i)->Tmin = 0;
  (result+i)->dTmax = 7150.0;
  (result+i)->dphi_max = 0.5;
  (result+i)->Tplus = argInit_real_T();
  (result+i)->T_ = argInit_real_T();
  (result+i)->phiPlus = argInit_real_T();
  (result+i)->phi_ = argInit_real_T();
  (result+i)->phi_min = -0.1;
  (result+i)->phi_max = -0.2;
  argInit_1x2_real_T((result+i)->weight);
  argInit_3x1_real_T((result+i)->weight_s);
}
static void Init_struct0_T(struct0_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
	int i;
	for (i = 0; i<8; i++)
	{
  (result+i)->label = i + 1;
  (result+i)->type = 3;
  (result+i)->enable = 1;
  (result+i)->constr_angle = 0;
  (result+i)->dt = 0.5;
  (result+i)->x0 = 0.0;
  (result+i)->y0 = 0.0;
  (result+i)->x = 3;
  (result+i)->y = 5;
  (result+i)->T = 536.5 + i * 500;
  (result+i)->phi = 0.0;
  (result+i)->Tmax = 71500;
  (result+i)->Tmin = 0;
  (result+i)->dTmax = 7150.0;
  (result+i)->dphi_max = 0.5;
  (result+i)->Tplus = argInit_real_T();
  (result+i)->T_ = argInit_real_T();
  (result+i)->phiPlus = argInit_real_T();
  (result+i)->phi_ = argInit_real_T();
  (result+i)->phi_min = -0.1;
  (result+i)->phi_max = -0.2;
  argInit_1x2_real_T((result+i)->weight);
  argInit_3x1_real_T((result+i)->weight_s);
	}
  (result+0)->x = 3+3;
  (result+0)->y = 5;
  (result+1)->x = 3+3;
  (result+1)->y = -5;
  (result+2)->x = 10+3;
  (result+2)->y = 5;
  (result+3)->x = 10+3;
  (result+3)->y = -5;
  (result+4)->x = 7+2;
  (result+4)->y = -7;
  (result+5)->x = 7+2;
  (result+5)->y = 7;
  (result+6)->x = 12-2;
  (result+6)->y = -5;
  (result+7)->x = 12-2;
  (result+7)->y = 5;

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
  result.Tx = 5670.0;
  result.Ty = 19700.0;
  result.Tm = -20000.0;
  return result;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_qpsolver(void)
{

  struct0_T thruster_p[8];
  double dv1[45];
  double rudder_table[45] = 
{0.0,	0.980000000000000,	0.0,	0.0,	0.980000000000000,
0.0872664625997165,	0.965000000000000,	0.0800000000000000,	0.0827124154481730,	0.968310384122777,
0.174532925199433,	0.956000000000000,	0.163000000000000,	0.168878105792867,	0.969796370378854,
0.261799387799149,	0.920000000000000,	0.200000000000000,	0.214060683563822,	0.941488183675186,
0.349065850398866,	0.895000000000000,	0.287000000000000,	0.310310945764607,	0.939890419144700,
0.436332312998582,	0.840000000000000,	0.300000000000000,	0.343023940420703,	0.891964124839110,
0.523598775598299,	0.811000000000000,	0.341000000000000,	0.398026222516986,	0.879773834573409,
0.610865238198015,	0.770000000000000,	0.345000000000000,	0.421232744089294,	0.843756481456587,
0.698131700797732,	0.743000000000000,	0.347000000000000,	0.436921840812821,	0.820035365091043};
  struct1_T T_reqr;
  struct2_T alloc_out[8];
  int alloc_out_size[2];
  int i;
  int j;
  double status;
  double temp_i;
  double temp_T[8]={363.212240412854,
314.516432393549,
368.190183808438,
309.537637353717,
0.0,
336.374513061158,
310.599456662282,
0.0};
  double temp_phi[8] = {-0.0478344712504149,
0.0500000000000000,
0.0500000000000000,
-0.0499999999999998,
0.0,
-0.0499999999999998,
0.0500000000000000,
0.0};
  FILE *labelt;
    HINSTANCE   * hDLL = LoadLibrary("AsynchDLL.dll");
alloation_out* dllfcn;
thruster_data_in* pInput = (thruster_data_in*) calloc(1,sizeof(thruster_data_in));//可以当一维数组的指针例子
alloation_out* all_out = (alloation_out*) calloc(8,sizeof(alloation_out));


  /* Initialize function 'Qpsolver' input arguments. */
  /* Initialize function input argument 'thruster_data'. */
  /* Initialize function input argument 'T_r'. */
  /* Initialize function input argument 'rudder_table'. */
  /* Call the entry-point 'Qpsolver'. */
  Init_struct0_T(thruster_p);
  argInit_9x5_real_T(dv1);

reald_file(thruster_p, rudder_table);

for (i = 0; i < 16; i++)
{	j = i%8;
	temp_i = i / 8;
	(pInput )->dt = 0.5;
	(pInput )->enable = 1;
	(pInput )->label = j + 1;
	(pInput )->T = 3000.0;
	//(pInput )->T = temp_T[i];
	//(pInput )->phi = temp_phi[i];
	
	(pInput )->phi =0.0;
	(pInput )->Trx = 10000.0 * cos(temp_i / 5.0);
	(pInput )->Try = 10000.0 * sin(temp_i / 5.0);
	(pInput )->Trm = 0.0;
	if (j == 0) (pInput )->phi =pi / 2.0;
	if (i >= 8)
	{
		(pInput )->T = (all_out+j)->T;
		(pInput )->phi =(all_out+j)->phi;
	}
	FCN_ALLOCATION(pInput,&all_out[0]);
//alloc_out[j] = all_out;

}
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
