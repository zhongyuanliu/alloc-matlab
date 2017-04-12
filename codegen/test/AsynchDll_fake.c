#include "qpmn.h"
#include "stdafx.h"
#include "string.h"
#include "qpsolver.h"
#include "qpsolver_emxAPI.h"
#include "qpsolver_initialize.h"
#include <stdio.h>
struct0_T thruster_p[numofThusters];
//static struct0_T *thruster_p = (struct0_T*) calloc(numofThusters,sizeof(struct0_T));
//static struct2_T *alloc_out=(struct2_T*) calloc(numofThusters,sizeof(struct2_T));

// = {{}};//init????????????????????
//static thruster_data thruster_in[numofThusters];

//static struct0_T *thruster_in_p = (struct0_T*) calloc(numofThusters,sizeof(struct0_T));
short counter8[numofThusters]={0,0,0,0,0,0,0,0};
short sum_counter8 = 0;
short counter_out8[numofThusters]={0,0,0,0,0,0,0,0};
short sum_counter_out8 = 0;
short count_out = 0;
int N_enabled_thr = 0;//count number of enabled thrusters
double flag = 0.0;
double flag1 = 0.0;
double flag2 = 0.0;
double flag3 = 0.0;
double flag4 = 0.0;
double flag5 = 0.0;
double flag6 = 0.0;
double flag7 = 0.0;
double flag8 = 0.0;

short output_enable = 0;
short init_thr = 0;
short init_allc_out = 0;
short init_read = 0;
short label_t = 1;
double status = 0.0;
emxArray_real_T *solution;
struct0_T *thruster_in;
int alloc_out_size[2];
double rudder_table[45];
/*
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
*/

extern	void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void reald_file(struct0_T *thruster_in, double *rudder_table);
static struct2_T alloc_out[numofThusters];
void FCN_ALLOCATION(thruster_data_in* pInput,alloation_out (*pOutput) )
{
	int i;
	
	//alloation_out* pOutput = (alloation_out*) calloc(1,sizeof(alloation_out));//可以当一维数组的指针例子
	
	struct1_T T_reqr;
	
//	 FILE* test_file;

	FILE* labelt;
    char* fpfile = "test_input.txt";
	char* solout = "solout.txt";

	if (init_read == 0)//read thruster config data from file
	{
		qpsolver_initialize();
		thruster_in = (struct0_T*) calloc(numofThusters,sizeof(struct0_T));
		reald_file(thruster_in, &rudder_table[0]);
		init_read = 1;
	
	}
		
	emxInit_real_T(&solution, 1);

			if(label_t >= 9)
		{
			label_t = 1;
		}
			/*
			labelt = fopen("tt.txt","w");
			fprintf(
				labelt,"%E\t\n", 
			solution->size);
			fclose(labelt);
*/

	//    test_file = fopen(fpfile,"w");
		
 /*
		fprintf(test_file,"%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", 
		
			(double)pInput->label,
			(double)pInput->enable,
			pInput->TSP,
			pInput->ASP,
			pInput->Tx,
			pInput->Ty,
			pInput->T,
			pInput->phi,
			pInput->Tm
			);
	 fclose(test_file);

	 */	
/*
	fprintf(test_file,"%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n", 
		
			(double)pInput->label,
			(double)pInput->type,
			(double)pInput->enable,
			pInput->dt,
			pInput->y0,
			pInput->x,
			pInput->y,
			pInput->T,
			pInput->phi,
			pInput->Tmax,
			pInput->Tmin,
			pInput->dTmax,
			pInput->dphi_max,
			pInput->phi_min,
			pInput->phi_max,
			pInput->Trx,
			pInput->Try,
			pInput->Trm);
	 fclose(test_file);
	

			pOutput->label = 0;
			pOutput->enable = 0;
			pOutput->TSP= 0.0;
			pOutput->ASP= 0.0;
			pOutput->Tx= 0.0;//current force
			pOutput->Ty= 0.0;//current angle
			pOutput->T= 0.0;
			pOutput->phi= 0.0;
			pOutput->Tm= 0.0;
*/

/*
	pOutput->label = alloc_out[count_out].label;
	pOutput->enable = alloc_out[count_out].enable;
	pOutput->Tm = alloc_out[count_out].Tm;
	pOutput->Tx = alloc_out[count_out].Tx;
	pOutput->Ty = alloc_out[count_out].Ty;
	memcpy(&pOutput->T , &alloc_out[count_out].T, sizeof(alloc_out[count_out].T));	
	memcpy(&pOutput->phi , &alloc_out[count_out].phi, sizeof(alloc_out[count_out].phi));	
	
	pOutput->TSP = (double)status;
	pOutput->ASP = (double)flag;//test
		count_out++;
	if (count_out>=8) count_out = 0;
*/



	if(init_thr == 0)//init
	{

		for(i = 0; i<numofThusters; i++)
		{
		//	init_thruster(&thruster[i]);
		//	init_thruster(&thruster_in[i]);
			
		}
		init_thr = 1;
		
	}

		if(init_allc_out == 0)//init
	{
		for(i = 0; i<numofThusters; i++)
		{
		//	init_alloation_out(&alloc_out[i]);
			
		}
		init_allc_out = 1;
		
	}



	if (label_t == pInput->label)//change to int from doulbe!!!!!!!!!!next week
	{
			
			memcpy(&thruster_in[label_t-1].label, &pInput->label, sizeof(pInput->label));
			//memcpy(&thruster_in[label_t-1].type, &pInput->type, sizeof(pInput->type));
			memcpy(&thruster_in[label_t-1].enable, &pInput->enable, sizeof(pInput->enable));
			memcpy(&thruster_in[label_t-1].dt, &pInput->dt, sizeof(pInput->dt));
			//memcpy(&thruster_in[label_t-1].x0, &pInput->x0, sizeof(pInput->x0));
			//memcpy(&thruster_in[label_t-1].y0, &pInput->y0, sizeof(pInput->y0));
			//memcpy(&thruster_in[label_t-1].x, &pInput->x, sizeof(pInput->x));
			//memcpy(&thruster_in[label_t-1].y, &pInput->y, sizeof(pInput->y));
			memcpy(&thruster_in[label_t-1].T, &pInput->T, sizeof(pInput->T));
			memcpy(&thruster_in[label_t-1].phi, &pInput->phi, sizeof(pInput->phi));
			//memcpy(&thruster_in[label_t-1].Tmax, &pInput->Tmax, sizeof(pInput->Tmax));
			//memcpy(&thruster_in[label_t-1].Tmin, &pInput->Tmin, sizeof(pInput->Tmin));
			//memcpy(&thruster_in[label_t-1].dTmax, &pInput->dTmax, sizeof(pInput->dTmax));
			//memcpy(&thruster_in[label_t-1].dphi_max, &pInput->dphi_max, sizeof(pInput->dphi_max));
			//memcpy(&thruster_in[label_t-1].phi_min, &pInput->phi_min, sizeof(pInput->phi_min));
			//memcpy(&thruster_in[label_t-1].phi_max, &pInput->phi_max, sizeof(pInput->phi_max));

		//	thruster_in[label_t-1].weight[0] = 1.0;
		//	thruster_in[label_t-1].weight[1] = 1.0;
		//	thruster_in[label_t-1].weight_s[0] = 1000000.0;
		//	thruster_in[label_t-1].weight_s[1] = 1000000.0;
		//	thruster_in[label_t-1].weight_s[2] = 100000000.0;



			memcpy(&T_reqr.Tx, &pInput->Trx, sizeof(pInput->Trx));
			memcpy(&T_reqr.Ty, &pInput->Try, sizeof(pInput->Try));
			memcpy(&T_reqr.Tm, &pInput->Trm, sizeof(pInput->Trm));




			label_t = label_t + 1;



	}



	if(label_t == 9)//if all 8 data are retrieved
	{

		N_enabled_thr = 0;
		for (i = 0; i<8; i++)
		{
			
			if(thruster_in[i].enable == 1)//copy enabled thruster for calculation
			{
				
				
				N_enabled_thr = N_enabled_thr * (N_enabled_thr<8);
				memcpy(&thruster_p[N_enabled_thr],&thruster_in[i],sizeof(thruster_in[i]));
				


				N_enabled_thr++;

				
			}
		}
		
		if (N_enabled_thr>=1)
		{

/*		void qpsolver(struct0_T thruster_data[8], const struct1_T *T_r, double
              N_enabled_thruster, const double rudder_table[45], double method,
              emxArray_real_T *solution, struct2_T alloc_out_data[], int
              alloc_out_size[2], double *status)
*/
			qpsolver(thruster_p, &T_reqr, N_enabled_thr, rudder_table, 0.0,1.0, solution, alloc_out,
                      &status);

	  	  	labelt = fopen("yes0.txt","w");
			fprintf(
				labelt,"%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t\n\
%s\t\n\
%E\t%E\t%E\t\n\
%E\t\n", 

			status,
			thruster_p[0].dt,
			rudder_table[6],
			(double)N_enabled_thr,

			"T",
			thruster_p[0].T,
			thruster_p[1].T,
			thruster_p[2].T,
			thruster_p[3].T,
			thruster_p[4].T,
			thruster_p[5].T,
			thruster_p[6].T,
			thruster_p[7].T,

			"phi",
			thruster_p[0].phi,
			thruster_p[1].phi,
			thruster_p[2].phi,
			thruster_p[3].phi,
			thruster_p[4].phi,
			thruster_p[5].phi,
			thruster_p[6].phi,
			thruster_p[7].phi,

			"Tmax",
			thruster_p[0].Tmax,
			thruster_p[1].Tmax,
			thruster_p[2].Tmax,
			thruster_p[3].Tmax,
			thruster_p[4].Tmax,
			thruster_p[5].Tmax,
			thruster_p[6].Tmax,
			thruster_p[7].Tmax,

			"Tplus",
			thruster_p[0].Tplus,
			thruster_p[1].Tplus,
			thruster_p[2].Tplus,
			thruster_p[3].Tplus,
			thruster_p[4].Tplus,
			thruster_p[5].Tplus,
			thruster_p[6].Tplus,
			thruster_p[7].Tplus,

			"T_",
			thruster_p[0].T_,
			thruster_p[1].T_,
			thruster_p[2].T_,
			thruster_p[3].T_,
			thruster_p[4].T_,
			thruster_p[5].T_,
			thruster_p[6].T_,
			thruster_p[7].T_,

			"phiPlus",
			thruster_p[0].phiPlus,
			thruster_p[1].phiPlus,
			thruster_p[2].phiPlus,
			thruster_p[3].phiPlus,
			thruster_p[4].phiPlus,
			thruster_p[5].phiPlus,
			thruster_p[6].phiPlus,
			thruster_p[7].phiPlus,

			"phi_",
			thruster_p[0].phi_,
			thruster_p[1].phi_,
			thruster_p[2].phi_,
			thruster_p[3].phi_,
			thruster_p[4].phi_,
			thruster_p[5].phi_,
			thruster_p[6].phi_,
			thruster_p[7].phi_,

			"phi_min",
			thruster_p[0].phi_min,
			thruster_p[1].phi_min,
			thruster_p[2].phi_min,
			thruster_p[3].phi_min,
			thruster_p[4].phi_min,
			thruster_p[5].phi_min,
			thruster_p[6].phi_min,
			thruster_p[7].phi_min,

			"phi_max",
			thruster_p[0].phi_max,
			thruster_p[1].phi_max,
			thruster_p[2].phi_max,
			thruster_p[3].phi_max,
			thruster_p[4].phi_max,
			thruster_p[5].phi_max,
			thruster_p[6].phi_max,
			thruster_p[7].phi_max,

			"alloc_out-T",
			alloc_out[0].T,
			alloc_out[1].T,
			alloc_out[2].T,
			alloc_out[3].T,
			alloc_out[4].T,
			alloc_out[5].T,
			alloc_out[6].T,
			alloc_out[7].T,

			"alloc_out-phi",
			alloc_out[0].phi,
			alloc_out[1].phi,
			alloc_out[2].phi,
			alloc_out[3].phi,
			alloc_out[4].phi,
			alloc_out[5].phi,
			alloc_out[6].phi,
			alloc_out[7].phi,


			"solution1-8",
			solution->data[0],
			solution->data[1],
			solution->data[2],
			solution->data[3],
			solution->data[4],
			solution->data[5],
			solution->data[6],
			solution->data[7],

			"solution9-16",
			solution->data[8],
			solution->data[9],
			solution->data[10],
			solution->data[11],
			solution->data[12],
			solution->data[13],
			solution->data[14],
			solution->data[15],

			"T_reqr",
			T_reqr.Tx,
			T_reqr.Ty,
			T_reqr.Tm,

			solution->data[10],
			alloc_out[5].T);
			fclose(labelt);
			//
				for(i = 0; i<8; i++)
				{
					(pOutput+i)->label = alloc_out[i].label;
					//memcpy(&((pOutput+i)->label) , &alloc_out[i].label, sizeof(alloc_out[i].label));	
					(pOutput+i)->enable = alloc_out[i].enable;
					(pOutput+i)->Tm = alloc_out[i].Tm;
					(pOutput+i)->Tx = alloc_out[i].Tx;
					(pOutput+i)->Ty = alloc_out[i].Ty;
					memcpy(&(pOutput+i)->T , &alloc_out[i].T, sizeof(alloc_out[i].T));	
					memcpy(&(pOutput+i)->phi , &alloc_out[i].phi, sizeof(alloc_out[i].phi));	
					
					(pOutput+i)->TSP = (double)status;
					(pOutput+i)->ASP = (double)flag;//test
				}
	//		output_enable = 1;
		//	void qpsolver(struct0_T thruster_data[8], const struct1_T *T_r, double
         //     N_enabled_thruster, const double rudder_table[45], double method,
          //    emxArray_real_T *solution, struct2_T *alloc_out, double *status)

		}

	}
	emxDestroyArray_real_T(solution);


	/*
if (output_enable == 1)
{
	
	pOutput->label = alloc_out[count_out].label;
	pOutput->enable = alloc_out[count_out].enable;
	pOutput->Tm = alloc_out[count_out].Tm;
	pOutput->Tx = alloc_out[count_out].Tx;
	pOutput->Ty = alloc_out[count_out].Ty;
	memcpy(&pOutput->T , &alloc_out[count_out].T, sizeof(alloc_out[count_out].T));	
	memcpy(&pOutput->phi , &alloc_out[count_out].phi, sizeof(alloc_out[count_out].phi));	
	
	pOutput->TSP = flag1;
	pOutput->ASP = (double)flag;//test


	count_out++;
	if (count_out>=8)
	{count_out = 0;}

	output_enable = 0;
}
*/


}