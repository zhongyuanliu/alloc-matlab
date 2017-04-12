#include <stdio.h>
#include "qpmn.h"
#include "string.h"
#include "qpsolver.h"
//#include "stdlib.h"
#define format_sign ','
#define MAXCPL   200   //每行最大字符数
#define MAXCITY  100  //每组数据中DATA最多项数，DIMENSION的最大值
#define MAXNAMEL 32   //NAME最大长度

void reald_file(struct0_T *thruster_in, double *rudder_table)
{
	int length = 10;
	int card = 0;
	int c_i = 0;

	double fdata;
	FILE *f;
	char ln[MAXCPL];
	char *lnp = ln;
	f = fopen("thruster_config.txt","r");
	if(f == NULL)
	{
		return;
	}
	
	while(1)
	{
		if (NULL==fgets(ln,MAXCPL,f)) 
			{
				
				break;
			}
		lnp = ln;

		//printf("%s",ln);
		
		switch (card)
		{
		case 1://Xcg
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].x0 = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;
		case 2://Ycg
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
				thruster_in[c_i].y0 = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;
		case 3://Type
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
				thruster_in[c_i].type = (int)fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;
		case 4://X
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
				thruster_in[c_i].x = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

		case 5://Y
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].y = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

		case 6://Wx
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].weight[0] = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 7://Wy
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].weight[1] = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 8://Ws1
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].weight_s[0] = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 9://Ws2
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].weight_s[1] = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 10://Ws3
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].weight_s[2] = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 11://Tmin
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].Tmin = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 12://Tmax
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].Tmax = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 13://dTmax
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].dTmax = fdata;
			c_i++;
			
			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 14://dphiMax
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].dphi_max = fdata;
			c_i++;

			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 15://phi_min
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].phi_min = fdata;
			c_i++;

			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;

			case 16://phi_max
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			thruster_in[c_i].phi_max = fdata;
			c_i++;

			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			
			}
			break;


		case 17:
			c_i = 0;
			while( *lnp == format_sign ) lnp++ ;
			while(sscanf(lnp, "%lf",&fdata)==1)
			{
			//atof
			*(rudder_table + c_i) = fdata;
			c_i++;

			while( *lnp != format_sign && *lnp != '\n') 
				{
					lnp++ ; //跳过已读过的数		
				}

			while( *lnp == format_sign )
				{
					lnp++ ;
				}
			if (*lnp == '\n')
			{
				fgets(ln,MAXCPL,f);
				lnp = ln;
			}
			
			}
			break;
		}

		if (strcmp(ln,"Xcg\n") == 0) card = 1;
		
		else if (strcmp(ln,"Ycg\n") == 0) card = 2;
		
		else if (strcmp(ln,"Type\n") == 0) card = 3;
		
		else if (strcmp(ln,"X\n") == 0) card = 4;
		
		else if (strcmp(ln,"Y\n") == 0) card = 5;
		
		else if (strcmp(ln,"Wx\n") == 0) card = 6;
		
		else if (strcmp(ln,"Wy\n") == 0) card = 7;
		
		else if (strcmp(ln,"Ws1\n") == 0) card = 8;
		
		else if (strcmp(ln,"Ws2\n") == 0) card = 9;
		
		else if (strcmp(ln,"Ws3\n") == 0) card = 10;
		
		else if (strcmp(ln,"Tmin\n") == 0) card = 11;
		
		else if (strcmp(ln,"Tmax\n") == 0) card = 12;
		
		else if (strcmp(ln,"dTmax\n") == 0) card = 13;
		
		else if (strcmp(ln,"dphiMax\n") == 0) card = 14;
		
		else if (strcmp(ln,"phi_min\n") == 0) card = 15;
		
		else if (strcmp(ln,"phi_max\n") == 0) card = 16;
		
		else if (strcmp(ln,"rudder_table\n") == 0) card = 17;
		
		else if (strcmp(ln,"END\n") == 0) break;

		


	}

}
