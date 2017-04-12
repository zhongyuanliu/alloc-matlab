#define size_template 1000
#define pi 3.141592653589793
#define numofThusters 8//max number, it does not mean there are 10 available thrusters
#define max_var 100//max number of variables

#pragma pack(1)
typedef struct
{
	int label;
    int type;//type of thruster: tunnel:0,azi:1
	int enable;
	double x0;//cog
	double y0;//cog
    double x;
    double y;
	double T;//current force
	double phi;//current angle
    double Tmax;
    double Tmin;
	double dTmax;//max dT
	double dphi_max;//max dphi
	double Tplus;//=min{T + dt * dTmax, Tmax}
	double T_;//={T - dt * dTmax, Tmin}
	double phiPlus;//={phi + dt * dphi_max, phi_max}
	double phi_;//={phi - dt * dphi_max, phi_max}
    double phi_min; //forbidden zone
    double phi_max; //forbidden zone

} thruster_data;
typedef struct
{
	int label;//NEEDED
	int enable;//NEEDED
	double dt;//NEEDED
	double T;//NEEDED//current force
	double phi;//NEEDED//current angle
	double Trx;//NEEDED//required force
	double Try;//NEEDED//required force
	double Trm;//NEEDED//required force

} thruster_data_in;
typedef struct
{
	int label;
	int enable;
    double TSP;
    double ASP;
	double Tx;//current force
	double Ty;//current angle
    double T;
    double phi;
	double Tm;
	//char outend;
} alloation_out;

typedef struct
{
    int N;//number of rows
    char typeofCons[size_template][64];
    char cons_name[size_template][64]; //constraint names
} input_ROWS;
typedef struct
{
    int N;//number of rows
    char var[size_template][64];//variable names
    char cons_name[size_template][64]; //constraint names
    double values[size_template];

} input_COLUMNS;
typedef struct
{
    int N;//number of rows
    char name[size_template][64];
    char cons_name[size_template][64]; //constraint names
    double values[size_template];
} input_RHS;
typedef struct
{
    int N;//number of rows
    char name[size_template][64];
    char var[size_template][64]; //variable names
    double values[size_template];
} input_BOUNDS;
typedef struct
{
    int N;//number of rows
    char var1[size_template][64];
    char var2[size_template][64]; //variable names
    double values[size_template];
} input_QUADS;


typedef struct
{
    int N;//number of rows
	int sol_code;/*OPT_UNSOL=-1,OPT_OK=0,
    OPT_FOUND=1,INF_PROB=2,UNB_PROB=3,EXC_ITER=4,
    FEA_PROB=5, UNK_SOL=6, UNK_PROB=7,NUM_PROB=8*/
    char var[max_var][64];
	double UP[max_var];
	double LO[max_var];
    double values[max_var];
	double pobj[200];//max number of iter
	double dobj[200];//max number of iter
	double elapsed_t;//elapsed time
} quad_solution;

typedef struct
{
	double Tx;
	double Ty;
    double Tm;
} T_r;//required force

#pragma pack()

//void init_thruster(thruster_data* thruster);
//void init_alloation_out(alloation_out *alloc_out);
//int quad_run(thruster_data* thruster_in, T_r* T_reqr , int N_enabled_thruster, alloation_out* alloc_out);
//extern T_r T_req;//required force
//extern thruster_data thruster[numofThusters];
//extern int quad_run();