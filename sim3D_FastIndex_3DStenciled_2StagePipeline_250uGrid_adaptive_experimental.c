/*This is now a heart with a small ischemic region*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define CANINE 0
#define HUMAN  1

#define PIPELINE 1
#define NORMAL   0

#define LONG_QT 0
#define ISCHEMIA 0

#define QTR_mm_SPACING  1
#define HALF_mm_SPACING 0
#define HALF_mm_SPACING_ALL_TRANSVERSE 0
#define THIRD_mm_SPACING 0

#define BZ_IMPORT

#define CHUNK 5
#define NUMLOOPS 30000

#ifdef BZ_IMPORT
#include "bzip2-1.0.2/bzlib.h"

#ifdef _WIN32
#define BZ2_LIBNAME "bzip2-1.0.2/libbz2.dll"
#include <windows.h>
static int BZ2DLLLoaded = 0;
static HINSTANCE BZ2DLLhLib;
int BZ2DLLLoadLibrary(void)
{
   HINSTANCE hLib;

   if(BZ2DLLLoaded==1){return 0;}
   hLib=LoadLibrary(BZ2_LIBNAME);
   if(hLib == NULL){
      fprintf(stderr,"Can't load %s\n",BZ2_LIBNAME);
      return -1;
   }
   BZ2_bzlibVersion=GetProcAddress(hLib,"BZ2_bzlibVersion");
   BZ2_bzopen=GetProcAddress(hLib,"BZ2_bzopen");
   BZ2_bzdopen=GetProcAddress(hLib,"BZ2_bzdopen");
   BZ2_bzread=GetProcAddress(hLib,"BZ2_bzread");
   BZ2_bzwrite=GetProcAddress(hLib,"BZ2_bzwrite");
   BZ2_bzflush=GetProcAddress(hLib,"BZ2_bzflush");
   BZ2_bzclose=GetProcAddress(hLib,"BZ2_bzclose");
   BZ2_bzerror=GetProcAddress(hLib,"BZ2_bzerror");

   if (!BZ2_bzlibVersion || !BZ2_bzopen || !BZ2_bzdopen
       || !BZ2_bzread || !BZ2_bzwrite || !BZ2_bzflush
       || !BZ2_bzclose || !BZ2_bzerror) {
      fprintf(stderr,"GetProcAddress failed.\n");
      return -1;
   }
   BZ2DLLLoaded=1;
   BZ2DLLhLib=hLib;
   return 0;

}
int BZ2DLLFreeLibrary(void)
{
   if(BZ2DLLLoaded==0){return 0;}
   FreeLibrary(BZ2DLLhLib);
   BZ2DLLLoaded=0;
}
#endif /* WIN32 */
#endif

#define True				1
#define False				0

/****************************************/
#if LONG_QT
/* Long QT (APD Restitution parameters) - Endocardial */
#define END_MAX	400.0		/* Maximum (asymptotic APD) ms */
#define END_MIN	150.0		/* Minimum APD at minimum DI */
#define END_MID	250.0		/* APD where line and exponent meet */
#define END_DMN	80			/* Minimum diastolic interval */
#define END_DMD 160			/* diastolic iterval where line and exponent meet */
#define END_DMX 350 		/* Maximum diastolic interval */

/* Long QT (APD Restitution parameters) - Epicardial */
#define EPI_MAX	325.0		/* Maximum (asymptotic APD) ms */
#define EPI_MIN	100.0		/* Minimum APD at minimum DI */
#define EPI_MID	200.0		/* APD where line and exponent meet */
#define EPI_DMN	80			/* Minimum diastolic interval */
#define EPI_DMD 160			/* diastolic iterval where line and exponent meet */
#define EPI_DMX 350			/* Maximum diastolic interval */

/* Long QT (APD Restitution parameters) - Mid Myocardial */
#define MID_MAX	450.0		/* Maximum (asymptotic APD) ms */
#define MID_MIN	200.0		/* Minimum APD at minimum DI */
#define MID_MID	300.0		/* APD where line and exponent meet */
#define MID_DMN	80			/* Minimum diastolic interval */
#define MID_DMD 160			/* diastolic iterval where line and exponent meet */
#define MID_DMX 350			/* Maximum diastolic interval */

#elif ISCHEMIA

/* Ischemia (APD Restitution parameters) - Endocardial */
#define END_MAX	375.0		/* Maximum (asymptotic APD) ms */
#define END_MIN	125.0		/* Minimum APD at minimum DI */
#define END_MID	275.0		/* APD where line and exponent meet */
#define END_DMN	145			/* Minimum diastolic interval */
#define END_DMD 245			/* diastolic iterval where line and exponent meet */
#define END_DMX 250 		/* Maximum diastolic interval */

/* Ischemia (APD Restitution parameters) - Epicardial */
#define EPI_MAX	275.0		/* Maximum (asymptotic APD) ms */
#define EPI_MIN	125.0		/* Minimum APD at minimum DI */
#define EPI_MID	175.0		/* APD where line and exponent meet */
#define EPI_DMN	145			/* Minimum diastolic interval */
#define EPI_DMD 245			/* diastolic iterval where line and exponent meet */
#define EPI_DMX 250			/* Maximum diastolic interval */

/* Ischemia (APD Restitution parameters) - Mid Myocardial */
#define MID_MAX	400.0		/* Maximum (asymptotic APD) ms */
#define MID_MIN	125.0		/* Minimum APD at minimum DI */
#define MID_MID	300.0		/* APD where line and exponent meet */
#define MID_DMN	145			/* Minimum diastolic interval */
#define MID_DMD 245			/* diastolic iterval where line and exponent meet */
#define MID_DMX 250			/* Maximum diastolic interval */

#else
/* Normal (APD Restitution parameters) - Endocardial */
#define END_MAX	375.0		/* Maximum (asymptotic APD) ms */
#define END_MIN	125.0		/* Minimum APD at minimum DI */
#define END_MID	275.0		/* APD where line and exponent meet */
#define END_DMN	25			/* Minimum diastolic interval */
#define END_DMD 125			/* diastolic iterval where line and exponent meet */
#define END_DMX 250 		/* Maximum diastolic interval */

/* Normal (APD Restitution parameters) - Epicardial */
#define EPI_MAX	275.0		/* Maximum (asymptotic APD) ms */
#define EPI_MIN	125.0		/* Minimum APD at minimum DI */
#define EPI_MID	175.0		/* APD where line and exponent meet */
#define EPI_DMN	25			/* Minimum diastolic interval */
#define EPI_DMD 125			/* diastolic iterval where line and exponent meet */
#define EPI_DMX 250			/* Maximum diastolic interval */

/* Normal (APD Restitution parameters) - Mid Myocardial */
#define MID_MAX	400.0		/* Maximum (asymptotic APD) ms */
#define MID_MIN	125.0		/* Minimum APD at minimum DI */
#define MID_MID	300.0		/* APD where line and exponent meet */
#define MID_DMN	25			/* Minimum diastolic interval */
#define MID_DMD 125			/* diastolic iterval where line and exponent meet */
#define MID_DMX 250			/* Maximum diastolic interval */
#endif

/* Ischemia (APD Restitution parameters) - Endocardial */
#define ISCHEMIC_END_MAX 375.0		/* Maximum (asymptotic APD) ms */
#define ISCHEMIC_END_MIN 125.0		/* Minimum APD at minimum DI */
#define ISCHEMIC_END_MID 275.0		/* APD where line and exponent meet */
#define ISCHEMIC_END_DMN	170			/* Minimum diastolic interval */
#define ISCHEMIC_END_DMD 270			/* diastolic iterval where line and exponent meet */
#define ISCHEMIC_END_DMX 275 		/* Maximum diastolic interval */
//#define ISCHEMIC_END_DMN 200		/* Minimum diastolic interval */
//#define ISCHEMIC_END_DMD 300		/* diastolic iterval where line and exponent meet */
//#define ISCHEMIC_END_DMX 305 		/* Maximum diastolic interval */

/* Ischemia (APD Restitution parameters) - Epicardial */
#define ISCHEMIC_EPI_MAX	275.0		/* Maximum (asymptotic APD) ms */
#define ISCHEMIC_EPI_MIN	125.0		/* Minimum APD at minimum DI */
#define ISCHEMIC_EPI_MID	175.0		/* APD where line and exponent meet */
#define ISCHEMIC_EPI_DMN	170			/* Minimum diastolic interval */
#define ISCHEMIC_EPI_DMD 270			/* diastolic iterval where line and exponent meet */
#define ISCHEMIC_EPI_DMX 275			/* Maximum diastolic interval */
//#define ISCHEMIC_EPI_DMN 200			/* Minimum diastolic interval */
//#define ISCHEMIC_EPI_DMD 300			/* diastolic iterval where line and exponent meet */
//#define ISCHEMIC_EPI_DMX 305			/* Maximum diastolic interval */

/* Ischemia (APD Restitution parameters) - Mid Myocardial */
#define ISCHEMIC_MID_MAX	400.0		/* Maximum (asymptotic APD) ms */
#define ISCHEMIC_MID_MIN	125.0		/* Minimum APD at minimum DI */
#define ISCHEMIC_MID_MID	300.0		/* APD where line and exponent meet */
#define ISCHEMIC_MID_DMN	170			/* Minimum diastolic interval */
#define ISCHEMIC_MID_DMD 270			/* diastolic iterval where line and exponent meet */
#define ISCHEMIC_MID_DMX 275			/* Maximum diastolic interval */
//#define ISCHEMIC_MID_DMN	200			/* Minimum diastolic interval */
//#define ISCHEMIC_MID_DMD 300			/* diastolic iterval where line and exponent meet */
//#define ISCHEMIC_MID_DMX 305			/* Maximum diastolic interval */


#define ENDOCARDIAL_CELL	0
#define MIDMYOCARDIAL_CELL	1
#define EPICARDIAL_CELL		2

/* Set up PDE restitution curve parameters endo */
static float ischemic_apd_end_slp;
static float ischemic_apd_end_exp;

/* Set up PDE restitution curve parameters epi */
static float ischemic_apd_epi_slp;
static float ischemic_apd_epi_exp;

/* Set up PDE restitution curve parameters mid */
static float ischemic_apd_mid_slp;
static float ischemic_apd_mid_exp;

static float apd_end_slp;	/* Slope for endo restitution curve (normal) */
static float apd_end_exp;	/* Exponent for endo restitution curve (normal) */

static float apd_epi_slp;	/* Slope for epi restitution curve (normal) */
static float apd_epi_exp;	/* Exponent for epi restitution curve (normal) */

static float apd_mid_slp;	/* Slope for mid restitution curve (normal) */
static float apd_mid_exp;	/* Exponent for mid restitution curve (normal) */

/*****************************************/


/****************************************/
#if QTR_mm_SPACING
#define SPACING		0.025	/* Grid step in cm*/
#define Z_MAX		392		/* Lowest z level heart geometry is still valid*/

#define R_D		2.262	/* Dimensionless fast axis neighborhood size	  */
#define K		2.811	/* Threshold in element units for anisotropy   */

#define X_GRID	442	/* Number of elements in X for 12 cm (0.05 space)  */
#define Y_GRID	442	/* Number of elements in Y for 12 cm (0.05 space)  */
#define Z_GRID	442	/* Number of elements in Z for 12 cm (0.05 space)  */

/* Actual time step and grid step for CA solutions */		
#define DT     0.002		/* Time step (sec)  2.0 ms				*/
#define DX	   0.025		/* Space step (cm)  250 micrometers	*/ //changed

/* ADI method parameters (time step is half DT) */
#define DTP    0.0005	/* Time step (sec)  0.5 ms			*/	//changed
#define DIFF   0.1		/* Isotropic diffusion constant 		*/ //changed
#define FACT   0.08		/* DIFF*DTP/(DX**2) 						*/ //changed

#define PM_MAX 1		/* ((DT/DTP/2) - 1) # of double steps PDE must take*/

/****************************************/
#elif HALF_mm_SPACING
#define SPACING		0.05	/* Grid step in cm */
#define Z_MAX			196	/* Lowest z level heart geometry is still valid*/

#define R_D		3.60	/* Dimensionless fast axis neighborhood size 	  */
#define K		4.95	/* Threshold in element units for anisotropy=     */

#define X_GRID	221	/* Number of elements in X for 12 cm (0.05 space)  */
#define Y_GRID	221	/* Number of elements in Y for 12 cm (0.05 space)  */
#define Z_GRID	221	/* Number of elements in Z for 12 cm (0.05 space)  */

/* Actual time step and grid step for CA solutions */		
#define DT     0.002		/* Time step (sec)  2.0 ms				*/
#define DX	   0.05		/* Space step (cm)  500 micrometers	*/

/* ADI method parameters (time step is half DT) */
#define DTP    0.00025	/* Time step (sec)  0.5 ms				*/
#define DIFF   0.7		/* Isotropic diffusion constant 		*/ 
#define FACT   0.07		/* DIFF*DTP/(DX**2) 						*/ 

#define PM_MAX 3		/* ((DT/DTP/2) - 1) # of double steps PDE must take*/

/****************************************/
#elif HALF_mm_SPACING_ALL_TRANSVERSE
#define SPACING		0.05	/* Grid step in cm */
#define Z_MAX		196	/* Lowest z level heart geometry is still valid*/

#define R_D		1.131	/* Dimensionless fast axis neighborhood size 	  */
#define K		0.352	/* Threshold in element units for anisotropy=     */

#define X_GRID	221	/* Number of elements in X for 12 cm (0.05 space)  */
#define Y_GRID	221	/* Number of elements in Y for 12 cm (0.05 space)  */
#define Z_GRID	221	/* Number of elements in Z for 12 cm (0.05 space)  */

/* Actual time step and grid step for CA solutions */		
#define DT     0.002		/* Time step (sec)  2.0 ms				*/
#define DX	   0.05		/* Space step (cm)  500 micrometers	*/

/* ADI method parameters (time step is half DT) */
#define DTP    0.002	/* Time step (sec)  0.5 ms				*/
#define DIFF   0.1		/* Isotropic diffusion constant 		*/ 
#define FACT   0.08		/* DIFF*DTP/(DX**2) 						*/ 

#define PM_MAX 0		/* ((DT/DTP/2) - 1) # of double steps PDE must take*/

/****************************************/
#elif THIRD_mm_SPACING
#define SPACING		0.035	/* Grid step in cm */
#define Z_MAX		279	/* Lowest z level heart geometry is still valid*/

#define R_D		1.61	/* Dimensionless fast axis neighborhood size 	  */
#define K		1.01	/* Threshold in element units for anisotropy=     */

#define X_GRID	315	/* Number of elements in X for 12 cm (0.035 space)  */
#define Y_GRID	315	/* Number of elements in Y for 12 cm (0.035 space)  */
#define Z_GRID	315	/* Number of elements in Z for 12 cm (0.035 space)  */

/* Actual time step and grid step for CA solutions */		
#define DT     0.001		/* Time step (sec)  2.0 ms				*/
#define DX	   0.035		/* Space step (cm)  500 micrometers	*/

/* ADI method parameters (time step is half DT) */
#define DTP    0.001	/* Time step (sec)  0.5 ms				*/
#define DIFF   0.1		/* Isotropic diffusion constant 		*/ 
#define FACT   0.0816	/* DIFF*DTP/(DX**2) 						*/ 

#define PM_MAX 0		/* ((DT/DTP/2) - 1) # of double steps PDE must take*/
#endif
/****************************************/

/* Action potential parameters */
#define UMAX	0.99	/* Maximum dimensionless potential value	*/
#define UMIN	0.01	/* Minimum dimensionless potential value	*/
#define UAPD	0.05  	/* 95% repolarization (end of APD)			*/

/* Normal (APD Restitution parameters) */
#define DI_INIT  500	/* Initial diastolic interval 			*/

#define S2		2000.0	/* Parameter for Excited Current Source	*/

#define N_EXCS	1		/* Excite for 3 steps (trigger waves for now)	  */
#define H_MASK	7		/* Half of voxel width for calculating offsets	  */

#define T_GRID	1		/* Number of theta divisions */
#define P_GRID	1		/* Number of phi divisions   */
#define N_TABS	100		/* Number of random voxel tables at theta, phi */

#define ANISOTROPY		1.0		/* Global Tissue anisotropy ratio for fibers  */

#define PI 3.1415927

typedef unsigned char Boolean;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define NINT(a) ((int)((a)+0.5))

#define BOUNDARY_TRUE	|= (unsigned char)1
#define BOUNDARY_FALSE	&= (unsigned char)254
#define EXC_TRUE 	|= (unsigned char)2
#define EXC_FALSE 	&= (unsigned char)253
#define REF_TRUE	|= (unsigned char)4
#define REF_FALSE	&= (unsigned char)251
#define REP_TRUE	|= (unsigned char)8
#define REP_FALSE	&= (unsigned char)247
#define INT_TRUE  |= (unsigned char)16
#define INT_FALSE &= (unsigned char)239
#define EPI_TRUE  |= (unsigned char)32
#define EPI_FALSE &= (unsigned char)223
#define END_TRUE  |= (unsigned char)64
#define END_FALSE &= (unsigned char)191
#define MID_TRUE  |= (unsigned char)128
#define MID_FALSE &= (unsigned char)127

#define BOUNDARY_TEST	& (unsigned char)1
#define EXC_TEST 	& (unsigned char)2
#define REF_TEST	& (unsigned char)4
#define REP_TEST	& (unsigned char)8
#define INT_TEST  & (unsigned char)16
#define EPI_TEST  & (unsigned char)32
#define END_TEST  & (unsigned char)64
#define MID_TEST  & (unsigned char)128

static unsigned short bitMasks[16] = 
	{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768};
	
typedef struct voxel {
	double x;
	double y; 
	double z;
		
} VOXEL;

/*Jon's code for fast loop indexing */
typedef struct spacial_element {
	unsigned short int	start;		/* index of first element inside heart  */
	unsigned short int	end;		/* index of last element inside heart */

} SPACIAL_ELEMENT;

SPACIAL_ELEMENT zBounds;
SPACIAL_ELEMENT *xBounds;
SPACIAL_ELEMENT **yBounds;


void initialize_fast_index(void);
/*end Jon's code*/

/* 3D array with bits defining local neighborhood  */
typedef struct voxel_table {
	unsigned short bits[16][16];	
	
} VOXEL_TABLE;

typedef struct cam_pde_element {
	unsigned char  src;		/* source strength					*/
	unsigned char  snk;		/* sinking value					*/		
	unsigned short tau;		/* diastolic interval				*/
	unsigned short nsr;		/* front sourcing steps				*/
	unsigned char  the_bin;	/* theta bin for voxel table 		*/
	unsigned char  phi_bin;	/* phi bin for voxel table			*/
	unsigned char  flags;	/* various flags condensed into one variable
						bit1: (boundary)	Is volume element on boundry	
						bit2: (exc) Is element EXCITING
						bit3: (ref)	Is element REFRACTORY 
						bit4: (rep)	Is element REPOLARIZED
						bit5: (int) Is element internal Boundary
						bit6: (epi) Is element epicardium
						bit7: (end) Is element endocardium
						bit8: (mid) Is element midmyocardium
						*/
	unsigned char ischemic;
	float u;			/* pseudo potential value				*/
	float tmp;			/* temporary potential value			*/
	float bet;			/* I-V curve parameter for desired apd	*/
} CAM_PDE_ELEMENT;

static VOXEL rvoxels[16][16][16];

static VOXEL_TABLE ***tableArray;

static CAM_PDE_ELEMENT	****grid;

typedef struct ecg {
	unsigned short time;	/* time step 							*/
	unsigned char i, j, k;	/* indices 								*/
	float source;			/* sum of 3 second spatial derivatives	*/

} ECG_ELEMENT;

static int iSeed = -9;
static float epsilon = 1.e-8;

/* Function Declarations */

#if CANINE
void initialize_geometry_canine(void);
void simple_pacing_from_bottom_canine(void);

#elif HUMAN
void initialize_geometry(void);
void initialize_ischemic_region(void);
void purkinje_fibre_stim_human(int);
void purkinje_fibre_stim_human2(int);
void purkinje_fibre_stim_human3(void);
void simple_pacing_from_bottom_human(void);
void simple_pacing_from_top_human(void);
#endif

static void set_boundary_conditions (void);
static void differentiate_heart_tissue(void);
static void differentiate_heart_tissue2(void);
void initialize_apds (void);
static int initialize_pde (void );

static int voxel_get_table_entry (VOXEL_TABLE *table, int i, int j, int k);
static void voxel_init_random_table_array (void);
static VOXEL_TABLE *voxel_get_random_table (unsigned char i, unsigned char j);
static int voxel_get_table (float the, float phi);
static void voxel_build_random_table (VOXEL_TABLE *table, double speed_ratio, double the, double phi);
static void cam_initialize ( void );

static float ran1(int *idum);
static double atanch_andy (double val);
static void compute_apd_parameters (float apd, float umax, float *tsh, float *bet);
static float apd_from_restitution (unsigned short tau, int cellType);

static time_t wrapper_start;
char timerName[256];

void write_2D_data(int);
void write_Slice_Geometry(void);

void write_3D_data(int);
void write_Compressed_3D_data(int);

void startTime(char*);
void stopTime(void);

#define LEAD_X 	0.0
#define LEAD_Z	-5.0
#define LEAD_Y	-9.0

#define V3_LEAD_X 	-9.0
#define V3_LEAD_Z	-5.0
#define V3_LEAD_Y	0.0

int Read_ECG ( void ) 
{
	long cnt, time;
	double dx2, dy2, dz2;
	double den, x, y, z, ecg;
	ECG_ELEMENT ecgElement;
	FILE *fp, *fo;
	
	char fname[256];
	
	fo = fopen("lead1_new.dat", "w");
	for (time=0; time < 500; time++) {
 	 	sprintf(fname, ":ECGFiles:ecg%d", time);
 		fp = fopen(fname, "rb");

		cnt = 0;
		ecg = 0.0;
		while (!feof(fp)) {
			fread(&ecgElement, sizeof(ECG_ELEMENT), 1, fp);
		
			z = -ecgElement.i*SPACING;
			x = -5.0 + ecgElement.j*SPACING;
			y = -5.0 + ecgElement.k*SPACING;
		
			dx2 = (x - LEAD_X)*(x - LEAD_X);
			dy2 = (y - LEAD_Y)*(y - LEAD_Y);
			dz2 = (z - LEAD_Z)*(z - LEAD_Z);
		
			den = sqrt(dx2 + dy2 + dz2);
			ecg += ecgElement.source/den;
			
			if ((++cnt % 100000 == 0))
				printf("count = %d\n", cnt);
			
		}
		fclose(fp);
		fprintf(fo, "%d\t%f\n", 2*time, ecg);
		printf("time = %d\tecg = %lf\n", 2*time, ecg);
	}
	fclose(fo);
	exit(1);
}

CAM_PDE_ELEMENT ***grid_i;
CAM_PDE_ELEMENT **grid_i_j;
CAM_PDE_ELEMENT *grid_i_j_k;

#if PIPELINE
CAM_PDE_ELEMENT ***grid_i2;
CAM_PDE_ELEMENT **grid_i2_j2;
CAM_PDE_ELEMENT *grid_i2_j2_k2;
#endif

int numACT = 0;
int numEXC = 0;

int main( void )
{

#if NORMAL
	float ux_m, uy_m, uz_m;
#elif PIPELINE
	int pipeline;
	int i2, j2, k2;
	float uxyz_m;
	int yStart, yEnd;
	int xStart, xEnd;
#endif

	float ux_p, uy_p, uz_p;
	time_t start, stop;
	unsigned short bitz;
	float diff;
	int tm, pm;
	int i, j, k;
	int ir, jr, kr;
	int iv, jv, kv;
	float apd, tsh, bet;
	float i_0, u_0;
	double ecgLead = 0.0;
	double V3_ecgLead = 0.0;
	double dx2, dy2, dz2, x, y, z;
	double V3_dx2, V3_dy2, V3_dz2;
		
	SPACIAL_ELEMENT* yBounds_i;
	SPACIAL_ELEMENT yBounds_i_j;
	FILE *ecg;
	CAM_PDE_ELEMENT *src_p, *snk_p;
	VOXEL_TABLE *table;
	
	/* Jon Labin */
	char BELL = 7;
	//get return value anc check for errors
/*	i = init_compression();
	if(i){
		if(i = BZ_CONFIG_ERROR)
			printf("the library has been mis-compiled");
		if(i = BZ_PARAM_ERROR)
			printf("strm is NULL\nor blockSize < 1 or blockSize > 9\nor verbosity < 0 or verbosity > 4\nor workFactor < 0 or workFactor > 250");
		if(i = BZ_MEM_ERROR)
			printf("not enough memory is available");
		exit(1);
	}

*/
#ifdef BZ_IMPORT
#ifdef _WIN32
   BZ2DLLLoadLibrary();
#endif
#endif

 /* Output Simulation Characteristics*/
#if HUMAN
	printf("3D HUMAN Heart Model\n");
#elif CANINE
	printf("3D CANINE Heart Model\n");
#endif

#if QTR_MICRON_SPACING
	printf("250 micron (0.25 cm) spacing between cellular automata\n");
#elif HALF_MICRON_SPACING
	printf("500 micron (0.5 cm) spacing between cellular automata\n");
#endif

	printf("Using Fast Indexing and Bit Masked Compact Booleans\n");
#if PIPELINE
	printf("Plus using 3D stenciling to speed up the PDE Solver\n");
	printf("Plus using a pipelined PDE Solver\n");
#endif
	printf("\nSimulation will run %i iterations\n\n",NUMLOOPS);

#if CANINE
#if QTR_MICRON_SPACING
	printf("Canine heart model can NOT be run with 250 micron grid spacing\n");
	printf("Simulation will exit!\n");
	exit(1);
#endif
#endif



 /* New Initialization Sequence*/
	
#if HUMAN
	initialize_geometry();
	
#elif CANINE
	initialize_geometry_canine();
#endif

	
	initialize_fast_index();
	set_boundary_conditions();
	differentiate_heart_tissue();
	initialize_apds();
	
	printf("Initialize Ischemic Region\n");
	initialize_ischemic_region();

	write_Compressed_3D_data(0);
	printf("Compressed Data WROTE!");

	voxel_init_random_table_array ();
	

 /* Open ECG file */
 	ecg = fopen("ECG1_Full.dat", "w");
	
	printf("%c", BELL);
 
	start = time(NULL);
	printf("Benchmarking Clock Started\n");

	printf("\n\n\n");
	printf("##############################################\n");
	printf("##############################################\n");
	printf("#########Beginning Cardiac Simulation#########\n");
	printf("##############################################\n");
	printf("##############################################\n");
	printf("\n\n\n");

	//write_Slice_Geometry();
	
/* Loop over the total number of time steps */
	for (tm = 0; tm < NUMLOOPS; tm++) {
		printf("CA iteration = %d\n", tm);
		fprintf(ecg, "%lf\t%lf\t%lf \t%i \t%i\n", (tm*2.0) , ecgLead, V3_ecgLead,numEXC,numACT);
		fflush(ecg);
		printf("time = %lf (ms) , ECG = %lf , ECG = %lf , numEXC = %i , numACT =%i\n", tm*2.0 , ecgLead, V3_ecgLead,numEXC,numACT);		
		ecgLead = 0.0;
		V3_ecgLead = 0.0;

	 /* Pacing the Heart at 60 bpm (tm%500)*/

		/* For 50000 iteration run
		if ((tm/1000) < 10)
			purkinje_fibre_stim_human2(tm % 250);
		else if ((tm/1000) < 20)
			purkinje_fibre_stim_human2(tm % 275);
		else if ((tm/1000) < 30)
			purkinje_fibre_stim_human2(tm % 300);
		else if ((tm/1000) < 40)
			purkinje_fibre_stim_human2(tm % 333);
		else if ((tm/1000) < 50)
			purkinje_fibre_stim_human2(tm % 500);
		*/

		/*
		if ( (tm < 5000) && ((tm % 353) == 0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 5000) && (tm < 10000) && ((tm % 333)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 10000) && (tm < 15000) && ((tm % 316)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 15000) && (tm < 20000) && ((tm % 300)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 20000) && (tm < 25000) && ((tm % 286)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		}
		*/

		
		if ( (tm < 5000) && ((tm % 500) == 0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 5000)&& ((tm % 261)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} 
		
		/*else if ( (tm >= 10000) && (tm < 15000) && ((tm % 300)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 15000) && (tm < 20000) && ((tm % 286)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 20000) && (tm < 25000) && ((tm % 273)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 25000) && (tm < 30000) && ((tm % 261)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 30000) && (tm < 35000) && ((tm % 250)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 35000) && (tm < 40000) && ((tm % 240)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		} else if ( (tm >= 40000) && (tm < 45000) && ((tm % 231)==0)) {
			purkinje_fibre_stim_human3();
			printf("STIMULUS APPLIED, tm = %i\n",tm);
		}
		*/
		if (numEXC) {
	//	if (1) {
			numEXC = 0;


	 /* EXCITING elements source to neighbors */
		for (i=zBounds.start; i<=zBounds.end; i++) {
		grid_i = grid[i];
		yBounds_i = yBounds[i];
		for (j=xBounds[i].start; j<=xBounds[i].end; j++) {
		grid_i_j = grid_i[j];
		yBounds_i_j = yBounds_i[j];
		for (k=yBounds_i_j.start; k<=yBounds_i_j.end; k++) {
			if (grid_i_j[k]) {
					//continue; /* Loop up */
		
		 /* Get element pointer */
			src_p = grid_i_j[k];
			
		 /* Ignore non-EXCITING */
		 	if (src_p->flags EXC_TEST) {

				numEXC++;
						
		 	/* Get a random neighborhood table for this element */
				table = voxel_get_random_table (src_p->the_bin, src_p->phi_bin);
					 	
		  	/* Loop over voxel table and source to neighbors as defined in table */
		  		for (ir=0; ir < 15; ir++) {
		  		for (jr=0; jr < 15; jr++) {
					bitz = table->bits[ir][jr];
		  		for (kr=0; kr < 15; kr++) {

				/* See if element is a neighbor given local fiber direction */
				//	if (!voxel_get_table_entry(table, ir, jr, kr))
					if (!(bitMasks[kr] & bitz))
						continue;
				
				/* Compute array grid index from voxel index */
					iv = i- H_MASK+ir;
					jv = j- H_MASK+jr;
					kv = k- H_MASK+kr;
				
				/* Check array bounds on neighbors for off array status */
			 		if ((iv < 0) || (iv > (Z_GRID-1)) ||
						(jv < 0) || (jv > (X_GRID-1)) ||
						(kv < 0) || (kv > (Y_GRID-1))) {
						continue;
					}

				/* Don't source if neighbor not in the cardiac muscle */
			 		if (!grid[iv][jv][kv]) {
						continue;
					}
	
				 /* Get neighboring CAM element */
			 		snk_p = grid[iv][jv][kv];
				
				 /* Add source unit if element is not refractory */
					if (!(snk_p->flags REF_TEST)) {
			 			snk_p->snk += src_p->src;	/* send source to non front elements */
					}
				}}}
				
			 /* Update phase for element (after sourcing) */
		 		src_p->nsr++;
			}
			}
		 
		}}}  /* end loop over EXCITING elements */

		}

		if ( (numEXC>0) || (numACT >0)) {

	 /* Check for excitation of elements */
	//	for (i=0; i<Z_GRID; i++) {
		for (i=zBounds.start; i<=zBounds.end; i++) {
		grid_i = grid[i];
		yBounds_i = yBounds[i];
		for (j=xBounds[i].start; j<=xBounds[i].end; j++) {
		grid_i_j = grid_i[j];
		yBounds_i_j = yBounds_i[j];
		for (k=yBounds_i_j.start; k<=yBounds_i_j.end; k++) {

			if (grid_i_j[k]) {
		
			snk_p = grid_i_j[k];

		 /* Ignore exciting elements and refractory elements */
		 	if (snk_p->flags EXC_TEST || snk_p->flags REF_TEST) {
				/* clear all sinks! */
				snk_p->snk = 0;
				continue;
			} 

																		
		 /* Check for switch to EXCITING state for RESTING elements */
			if (snk_p->snk >= K) {
				snk_p->flags EXC_TRUE; 		/* EXCITING state 		*/
				snk_p->flags REF_TRUE;
				snk_p->flags REP_FALSE;		/* Not repolarized 		*/		
				snk_p->nsr = 0;			/* sourcing step zero 	*/
				numEXC++;
			} 
			/* clear all sinks! */
			snk_p->snk = 0;
			}
		}}}  /* end loop over excitation check */

		}


		if ( (numEXC>0) || (numACT >0)) {
	//if (1) {
			numACT=0;
		
	 /* Now, unfortunately for explicit solutions, we will have to do 8 pde steps for every CA step */
	 /* Now doing two pde steps per loop iteration to remove variable swapping */
#if NORMAL
	 	for (pm=0; pm <= PM_MAX; pm++) {
			printf("PDE iteration: %d\n", pm);
				
			for (i=zBounds.start; i<=zBounds.end; i++) {
			grid_i = grid[i];
			yBounds_i = yBounds[i];

			/*ECG CODE	Calculate ECG on last iteration of PDE loop */
			if (pm == PM_MAX) {
				z = -i*SPACING;
				dz2 = (z - LEAD_Z)*(z - LEAD_Z);
				V3_dz2 = (z - V3_LEAD_Z)*(z - V3_LEAD_Z);
			}
			for (j=xBounds[i].start; j<=xBounds[i].end; j++) {

			grid_i_j = grid_i[j];
			yBounds_i_j = yBounds_i[j];

			/*ECG CODE*/
			if (pm == PM_MAX) {
				x = -5.5 + j*SPACING;
				dx2 = (x - LEAD_X)*(x - LEAD_X);
				V3_dx2 = (x - V3_LEAD_X)*(x - V3_LEAD_X);
			}


			for (k=yBounds_i_j.start; k<=yBounds_i_j.end; k++) {
				//printf("(%i,%i,%i)\n",i,j,k);

			/* Make sure we are in the heart */
				if (!grid_i_j[k]) {
					continue; /* Loop up */
				}

				grid_i_j_k = grid_i_j[k];

			 /* Local potential and restitution parameter */
				u_0 = grid_i_j_k->u;
				bet = grid_i_j_k->bet;
						
			 /* Can turn this code into a "stencils" array and save lots of time! */
			 /* Get neighboring potentials, Implement no flux as required */
				if (!(grid_i_j_k->flags BOUNDARY_TEST)) {
					uz_p = grid[i+1][j][k]->u - u_0;
					uz_m = grid[i-1][j][k]->u - u_0;
					ux_p = grid_i[j+1][k]->u  - u_0;
					ux_m = grid_i[j-1][k]->u  - u_0;
					uy_p = grid_i_j[k+1]->u   - u_0;
					uy_m = grid_i_j[k-1]->u   - u_0;
				
				} else {
					if (i==0) {
						uz_m = 0.0;
					} else {
						uz_m = (grid[i-1][j][k]) ? (grid[i-1][j][k]->u - u_0): 0.0;
					}
					uz_p = (grid[i+1][j][k]) ? (grid[i+1][j][k]->u - u_0) : 0.0;
					ux_m = (grid_i[j-1][k]) ? (grid_i[j-1][k]->u - u_0) : 0.0;
					ux_p = (grid_i[j+1][k]) ? (grid_i[j+1][k]->u - u_0) : 0.0;
					uy_m = (grid_i_j[k-1]) ? (grid_i_j[k-1]->u - u_0) : 0.0;
					uy_p = (grid_i_j[k+1]) ? (grid_i_j[k+1]->u - u_0): 0.0;
				}
			
			 /* Compute 2nd differences for (isotropic) laplacian */
			//	d2x = ux_p + ux_m - 2*u_0;	
			//	d2y = uy_p + uy_m - 2*u_0;	
			//	d2z = uz_p + uz_m - 2*u_0;	
				
			 /* Set current depending on whether exc or rep */
				if (grid_i_j_k->flags EXC_TEST) {
					i_0 = S2*(1.0 - u_0);
				} else {
					i_0 = 2*bet*u_0*(u_0 - 1.0);
				}
				
			 /* Excplicit integration */
				//gpde_i_j_k->tmp = u_0 + FACT*gpde_i_j_k->dewnsud + DTP*i_0;
				grid_i_j_k->tmp = MAX((MIN((u_0 + FACT* (uz_m + uz_p + ux_m + ux_p + uy_m + uy_p) + DTP*i_0), UMAX)),UMIN);

				
				/*ECG CALCULATIONS*/
				if (pm == PM_MAX) {
					y = -5.5 + k*SPACING;
					dy2 = (y - LEAD_Y)*(y - LEAD_Y);
					V3_dy2 = (y - V3_LEAD_Y)*(y - V3_LEAD_Y);
					ecgLead += (ux_p + ux_m + uy_p + uy_m + uz_p + uz_m)/sqrt(dx2 + dy2 + dz2);
					V3_ecgLead += (ux_p + ux_m + uy_p + uy_m + uz_p + uz_m)/sqrt(V3_dx2 + V3_dy2 + V3_dz2);
					 /* Check for EXC -> REF */
				 	
				}
		
			}}}

			for (i=zBounds.start; i<=zBounds.end; i++) {
				grid_i = grid[i];
				yBounds_i = yBounds[i];
			for (j=xBounds[i].start; j<=xBounds[i].end; j++) {
				grid_i_j = grid_i[j];
				yBounds_i_j = yBounds_i[j];
			for (k=yBounds_i_j.start; k<=yBounds_i_j.end; k++) {

				grid_i_j_k = grid_i_j[k];

			/* Make sure we are in the heart */
				if (!grid_i_j_k) {
					continue; /* Loop up */
				}
				grid_i_j_k->u = grid_i_j_k->tmp;
				
				if (grid_i_j_k->u > UMIN)
					numACT++;


				grid_i_j_k->tmp = 0.0;
				if (grid_i_j_k->flags EXC_TEST) {
						
				 /* Check whether we are still in the front region for element */
					if (grid_i_j_k->nsr == N_EXCS) {
						grid_i_j_k->flags EXC_FALSE;				/* Cease CA and PDE sourcing */

						/* Compute new I-V parameters from peak potential for apd */
						if (grid_i_j_k->ischemic) {
							if(grid_i_j_k->flags EPI_TEST) {
								apd = ischemic_apd_from_restitution (grid_i_j_k->tau, EPICARDIAL_CELL); 
							} else if (grid_i_j_k->flags END_TEST) {
								apd = ischemic_apd_from_restitution (grid_i_j_k->tau, ENDOCARDIAL_CELL); 
							} else {
								apd = ischemic_apd_from_restitution (grid_i_j_k->tau, MIDMYOCARDIAL_CELL); 
							}
						} else {
							if(grid_i_j_k->flags EPI_TEST) {
								apd = apd_from_restitution (grid_i_j_k->tau, EPICARDIAL_CELL); 
							} else if (grid_i_j_k->flags END_TEST) {
								apd = apd_from_restitution (grid_i_j_k->tau, ENDOCARDIAL_CELL); 
							} else {
								apd = apd_from_restitution (grid_i_j_k->tau, MIDMYOCARDIAL_CELL); 
							}
						}

					   compute_apd_parameters (apd*0.001, grid_i_j_k->u, &tsh, 
								&(grid_i_j_k->bet));
			 
				 /* Assign new source strength at next autochrone based on sourcing table */
					} else {	
						/* This is a NOOP for now, will change with ->nsr; recovery effects */
						grid_i_j_k->src = 1;
					}		
				//continue;	/* Loop up */
				} else {
				/* Update diastolic interval for repolarized elements */
					if (grid_i_j_k->flags REP_TEST) {
						grid_i_j_k->tau += 2;		/* CA time steps in 2 milliseconds, restitution curves in ms */
					 /* Check for REF-> Excitable */
						if(grid_i_j_k->flags REF_TEST) {
							if(grid_i_j_k->flags EPI_TEST) {
								if( grid_i_j_k->tau >= EPI_DMN) 
									grid_i_j_k->flags REF_FALSE;
							} else if (grid_i_j_k->flags END_TEST) {
								if ( grid_i_j_k->tau >= END_DMN)
									grid_i_j_k->flags REF_FALSE;
							} else {
								if( grid_i_j_k->tau >= MID_DMN) 
									grid_i_j_k->flags REF_FALSE; 
							}
						}
					} else {
					/* Check for repolarization via PDEs */
						if (grid_i_j_k->flags REF_TEST &&  (grid_i_j_k->u <= UAPD) ) {
							grid_i_j_k->flags REP_TRUE;		/* end of action potential */
							grid_i_j_k->tau = 0;			/* start diastolic interval*/
						}
					}
				}
				
			}}}


		}
#elif PIPELINE
	 /* Now, unfortunately for explicit solutions, we will have to do 8 pde steps for every CA step */
	 /* Now doing two pde steps per loop iteration to remove variable swapping */
	 	for (pm=0; pm <= (PM_MAX); pm++) {
			printf("PDE iteration: %d\n", pm);

			for (i=0; i<= zBounds.end + 2; i++) {
				//printf("i = %i", i);
				i2 = i-1;
				grid_i = grid[i];
				if (i > 0) {
					pipeline = 1;
					grid_i2= grid[i2];
					xStart = MIN (xBounds[i].start,xBounds[i-1].start+1);
					xStart = MAX (xStart, 1);
					xEnd = MAX(xBounds[i].end,xBounds[i-1].end+1);
				/* ECG CODE	Calculate ECG on last iteration of PDE loop */
					if (pm == PM_MAX) {
						z = -i2*SPACING;
						dz2 = (z - LEAD_Z)*(z - LEAD_Z);
						V3_dz2 = (z - V3_LEAD_Z)*(z - V3_LEAD_Z);
					}
				} else {
					xStart = xBounds[i].start;
					xEnd = xBounds[i].end;
					pipeline = 0;
				}
				yBounds_i = yBounds[i];

			for (j= xStart; j<=xEnd; j++) {
				j2 = j-1;
				grid_i_j = grid_i[j];
				yBounds_i_j = yBounds_i[j];
				if (pipeline) {
					grid_i2_j2 = grid_i2[j2];
					yStart = MIN (yBounds_i_j.start, yBounds[i-1][j-1].start+1);
					yStart = MAX (yStart, 1);
					yEnd = MAX(yBounds_i_j.end,yBounds[i-1][j-1].end+1);
					
					 /* ECG CODE*/
					if (pm == PM_MAX) {
						x = -5.5 + j2*SPACING;
						dx2 = (x - LEAD_X)*(x - LEAD_X);
						V3_dx2 = (x - V3_LEAD_X)*(x - V3_LEAD_X);
					}
				} else {
					yStart = yBounds_i_j.start;
					yEnd = yBounds_i_j.end;
				}
				
			for (k=yStart; k<=yEnd; k++) {
				k2 = k-1;

			/* Make sure we are in the heart */
				if (grid_i_j[k]) {
				
					grid_i_j_k = grid_i_j[k];	

				 /* Local potential and restitution parameter */
					u_0 = grid_i_j_k->u;
				//	bet = grid_i_j_k->bet;
						
				 /* Get neighboring potentials, Implement no flux as required */
					uxyz_m = 0.0;
					if (!(grid_i_j_k->flags BOUNDARY_TEST)) {
						uz_p = grid[i+1][j][k]->u - u_0;
						grid[i+1][j][k]->tmp += -uz_p;
						ux_p = grid_i[j+1][k]->u - u_0;
						grid_i[j+1][k]->tmp += -ux_p;
						uy_p = grid_i_j[k+1]->u - u_0;
						grid_i_j[k+1]->tmp += -uy_p;	
					} else {
						if (grid_i_j[k+1]) {
							uy_p = grid_i_j[k+1]->u - u_0;
							grid_i_j[k+1]->tmp += -uy_p;
						} else
							uy_p = 0.0;
						if (grid_i[j+1][k]) {
							ux_p = grid_i[j+1][k]->u - u_0;
							grid_i[j+1][k]->tmp += -ux_p;
						} else
							ux_p = 0.0;
						if (grid[i+1][j][k]) {
							uz_p = grid[i+1][j][k]->u - u_0;
							grid[i+1][j][k]->tmp += -uz_p;
						} else
							uz_p = 0.0;
					}
			
					uxyz_m = grid_i_j_k->tmp;

				 /* Compute 2nd differences for (isotropic) laplacian */
				//	d2x = ux_p + ux_m - 2*u_0;	
				//	d2y = uy_p + uy_m - 2*u_0;	
				//	d2z = uz_p + uz_m - 2*u_0;	
				
				 /* Set current depending on whether exc or rep */
					if (grid_i_j_k->flags EXC_TEST)
						i_0 = S2*(1.0 - u_0);
					else
						i_0 = 2*grid_i_j_k->bet*u_0*(u_0 - 1.0);
				
				 /* Excplicit integration */
					//gpde_i_j_k->tmp = u_0 + FACT*gpde_i_j_k->dewnsud + DTP*i_0;
					grid_i_j_k->tmp = MAX((MIN((u_0 + FACT* (ux_p + uy_p + uz_p + uxyz_m) + DTP*i_0), UMAX)),UMIN);
					grid_i_j_k->u = 0.0;
				}

				if ( pipeline && (grid_i2_j2[k2])) {

					grid_i2_j2_k2 = grid_i2_j2[k2];	

				 /* Local potential and restitution parameter */
					u_0 = grid_i2_j2_k2->tmp;
				 //bet = grid_i_j_k->bet;
						
				 /* Get neighboring potentials, Implement no flux as required */
					uxyz_m = 0.0;
					if (!(grid_i2_j2_k2->flags BOUNDARY_TEST)) {
						//uz_p = grid[i2+1][j2][k2]->tmp - u_0;
						uz_p = grid_i[j2][k2]->tmp - u_0;
						//grid[i2+1][j2][k2]->u += -uz_p;
						grid_i[j2][k2]->u += -uz_p;
						//ux_p = grid_i2[j2+1][k2]->tmp - u_0;
						ux_p = grid_i2[j][k2]->tmp - u_0;
						//grid_i2[j2+1][k2]->u += -ux_p;
						grid_i2[j][k2]->u += -ux_p;
						//uy_p = grid_i2_j2[k2+1]->tmp - u_0;
						uy_p = grid_i2_j2[k]->tmp - u_0;
						//grid_i2_j2[k2+1]->u += -uy_p;	
						grid_i2_j2[k]->u += -uy_p;	
					} else {
						//k2+1 -> k
						if (grid_i2_j2[k]) {
							uy_p = grid_i2_j2[k]->tmp - u_0;
							grid_i2_j2[k]->u += -uy_p;
						} else
							uy_p = 0.0;
						//j2+1 -> j
						if (grid_i2[j][k2]) {
							ux_p = grid_i2[j][k2]->tmp - u_0;
							grid_i2[j][k2]->u += -ux_p;
						} else
							ux_p = 0.0;
						//grid[i2+1] -> grid_i
						if (grid_i[j2][k2]) {
							uz_p = grid_i[j2][k2]->tmp - u_0;
							grid_i[j2][k2]->u += -uz_p;
						} else
							uz_p = 0.0;
					}
			
					uxyz_m = grid_i2_j2_k2->u;

				 /* Compute 2nd differences for (isotropic) laplacian */
				//	d2x = ux_p + ux_m - 2*u_0;	
				//	d2y = uy_p + uy_m - 2*u_0;	
				//	d2z = uz_p + uz_m - 2*u_0;	
				
				 /* Set current depending on whether exc or rep */
					if (grid_i2_j2_k2->flags EXC_TEST)
						i_0 = S2*(1.0 - u_0);
					else
						i_0 = 2*grid_i2_j2_k2->bet*u_0*(u_0 - 1.0);
				
				 /* Excplicit integration */
					//gpde_i_j_k->tmp = u_0 + FACT*gpde_i_j_k->dewnsud + DTP*i_0;
					grid_i2_j2_k2->u = MAX((MIN((u_0 + FACT* (ux_p + uy_p + uz_p + uxyz_m) + DTP*i_0), UMAX)),UMIN);
					grid_i2_j2_k2->tmp = 0.0;

					if (grid_i2_j2_k2->u > UMIN)
						numACT++;

				/* ECG CALCULATIONS*/
					if (pm ==PM_MAX) {
						y = -5.5 + k2*SPACING;
						dy2 = (y - LEAD_Y)*(y - LEAD_Y);
						V3_dy2 = (y - V3_LEAD_Y)*(y - V3_LEAD_Y);
						ecgLead += (ux_p + uy_p + uz_p + uxyz_m)/sqrt(dx2 + dy2 + dz2);
						V3_ecgLead += (ux_p + uy_p + uz_p + uxyz_m)/sqrt(V3_dx2 + V3_dy2 + V3_dz2);

							/* Check for EXC -> REF */
				 	if (grid_i2_j2_k2->flags EXC_TEST) {
						
					 /* Check whether we are still in the front region for element */
				 		if (grid_i2_j2_k2->nsr == N_EXCS) {
							grid_i2_j2_k2->flags EXC_FALSE;				/* Cease CA and PDE sourcing */
			 
							u_0 = grid_i2_j2_k2->u;

							/* Compute new I-V parameters from peak potential for apd */
							if (grid_i2_j2_k2->ischemic) {
								if(grid_i2_j2_k2->flags EPI_TEST) {
									apd = ischemic_apd_from_restitution (grid_i2_j2_k2->tau, EPICARDIAL_CELL); 
								} else if (grid_i2_j2_k2->flags END_TEST) {
									apd = ischemic_apd_from_restitution (grid_i2_j2_k2->tau, ENDOCARDIAL_CELL); 
								} else {
									apd = ischemic_apd_from_restitution (grid_i2_j2_k2->tau, MIDMYOCARDIAL_CELL); 
								}
							} else {
								if(grid_i2_j2_k2->flags EPI_TEST) {
									apd = apd_from_restitution (grid_i2_j2_k2->tau, EPICARDIAL_CELL); 
								} else if (grid_i2_j2_k2->flags END_TEST) {
									apd = apd_from_restitution (grid_i2_j2_k2->tau, ENDOCARDIAL_CELL); 
								} else {
									apd = apd_from_restitution (grid_i2_j2_k2->tau, MIDMYOCARDIAL_CELL); 
								}
							}
						/*
							if(grid_i2_j2_k2->flags EPI_TEST) {
								apd = apd_from_restitution (grid_i2_j2_k2->tau, EPICARDIAL_CELL); 

							} else if (grid_i2_j2_k2->flags END_TEST) {
								apd = apd_from_restitution (grid_i2_j2_k2->tau, ENDOCARDIAL_CELL); 
								
							} else {
								apd = apd_from_restitution (grid_i2_j2_k2->tau, MIDMYOCARDIAL_CELL); 
							}
							*/
						   compute_apd_parameters (apd*0.001, grid_i2_j2_k2->u, &tsh, 
								&(grid_i2_j2_k2->bet));

					 /* Assign new source strength at next autochrone based on sourcing table */
						} else {	
							/* This is a NOOP for now, will change with ->nsr; recovery effects */
							grid_i2_j2_k2->src = 1;
						}		
					//continue;	/* Loop up */
					} else {
					/* Update diastolic interval for repolarized elements */
						if (grid_i2_j2_k2->flags REP_TEST) {
							grid_i2_j2_k2->tau += 2;		/* CA time steps in 2 milliseconds, restitution curves in ms */
						 /* Check for REF-> Excitable */
							if(grid_i2_j2_k2->flags REF_TEST) {
								if(grid_i2_j2_k2->flags EPI_TEST) {
									if( grid_i2_j2_k2->tau >= EPI_DMN) 
										grid_i2_j2_k2->flags REF_FALSE;
								} else if (grid_i2_j2_k2->flags END_TEST) {
									if ( grid_i2_j2_k2->tau >= END_DMN)
										grid_i2_j2_k2->flags REF_FALSE;
								} else {
									if( grid_i2_j2_k2->tau >= MID_DMN) 
										grid_i2_j2_k2->flags REF_FALSE; 
								}
							}
							/*
							if (grid_i2_j2_k2->flags REF_TEST && (grid_i2_j2_k2->tau >= AP_DMN)) {
								grid_i2_j2_k2->flags REF_FALSE;
							}
							*/
						} else {
						/* Check for repolarization via PDEs */
							if (grid_i2_j2_k2->flags REF_TEST &&  (grid_i2_j2_k2->u <= UAPD) ) {
								grid_i2_j2_k2->flags REP_TRUE;		/* end of action potential */
								grid_i2_j2_k2->tau = 0;			/* start diastolic interval*/
							}
						}
					}
	
					}
				}
		
			}	//End k loop
			}	//End j loop
			}	//End i loop
		} // End pm 0 -> PM_MAX loop
#endif
		}

	//*			
	//if( (tm%4)==0)
	//	write_Compressed_3D_data(tm);
	//*/
	//write_3D_data(tm);
	//write_Compressed_3D_data(tm);
	//write_2D_data(tm);


	} /* End simulation main loop */

	
	stop = time(NULL); 
	diff = difftime(stop, start); 
	printf("Difference is %e\n", diff/60); 
	
	printf("%c", BELL);
	fclose(ecg);
	
	//HC_Pause();
#ifdef _WIN32
#ifdef BZ_IMPORT
   BZ2DLLFreeLibrary();
#endif
#endif
	return (0);
}

static double atanch_andy (double val) {
	return (0.5*log((1+val)/(1-val)));
}
		
static void compute_apd_parameters (float apd, float umax, float *tsh, float *bet) {
	float v0 = UAPD;
	float vm;
	
	vm = umax;
	*bet = (float) ((1/apd) * (atanch_andy(2*vm-1) - atanch_andy(2*v0-1)));
	*tsh = (float) (apd * (1/(1 - atanch_andy(2*v0-1)/atanch_andy(2*vm-1))));
}

/*
static void compute_apd_parameters (float apd, float umax, float *tsh, float *bet) {
	//float u0 = UAPD;
	//float um;
	
	//um   = umax;
	*bet = (float) ((1/apd) * ( (0.5*log(umax/(1-umax))) + 1.472219489583));
	*tsh = (float) (apd * (1/(1 +(1.472219489583)/(0.5*log(umax/(1-umax)))   )));
}
*/
/* ORIGINAL
static float apd_from_restitution (unsigned short tau) {
	double alpha = 1.2;
	double scale = 0.3;
	double b, apd;
	double t;
	
	t = 0.001*tau;
	b = scale*alpha;
	apd = 0.5*(b*b + sqrt(b*b*b*b + 4*b*b*t));
	return ((float) apd*1000);
}
*/
/*
static float apd_from_restitution (unsigned short tau) {
	//double alpha = 1.2;
	//double scale = 0.3;
	//double b;
	//double apd;
	//double t;
	
	//t = 0.001*tau;
	//b = scale*alpha;
	return (float) 500*(0.1296 + sqrt(0.01679616 + 0.5184*0.001*tau));

}
*/
static VOXEL_TABLE *voxel_get_random_table (unsigned char i, unsigned char j) {
	int k = (int) (fabs(N_TABS*ran1(&iSeed)-epsilon)); 
	return (&(tableArray[i][j][0]));
}	


static int voxel_get_table_entry (VOXEL_TABLE *table, int i, int j, int k) {
	unsigned short mask = bitMasks[k];
	unsigned short data = table->bits[i][j];
	return ((int) ((mask & data) ? 1: 0));
}	

void voxel_init_random_table_array (void) {
	int i, j, k;
	int the_d, phi_d;
	double the_r, phi_r;
	FILE *fp;
	VOXEL_TABLE *table;
	
	printf("Allocating neighborhood tables ...\n");
		
 /* Allocate voxel table Array - calloc should init to all bits to zeros */
	tableArray = (VOXEL_TABLE ***) calloc(T_GRID+1, sizeof(VOXEL_TABLE**));
	if (tableArray == NULL) {
		printf("Insufficient Dynamic Memory\n");
		exit(1);
	}
	for (i=0; i<=T_GRID; i++) {
		tableArray[i] = (VOXEL_TABLE **) calloc(P_GRID+1, sizeof(VOXEL_TABLE *));
		if (tableArray[i] == NULL) {
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}
	for (i=0; i<=T_GRID; i++) {
	for (j=0; j<=P_GRID; j++) {
		tableArray[i][j] = (VOXEL_TABLE *) calloc(N_TABS, sizeof(VOXEL_TABLE));
		if (tableArray[i][j] == NULL) {	
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}}			
	
 /* Fill voxel tables for 3D fiber orientations in increments of 2 degrees in theta, phi*/
	for (k=0; k<N_TABS; k++) {
		voxel_build_random_table(&(tableArray[0][0][k]), ANISOTROPY, 0.0, 0.0);
	}	 
}


void voxel_build_random_table (VOXEL_TABLE *table, double speed_ratio, double the, double phi)
{
	int i, j, k;
	double t, p;
	double v1_x, v1_y, v1_z, v1_n;
	double v2_x, v2_y, v2_z, v2_n;
	#if 0
	double v3_x, v3_y, v3_z, v3_n;
	#endif
	double rot[3][3];
	
	double a, b, r;
	double x, y, z;
	double xp, yp, zp;
	double c, s, u, T;
	
	unsigned short mask = 0;
	double dh;
	double min = -1.0; 
	double max =  1.0;
	
	t = the;
	p = phi;
	
 /* Fill voxel array with randomized seed point locations */
 	dh = (max-min)/15.0;
 	for (i=0; i<15; i++) {
 	for (j=0; j<15; j++) {
 	for (k=0; k<15; k++) {
 		rvoxels[i][j][k].x = min + dh*(i + ran1(&iSeed));
		rvoxels[i][j][k].y = min + dh*(j + ran1(&iSeed));
		rvoxels[i][j][k].z = min + dh*(k + ran1(&iSeed));
	}}}
	
 /* Prinicipal axis vector for fiber */
	v1_n = 1.0;
	v1_x = sin(t)*cos(p);
	v1_y = sin(t)*sin(p);
	v1_z = cos(t);
	
 /* Gram Shmidt Orthogonalization (v1 x X) */
	v2_x = 0.0;
	v2_y = -cos(t);
	v2_z = sin(t)*sin(p);
	v2_n = sqrt(v2_y*v2_y + v2_z*v2_z);
	v2_y /= v2_n;
	v2_z /= v2_n;
	
 /* (v1 x v2) */
 	#if 0	/* This 3rd axis vector is not needed to compute transformation */
	v3_x = -cos(t)*cos(t) - sin(t)*sin(t)*sin(p)*cos(p);
	v3_y = sin(t)*sin(t)*cos(p)*cos(p);
	v3_z = cos(t)*sin(t)*cos(p);
	v3_n = sqrt(v3_x*v3_x + v3_y*v3_y + v3_z*v3_z);
	v3_x /= v3_n;
	v3_y /= v3_n;
	v3_z /= v3_n;
	#endif
	
 /* Compute rotation matrix for points (angle T about v2) */	
	T = acos(v1_x);
	c = cos(T);
	s = sin(T);
	u = 1-c;
	rot[0][0] = u*v2_x*v2_x + c;
	rot[0][1] = u*v2_x*v2_y + s*v2_z;
	rot[0][2] = u*v2_x*v2_z - s*v2_y;
	rot[1][0] = u*v2_x*v2_y - s*v2_z;
	rot[1][1] = u*v2_y*v2_y + c;
	rot[1][2] = u*v2_y*v2_z + s*v2_x;
	rot[2][0] = u*v2_x*v2_z + s*v2_y;
	rot[2][1] = u*v2_y*v2_z - s*v2_x;
	rot[2][2] = u*v2_z*v2_z + c;
 	
 /* Set ellipsoid parameters */
 	a = R_D/7.5;
 	b = speed_ratio*a;	
 
 /* Loop over voxels and calculate inside/outside rotated ellipsoid */
 	for (i=0; i<15; i++) {
 	for (j=0; j<15; j++) {
 	for (k=0; k<15; k++) {
 		x = rvoxels[i][j][k].x;
 		y = rvoxels[i][j][k].y;
 		z = rvoxels[i][j][k].z;
 		xp = rot[0][0]*x + rot[0][1]*y + rot[0][2]*z;
  		yp = rot[1][0]*x + rot[1][1]*y + rot[1][2]*z;
  		zp = rot[2][0]*x + rot[2][1]*y + rot[2][2]*z;
 		r = xp*xp/(a*a) + (yp*yp + zp*zp)/(b*b);
 		if (r <= 1.0) {
  			table->bits[i][j] |= bitMasks[k];	/* Set correct bit in voxel table */
 		}
 	}}}
}

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349
            
static float ran1(int *idum)
{
	static long ix1,ix2,ix3;
	static float r[98];
	float temp;
	static int iff=0;
	int j;
            
	if (*idum < 0 || iff == 0) {
		iff=1;
		ix1=(IC1-(*idum)) % M1;
		ix1=(IA1*ix1+IC1) % M1;
		ix2=ix1 % M2;
		ix1=(IA1*ix1+IC1) % M1;
		ix3=ix1 % M3;
		for (j=1;j<=97;j++) {
			ix1=(IA1*ix1+IC1) % M1;
			ix2=(IA2*ix2+IC2) % M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		}
		*idum=1;
	}
	ix1=(IA1*ix1+IC1) % M1;
	ix2=(IA2*ix2+IC2) % M2;
	ix3=(IA3*ix3+IC3) % M3;
	j=1 + ((97*ix3)/M3);
	if (j > 97 || j < 1) fprintf (stderr, "RAN1: This cannot happen.");
	temp=r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
}

#define A_0		1.1
#define ALPHA	0.1
#define B_0		0.7
#define BETA	0.12
#define C_0		0.7
#define GAMMA	0.1

#define R_A		9.8
#define R_B		5.4

double x_ctr_LV, x_ctr_RV;
double r_LV, r_A, r_B, r_C;
double r_RV, rR_B, rR_C;
double lv_thick, sp_thick, rv_thick;

int Inside_Outer_Ellipse (double x, double y)
{
	return ((x*x/(r_B*r_B) + y*y/(r_C*r_C)) <= 1.0);
}

int Inside_Left_Ventricle (double x, double y)
{
	return (((x-x_ctr_LV)*(x-x_ctr_LV) + y*y) <= r_LV*r_LV);
}

int Inside_Right_Ventricle (double x, double y)
{
	return ((x < x_ctr_LV) && ((x*x/(rR_B*rR_B) + y*y/(rR_C*rR_C)) <= 1.0) &&
			!((((x-x_ctr_RV)*(x-x_ctr_RV)/(r_RV*r_RV) + y*y/(2*r_RV*r_RV)) <= 1.0)));
}

void SetNewZValue (double z) 
{
	r_B  = R_B*sqrt(1 - z*z/(R_A*R_A));	
	r_C  = 0.80*r_B;
	sp_thick = B_0*exp(fabs(z)*BETA);
	rv_thick = C_0*exp(fabs(z)*GAMMA);
	lv_thick = A_0*exp(fabs(z)*ALPHA);
	r_LV = MAX((0.75*r_B - lv_thick), 0.0);
	r_RV = r_LV + sp_thick;
	rR_B = MAX((r_B-rv_thick), 0.0);
	rR_C = MAX((r_C-rv_thick), 0.0);
	x_ctr_LV = 0.25*r_B;
	x_ctr_RV = x_ctr_LV;
}

void CalculateVolumes (void)
{
	double x, y, z;
	long nRV, nLV, nMu;
	int InsideOuterEllipse;
	int InsideLeftVentricle;
	int InsideRightVentricle;
	
	nRV = 0;
	nLV = 0;
	nMu = 0;
	for (z=0.0; z <= 9.5; z += SPACING) {
		SetNewZValue (z);
	
		for (x = -5.5; x <= 5.5; x += SPACING) {
		for (y = -5.5; y <= 5.5; y += SPACING) {

			InsideOuterEllipse   = Inside_Outer_Ellipse   (x, y);
			InsideLeftVentricle  = Inside_Left_Ventricle  (x, y);
			InsideRightVentricle = Inside_Right_Ventricle (x, y); 
			
			if (InsideLeftVentricle)  nLV++;
			if (InsideRightVentricle) nRV++;
			if (InsideOuterEllipse && !(InsideLeftVentricle) && !(InsideRightVentricle)) nMu++;

		}}
	}
	printf("nLV = %ld, nRV = %ld, nMu = %ld\n", nLV, nRV, nMu);
}

/* Compute normal cell apd from prior diastolic interval (restitution function) */
static float apd_from_restitution (unsigned short tau, int cellType)
{
	float apd; 
	unsigned short di;
	
	if (cellType == ENDOCARDIAL_CELL) {
		di = MIN(tau, END_DMX);
		if (di < END_DMN) {
			apd = END_MIN;
		} else if (di < END_DMD) {
			apd = END_MIN + apd_end_slp*(di-END_DMN);
		} else {
			apd = END_MID + (END_MAX-END_MID)*(1 - exp(-apd_end_exp*(di - END_DMD)));
		}
	
	} else if (cellType == EPICARDIAL_CELL ) {

		di = MIN(tau, EPI_DMX);
		if (di < EPI_DMN) {
			apd = EPI_MIN;
		} else if (di < EPI_DMD) {
			apd = EPI_MIN + apd_epi_slp*(di-EPI_DMN);
		} else {
			apd = EPI_MID + (EPI_MAX-EPI_MID)*(1 - exp(-apd_epi_exp*(di - EPI_DMD)));
		}
	
	} else { /* Mid myocardial cell */
		
		di = MIN(tau, MID_DMX);
		if (di < MID_DMN) {
			apd = MID_MIN;
		} else if (di < MID_DMD) {
			apd = MID_MIN + apd_mid_slp*(di-MID_DMN);
		} else {
			apd = MID_MID + (MID_MAX-MID_MID)*(1 - exp(-apd_mid_exp*(di - MID_DMD)));
		}
	}
	return (apd);
}

/* Compute normal cell apd from prior diastolic interval (restitution function) */
static float ischemic_apd_from_restitution (unsigned short tau, int cellType)
{
	float apd; 
	unsigned short di;
	
	if (cellType == ENDOCARDIAL_CELL) {
		di = MIN(tau, ISCHEMIC_END_DMX);
		if (di < ISCHEMIC_END_DMN) {
			apd = ISCHEMIC_END_MIN;
		} else if (di < ISCHEMIC_END_DMD) {
			apd = ISCHEMIC_END_MIN + apd_end_slp*(di-ISCHEMIC_END_DMN);
		} else {
			apd = ISCHEMIC_END_MID + (ISCHEMIC_END_MAX-ISCHEMIC_END_MID)*(1 - exp(-apd_end_exp*(di - ISCHEMIC_END_DMD)));
		}
	
	} else if (cellType == EPICARDIAL_CELL ) {

		di = MIN(tau, ISCHEMIC_EPI_DMX);
		if (di < ISCHEMIC_EPI_DMN) {
			apd = ISCHEMIC_EPI_MIN;
		} else if (di < ISCHEMIC_EPI_DMD) {
			apd = ISCHEMIC_EPI_MIN + apd_epi_slp*(di-ISCHEMIC_EPI_DMN);
		} else {
			apd = ISCHEMIC_EPI_MID + (ISCHEMIC_EPI_MAX-ISCHEMIC_EPI_MID)*(1 - exp(-apd_epi_exp*(di - ISCHEMIC_EPI_DMD)));
		}
	
	} else { /* Mid myocardial cell */
		
		di = MIN(tau, ISCHEMIC_MID_DMX);
		if (di < ISCHEMIC_MID_DMN) {
			apd = ISCHEMIC_MID_MIN;
		} else if (di < ISCHEMIC_MID_DMD) {
			apd = ISCHEMIC_MID_MIN + apd_mid_slp*(di-ISCHEMIC_MID_DMN);
		} else {
			apd = ISCHEMIC_MID_MID + (ISCHEMIC_MID_MAX-ISCHEMIC_MID_MID)*(1 - exp(-apd_mid_exp*(di - ISCHEMIC_MID_DMD)));
		}
	}
	return (apd);
}
////////////////////////////////////////////
//////////INITIALIZE HEART GEOMETRY/////////
////////////////////////////////////////////

#if HUMAN

void initialize_apds ( void )
{	int i, j, k;
	float bet, apd, tsh;

/* Set up PDE restitution curve parameters endo */
	apd_end_slp = (END_MID - END_MIN)/(END_DMD - END_DMN);
	apd_end_exp = apd_end_slp/(END_MAX - END_MID);

/* Set up PDE restitution curve parameters epi */
	apd_epi_slp = (EPI_MID - EPI_MIN)/(EPI_DMD - EPI_DMN);
	apd_epi_exp = apd_epi_slp/(EPI_MAX - EPI_MID);

/* Set up PDE restitution curve parameters mid */
	apd_mid_slp = (MID_MID - MID_MIN)/(MID_DMD - MID_DMN);
	apd_mid_exp = apd_mid_slp/(MID_MAX - MID_MID);

  /* Loop over all labelled cells */
  /* Get APDs from restitution based on cell type */
  /* Assign repolarization current source parameters cell by cell */
	for (i=0; i<Z_GRID; i++) {
		grid_i = grid[i];
		for (j=0; j<X_GRID; j++) {
			grid_i_j= grid_i[j];
			for (k=0; k<Y_GRID; k++) {
				grid_i_j_k = grid_i_j[k];
				if (grid_i_j_k) {
					if(grid_i_j_k->flags EPI_TEST) {
							apd = apd_from_restitution (225, EPICARDIAL_CELL); 
					} else if (grid_i_j_k->flags END_TEST) {
							apd = apd_from_restitution (225, ENDOCARDIAL_CELL); 
					} else {
							apd = apd_from_restitution (225, MIDMYOCARDIAL_CELL); 
					}
					compute_apd_parameters (apd*0.001, UMAX, &tsh, &bet);
					grid_i_j_k->bet = bet;
				}
	}}}
}

void initialize_ischemic_region(void) {
	int i, j, k;
	float bet, apd, tsh;
	int nIschemic = 0;

/* Set up PDE restitution curve parameters endo */
	ischemic_apd_end_slp = (ISCHEMIC_END_MID - ISCHEMIC_END_MIN)/(ISCHEMIC_END_DMD - ISCHEMIC_END_DMN);
	ischemic_apd_end_exp = ischemic_apd_end_slp/(ISCHEMIC_END_MAX - ISCHEMIC_END_MID);

/* Set up PDE restitution curve parameters epi */
	ischemic_apd_epi_slp = (ISCHEMIC_EPI_MID - ISCHEMIC_EPI_MIN)/(ISCHEMIC_EPI_DMD - ISCHEMIC_EPI_DMN);
	ischemic_apd_epi_exp = ischemic_apd_epi_slp/(ISCHEMIC_EPI_MAX - ISCHEMIC_EPI_MID);

/* Set up PDE restitution curve parameters mid */
	ischemic_apd_mid_slp = (ISCHEMIC_MID_MID - ISCHEMIC_MID_MIN)/(ISCHEMIC_MID_DMD - ISCHEMIC_MID_DMN);
	ischemic_apd_mid_exp = ischemic_apd_mid_slp/(ISCHEMIC_MID_MAX - ISCHEMIC_MID_MID);


	/* Initialize the state of all CAM_PDE elements */
	for (i=(Z_GRID/2-50); i<(Z_GRID/2+50); i++) {
		grid_i=grid[i];
	for (k=(X_GRID/2-30); k<(X_GRID/2+30); k++) {
		//grid_i_j=grid_i[j];
	for (j=(Y_GRID/3*2); j<(Y_GRID); j++) {
		//grid_i_j_k = grid_i_j[k];
		grid_i_j_k = grid_i[j][k];
		if (grid_i_j_k) {
			grid_i_j_k->ischemic = 1;
			if(grid_i_j_k->flags EPI_TEST) {
				apd = ischemic_apd_from_restitution (225, EPICARDIAL_CELL); 
			} else if (grid_i_j_k->flags END_TEST) {
				apd = ischemic_apd_from_restitution (225, ENDOCARDIAL_CELL); 
			} else {
				apd = ischemic_apd_from_restitution (225, MIDMYOCARDIAL_CELL); 
			}
			compute_apd_parameters (apd*0.001, UMAX, &tsh, &bet);
			grid_i_j_k->bet = bet;
			nIschemic++;
		}
	}}}

	printf("Number of Ischemic Cells = %i\n",nIschemic);
}

void initialize_geometry ( void ) {
	int i, j, k;
	int nCellularAutomata = 0;
	double R, x, y, z;
	double rv_thick;
	int InsideOuterEllipse;
	int InsideLeftVentricle;
	int InsideRightVentricle;
	float apd, tsh, bet;
	CAM_PDE_ELEMENT *ptr;
	
	//printf("Calculating Volumes...\n");
 	CalculateVolumes();

	printf("Initializing Geometry ...\n");
			
	/* Allocate voxel table Array - calloc should init to all bits to zeros */
	grid = (CAM_PDE_ELEMENT ****) calloc(Z_GRID, sizeof(CAM_PDE_ELEMENT ***));
	if (grid  == NULL) {
		printf("Insufficient Dynamic Memory\n");
		exit(1);
	}
	for (i=0; i< Z_GRID; i++) {
		grid[i] = (CAM_PDE_ELEMENT ***) calloc(X_GRID, sizeof(CAM_PDE_ELEMENT **));
		if (grid[i] == NULL) {
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}
	for (i=0; i< Z_GRID; i++) {
	for (j=0; j< X_GRID; j++) {
		grid[i][j] = (CAM_PDE_ELEMENT **) calloc(Y_GRID, sizeof(CAM_PDE_ELEMENT*));
		if (grid[i][j] == NULL) {	
			printf("Insuffifcient Dynamic Memory\n");
			exit(1);
		}
	}}	

	/* Compute inside/outside heart field for elements */
	for (i=0; i<Z_GRID; i++) {
		z = -i*SPACING;
		SetNewZValue (z);
		grid_i = grid[i];
		
		for (j=0; j<X_GRID; j++) {
			x = -5.5 + j*SPACING;
			grid_i_j= grid_i[j];

			for (k=0; k<Y_GRID; k++) {
				y = -5.5 + k*SPACING;
				
				InsideOuterEllipse   = Inside_Outer_Ellipse   (x, y);
				InsideLeftVentricle  = Inside_Left_Ventricle  (x, y);
				InsideRightVentricle = Inside_Right_Ventricle (x, y); 
		
				if ((i<Z_MAX) && InsideOuterEllipse && !(InsideLeftVentricle) && !(InsideRightVentricle)) {
					grid_i_j[k] = (CAM_PDE_ELEMENT *)1;
					nCellularAutomata++;
				} else {
					grid_i_j[k] = NULL;
				}
	}}}

	printf("Number of Cellular Automata . . . %i \n",nCellularAutomata);

	/*insert 3d array to allocate space for cam_pde*/
	ptr = (CAM_PDE_ELEMENT *) calloc(nCellularAutomata, sizeof(CAM_PDE_ELEMENT));
	if(!ptr) {
		printf("not enough memory for CAM_PDE_ELEMENT array\n");
		exit(1);
	}
	
	/* Get IV-Curve characteristics and add to elements */
	/* Moved to a new function
	apd = apd_from_restitution ((unsigned short) (DI_INIT)); 		
	apd_parameters (apd*0.001, &tsh, &bet);	
	printf("Initial APD = %f\n", apd);
	*/
	printf("Initializing array element\n");

	/* Initialize the state of all CAM_PDE elements */
	for (i=0; i<Z_GRID; i++) {
		grid_i=grid[i];
	for (j=0; j<X_GRID; j++) {
		grid_i_j=grid_i[j];
	for (k=0; k<Y_GRID; k++) {
		if (grid_i_j[k]) {
			ptr->ischemic = 0;
			ptr->flags = 0;
			ptr->src = 1;
			ptr->nsr = 0;
			ptr->snk = 0;
			ptr->tau = DI_INIT;
				
			ptr->the_bin = 0;
			ptr->phi_bin = 0;
			
			ptr->u   = UMIN;	/* Potential value 		*/
			ptr->tmp = 0.0;		/* Temp potential value	 	*/
			ptr->bet = bet;		/* APD Restitution parameter 	*/
	
			grid_i_j[k] = ptr;
			ptr++;
		}
	}}}
	 
}

#elif CANINE
void initialize_geometry_canine ( void )
{
	int i, j, k;
	int valid;
	int nCellularAutomata = 0;
	FILE *fp_canine;
	float apd, tsh, bet;
	CAM_PDE_ELEMENT *ptr;
	
	printf("Initializing Geometry ...\n");
			
	/* Allocate voxel table Array - calloc should init to all bits to zeros */
	grid = (CAM_PDE_ELEMENT ****) calloc(Z_GRID, sizeof(CAM_PDE_ELEMENT ***));
	if (grid  == NULL) {
		printf("Insufficient Dynamic Memory\n");
		exit(1);
	}
	for (i=0; i< Z_GRID; i++) {
		grid[i] = (CAM_PDE_ELEMENT ***) calloc(X_GRID, sizeof(CAM_PDE_ELEMENT **));
		if (grid[i] == NULL) {
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}
	for (i=0; i< Z_GRID; i++) {
	for (j=0; j< X_GRID; j++) {
		grid[i][j] = (CAM_PDE_ELEMENT **) calloc(Y_GRID, sizeof(CAM_PDE_ELEMENT*));
		if (grid[i][j] == NULL) {	
			printf("Insuffifcient Dynamic Memory\n");
			exit(1);
		}
	}}	

	if ( (fp_canine = fopen("3d_dog_heart_geometry.dat","rb")) == NULL)
		printf("Error reading canine geometry file\n");
	else {
		/* Read in heart geometry from external file*/
		
		for (i=0; i<Z_GRID; i++) {
			grid_i = grid[i];
			for (j=0; j<X_GRID; j++) {
				grid_i_j= grid_i[j];
				for (k=0; k<Y_GRID; k++) {
					fread(&valid, sizeof(int),1,fp_canine);
					if (valid) {
						grid_i_j[k] = (CAM_PDE_ELEMENT *)1;
						nCellularAutomata++;
					} else {
						grid_i_j[k] = NULL;
					}
		}}}	
	}
	printf("Number of Cellular Automata . . . %i \n",nCellularAutomata);

	/*insert 3d array to allocate space for cam_pde*/
	ptr = (CAM_PDE_ELEMENT *) calloc(nCellularAutomata, sizeof(CAM_PDE_ELEMENT));
	if(!ptr) {
		printf("not enough memory for CAM_PDE_ELEMENT array\n");
		exit(1);
	}
	 /* Get IV-Curve characteristics and add to elements */
	apd = apd_from_restitution ((unsigned short) (DI_INIT)); 		
	apd_parameters (apd*0.001, &tsh, &bet);	
	printf("Initial APD = %f\n", apd);
	printf("Initializing array element\n");

	/* Initialize the state of all CAM_PDE elements */
	for (i=0; i<Z_GRID; i++) {
		grid_i=grid[i];
	for (j=0; j<X_GRID; j++) {
		grid_i_j=grid_i[j];
	for (k=0; k<Y_GRID; k++) {
		if (grid_i_j[k]) {
			ptr->flags = 0;
			ptr->src = 1;
			ptr->nsr = 0;
			ptr->snk = 0;
			ptr->tau = DI_INIT;
				
			ptr->the_bin = 0;
			ptr->phi_bin = 0;
			
			ptr->u   = UMIN;	/* Potential value 		*/
			ptr->tmp = 0.0;	/* Temp potential value	 	*/
			ptr->bet = bet;		/* APD Restitution parameter 	*/
	
			grid_i_j[k] = ptr;
			ptr++;
		}
	}}}
	 
}

#endif

void set_boundary_conditions(void) {

	int i,j,k;
	int lastBoundary;
	int nBoundary = 0;

 /* Loop over CA grid elements and set boundary conditions */
	for (i=0; i<Z_GRID-1; i++) {
		grid_i = grid[i];
	for (j=1; j<X_GRID-1; j++) {
		grid_i_j = grid_i[j];
	for (k=1; k<Y_GRID-1; k++) {

		grid_i_j_k = grid_i_j[k];
		
	 	/* We are far enough away from array boundaries using our */
	 	/* grid, so don't bother checking on neighbor indices */
		if (grid_i_j_k) {
		
			grid_i_j_k->flags BOUNDARY_FALSE;
		 /* Check whether up (i-1) neighbor is a boundary element */
		 /* This is temporary and assumes first sheet has all boundary elements */
		 	if (i==0)
				grid_i_j_k->flags BOUNDARY_TRUE;
			else {
				if (!grid[i-1][j][k])
					grid_i_j_k->flags BOUNDARY_TRUE;
			}
			
		 /* Check whether dn (i+1) neighbor is a boundary element */
			if (!grid[i+1][j][k]) 
				grid_i_j_k->flags BOUNDARY_TRUE;

		 /* Check whether north (j-1) neighbor is a boundary element */
			if (!grid[i][j-1][k])
				grid_i_j_k->flags BOUNDARY_TRUE;

		 /* Check whether south (j+1) neighbor is a boundary element */
			if (!grid[i][j+1][k]) 
				grid_i_j_k->flags BOUNDARY_TRUE;

		 /* Check whether west (k-1) neighbor is a boundary element */
			if (!grid[i][j][k-1]) 
				grid_i_j_k->flags BOUNDARY_TRUE;
	
		 /* Check whether east (k+1) neighbor is a boundary element */
			if (!grid[i][j][k+1]) 
				grid_i_j_k->flags BOUNDARY_TRUE;

			if (grid_i_j_k->flags BOUNDARY_TEST){
				nBoundary++;
				grid_i_j_k->flags INT_TRUE;
			}
		}
	}}}
	printf("Number of boundary elements: %d\n", nBoundary);

	printf("Trying the negative z direction\n");
	for (k=1; k<Y_GRID-1; k++) {
	for (j=1; j<X_GRID-1; j++) {
		nBoundary = 0;
		for (i=(Z_GRID-1); i>=0; i--) {
		grid_i_j_k = grid[i][j][k];
		if (grid_i_j_k) {
			if (grid_i_j_k->flags BOUNDARY_TEST){
				nBoundary++;
				if(nBoundary == 1) {
					grid_i_j_k->flags INT_FALSE; 
					continue;
				} 
			} else 
				grid_i_j_k->flags INT_FALSE;
		}
	}}}


	printf("Determine internal ventricle surfaces\n");
	for (i=0; i<Z_GRID-1; i++) {
		grid_i = grid[i];
	for (j=1; j<X_GRID-1; j++) {
		grid_i_j = grid_i[j];
		nBoundary = 0;
		lastBoundary = -1;
	for (k=1; k<Y_GRID-1; k++) {

		grid_i_j_k = grid_i_j[k];
		if (grid_i_j_k) {
			if (grid_i_j_k->flags BOUNDARY_TEST){
			//if ( !(grid_i[j+1][k]) || !(grid_i[j-1][k]) || !(grid_i[j][k+1]) || !(grid_i[j][k-1])){
				nBoundary++;
				if(nBoundary == 1) {
					grid_i_j_k->flags INT_FALSE; 
					lastBoundary = k;
				} else
					lastBoundary = k;
			} else 
				grid_i_j_k->flags INT_FALSE;
		}
		
	}
	if (lastBoundary != -1)
		grid_i_j[lastBoundary]->flags INT_FALSE;
	}
	}

	printf("Trying the other direction\n");
	for (i=0; i<Z_GRID-1; i++) {
	for (j=1; j<X_GRID-1; j++) {
		nBoundary = 0;
		lastBoundary = -1;
	for (k=1; k<Y_GRID-1; k++) {
		grid_i_j_k = grid[i][k][j];
		if (grid_i_j_k) {
			if (grid_i_j_k->flags BOUNDARY_TEST){
				nBoundary++;
				if(nBoundary == 1) {
					grid_i_j_k->flags INT_FALSE; 
					lastBoundary = k;
				} else {
					lastBoundary = k;
				}
			} else {
				grid_i_j_k->flags INT_FALSE;
			}
		}
	}
	if (lastBoundary != -1)
		grid[i][lastBoundary][j]->flags INT_FALSE;
	}
	}

}	//set_boundary_conditions();

void differentiate_heart_tissue(void) {
	int i, j, k;
	int m,n;
	int outerWallStart = 0;
	int outerWallEnd = 0;
	int rvWallStart = 0;
	int rvWallEnd = 0;
	int lvWallStart = 0;
	int lvWallEnd = 0;
	int nBoundary = 0;
	int rvWallThickness = 0;
	int lvWallThickness = 0;
	int jThickness = 0;
	int septumCenter = 0;

	int septum_width = 0;
	int septum_wall = 0;
	int external_inner_surface_transition = 0;
	int inner_surface_external_transition = 0;

	int nMid =0;
	int nEpi =0;
	int nEnd =0;

	int radius = 0;

	printf("Set all heart tissue that is not a boundary to be midmyocardium tissue\n");
	for (i=0; i<Z_GRID-1; i++) {
		grid_i = grid[i];
	for (j=1; j<X_GRID-1; j++) {
		grid_i_j = grid_i[j];
	for (k=1; k<Y_GRID-1; k++) {

		grid_i_j_k = grid_i_j[k];
		if (grid_i_j_k)
			grid_i_j_k->flags MID_TRUE;	
	}}}

	//printf("Loop through the slices and differentiate the heart tissue types\n");
	for (i=0; i<Z_GRID-1; i++) {

		grid_i = grid[i];

		outerWallStart = 0;
		outerWallEnd = 0;
		rvWallStart = 0;
		rvWallEnd = 0;
		lvWallStart = 0;
		lvWallEnd = 0;
		nBoundary = 0;
		rvWallThickness = 0;
		lvWallThickness = 0;
		jThickness = 0;
		septumCenter = 0;


		//printf("Compute statistics of heart wall sizes for each slice\n");
		for (j=1; j<X_GRID-1; j++) {

			//grid_i_j_k = grid_i[110][j];
			grid_i_j_k = grid_i[j][Y_GRID/2];
			
			/*Am I in the heart?*/
			if (grid_i_j_k) {

				if (grid_i_j_k->flags BOUNDARY_TEST)
					nBoundary++;

				if ((grid_i_j_k->flags BOUNDARY_TEST) && (!(grid_i_j_k->flags INT_TEST))) {
					if (!outerWallStart)
						outerWallStart = j;
					else {
						outerWallEnd = j;
						//printf("outerWallEnd = %i\n",outerWallEnd);
						jThickness = outerWallEnd - outerWallStart;
						lvWallThickness = outerWallEnd - lvWallEnd;
					}
				}

				if (grid_i_j_k->flags INT_TEST) {
					if (!rvWallStart) 
						rvWallStart = j;
					else if (!rvWallEnd) {
						rvWallEnd = j;
						rvWallThickness = rvWallStart - outerWallStart;
					} else if (!lvWallStart) {
						lvWallStart = j;
						septumCenter = (lvWallStart - rvWallEnd)/2;
					} else // if (~lvWallEnd)
						lvWallEnd = j;
				}
			}  

		} /*j loop*/
		/*END compute statistics of heart wall sizes for each slice*/

		//printf("rvWallThickness = %i\n",rvWallThickness);
		//printf("lvWallThickness = %i\n",lvWallThickness);
		//printf("jThickness = %i\n",jThickness);

		//printf("Determine epicardium\n");
		radius =0;
		for (j=1; j<X_GRID-1; j++) {
			if (nBoundary == 6) {
                radius = rvWallThickness + (lvWallThickness - rvWallThickness)*(j-outerWallStart)/jThickness;
                radius = floor(radius/3);
			} else
				radius = 4;
			
			for (k=1; k<Y_GRID-1; k++) {

				grid_i_j_k = grid_i[j][k];

				if (grid_i_j_k) {
					if ((grid_i_j_k->flags BOUNDARY_TEST) && (!(grid_i_j_k->flags INT_TEST))) {
						for (m = (j-radius); m<=(j+radius); m++) {
						for (n = (k-radius); n<=(k+radius); n++) {
							if ( (m>0) && (n>0) && (m<X_GRID) && (n<Y_GRID) ) {
								if (grid_i[m][n]) {
									if (grid_i[m][n]->flags MID_TEST) {
										if (sqrt( (m-j)*(m-j) + (n-k)*(n-k)) < radius) {
											grid_i[m][n]->flags MID_FALSE;
											grid_i[m][n]->flags EPI_TRUE;
										}
									}
								}
							}
						}} /*m,n loop*/
					}
				}
		}} /*j,k loop*/
		/*End determine epicardium */

		//printf("Determine endocardium\n");
		radius =0;
		for (j=1; j<X_GRID-1; j++) {
			if (nBoundary == 6) {
                radius = rvWallThickness + (lvWallThickness - rvWallThickness)*(j-outerWallStart)/jThickness;
                radius = floor(radius/2);
				//changed from floor(radius/3) -> floor(radius/2) 1/8/04
			} else
				radius = 4;
			
			for (k=1; k<Y_GRID-1; k++) {

				grid_i_j_k = grid_i[j][k];

				if (grid_i_j_k) {
					if (grid_i_j_k->flags INT_TEST) {
						for (m = (j-radius); m<=(j+radius); m++) {
						for (n = (k-radius); n<=(k+radius); n++) {
							if ( (m>0) && (n>0) && (m<X_GRID) && (n<Y_GRID) ) {
								if (grid_i[m][n]) {
									if (grid_i[m][n]->flags MID_TEST) {
										if (sqrt( (m-j)*(m-j) + (n-k)*(n-k)) < radius) {
											grid_i[m][n]->flags MID_FALSE;
											grid_i[m][n]->flags END_TRUE;
										}
									}
								}
							}
						}} /*m,n loop*/
					}
				}
		}} /*j,k loop*/
		/*End determine endocardium */

		//printf("Set septum to be midmyocardium cells\n");
#if 0
		for (k=1; k < Y_GRID; k++) {
			septum_width = 0;
			septum_wall = 0;
			external_inner_surface_transition = 0;
			inner_surface_external_transition = 0;

			for (j = 1; j < X_GRID; j++) {
				if (grid_i[j][k]) {
					if ( (grid_i[j][k]->flags INT_TEST) && (!(grid_i[j-1][k]) ) )
						external_inner_surface_transition = j;
					
					if (external_inner_surface_transition)
						septum_wall++;
            
					if ( (grid_i[j][k]->flags INT_TEST) && (!grid_i[j+1][k]) && external_inner_surface_transition) {
						inner_surface_external_transition = j;
						j = 500;
					}
				}


			} /*j loop*/
			
			
			if (inner_surface_external_transition && external_inner_surface_transition){
				for (m = external_inner_surface_transition; m<=inner_surface_external_transition; m++) {
					if (grid[i][m][k]) {
						if (grid_i[m][k]->flags END_TEST) {
							grid_i[m][k]->flags MID_TRUE;
							grid_i[m][k]->flags END_FALSE;
						}
					}
				}
			}
			  
		} /*k loop*/
#endif
		/* End set septum to be midmyocardium cells */
	

		//write_2D_data(i);

	} /*i loop*/

	for (i=0; i<Z_GRID-1; i++) {
		grid_i = grid[i];
	for (j=1; j<X_GRID-1; j++) {
		grid_i_j = grid_i[j];
	for (k=1; k<Y_GRID-1; k++) {
		grid_i_j_k = grid_i_j[k];
		
		if (grid_i_j_k) {
			if ((grid_i_j_k->flags MID_TEST))
				nMid++;
			else if (grid_i_j_k->flags END_TEST)
				nEnd++;
			else if (grid_i_j_k->flags EPI_TEST)
				nEpi++;
		}
	
	}}}

	printf("nMid = %i\n",nMid);
	printf("nEnd = %i\n",nEnd);
	printf("nEpi = %i\n",nEpi);
	printf("Total= %i\n",(nEpi+nEnd+nMid));
}

void differentiate_heart_tissue2(void) {
	int i, j, k;
	int m,n;
	int outerWallStart = 0;
	int outerWallEnd = 0;
	int rvWallStart = 0;
	int rvWallEnd = 0;
	int lvWallStart = 0;
	int lvWallEnd = 0;
	int nBoundary = 0;
	int rvWallThickness = 0;
	int lvWallThickness = 0;
	int jThickness = 0;
	int septumCenter = 0;

	int septum_width = 0;
	int septum_wall = 0;
	int external_inner_surface_transition = 0;
	int inner_surface_external_transition = 0;

	int epi_mid_boundary = 0;

	int nMid =0;
	int nEpi =0;
	int nEnd =0;

	int radius = 0;

	printf("Set all heart tissue that is not a boundary to be midmyocardium tissue\n");
	for (i=0; i<Z_GRID-1; i++) {
		grid_i = grid[i];
	for (j=1; j<X_GRID-1; j++) {
		grid_i_j = grid_i[j];
	for (k=1; k<Y_GRID-1; k++) {
		grid_i_j_k = grid_i_j[k];
		
		if (grid_i_j_k) {

			grid_i_j_k->flags MID_FALSE;
			if ((grid_i_j_k->flags BOUNDARY_TEST) && (!(grid_i_j_k->flags INT_TEST)))
				grid_i_j_k->flags EPI_TRUE;
			else
				grid_i_j_k->flags END_TRUE;		
		}
	
	}}}

	//printf("Loop through the slices and differentiate the heart tissue types\n");
	for (i=0; i<Z_GRID-1; i++) {
		grid_i = grid[i];

		outerWallStart = 0;
		outerWallEnd = 0;
		rvWallStart = 0;
		rvWallEnd = 0;
		lvWallStart = 0;
		lvWallEnd = 0;
		nBoundary = 0;
		rvWallThickness = 0;
		lvWallThickness = 0;
		jThickness = 0;
		septumCenter = 0;


		//printf("Compute statistics of heart wall sizes for each slice\n");
		for (j=1; j<X_GRID-1; j++) {

			grid_i_j_k = grid_i[j][Y_GRID/2];
			
			/*Am I in the heart?*/
			if (grid_i_j_k) {

				if (grid_i_j_k->flags BOUNDARY_TEST)
					nBoundary++;

				if ((grid_i_j_k->flags BOUNDARY_TEST) && (!(grid_i_j_k->flags INT_TEST))) {
					if (!outerWallStart) {
						outerWallStart = j;
					} else {
						outerWallEnd = j;
						jThickness = outerWallEnd - outerWallStart;
						lvWallThickness = outerWallEnd - lvWallEnd;
					}
				}

				if (grid_i_j_k->flags INT_TEST) {
					if (!rvWallStart) {
						rvWallStart = j;
					} else if (!rvWallEnd) {
						rvWallEnd = j;
						rvWallThickness = rvWallStart - outerWallStart;
					} else if (!lvWallStart) {
						lvWallStart = j;
						septumCenter = (lvWallStart - rvWallEnd)/2;
					} else // if (~lvWallEnd)
						lvWallEnd = j;
				}
			}  

		} /*j loop*/
		/*END compute statistics of heart wall sizes for each slice*/

		//printf("Determine epicardium\n");
		radius =0;
		for (j=1; j<X_GRID-1; j++) {
			//printf("%i\n",j);
			if (nBoundary == 6) {
                radius = rvWallThickness + (lvWallThickness - rvWallThickness)*(j-outerWallStart)/jThickness;
                radius = ceil(radius/10);
			} else
				radius = 1/(SPACING*10); //creates a 1 mm thick epicardium

			
			for (k=1; k<Y_GRID-1; k++) {

				grid_i_j_k = grid_i[j][k];

				if (grid_i_j_k) {
					if ((grid_i_j_k->flags BOUNDARY_TEST) && (!(grid_i_j_k->flags INT_TEST))) {
						for (m = (j-radius); m<=(j+radius); m++) {
						for (n = (k-radius); n<=(k+radius); n++) {
							if ( (m>0) && (n>0) && (m<X_GRID) && (n<Y_GRID) ) {
								if (grid_i[m][n]) {
									if (grid_i[m][n]->flags END_TEST) {
										if (sqrt( (m-j)*(m-j) + (n-k)*(n-k)) < radius) {
											grid_i[m][n]->flags END_FALSE;
											grid_i[m][n]->flags EPI_TRUE;
										}
									}
								}
							}
						}} /*m,n loop*/
					}
				}
		}} /*j,k loop*/
		/*End determine epicardium */

		//printf("Determine mid_myocardium\n");
		radius =0;
		for (j=1; j<X_GRID-1; j++) {
			if (nBoundary == 6) {
                radius = rvWallThickness + (lvWallThickness - rvWallThickness)*(j-outerWallStart)/jThickness;
                radius = ceil(radius * 3/10); //60% of ventricle wall

			} else
				radius = 3/(SPACING*10); //creates a 3 mm thick midmyocardium
			
			for (k=1; k<Y_GRID-1; k++) {

				grid_i_j_k = grid_i[j][k];

				if (grid_i_j_k) {
					if (grid_i_j_k->flags EPI_TEST) {
								 /* Check whether north (j-1) neighbor is a boundary element */
						epi_mid_boundary = FALSE;

						if (grid[i][j-1][k]) 
							if (grid[i][j-1][k]->flags END_TEST)
								epi_mid_boundary = TRUE;
						if (grid[i][j+1][k]) 
							if (grid[i][j+1][k]->flags END_TEST)
								epi_mid_boundary = TRUE;

						if (grid[i][j][k-1]) 
							if (grid[i][j][k-1]->flags END_TEST)
								epi_mid_boundary = TRUE;
						if (grid[i][j][k+1]) 
							if (grid[i][j][k+1]->flags END_TEST)
								epi_mid_boundary = TRUE;

						if (epi_mid_boundary) {
							for (m = (j-radius); m<=(j+radius); m++) {
							for (n = (k-radius); n<=(k+radius); n++) {
								if ( (m>0) && (n>0) && (m<X_GRID) && (n<Y_GRID) ) {
									if (grid_i[m][n]) {
										if (grid_i[m][n]->flags END_TEST) {
					//						if (sqrt( (m-j)*(m-j) + (n-k)*(n-k)) < radius) {
												grid_i[m][n]->flags MID_TRUE;
												grid_i[m][n]->flags END_FALSE;
					//						}
										}
									}
								}
							}} /*m,n loop*/
						}
					}
				}
		}} /*j,k loop*/
		/*End determine endocardium */



	} /*i loop*/

	for (i=0; i<Z_GRID-1; i++) {
		grid_i = grid[i];
	for (j=1; j<X_GRID-1; j++) {
		grid_i_j = grid_i[j];
	for (k=1; k<Y_GRID-1; k++) {
		grid_i_j_k = grid_i_j[k];
		
		if (grid_i_j_k) {
			if ((grid_i_j_k->flags MID_TEST))
				nMid++;
			else if (grid_i_j_k->flags END_TEST)
				nEnd++;
			else if (grid_i_j_k->flags EPI_TEST)
				nEpi++;
		}
	
	}}}

	printf("nMid = %i\n",nMid);
	printf("nEnd = %i\n",nEnd);
	printf("nEpi = %i\n",nEpi);
	printf("Total= %i\n",(nEpi+nEnd+nMid));
}


void startTime( char* n){
  printf("starting %s\n" , n);
  wrapper_start = time(NULL);	
  strcpy(&timerName, n);
}

void stopTime( void ){
  time_t stop = time(NULL);
  double diff = difftime(stop, wrapper_start);
  printf("%s Took:%e seconds\n", timerName, diff);
}

void initialize_fast_index(void) {

	int i, j, k, startedY, startedX, startedZ;
	printf("Initializing fast indexing\n");
	
	yBounds = (SPACIAL_ELEMENT **) calloc(Z_GRID, sizeof(SPACIAL_ELEMENT *));
	if (yBounds == NULL) {
		printf("Insufficient Dynamic Memory\n");
		exit(1);
	}
	for (i=0; i< Z_GRID; i++) {
		yBounds[i] = (SPACIAL_ELEMENT *) calloc(X_GRID, sizeof(SPACIAL_ELEMENT));
		if (yBounds[i] == NULL) {
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}

	xBounds = (SPACIAL_ELEMENT *) calloc(Z_GRID, sizeof(SPACIAL_ELEMENT));
	if(xBounds == NULL){
		printf("Insufficient Dynamic Memory\n");
		exit(1);
	}

	startedZ = 0;
	zBounds.start = 0;
	zBounds.end = 0;
	for(i=0; i<Z_GRID; i++){

		startedX = 0;
		xBounds[i].start = 0;
		xBounds[i].end = 0;
		for(j=0; j<X_GRID; j++){
			startedY = 0;
			yBounds[i][j].start = 0;
			yBounds[i][j].end = 0;
			for(k=0; k<Y_GRID; k++){
				if(grid[i][j][k]){
					/*set up yBounds */
					if(!startedY){
						startedY = 1;
						yBounds[i][j].start = k;
					}
					yBounds[i][j].end = k;
					
					if(!startedX){
						startedX = 1;
						xBounds[i].start = j;
					}
					xBounds[i].end = j;
					if(!startedZ){
						startedZ = 1;
						zBounds.start = i;
					}
					zBounds.end = i;
				}
	}
	}
	}

	printf("zbounds.start = %i\n", zBounds.start);
	printf("zbounds.end = %i\n", zBounds.end);
}

////////////////////////////////////////////
//PACING and PURKINJE SIMULATION CODE///////
////////////////////////////////////////////

#if HUMAN
void purkinje_fibre_stim_human(int t) {

	int zstart, zend;
	int i, j, k;
	SPACIAL_ELEMENT* yBounds_i;
	SPACIAL_ELEMENT yBounds_i_j;

	//simulate the purkinje fibre network of a human

	if(t > 7) return;

	//zend = 166 - (t*20);
	zend = 340 - (t*40);
	//zstart = 166 - ((t+1)*20);
	zstart = 340 - ((t+1)*40);

	for (i=zstart; i<zend; i++) {
	grid_i = grid[i];
	yBounds_i = yBounds[i];
	for (j=xBounds[i].start; j<=xBounds[i].end; j++) {
	grid_i_j = grid_i[j];
	yBounds_i_j = yBounds_i[j];
	for (k=yBounds_i_j.start; k<=yBounds_i_j.end; k++) {
		grid_i_j_k = grid_i_j[k];
		if (grid_i_j_k) {
			if (grid_i_j_k->flags INT_TEST){
				grid_i_j_k->flags EXC_TRUE;
				grid_i_j_k->flags REF_TRUE;
				grid_i_j_k->flags REP_FALSE;
				grid_i_j_k->src = 1;
				grid_i_j_k->snk = 0;
				grid_i_j_k->nsr = 0;

				printf("grid(%i,%i,%i).u = %f\n",i,j,k,grid_i_j_k->u);
			}
		}
	}}}
}

void purkinje_fibre_stim_human3(void) {

int i,j,k;
int nPaced = 0;
int radius = 5;
int m,n;
SPACIAL_ELEMENT* yBounds_i;
SPACIAL_ELEMENT yBounds_i_j;

		for (i=0; i<= zBounds.end; i++) {
		grid_i = grid[i];
		yBounds_i = yBounds[i];
		for (j=xBounds[i].start; j<=xBounds[i].end; j++) {
		grid_i_j = grid_i[j];
		yBounds_i_j = yBounds_i[j];
		for (k=yBounds_i_j.start; k<=yBounds_i_j.end; k++) {
			grid_i_j_k = grid_i_j[k];
			/*
			if (grid_i_j_k) {
				if (grid_i_j_k->flags END_TEST){
					if(!(grid_i_j_k->flags EXC_TEST)) {
						grid_i_j_k->flags EXC_TRUE;
						grid_i_j_k->flags REF_TRUE;
						grid_i_j_k->flags REP_FALSE;
						grid_i_j_k->src = 1;
						grid_i_j_k->snk = 0;
						grid_i_j_k->nsr = 0;
						//printf("grid(%i,%i,%i).u = %f\n",i,j,k,grid_i_j_k->u);
						nPaced++;
					}
				}
			}
			*/

			if (grid_i_j_k) {
					if (grid_i_j_k->flags INT_TEST) {
						for (m = (j-radius); m<=(j+radius); m++) {
						for (n = (k-radius); n<=(k+radius); n++) {
							if ( (m>0) && (n>0) && (m<X_GRID) && (n<Y_GRID) ) {
								if (grid_i[m][n]) {
									if ( !(grid_i[m][n]->flags EXC_TEST)) {
										if (sqrt( (m-j)*(m-j) + (n-k)*(n-k)) < radius) {
											grid_i[m][n]->flags EXC_TRUE;
											grid_i[m][n]->flags REF_TRUE;
											grid_i[m][n]->flags REP_FALSE;
											grid_i[m][n]->src = 1;
											grid_i[m][n]->snk = 0;
											grid_i[m][n]->nsr = 0;
						
											nPaced++;
										}
									}
								}
							}
						}} /*m,n loop*/
					}
				}
		}}}
		printf("Number Paced = %i\n",nPaced);
		numEXC += nPaced;
}

void purkinje_fibre_stim_human2(int t) {

	int zstart, zend;
	int i, j, k;
	int m,n;
	int nPaced = 0;

	SPACIAL_ELEMENT* yBounds_i;
	SPACIAL_ELEMENT yBounds_i_j;

	int septum_width = 0;
	int septum_wall = 0;
	int external_inner_surface_transition = 0;
	int inner_surface_external_transition = 0;

	//simulate the purkinje fibre network of a human

	if(t >= 15) 
		return;
	else if (t>=5) {
		printf("Activate the interior endocardium, t = %i\n",t);
		
		zstart = 320 - (((t-5)+1)*32) ;
		zend = 320 - ((t-5)*32) - 1;
		zend = zend/ (.035/.025);
		zstart = zstart/ (.035/.025);

		printf("zstart = %i\n",zstart);
		printf("zend = %i\n",zend);

		for (i=zstart; i<zend; i++) {
		grid_i = grid[i];
		yBounds_i = yBounds[i];
		for (j=xBounds[i].start; j<=xBounds[i].end; j++) {
		grid_i_j = grid_i[j];
		yBounds_i_j = yBounds_i[j];
		for (k=yBounds_i_j.start; k<=yBounds_i_j.end; k++) {
			grid_i_j_k = grid_i_j[k];
			if (grid_i_j_k) {
				if (grid_i_j_k->flags END_TEST){
					if(!(grid_i_j_k->flags EXC_TEST)) {
						grid_i_j_k->flags EXC_TRUE;
						grid_i_j_k->flags REF_TRUE;
						grid_i_j_k->flags REP_FALSE;
						grid_i_j_k->src = 1;
						grid_i_j_k->snk = 0;
						grid_i_j_k->nsr = 0;
						//printf("grid(%i,%i,%i).u = %f\n",i,j,k,grid_i_j_k->u);
						nPaced++;
					}
				}
			}
		}}}

		zend = 300 + (((t-5)+1)*10) ;
		zstart = 300 + ((t-5)*10) - 1;
		zend = zend/(.035/.025);
		zstart = zstart/(.035/.025);

		printf("zstart = %i\n",zstart);
		printf("zend = %i\n",zend);

		for (i=zstart; i<zend; i++) {
		grid_i = grid[i];
		yBounds_i = yBounds[i];
		for (j=xBounds[i].start+15; j<=xBounds[i].end-15; j=j+5) {
		grid_i_j = grid_i[j];
		yBounds_i_j = yBounds_i[j];
		for (k=yBounds_i_j.start+15; k<=yBounds_i_j.end-15; k=k+5) {
			grid_i_j_k = grid_i_j[k];
			if (grid_i_j_k) {
				if(!(grid_i_j_k->flags EXC_TEST)) {
					grid_i_j_k->flags EXC_TRUE;
					grid_i_j_k->flags REF_TRUE;
					grid_i_j_k->flags REP_FALSE;
					grid_i_j_k->src = 1;
					grid_i_j_k->snk = 0;
					grid_i_j_k->nsr = 0;
					//printf("grid(%i,%i,%i).u = %f\n",i,j,k,grid_i_j_k->u);
					nPaced++;
				}

			}
		}}}

		printf("Number Paced = %i\n",nPaced);
	} else {
		
		printf("Activate the bundle of His, t = %i\n",t);
		zstart = t * 60;
		zend = ((t+1)*60)-1;
		zend = zend/(.035/.025);
		zstart = zstart/(.035/.025);

		printf("zstart = %i\n",zstart);
		printf("zend = %i\n",zend);

		for (i=zstart; i<zend; i++) {
			grid_i = grid[i];
			for (k=1; k < Y_GRID; k++) {
				septum_width = 0;
			 	septum_wall = 0;
				external_inner_surface_transition = 0;
				inner_surface_external_transition = 0;

				for (j = 1; j < X_GRID; j++) {
					if (grid_i[j][k]) {
						if ( (grid_i[j][k]->flags INT_TEST) && (!(grid_i[j-1][k]) ) )
							external_inner_surface_transition = j;
				
						if (external_inner_surface_transition)
							septum_wall++;
           
						if ( (grid_i[j][k]->flags INT_TEST) && (!grid_i[j+1][k]) && external_inner_surface_transition) {
							inner_surface_external_transition = j;
							j = 500;
						}
					}
				} //j loop

			//	printf("ext = %i\n",external_inner_surface_transition);
			//	printf("in = %i\n",inner_surface_external_transition);
				
				if (inner_surface_external_transition && external_inner_surface_transition){
					for (m = external_inner_surface_transition; m<=inner_surface_external_transition; m++) {
						if (grid_i[m][k]) {
							if (!(grid_i[m][k]->flags EXC_TEST)) {
								grid_i[m][k]->flags EXC_TRUE;
								grid_i[m][k]->flags REF_TRUE;
								grid_i[m][k]->flags REP_FALSE;
								grid_i[m][k]->src = 1;
								grid_i[m][k]->snk = 0;
								grid_i[m][k]->nsr = 0;
								nPaced++;
							}
						}
					}
				}
			} //k loop
		} //i loop

		printf("Number Paced = %i\n",nPaced);
	}
}

void simple_pacing_from_top_human(void) {
	int i,j,k;
	int nPaced = 0;

	for (i=8-5; i<= 8+5; i++) {
		grid_i=grid[i];
	for (j=140/2-5; j<= 140/2+5; j++) {
		grid_i_j=grid_i[j];
 	for (k=Y_GRID/2-5; k<= Y_GRID/2+5; k++) {
		grid_i_j_k=grid_i_j[k];
		if (grid_i_j_k) {
			grid_i_j_k->flags EXC_TRUE;
			grid_i_j_k->flags REF_TRUE;
			grid_i_j_k->flags REP_FALSE;
			grid_i_j_k->src = 1;
			grid_i_j_k->snk = 0;
			grid_i_j_k->nsr = 0;
			nPaced++;
		}
 	}}}
	printf("Number of CA excited = %i\n",nPaced);
}

void simple_pacing_from_bottom_human(void) {
	int i,j,k;

	for (i=190-3; i<= 190+3; i++) {
		grid_i=grid[i];
 	for (j=70-3; j<= 70+3; j++) {
		grid_i_j=grid_i[j];
 	for (k=Y_GRID/2-3; k<= Y_GRID/2+3; k++) {
		grid_i_j_k=grid_i_j[k];
		if (grid_i_j_k) {
			grid_i_j_k->flags EXC_TRUE;
			grid_i_j_k->flags REF_TRUE;
			grid_i_j_k->flags REP_FALSE;
			grid_i_j_k->src = 1;
			grid_i_j_k->snk = 0;
			grid_i_j_k->nsr = 0;
		}
 	}}}
}

#elif CANINE
void simple_pacing_from_bottom_canine(void) {
	int i,j,k;

	for (i=110-3; i<= 110+3; i++) {
		grid_i=grid[i];
 	for (j=70-3; j<= 70+3; j++) {
		grid_i_j=grid_i[j];
 	for (k=Y_GRID/2-3; k<= Y_GRID/2+3; k++) {
		grid_i_j_k=grid_i_j[k];
		if (grid_i_j_k) {
			grid_i_j_k->flags EXC_TRUE;
			grid_i_j_k->flags REF_TRUE;
			grid_i_j_k->flags REP_FALSE;
			grid_i_j_k->src = 1;
			grid_i_j_k->snk = 0;
			grid_i_j_k->nsr = 0;
		}
 	}}}
}
#endif

////////////////////////////////////////////
//DATA OUTPUT AND VISUALIZATION CODE////////
////////////////////////////////////////////
void write_2D_data(int tm) {
		int u;
		int i,j,k;
		FILE *fp;
		char fname[256];
		unsigned char image[Y_GRID*X_GRID];
		printf("Writing 2d Data file\n");
	 	sprintf(fname, "2D_DataFile%i.bmp", tm);

		i = Z_GRID/2;
		grid_i=grid[i];
		for (j=0; j<X_GRID; j++) {
			grid_i_j=grid_i[j];
			for (k=0; k<Y_GRID; k++)  {
				if (grid_i_j[k]) {
					u = grid_i_j[k]->u*255;		
				} else
					u =0;
				image[(j*X_GRID)+k] = u;
		}}
		writeGrayScaleDataToBmpFile(fname, Y_GRID, X_GRID,image);

}

void write_Slice_Geometry(void) {
		int R,G,B;
		int i,j,k;
		FILE *fp;
		char fname[256];
		unsigned char image[Y_GRID*X_GRID*3];
		printf("Writing 2D, RGB Data file\n");
	 	

		for (i=0; i<Z_MAX; i++) {
			sprintf(fname, "Slice_%i_Tissue_Types.bmp", i);
			grid_i = grid[i];
			for (j=0; j<X_GRID; j++) {
				grid_i_j=grid_i[j];
				for (k=0; k<Y_GRID; k++)  {
					grid_i_j_k = grid_i_j[k];
					R = 0;
					G = 0;
					B = 0;
					if (grid_i_j_k) {

						if (grid_i_j_k->flags EPI_TEST )
							R = 255;
						if (grid_i_j_k->flags END_TEST )
							G = 255;
						if (grid_i_j_k->flags MID_TEST )
							B = 255;
					} 
					image[3*((j*X_GRID)+k)] = R;
					image[3*((j*X_GRID)+k)+1] = G;
					image[3*((j*X_GRID)+k)+2] = B;
			}}
			write24BitBmpFile(fname, Y_GRID, X_GRID, image);
		}

}

void write_3D_data(int tm) {
	unsigned short int value;
	float u;
	int i,j,k;
	FILE *fp;
	char fname[256];

	printf("Writing 3D Data file\n");
	sprintf(fname, "3DDataFile%d", tm);		
	fp = fopen(fname, "wb");

	for (i=0; i<Z_GRID; i++) {
		grid_i = grid[i];
		for (j=0; j<X_GRID; j++) {
			grid_i_j = grid_i[j];
			for (k=0; k<Y_GRID; k++) {
				if (grid_i_j[k])
					u = grid_i_j[k]->u * 255;
				else 
					u = 0.0;
				value = (unsigned short int)u;
				fwrite(&value, sizeof(unsigned short int), 1, fp);
	}}}
	fclose(fp);
}

#ifdef BZ_IMPORT
void write_Compressed_3D_data(int tm) {
	BZFILE* b;
	int     level = 1;
	int     nBuf = 0;
	unsigned short int    buf[Y_GRID];
	int     bzerror;
	int     nWritten;
	char mode[10];
	float u;
	int i,j,k;
	FILE *f;
	char fname[256];
	time_t compress_start, compress_stop;
	double compress_diff;

	compress_start = time(NULL);
	printf("Writing Compressed 3D Data file\n");
	sprintf(fname, "3DDataFile%d", tm);

#ifdef _WIN32
	mode[0]='w';
    mode[1] = '0' + level;
    mode[2] = '\0';
	b = BZ2_bzopen(fname, mode);
#else
	f = fopen(fname, "wb");

	if (!f) {
	   /* handle error */
	   printf("Error opening file%s", fname);
	   exit(1);
	}
	b = BZ2_bzWriteOpen ( &bzerror, f, level, 1, 30);
	if (bzerror != BZ_OK) {
	   BZ2_bzWriteClose (&bzerror, b, 0, NULL, NULL);
	   /* handle error */
	   printf("Error calling BZ2_bzWriteOpen");
	   exit(1);
	}
#endif

	u = 0.0;
	for (i=0; i<Z_GRID; (i=i+2)) { 
	grid_i = grid[i];
	for (j=0; j<X_GRID; (j=j+2)) {
		grid_i_j = grid_i[j];
		nBuf = 0;
		for (k=0; k<Y_GRID; (k=k+2)) { 
			if (grid_i_j[k]) {
				if (i>=420)
					u = 0.0;
				else {
					/*
					u = grid_i_j[k]->u;
					*/
					if (grid_i_j[k]->ischemic)
						u = 0.75;
					else
						u = 0.25;
				}
			} else {
				u = 0.0;
			}
			u = 255 * u;
			buf[nBuf] = (unsigned short int)u;
			nBuf += 1;
		}
		//printf("BZ2_bzWrite\n");
#ifdef _WIN32
		BZ2_bzwrite(b, buf, nBuf * 2);
#else
		BZ2_bzWrite ( &bzerror, b, buf, nBuf);
		if (bzerror == BZ_IO_ERROR) { 
			BZ2_bzWriteClose ( &bzerror, b, 0, NULL, NULL);
			/* handle error */
			printf("Error calling BZ2_bzWrite");
			exit(1);
		}
#endif
	}}

#ifdef _WIN32
	BZ2_bzclose(b);
#else
	BZ2_bzWriteClose ( &bzerror, b, 0, NULL, NULL);
	if (bzerror == BZ_IO_ERROR) {
   		/* handle error */
		printf("Error closing file");
		exit(1);
	}
#endif	
	compress_stop = time(NULL); 
	compress_diff = difftime(compress_stop, compress_start); 
	printf("Compressing and writing took (in seconds) %e\n", compress_diff); 
}
#endif
