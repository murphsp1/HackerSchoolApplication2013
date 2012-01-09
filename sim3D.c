#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "nr.h"
//#include "nrutil.h"

typedef	struct {float x, y, z;}	Point;
#define	HC_POINT Point

//#include "hc.h"

#define True				1
#define False				0

#define R_D					3.60		/* Dimensionless fast axis neighborhood size 	  */
#define K					4.95		/* Threshold in element units for anisotropy=     */

#define N_EXCS				1			/* Excite for 3 steps (trigger waves for now)	  */
#define H_MASK				7			/* Half of voxel width for calculating offsets	  */

#define T_GRID	 			1			/* Number of theta divisions */
#define P_GRID	 			1			/* Number of phi divisions   */
#define N_TABS				100			/* Number of random voxel tables at theta, phi    */

#define ANISOTROPY			1.0			/* Global Tissue anisotropy ratio for fibers	  */

#define X_GRID				221			/* Number of elements in X for 12 cm (0.05 space)  */
#define Y_GRID				221			/* Number of elements in Y for 12 cm (0.05 space)  */
#define Z_GRID				221			/* Number of elements in Z for 12 cm (0.05 space)  */
#define SPACING				0.05		/* Grid step in cm								   */

#define PI 3.1415927

typedef unsigned char Boolean;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define NINT(a) ((int)((a)+0.5))

static unsigned short bitMasks[16] = 
	{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768};
	
typedef struct voxel {
	double x;
	double y; 
	double z;
		
} VOXEL;

typedef struct geometry_element {
	unsigned char	inside;			/* Is volume element inside heart?	*/
	unsigned char	boundary;		/* Is volume element on boundary?	*/

} GEOMETRY_ELEMENT;

/* 3D array with bits defining local neighborhood  */
typedef struct voxel_table {
	unsigned short bits[16][16];	
	
} VOXEL_TABLE;

typedef struct cam_element {
	unsigned char  exc;				/* is element EXCITING 			*/
	unsigned char  ref;				/* is element REFRACTORY		*/
	unsigned char  src;				/* source strength			    */
	unsigned char  snk;				/* sinking value				*/
	unsigned short tau;				/* diastolic interval			*/
	unsigned char  rep;				/* is element repolarized		*/
	unsigned short nsr;				/* front sourcing steps			*/
	unsigned char  the_bin;			/* theta bin for voxel table 	*/
	unsigned char  phi_bin;			/* phi bin for voxel table		*/
	
} CAM_ELEMENT;


static VOXEL rvoxels[16][16][16];

static VOXEL_TABLE ***tableArray;

static CAM_ELEMENT ***grid;

static GEOMETRY_ELEMENT	***geometry;

/* Actual time step and grid step for CA solutions */		
#define DT     0.002			/* Time step (sec)  2.0 ms						*/
#define DX	   0.05				/* Space step (cm)  500 micrometers				*/

/* ADI method parameters (time step is half DT) */
#define DTP    0.00025			/* Time step (sec)  0.25 ms						*/
#define DIFF   0.7				/* Isotropic diffusion constant 				*/
#define FACT   0.07				/* DIFF*DTP/(DX**2) 							*/

/* Action potential parameters */
#define UMAX	0.99			/* Maximum dimensionless potential value		*/
#define UMIN    0.01			/* Minimum dimensionless potential value		*/
#define UAPD	0.05  			/* 95% repolarization (end of APD)				*/

/* Normal (APD Restitution parameters) */
#define DI_INIT  640			/* Initial diastolic interval 					*/
#define AP_DMN	 20				/* Mininum Diastolic interval					*/

#define S2		2000.0			/* Parameter for Excited Current Source			*/

/* Basic PDE element structure */
typedef struct pde {
	float u;					/* pseudo potential value						*/
	float tmp;					/* temporary potential value					*/
	float dew;					/* second potential difference (east-west)		*/
	float dns;					/* second potential difference (north-south)	*/
	float dud;					/* second potential difference (up-down)		*/
	float bet;					/* I-V curve parameter for desired apd			*/
	
} PDE_ELEMENT;

typedef struct ecg {
	unsigned short time;		/* time step 									*/
	unsigned char i, j, k;		/* indices 										*/
	float source;				/* sum of 3 second spatial derivatives			*/

} ECG_ELEMENT;

static PDE_ELEMENT ***gpde;
static float ***u_xyz;

static float apd_slp = 0.0;		/* Slope for restitution curve (normal)			*/
static float apd_exp = 0.0;		/* Exponent for restitution curve (normal)		*/

static int iSeed = -9;
static float epsilon = 1.e-8;

/* Function Declarations */
static int initialize_geometry (void );
static int initialize_display ( void );
static int initialize_pde (void );


static int voxel_get_table_entry (VOXEL_TABLE *table, int i, int j, int k);
static int voxel_init_random_table_array (void);
static VOXEL_TABLE *voxel_get_random_table (unsigned char i, unsigned char j);
static int voxel_get_table (float the, float phi);
static int voxel_build_random_table (VOXEL_TABLE *table, double speed_ratio, double the, double phi);
static void cam_initialize ( void );

static float ran1(int *idum);
static double atanch (double val);
static int compute_apd_parameters (float apd, float vmin, float *tsh, float *bet);
static float apd_from_restitution (unsigned short tau);


#define LEAD_X 		0.0
#define LEAD_Z	   -5.0
#define LEAD_Y	   -9.0

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
			
			if ((++cnt % 100000 == 0)) {
				printf("count = %d\n", cnt);
			}
		}
		fclose(fp);
		fprintf(fo, "%d\t%f\n", 2*time, ecg);
		printf("time = %d\tecg = %lf\n", 2*time, ecg);
	}
	fclose(fo);
	exit(1);
}


double ComputeECG ( void );
double ComputeECG ( void ) 
{
	int i, j, k;
	double dx2, dy2, dz2;
	double den, x, y, z, ecg, source;
	
	ecg = 0.0;
	for (i=0; i<Z_GRID; i++) {
		z = -i*SPACING;
		dz2 = (z - LEAD_Z)*(z - LEAD_Z);
		
		for (j=0; j<X_GRID; j++) {
			x = -5.5 + j*0.05;
			dx2 = (x - LEAD_X)*(x - LEAD_X);
			
			for (k=0; k<Y_GRID; k++) {
				y = -5.5 + k*0.05;
				dy2 = (y - LEAD_Y)*(y - LEAD_Y);
				source = gpde[i][j][k].dew + gpde[i][j][k].dns + gpde[i][j][k].dud;
		
				den = sqrt(dx2 + dy2 + dz2);
				ecg += source/den;
			}
		}
	}
	return (ecg);
}
	

int main( void )
{
	char fname[256];
	int tm, pm;
	int i, j, k;
	int ir, jr, kr;
	int iv, jv, kv;
	float apd, tsh, bet;
	float d2x, d2y, d2z, i_0, u_0;
	float ux_p, ux_m, uy_p, uy_m, uz_p, uz_m;
	double ecgLead;
	
	FILE *ecg;
	ECG_ELEMENT ecgElement;

	CAM_ELEMENT *src_p, *snk_p;
	PDE_ELEMENT *pde_p;
	VOXEL_TABLE *table;
	
 /* Initialize Geometry */
 	initialize_geometry ();		/* Must do this first */

 /* Build voxel tables */
	voxel_init_random_table_array ();
	
 /* Initialize CAM grid	*/
 	cam_initialize ();

 /* Initialize PDE system */
 	initialize_pde ();

 /* Initialize Display */
	//initialize_display ();
	
 /* Update the display */
	//update_display ();
	
 /* Open ECG file */
 	ecg = fopen("ECG1_Full.dat", "w");
	
 /* Loop over the total number of time steps */
	for (tm = 0; tm < 5000; tm++) {
		printf("CA iteration = %d\n", tm);
	
	 /* Compute ECG */
	 	printf("Computing ECG ...\n");
	 	ecgLead = ComputeECG();
		fprintf(ecg, "%lf\t%lf\n", tm*2.0 , ecgLead);
		printf("time = %lf (ms) , ECG = %lf\n", tm*2.0 , ecgLead);		
	
 	 /* Pace at 60 bpm */
	 	if ((tm < 2000)) {
			if ((tm % 500) == 0) {
 				for (i=5-3; i<= 5+3; i++) {
 				for (j=70-3; j<= 70+3; j++) {
 				for (k=Y_GRID/2-3; k<= Y_GRID/2+3; k++) {
 					grid[i][j][k].exc = True;
					grid[i][j][k].ref = True;
					grid[i][j][k].rep = False;
		
 					grid[i][j][k].src = 1;
 					grid[i][j][k].snk = 0;
 					grid[i][j][k].nsr = 0;
 				}}}
			}
		}
		
	 /* Pace at 120 bpm */
	 	if ((tm >= 2000)) {
			if ((tm % 250) == 0) {
 				for (i=5-3; i<= 5+3; i++) {
 				for (j=70-3; j<= 70+3; j++) {
 				for (k=Y_GRID/2-3; k<= Y_GRID/2+3; k++) {
 					grid[i][j][k].exc = True;
					grid[i][j][k].ref = True;
					grid[i][j][k].rep = False;
		
 					grid[i][j][k].src = 1;
 					grid[i][j][k].snk = 0;
 					grid[i][j][k].nsr = 0;
 				}}}
			}
		}
		
		
	 /* EXCITING elements source to neighbors */
		for (i=0; i<Z_GRID; i++) {
		for (j=0; j<X_GRID; j++) {
		for (k=0; k<Y_GRID; k++) {
		
		 /* Get element pointer */
			src_p = &(grid[i][j][k]);
			
		 /* Ignore non-EXCITING */
		 	if (src_p->exc) {
						
		 	/* Get a random neighborhood table for this element */
				table = voxel_get_random_table (src_p->the_bin, src_p->phi_bin);
					 	
		  	/* Loop over voxel table and source to neighbors as defined in table */
		  		for (ir=0; ir < 15; ir++) {
		  		for (jr=0; jr < 15; jr++) {
		  		for (kr=0; kr < 15; kr++) {
				
			 	/* See if element is a neighbor given local fiber direction */
					if (!voxel_get_table_entry(table, ir, jr, kr)) {
						continue;
					}
					  	
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
			 		if (!geometry[iv][jv][kv].inside) {
						continue;
					}
	
				 /* Get neighboring CAM element */
			 		snk_p = &(grid[iv][jv][kv]);
				
				 /* Add source unit if element is not refractory */
					if (!snk_p->ref) {
			 			snk_p->snk += src_p->src;	/* send source to non front elements */
					}
				}}}
				
			 /* Update phase for element (after sourcing) */
		 		src_p->nsr++;
			}
		 
		}}}  /* end loop over EXCITING elements */

	 /* Check for excitation of elements */
		for (i=0; i<Z_GRID; i++) {
		for (j=0; j<X_GRID; j++) {
		for (k=0; k<Y_GRID; k++) {
			snk_p = &(grid[i][j][k]);
			
		 /* Ignore exciting elements and refractory elements */
		 	if (snk_p->exc || snk_p->ref) {
				continue;
			}
																		
		 /* Check for switch to EXCITING state for RESTING elements */
			if (snk_p->snk >= K) {
			
				snk_p->exc = True;				/* EXCITING state 		*/
				snk_p->ref = True;
				snk_p->rep = False;				/* Not repolarized 		*/
				
				snk_p->nsr = 0;					/* sourcing step zero 	*/
			} 
		}}}  /* end loop over excitation check */

	 /* Clear all sinks */
		for (i=0; i<Z_GRID; i++) {
		for (j=0; j<X_GRID; j++) {
		for (k=0; k<Y_GRID; k++) {
			grid[i][j][k].snk = 0;
			
		}}}  /* end loop over excitation check */
		
	 /* Now, unfortunately for explicit solutions, we will have to do 8 pde steps for every CA step */
	 	for (pm=0; pm < 8; pm++) {
			printf("PDE iteration: %d\n", pm);
		
			for (i=0; i<Z_GRID-1; i++) {
			for (j=1; j<X_GRID-1; j++) {
			for (k=1; k<Y_GRID-1; k++) {
			 
			 /* Get element pointer */
				src_p = &(grid[i][j][k]);
			
			 /* Local potential and restitution parameter */
				u_0 = gpde[i][j][k].u;
				bet = gpde[i][j][k].bet;
				
			 /* Make sure we are in the heart */
				if (!geometry[i][j][k].inside) {
					continue; /* Loop up */
				}
						
			 /* Can turn this code into a "stencils" array and save lots of time! */
			 /* Get neighboring potentials, Implement no flux as required */
				if (!geometry[i][j][k].boundary) {
					uz_p = gpde[i+1][j][k].u;
					uz_m = gpde[i-1][j][k].u;
					ux_p = gpde[i][j+1][k].u;
					ux_m = gpde[i][j-1][k].u;
					uy_p = gpde[i][j][k+1].u;
					uy_m = gpde[i][j][k-1].u;
				
				} else {
					if (i==0) {
						uz_m = u_0;
					} else {
						uz_m = (geometry[i-1][j][k].inside) ? gpde[i-1][j][k].u : u_0;
					}
					uz_p = (geometry[i+1][j][k].inside) ? gpde[i+1][j][k].u : u_0;
					ux_m = (geometry[i][j-1][k].inside) ? gpde[i][j-1][k].u : u_0;
					ux_p = (geometry[i][j+1][k].inside) ? gpde[i][j+1][k].u : u_0;
					uy_m = (geometry[i][j][k-1].inside) ? gpde[i][j][k-1].u : u_0;
					uy_p = (geometry[i][j][k+1].inside) ? gpde[i][j][k+1].u : u_0;
				}
			
			 /* Compute 2nd differences for (isotropic) laplacian */
				d2x = ux_p + ux_m - 2*u_0;
				d2y = uy_p + uy_m - 2*u_0;
				d2z = uz_p + uz_m - 2*u_0;
				
			 /* Set current depending on whether exc or rep */
				if (grid[i][j][k].exc) {
					i_0 = S2*(1.0 - u_0);
				} else {
					i_0 = 2*bet*u_0*(u_0 - 1.0);
				}
				
			 /* Excplicit integration */
				gpde[i][j][k].tmp = u_0 + FACT*(d2x + d2y + d2z) + DTP*i_0;
								
			 /* store 2nd differences for ECG computation */
				gpde[i][j][k].dew = d2x;
				gpde[i][j][k].dns = d2y;
				gpde[i][j][k].dud = d2z;

			}}}
			
		 /* Load temp potentials into actual potentials */
			for (i=0; i<Z_GRID; i++) {
			for (j=0; j<X_GRID; j++) {
			for (k=0; k<Y_GRID; k++) {
				gpde[i][j][k].u = MIN(gpde[i][j][k].tmp, UMAX);
				gpde[i][j][k].u = MAX(gpde[i][j][k].u, UMIN);
			}}}
		}

	 /* Check for Switching of EXC elements to REF AND REF to Excitable */
		for (i=0; i<Z_GRID; i++) {
		for (j=0; j<X_GRID; j++) {
		for (k=0; k<Y_GRID; k++) {
		
		 /* Get element pointer */
			src_p = &(grid[i][j][k]);
			
		 /* Check for EXC -> REF */
		 	if (src_p->exc) {
						
			 /* Check whether we are still in the front region for element */
		 		if (src_p->nsr == N_EXCS) {
					src_p->exc = False;				/* Cease CA and PDE sourcing */
			 
			 	 /* Compute new I-V parameters from peak potential for apd */
					apd = apd_from_restitution (grid[i][j][k].tau); 
				   	compute_apd_parameters (apd*0.001, gpde[i][j][k].u, &tsh, &(gpde[i][j][k].bet));

			 /* Assign new source strength at next autochrone based on sourcing table */
				} else {	
					/* This is a NOOP for now, will change with ->nsr; recovery effects */
					src_p->src = 1;
				}		
				continue;	/* Loop up */
			}
			
		 /* Update diastolic interval for repolarized elements */
		 	if (src_p->rep) {
				src_p->tau += 2;		/* CA time steps in 2 milliseconds, restitution curves in ms */
			}
				
		 /* Check for REF-> Excitable */
		 	if (src_p->ref && src_p->rep && (src_p->tau >= AP_DMN)) {
				src_p->ref = False;
			}

		}}}  /* end loop over EXCITING elements */
		
	 /* Check for repolarization via PDEs */
		for (i=0; i<Z_GRID; i++) {
		for (j=0; j<X_GRID; j++) {
		for (k=0; k<Y_GRID; k++) {

		 /* Get element pointer */
			src_p = &(grid[i][j][k]);
			
		 /* Ignore elements in wavefront region */
			if (src_p->exc) {
				continue;
			}
			
		 /* Check for repolarization at APD */
		    if (!src_p->rep) {
				if (src_p->ref && gpde[i][j][k].u <= UAPD) {
					src_p->rep = True;		/* end of action potential */
					src_p->tau = 0;			/* start diastolic interval*/
				}
			}
		}}}

				
	} /* End loop over time steps */
	fclose(ecg);
	

	return (0);
}

static double atanch (double val)
{
	return (0.5*log((1+val)/(1-val)));
}
		
static int apd_parameters (float apd, float *tsh, float *bet)
{
	float v0 = UAPD;
	float vm = UMAX;
	
	*bet = (float) ((1/apd) * (atanch(2*vm-1) - atanch(2*v0-1)));
	*tsh = (float) (apd * (1/(1 - atanch(2*v0-1)/atanch(2*vm-1))));
}

static int compute_apd_parameters (float apd, float umax, float *tsh, float *bet)
{
	float u0 = UAPD;
	float um;
	
	um   = umax;
	*bet = (float) ((1/apd) * (atanch(2*um-1) - atanch(2*u0-1)));
	*tsh = (float) (apd * (1/(1-atanch(2*u0-1)/atanch(2*um-1))));
}

static float apd_from_restitution (unsigned short tau)
{
	double alpha = 1.2;
	double scale = 0.3;
	double b, apd;
	double t;
	
	t = 0.001*tau;
	b = scale*alpha;
	apd = 0.5*(b*b + sqrt(b*b*b*b + 4*b*b*t));
	return ((float) apd*1000);
}

static int initialize_pde ( void )
{
	int i, j, k;
	float apd, tsh, bet;
	Boolean interior, boundary;
	
	printf("Initializing pde element memory ...\n");

 /* Allocate pde solution array */
	gpde = (PDE_ELEMENT ***) calloc(Z_GRID, sizeof(PDE_ELEMENT **));
	if (gpde == NULL) {
		printf("Insufficient Dynamic Memory\n");
		exit(1);
	}
	for (i=0; i< Z_GRID; i++) {
		gpde[i] = (PDE_ELEMENT **) calloc(X_GRID, sizeof(PDE_ELEMENT *));
		if (gpde[i] == NULL) {
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}
	for (i=0; i< Z_GRID; i++) {
	for (j=0; j< X_GRID; j++) {
		gpde[i][j] = (PDE_ELEMENT *) calloc(Y_GRID, sizeof(PDE_ELEMENT));
		if (gpde[i][j] == NULL) {	
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}}		
	
 /* Get IV-Curve characteristics and add to elements */
	apd = apd_from_restitution ((unsigned short) (DI_INIT)); 		
	apd_parameters (apd*0.001, &tsh, &bet);	
	printf("Initial APD = %f\n", apd);
	printf("Initializing array element\n");
	
 /* Initialize PDE elements to default values */
	for (i=0; i<Z_GRID; i++) {			
	for (j=0; j<X_GRID; j++) {    
	for (k=0; k<Y_GRID; k++) {    
		gpde[i][j][k].u   = UMIN;			/* Potential value 				*/
		gpde[i][j][k].tmp = 0.0;			/* Temp potential value		 	*/
		gpde[i][j][k].bet = bet;			/* APD Restitution parameter 	*/
		gpde[i][j][k].dns = 0.0;			/* For ECG computation 			*/
		gpde[i][j][k].dew = 0.0;			/* For ECG computation			*/
		gpde[i][j][k].dud = 0.0;			/* For ECG computaiton			*/
	}}}
 	return (0);
}


static void cam_initialize ( void )
{
	int i, j, k;

	printf("Initializing ca element memory ...\n");
	
 /* Allocate voxel table Array - calloc should init to all bits to zeros */
	grid = (CAM_ELEMENT ***) calloc(Z_GRID, sizeof(CAM_ELEMENT **));
	if (grid == NULL) {
		printf("Insufficient Dynamic Memory\n");
		exit(1);
	}
	for (i=0; i< Z_GRID; i++) {
		grid[i] = (CAM_ELEMENT **) calloc(X_GRID, sizeof(CAM_ELEMENT *));
		if (grid[i] == NULL) {
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}
	for (i=0; i< Z_GRID; i++) {
	for (j=0; j< X_GRID; j++) {
		grid[i][j] = (CAM_ELEMENT *) calloc(Y_GRID, sizeof(CAM_ELEMENT));
		if (grid[i][j] == NULL) {	
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}}			

 /* Initialize array elements to defaults */
    for (i=0; i<Z_GRID; i++) {
    for (j=0; j<X_GRID; j++) {
    for (k=0; k<Y_GRID; k++) {
		grid[i][j][k].exc = False;
		grid[i][j][k].ref = False;
		grid[i][j][k].src = 1;
		grid[i][j][k].nsr = 0;
		grid[i][j][k].snk = 0;
		grid[i][j][k].tau = DI_INIT;
	}}}
 	
 /* Assign fiber voxel array bin values to elements */
    for (i=0; i<Z_GRID; i++) {
    for (j=0; j<X_GRID; j++) {
    for (k=0; k<Y_GRID; k++) {
		grid[i][j][k].the_bin = 0;
		grid[i][j][k].phi_bin = 0;
	}}}

}

static VOXEL_TABLE *voxel_get_random_table (unsigned char i, unsigned char j)
{
	int k = (int) (fabs(N_TABS*ran1(&iSeed)-epsilon)); 
	return (&(tableArray[i][j][0]));
}	

static int voxel_get_table_entry (VOXEL_TABLE *table, int i, int j, int k)
{
	unsigned short mask = bitMasks[k];
	unsigned short data = table->bits[i][j];
	return ((int) ((mask & data) ? 1: 0));
}	

int voxel_init_random_table_array (void) 
{
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
	return (0);
}


int voxel_build_random_table (VOXEL_TABLE *table, double speed_ratio, double the, double phi)
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
	return 0;
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


static int CalculateMeshSize ( void );

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

int SetNewZValue (double z) 
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
	return (0);
}

int CalculateVolumes (void)
{
	double x, y, z;
	long nRV, nLV, nMu;
	int InsideOuterEllipse;
	int InsideLeftVentricle;
	int InsideRightVentricle;
	
	nRV = 0;
	nLV = 0;
	nMu = 0;
	for (z=0.0; z <= 9.5; z += 0.05) {
		SetNewZValue (z);
	
		for (x = -5.5; x <= 5.5; x += 0.05) {
		for (y = -5.5; y <= 5.5; y += 0.05) {

			InsideOuterEllipse   = Inside_Outer_Ellipse   (x, y);
			InsideLeftVentricle  = Inside_Left_Ventricle  (x, y);
			InsideRightVentricle = Inside_Right_Ventricle (x, y); 
			
			if (InsideLeftVentricle)  nLV++;
			if (InsideRightVentricle) nRV++;
			if (InsideOuterEllipse && !(InsideLeftVentricle) && !(InsideRightVentricle)) nMu++;

		}}
	}
	printf("nLV = %ld, nRV = %ld, nMu = %ld\n", nLV, nRV, nMu);
	return (0);
}

int initialize_geometry ( void )
{
	int i, j, k;
	int nBoundary;
	double R, x, y, z;
	double rv_thick;
	int InsideOuterEllipse;
	int InsideLeftVentricle;
	int InsideRightVentricle;
	
	printf("Calculating Volumes...\n");
 	CalculateVolumes();
	printf("Initializing Geometry ...\n");
	
 /* Allocate voxel table Array - calloc should init to all bits to zeros */
	geometry = (GEOMETRY_ELEMENT ***) calloc(Z_GRID, sizeof(GEOMETRY_ELEMENT **));
	if (geometry  == NULL) {
		printf("Insufficient Dynamic Memory\n");
		exit(1);
	}
	for (i=0; i< Z_GRID; i++) {
		geometry[i] = (GEOMETRY_ELEMENT **) calloc(X_GRID, sizeof(GEOMETRY_ELEMENT *));
		if (geometry[i] == NULL) {
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}
	for (i=0; i< Z_GRID; i++) {
	for (j=0; j< X_GRID; j++) {
		geometry[i][j] = (GEOMETRY_ELEMENT *) calloc(Y_GRID, sizeof(GEOMETRY_ELEMENT));
		if (geometry[i][j] == NULL) {	
			printf("Insufficient Dynamic Memory\n");
			exit(1);
		}
	}}			

 /* Compute inside/outside heart field for elements */
	for (i=0; i<Z_GRID; i++) {
		z = -i*SPACING;
		SetNewZValue (z);
		
		for (j=0; j<X_GRID; j++) {
			x = -5.5 + j*0.05;
			for (k=0; k<Y_GRID; k++) {
				y = -5.5 + k*0.05;
				
				InsideOuterEllipse   = Inside_Outer_Ellipse   (x, y);
				InsideLeftVentricle  = Inside_Left_Ventricle  (x, y);
				InsideRightVentricle = Inside_Right_Ventricle (x, y); 
		
				if (InsideOuterEllipse && !(InsideLeftVentricle) && !(InsideRightVentricle)) {
					geometry[i][j][k].inside = True;
				}
			}
		}
	}
 
 /* Loop over CA grid elements and set boundary conditions */
 	nBoundary = 0;
	for (i=0; i<Z_GRID-1; i++) {
	for (j=1; j<X_GRID-1; j++) {
	for (k=1; k<Y_GRID-1; k++) {
		geometry[i][j][k].boundary = False;
		
	 /* We are far enough away from array boundaries using our */
	 /* geometry, so don't bother checking on neighbor indices */
		if (geometry[i][j][k].inside) {
		
		 /* Check whether up (i-1) neighbor is a boundary element */
		 /* This is temporary and assumes first sheet has all boundary elements */
		 	if (i==0) {
				geometry[i][j][k].boundary = True;
			} else {
				if (!geometry[i-1][j][k].inside) {
					geometry[i][j][k].boundary = True;
				}
			}
			
		 /* Check whether dn (i+1) neighbor is a boundary element */
			if (!geometry[i+1][j][k].inside) {
				geometry[i][j][k].boundary = True;
			}

		 /* Check whether north (j-1) neighbor is a boundary element */
			if (!geometry[i][j-1][k].inside) {
				geometry[i][j][k].boundary = True;
			}

		 /* Check whether south (j+1) neighbor is a boundary element */
			if (!geometry[i][j+1][k].inside) {
				geometry[i][j][k].boundary = True;
			}

		 /* Check whether west (k-1) neighbor is a boundary element */
			if (!geometry[i][j][k-1].inside) {
				geometry[i][j][k].boundary = True;
			}

		 /* Check whether east (k+1) neighbor is a boundary element */
			if (!geometry[i][j][k+1].inside) {
				geometry[i][j][k].boundary = True;
			}
			if (geometry[i][j][k].boundary) nBoundary++;
		}
	}}}
	printf("Number of boundary elements: %d\n", nBoundary);
	return (0);
}	


#define COLORS  256
            
