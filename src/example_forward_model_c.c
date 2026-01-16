/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

/* Example C (not C++) driver program for simple forward model*/

#undef  GATDAEM1D_API_EXPORTS
#define GATDAEM1D_STATIC_LINKAGE
#include "gatdaem1d.h"

#define RPT report_error_exit(status,0) //Report last error
#define RPTEX report_error_exit(status,1) //Report last error and exit

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

//Simple array allocation
static double* malloc_array_double(size_t n)
{
	if (n == 0) {
		fprintf(stderr, "malloc_array_double: n == 0\n");
		return NULL;
	}

	double* p = malloc(n * sizeof * p);
	if (!p) {
		fprintf(stderr, "malloc_array_double: allocation failed (%zu elements)\n", n);
		return NULL;
	}

	return p;
};

//Simple array free
static void free_array_double(double** p)
{
	if (p && *p) {
		free(*p);
		*p = NULL; 
	}
};

void print_earth(const char* label, const gatdaem1d_earth* pE) {
	printf("%s\n",label);
	const int32_t nl = pE->nlayers;
	for (int i = 0; i < (nl - 1); i++) {
		printf("%d\t%lf\t%lf\n", i, pE->conductivity[i], pE->thickness[i]);
	}
	printf("%d\t%lf\n", nl - 1, pE->conductivity[nl - 1]);
};

void print_response(const char* label, const gatdaem1d_response* pR) {
	printf("%s\n",label);
	printf("Primary field\n");
	for (size_t i = 0; i < pR->nwindows; i++) {
		printf("P%zd\t%12e\t%12e\t%12e\n", i, pR->PX[i], pR->PY[i], pR->PZ[i]);
	}
	printf("Secondary field\n");
	for (size_t i = 0; i < pR->nwindows; i++) {
		printf("S%zd\t%12e\t%12e\t%12e\n", i, pR->SX[i], pR->SY[i], pR->SZ[i]);
	}
	printf("\n");
};

static void report_error_exit(gatdaem1d_status error_code, int exit_on_error) {
	if (error_code > 0) {
		const char* msg = gatdaem1d_get_last_error_string(error_code);
		fprintf(stderr, "%s", msg);
		if(exit_on_error != 0) exit(error_code);
	}
};

static gatdaem1d_system_geometry* allocate_system_geometry_struct() {
	gatdaem1d_system_geometry* p = malloc(sizeof(gatdaem1d_system_geometry));
	if (!p) {
		fprintf(stderr, "allocate_geometry_struct: allocation failed\n");
		return NULL;
	}
	return p;
};

static void free_system_geometry_struct(gatdaem1d_system_geometry** ppG) {
	if (ppG && *ppG) {
		free(*ppG);
		*ppG = NULL;
	}
};

static gatdaem1d_earth* allocate_earth_struct(const int32_t nlayers){
	gatdaem1d_earth* pE = malloc(sizeof(gatdaem1d_earth));
	if (!pE) {
		fprintf(stderr, "allocate_earth_struct: allocation failed\n");
		return NULL;
	}
	pE->nlayers = nlayers;
	pE->thickness = malloc_array_double(nlayers - 1);
	pE->conductivity = malloc_array_double(nlayers);
	pE->chargeability = malloc_array_double(nlayers);
	pE->timeconstant = malloc_array_double(nlayers);
	pE->frequencydependence = malloc_array_double(nlayers);
	pE->iptype = GATDAEM1D_IPTYPE_NONE;
	return pE;
};

static void free_earth_struct(gatdaem1d_earth** ppE){
	if (ppE && *ppE) {
		free_array_double(&(*ppE)->thickness);
		free_array_double(&(*ppE)->conductivity);
		free_array_double(&(*ppE)->chargeability);
		free_array_double(&(*ppE)->timeconstant);
		free_array_double(&(*ppE)->frequencydependence);
		free(*ppE);
		*ppE = NULL;
	}
};

static gatdaem1d_response* allocate_response_struct(const int32_t nwindows) {
	gatdaem1d_response* sR = malloc(sizeof(gatdaem1d_response));
	if (!sR) {
		fprintf(stderr, "allocate_response_struct: allocation failed\n");
		return NULL;
	}
	sR->nwindows = nwindows;
	sR->PX = malloc_array_double(nwindows);//X Primary field windows
	sR->PY = malloc_array_double(nwindows);//Y Primary field windows
	sR->PZ = malloc_array_double(nwindows);//Z Primary field windows
	sR->SX = malloc_array_double(nwindows);//X Secondary field windows
	sR->SY = malloc_array_double(nwindows);//Y Secondary field windows
	sR->SZ = malloc_array_double(nwindows);//Z Secondary field windows
	return sR;
};

static void free_response_struct(gatdaem1d_response** ppR) {	
	if (ppR  && *ppR) {
		free_array_double(&(*ppR)->PX);
		free_array_double(&(*ppR)->PY);
		free_array_double(&(*ppR)->PZ);
		free_array_double(&(*ppR)->SX);
		free_array_double(&(*ppR)->SY);
		free_array_double(&(*ppR)->SZ);
		free(*ppR);
		*ppR = NULL;
	}
};

int main(int argc, char* argv[])
{
	//The AEM system's .stm file
	//const char* systemfile = "Skytem-LM.stm";
	//const char* systemfile = argv[1];
	//const char* systemfile = "C:\\Users\\rossc\\Work\\code\\repos\\ga-aem\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm";	
	const char* systemfile = "C:\\Users\\rossc\\Work\\code\\repos\\ga-aem\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm";

	//An error status code
	gatdaem1d_status status = 0;
	
	gatdaem1d_system_geometry* pG = allocate_system_geometry_struct();;
	pG->tx_height = 35.0;
	pG->tx_roll = 0.0;		pG->tx_pitch = 0.0;		pG->tx_yaw  = 0.0;
	pG->txrx_dx = -12.62;	pG->txrx_dy  = 0.0;		pG->txrx_dz = +2.16;
	pG->rx_roll = 0.0;		pG->rx_pitch = 0.0;		pG->rx_yaw  = 0.0;

	//Create a AEM system handle
	gatdaem1d_system_handle SysHandle = 0;
	status = createhandle(systemfile, &SysHandle); RPTEX;
	//Get some basic system info
	int32_t nt, ns, nw;
	double la, pc, bf;
	status = nturns(SysHandle, &nt); RPT;
	status = looparea(SysHandle, &la); RPTEX;
	status = peakcurrent(SysHandle, &pc); RPTEX;
	status = basefrequency(SysHandle, &bf); RPTEX;
	status = nsamplesperwaveform(SysHandle, &ns); RPTEX;
	status = nwindows(SysHandle,&nw); RPTEX;

	//Get and print out window times
	double* wt_low = malloc_array_double(nw); RPTEX;
	double* wt_high = malloc_array_double(nw); RPTEX;
	windowtimes(SysHandle, wt_low, wt_high); RPTEX;
	for (int i = 0; i < nw; i++) {
		printf("%6d\t%12e\t%12e\n", i, wt_low[i], wt_high[i]);
	}
	free_array_double(&wt_low);
	free_array_double(&wt_high);
	
	//Get and print out waveform
	double* wv_time    = malloc_array_double(ns);
	double* wv_current = malloc_array_double(ns);
	double* wv_voltage = malloc_array_double(ns);
	waveform(SysHandle, wv_time,wv_current,wv_voltage); RPTEX;
	//for (int i = 0; i < ns; i++) {
	//	printf("%6d\t%12e\t%12e\t%12e\n", i, wv_time[i], wv_current[i], wv_voltage[i]);
	//}
	free_array_double(&wv_time);
	free_array_double(&wv_current);
	free_array_double(&wv_voltage);

	//The Earth model
	const int32_t nl = 3;
	gatdaem1d_earth* pE = allocate_earth_struct(nl);

	pE->iptype = GATDAEM1D_IPTYPE_NONE;

	pE->thickness[0] = 20.0;
	pE->thickness[1] = 20.0;

	pE->conductivity[0] = 0.010;
	pE->conductivity[1] = 0.100;
	pE->conductivity[2] = 0.001;
	
	pE->chargeability[0] = 0.0;
	pE->chargeability[1] = 0.3;
	pE->chargeability[2] = 0.0;
	//pE->chargeability = NULL;

	pE->timeconstant[0] = 0.0;
	pE->timeconstant[1] = 0.001;
	pE->timeconstant[2] = 0.0;

	pE->frequencydependence[0] = 0.0;
	pE->frequencydependence[1] = 0.5;
	pE->frequencydependence[2] = 0.0;

	print_earth("Earth",pE);
	
	gatdaem1d_response* pR = allocate_response_struct(nw);
	size_t nloops = 1;
	for (size_t i = 0; i < nloops; i++) {
		//For each forward model just change the geometry and earth as required

		//pE->iptype = GATDAEM1D_IPTYPE_NONE;
		pE->iptype = GATDAEM1D_IPTYPE_PELTON;
		status = forwardmodel(SysHandle, pG, pE, pR); RPTEX;
		print_response("Forward model: ", pR);

		pE->iptype = GATDAEM1D_IPTYPE_COLECOLE;
		status = forwardmodel(SysHandle, pG, pE, pR); RPTEX;
		print_response("Forward model IP Cole-Cole: ", pR);
		
		//Set the IP Type (pE->iptype) via enum GATDAEM1D_IPTYPE_NONE, GATDAEM1D_IPTYPE_COLECOLE or GATDAEM1D_IPTYPE_PELTON
		pE->iptype = GATDAEM1D_IPTYPE_PELTON;

		//Set the IP Type (pE->iptype) via string "NONE, "COLECOLE" or "PELTON"
		status = iptype_from_string("NONE", &pE->iptype); RPTEX;
		status = iptype_from_string("COLECOLE", &pE->iptype); RPTEX;
		status = iptype_from_string("PELTON", &pE->iptype); RPTEX;

		status = forwardmodel(SysHandle, pG, pE, pR); RPTEX;
		print_response("Forward model IP PELTON: ", pR);
	}

	//Example of computing individual derivatives
	static const char* const dtype_names[] = {"DC", "DT", "DH", "DR", "DX", "DY", "DZ", "DTX_HEIGHT", "DTX_ROLL", "DTX_PITCH", "DTX_YAW", "DRX_ROLL", "DRX_PITCH", "DRX_YAW"};
	size_t n_dtype_names = sizeof dtype_names / sizeof dtype_names[0];
	for (size_t i = 0; i < n_dtype_names; i++) {
		printf("\nComputing derivative type %s\n", dtype_names[i]);
		status = derivative(SysHandle, pG, pE, dtype_names[i], 0, pR); RPTEX;
		print_response(dtype_names[i], pR);
	};

	status = derivative(SysHandle, pG, pE, "DC", 0, pR); RPTEX;
	status = derivative(SysHandle, pG, pE, "DC", 1, pR); RPTEX;
	status = derivative(SysHandle, pG, pE, "DC", 2, pR); RPTEX;
	status = derivative(SysHandle, pG, pE, "DT", 0, pR); RPTEX;
	status = derivative(SysHandle, pG, pE, "DT", 1, pR); RPTEX;

	//Example of fm_dlogc FM and log-base-e layer conductivity derivatives function 
	const int32_t ncalc = nl + 1;//FM DC1 DC2 ... DCN
	const int32_t ncomp  = 3;//X Y Z
	const int32_t nfield = 2;//Primary and Secondary
	//Allocate linear-array for results with row-major [ncalc][ncomp][nvers][nwindow]
	int32_t ne = ncalc * ncomp * nfield * nw;
	double* fm_dlogc_results = malloc_array_double(ne);
	status = fm_dlogc(SysHandle, pG, pE, fm_dlogc_results); RPTEX;
	printf("\nfm_dlogc test\n");
	for (size_t calci = 0; calci < ncalc; calci++) {
		if (calci == 0) printf("Calculation %zd: Forward model\n", calci);
		else           printf("Calculation %zd: Layer %zd log conductivity derivative\n", calci, calci - 1);
		for (size_t wi = 0; wi < nw; wi++) {
			printf("W%02zd", wi);
			for (size_t compi = 0; compi < ncomp; compi++) {
				for (size_t fi = 0; fi < nfield; fi++) {
					size_t idx = calci * (ncomp * nfield * nw) + compi * (nfield * nw) + fi * nw + wi;
					printf("\t%12e", fm_dlogc_results[idx]);
				}
			}
			printf("\n");
		}
	};

	//Free the allocated structs and arrays
	free_earth_struct(&pE);
	free_response_struct(&pR);
	free_array_double(&fm_dlogc_results);

	//Release the AEM system handle
	status = deletehandle(&SysHandle); RPT;
};
