/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

/* Example C (not C++) driver program for simple forward model*/


#include "malloc.h"
#include "gatdaem1d.h"

int main(int argc, char* argv[])
{		
	//const char* systemfile = "Skytem-LM.stm";
	const char* systemfile = argv[1];

	void* SysHandle = createhandle(systemfile);	

	double tx_height = 35.0;
	double tx_roll   = 0.0;
	double tx_pitch  = 0.0; 
	double tx_yaw    = 0.0; 
	double txrx_dx   = -12.0; 
	double txrx_dy   = 0.0; 
	double txrx_dz   = +2.0; 
	double rx_roll   = 0.0; 
	double rx_pitch  = 0.0; 
	double rx_yaw    = 0.0; 
	
	int  nw = nwindows(SysHandle);

	int  nlayers = 3;
	double *conductivity = (double*)malloc(nlayers * sizeof(double));
	double *thickness = (double*)malloc((nlayers - 1) * sizeof(double));

	double PX;//X Primary field
	double PY;//Y Primary field
	double PZ;//Z Primary field
	double* SX = (double*)malloc(nw*sizeof(double));//X Secondary field windows
	double* SY = (double*)malloc(nw*sizeof(double));//Y Secondary field windows
	double* SZ = (double*)malloc(nw*sizeof(double));//Z Secondary field windows

	conductivity[0] = 0.05;
	conductivity[1] = 0.20;
	conductivity[2] = 0.01;

	thickness[0] = 30.0;
	thickness[1] = 20.0;
	
	
	int nloops = 10;
	for (int i = 0; i < nloops; i++) {
		//For each forward model just change the geometry and earth as required
		forwardmodel(SysHandle, tx_height,
						tx_roll, tx_pitch, tx_yaw,
						txrx_dx, txrx_dy, txrx_dz,
						rx_roll, rx_pitch, rx_yaw,
						nlayers, conductivity, thickness,
						&PX, &PY, &PZ, SX, SY, SZ);
	}
	
	for (int i = 0; i < nlayers-1; i++) {
		printf("%d\t%lf\t%lf\n", i, conductivity[i], thickness[i]);
	}
	printf("%d\t%lf\n", nlayers-1, conductivity[nlayers-1]);

	printf("nwindows %d\n", nw);
	printf("P \t%lg\t%lg\t%lg\n", PX, PY, PZ);
	for (int i = 0; i < nw; i++) {
		printf("S%d\t%lg\t%lg\t%lg\n", i, SX[i], SY[i], SZ[i]);
	}

	free(conductivity);
	free(thickness);
	free(SX);
	free(SY);
	free(SZ);
	deletehandle(SysHandle);
}