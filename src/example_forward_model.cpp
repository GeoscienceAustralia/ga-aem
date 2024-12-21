/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

/* Example driver program for simple forward model*/

#include <vector>
#include <cstring>
#include <iostream>

#include "general_utils.hpp"
#include "file_utils.hpp"
#include "random_utils.hpp"
#include "blocklanguage.hpp"
#include "vector_utils.hpp"
#include "lem.hpp"
#include "tdemsystem.hpp"
class cLogger glog; //The global instance of the log file manager

using namespace AEM;

int skytem_example_ip()
{
	//Load the AEM system specification files for the Skytem moments //only do this once				
	//cTDEmSystem S("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm");	
	cTDEmSystem S("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm");
	

	//Load the system geometry 
	TDEmGeometry G;
	G.tx_height = 40;
	G.tx_roll = 0;       G.tx_pitch = 0; G.tx_yaw  = 0;
	G.txrx_dx = -13.35;  G.txrx_dy  = 0; G.txrx_dz = +2.00;
	G.rx_roll = 0;       G.rx_pitch = 0; G.rx_yaw  = 0;

	//Create the  earth structure	
	Earth1D E(3);
	E.conductivity[0] = 0.010;
	E.conductivity[1] = 0.100;
	E.conductivity[2] = 0.001;
	E.thickness[0] = 20;
	E.thickness[1] = 40;

	E.chargeability[0] = 0.0;
	E.chargeability[1] = 1.0;	
	E.chargeability[2] = 0.0;

	E.timeconstant[0] = 0.0;
	E.timeconstant[1] = 0.001;
	E.timeconstant[2] = 0.0;

	E.frequencydependence[0] = 0.0;
	E.frequencydependence[1] = 1.0;
	E.frequencydependence[2] = 0.0;
	
	//Create a response object
	TDEmResponse R;
	//Run the forward model
	R = S.forward(E, G);
	for (size_t i = 0; i < R.size(); i++){
		double wct = S.window(i).centre_time();
		printf("%zu %10e %10e\n", i, wct, R.secondary(ZCOMP,i));
	}	
	return 0;
}

int skytem_example() {
	//Load the AEM system specification files for the Skytem moments
	//only do this once
	//LM = Low moment pulse
	cTDEmSystem LM("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm");
	//HM = high moment pulse	
	cTDEmSystem HM("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm");

	//Load the system geometry (same for both moments)
	//This changes every fiducial/station
	TDEmGeometry G;
	G.tx_height = 30;
	G.tx_roll = 0;       G.tx_pitch = 0; G.tx_yaw = 0;
	G.txrx_dx = -12.62;  G.txrx_dy = 0; G.txrx_dz = +2.16;
	G.rx_roll = 0;       G.rx_pitch = 0; G.rx_yaw = 0;

	//Create the earth structure
	//This changes every fiducial/station
	Earth1D E(3);
	E.conductivity[0] = 0.010;
	E.conductivity[1] = 0.100;
	E.conductivity[2] = 0.001;
	E.thickness[0] = 20;
	E.thickness[1] = 40;

	//bottom layer is infinite thickness and not set

	//Create a response object for each moment (they have different numbers of windwos)
	TDEmResponse LMR;
	TDEmResponse HMR;	
				
	//Run the forward model for each moment			
	LMR = LM.forward(E, G);
	HMR = HM.forward(E, G);
	
	for (size_t i = 0; i < LMR.size(); i++){
		printf("%zu %g\n", i, LMR.secondary(ZCOMP,i));
	}

	for (size_t i = 0; i < HMR.size(); i++){
		printf("%zu %g\n", i, HMR.secondary(ZCOMP,i));
	}
	
	return 0;
}

int skytem_computation_time() {
	//Load the AEM system specification files for the Skytem moments
	//only do this once
	//LM = Low moment pulse
	cTDEmSystem LM("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm");
	//HM = high moment pulse	
	cTDEmSystem HM("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm");

	//Load the system geometry (same for both moments)
	//This changes every fiducial/station
	TDEmGeometry G;
	G.tx_height = 30;
	G.tx_roll = 0;       G.tx_pitch = 0; G.tx_yaw = 0;
	G.txrx_dx = -12.62;  G.txrx_dy = 0; G.txrx_dz = +2.16;
	G.rx_roll = 0;       G.rx_pitch = 0; G.rx_yaw = 0;

	//Create the earth structure
	//This changes every fiducial/station
	Earth1D E(1);
	
	//Create a response object for each moment (they have different numbers of windwos)	
	TDEmResponse LMR;
	TDEmResponse HMR;
	double sum = 0.0;
	for (size_t j = 1; j <= 50; j++) {
		size_t nlayers = j;
		E = Earth1D(nlayers);
		double t1 = gettime();
		size_t nloops = 1000;
		for (size_t i = 0; i < nloops; i++) {
			//Run the forward model for each moment with random models
			for (size_t k = 0; k < nlayers; k++) E.conductivity[k] = urand(0.001, 2.0);
			for (size_t k = 0; k < nlayers - 1; k++) E.thickness[k] = urand(1.0, 10.0);
			LMR = LM.forward(E, G);
			HMR = HM.forward(E, G);
			sum += LMR.primary(XCOMP);//Just to make sure compiler does not optimize out the computation
			sum += HMR.primary(XCOMP);//Just to make sure compiler does not optimize out the computation
		}
		double t2 = gettime();
		std::cout << "Time: " << j << " " << t2 - t1 << std::endl;
	}	
	return 0;
}

int main(int argc, char* argv[]) {
	try {
		//skytem_example();
		//skytem_example_ip();
		skytem_computation_time();
	}
	catch (std::exception& e) {
		std::cout << e.what();
	}
	return 0;
}

