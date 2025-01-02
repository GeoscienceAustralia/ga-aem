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

#include "vector_utils.hpp"
#include "general_utils.hpp"
#include "file_utils.hpp"
#include "random_utils.hpp"
#include "blocklanguage.hpp"
#include "lem.hpp"
#include "tdemsystem.hpp"
class cLogger glog; //The global instance of the log file manager

using namespace AEM;

int skytem_example_ip() {
	//Load the AEM system specification files for the Skytem moments //only do this once				
	//TDEmSystem S("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm");	
	TDEmSystem S("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm");
	

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
	R = S.forward_model(E, G);
	for (size_t i = 0; i < R.nWindows(); i++){
		double wct = S.window(i).centre_time();
		printf("%zu %10e %10e\n", i, wct, R.secondary(ZCOMP,i));
	}	
	return 0;
}

int skytem_example() {
	//Load the AEM system specification files for the Skytem moments
	//only do this once
	//LM = Low moment pulse
	TDEmSystem LM("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm");
	//HM = high moment pulse	
	TDEmSystem HM("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm");

	//Load the system geometry (same for both moments)
	//This changes every fiducial/station
	TDEmGeometry G;
	G.tx_height = 30;
	G.tx_roll = 0;       G.tx_pitch = 0; G.tx_yaw = 0;
	G.txrx_dx = -12.62;  G.txrx_dy = 0; G.txrx_dz = +2.16;
	G.rx_roll = 1;       G.rx_pitch = 0; G.rx_yaw = 0;

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
	LMR = LM.forward_model(E, G);
	HMR = HM.forward_model(E, G);
	
	for (size_t i = 0; i < LMR.nWindows(); i++){
		printf("%zu %g\n", i, LMR.secondary(ZCOMP,i));
	}

	for (size_t i = 0; i < HMR.nWindows(); i++){
		printf("%zu %g\n", i, HMR.secondary(ZCOMP,i));
	}
	
	return 0;
}

int skytem_computation_time() {
	//Load the AEM system specification files for the Skytem moments
	//only do this once
	//LM = Low moment pulse
	TDEmSystem LM("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm");
	//HM = high moment pulse	
	TDEmSystem HM("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm");

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
			LMR = LM.forward_model(E, G);
			HMR = HM.forward_model(E, G);
			sum += LMR.primary(XCOMP);//Just to make sure compiler does not optimize out the computation
			sum += HMR.primary(XCOMP);//Just to make sure compiler does not optimize out the computation
		}
		double t2 = gettime();
		std::cout << "Time: " << j << " " << t2 - t1 << std::endl;
	}
	return 0;
};

static void test_derivatives() {
	fs::path stmpath = "../../examples/SkyTEM-BHMAR-2009/stmfiles/Skytem-LM.stm";
	TDEmSystem T(stmpath);
	TDEmGeometry G;
	G.tx_height = 30;
	G.tx_roll = 3;       G.tx_pitch = 10; G.tx_yaw = -5;
	G.txrx_dx = -12.62;  G.txrx_dy = 12; G.txrx_dz = +2.16;
	G.rx_roll = 4;       G.rx_pitch = -3; G.rx_yaw = 2;
	// The G.txrx_dy = 12 is so that we get some Y-component response

	std::vector<double> c = { 0.01, 0.1, 0.001 };
	std::vector<double> t = { 200, 30 };
	Earth1D E(c, t);
	TDEmResponse R = T.forward_model(E, G);

	TDEmResponse DA; // Analytic derivative
	TDEmResponse DN; // Numerical derivative
	TDEmResponse R1; // Numerical perturbation response
	TDEmResponse PCD; // Percentage difference

	std::cout << "Conductivity derivatives" << std::endl;
	for (size_t li = 0; li < c.size(); li++) {
		// Analytic derivative DA
		DA = T.derivative(CalculationType(CalculationType::Mode::DC, li));

		// Numerical derivative DN
		std::vector<double> c1 = c;
		double delta = c[li] * 0.0001;
		c1[li] = c[li] + delta;
		E = Earth1D(c1, t);
		TDEmResponse R1 = T.forward_model(E, G);
		DN = (R1 - R) / delta;

		// Difference
		TDEmResponse PCD = percent_difference(DN, DA);
		TDEmResponse::display_max_abs_percent_difference(PCD);
	}

	std::cout << "Thickness derivatives" << std::endl;
	for (size_t li = 0; li < c.size()-1; li++) {
		// Analytic derivative DA
		TDEmResponse DA = T.derivative(CalculationType(CalculationType::Mode::DT, li));

		// Numerical derivative DN
		std::vector<double> t1 = t;
		double delta = t[li] * 0.001;
		t1[li] = t[li] + delta;
		E = Earth1D(c, t1);
		TDEmResponse R1 = T.forward_model(E, G);
		TDEmResponse DN = (R1 - R) / delta;

		// Difference
		TDEmResponse PCD = percent_difference(DN, DA);
		TDEmResponse::display_max_abs_percent_difference(PCD);
	}

	TDEmGeometry G1; // Geometry perturbation
	double delta = 0.001; // 0.001 m
	E = Earth1D(c, t);

	std::cout << "DX derivative" << std::endl;
	DA = T.derivative(CalculationType(CalculationType::Mode::DX));
	G1 = G; G1.txrx_dx += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA); TDEmResponse::display_max_abs_percent_difference(PCD);

	std::cout << "DY derivative" << std::endl;
	DA = T.derivative(CalculationType(CalculationType::Mode::DY));
	G1 = G; G1.txrx_dy += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA); TDEmResponse::display_max_abs_percent_difference(PCD);
	//std::cout << DA;
	//std::cout << DN;

	std::cout << "DZ derivative" << std::endl;
	DA = T.derivative(CalculationType(CalculationType::Mode::DZ));
	G1 = G; G1.txrx_dz += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA); TDEmResponse::display_max_abs_percent_difference(PCD);

	std::cout << "DH derivative" << std::endl;
	DA = T.derivative(CalculationType(CalculationType::Mode::DH));
	G1 = G; // Tx moves independent of Rx
	G1.tx_height += delta; G1.txrx_dz -= delta;
	R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA); TDEmResponse::display_max_abs_percent_difference(PCD);

	std::cout << "Tx Height derivative" << std::endl;
	DA = T.derivative(CalculationType(CalculationType::Mode::DTX_HEIGHT));
	G1 = G; // Rx moves with Tx
	G1.tx_height += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA); TDEmResponse::display_max_abs_percent_difference(PCD);
	
	delta = 0.001; // 0.001 degree
	std::cout << "Rx roll derivative" << std::endl;
	DA = T.derivative(CalculationType(CalculationType::Mode::DRX_ROLL), G, R);
	G1 = G; G1.rx_roll += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA); TDEmResponse::display_max_abs_percent_difference(PCD);

	std::cout << "Rx pitch derivative" << std::endl;
	DA = T.derivative(CalculationType(CalculationType::Mode::DRX_PITCH), G, R);
	G1 = G; G1.rx_pitch += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA); TDEmResponse::display_max_abs_percent_difference(PCD);
	
	std::cout << "Rx yaw derivative" << std::endl;
	DA = T.derivative(CalculationType(CalculationType::Mode::DRX_YAW), G, R);
	G1 = G; G1.rx_yaw += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA); TDEmResponse::display_max_abs_percent_difference(PCD);
 };

int main(int argc, char* argv[]) {
	try {
		//skytem_example();
		//skytem_example_ip();
		//skytem_computation_time();
		test_derivatives();
	}
	catch (std::exception& e) {
		std::cout << e.what();
	}
	return 0;
}

