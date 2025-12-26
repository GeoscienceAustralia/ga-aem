/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

/* Example driver program for simple forward model*/
#include <exception>
#include <vector>
#include <cstring>
#include <iostream>
#include <ostream>

#include "logger.hpp"
#include "general_types.hpp"
#include "file_formats.hpp"
#include "outputmanager.hpp"
#include "vector_utils.hpp"
#include "general_utils.hpp"
#include "file_utils.hpp"
#include "random_utils.hpp"
#include "spectralaemsystem.hpp"
#include "tdemsystem.hpp"


class cLogger glog; //The global instance of the log file manager

using namespace AEM;


int skytem_example_ip() {
	//Load the AEM system specification files for the Skytem moments //only do this once				
	//TDEmSystem S("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm");	
	TDEmSystem S("..\\..\\examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm");
	

	//Load the system geometry 
	TDEmGeometry G;
	G.tx_height() = 40;
	G.tx_roll() = 0;       G.tx_pitch() = 0; G.tx_yaw() = 0;
	G.txrx_dx() = -13.35;  G.txrx_dy() = 0; G.txrx_dz() = +2.00;
	G.rx_roll() = 0;       G.rx_pitch() = 0; G.rx_yaw() = 0;

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
	TDEmResponse<double> R;
	//Run the forward model
	R = S.forward_model(E, G);
	for (size_t i = 0; i < R.nWindows(); i++){
		double wct = S.window(i).centre();
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
	G.tx_height() = 30;
	G.tx_roll() = 0;       G.tx_pitch() = 0; G.tx_yaw() = 0;
	G.txrx_dx() = -12.62;  G.txrx_dy() = 0; G.txrx_dz() = +2.16;
	G.rx_roll() = 1;       G.rx_pitch() = 0; G.rx_yaw() = 0;

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
	TDEmResponse<double> LMR;
	TDEmResponse<double> HMR;
				
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
	G.tx_height() = 30;
	G.tx_roll() = 0;       G.tx_pitch() = 0; G.tx_yaw() = 0;
	G.txrx_dx() = -12.62;  G.txrx_dy() = 0; G.txrx_dz() = +2.16;
	G.rx_roll() = 0;       G.rx_pitch() = 0; G.rx_yaw() = 0;

	//Create the earth structure
	//This changes every fiducial/station
	Earth1D E(1);

	//Create a response object for each moment (they have different numbers of windwos)	
	TDEmResponse<double> LMR;
	TDEmResponse<double> HMR;
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
			sum += LMR.primary(XCOMP,0);//Just to make sure compiler does not optimize out the computation
			sum += HMR.primary(XCOMP,0);//Just to make sure compiler does not optimize out the computation
		}
		double t2 = gettime();
		std::cout << "Time: " << j << " " << t2 - t1 << std::endl;
	}
	return 0;
};

void write_responses(const TDEmResponse<cdouble>& R, const TDEmResponse<cdouble>& R1, const TDEmResponse<cdouble>& DA, const TDEmResponse<cdouble>& DN, const TDEmResponse<cdouble>& PCD) {
	//using RT = TDEmResponse<cdouble>;
	std::ofstream ofr("C:/Users/rossc/Work/Tempest_Spectral/test/r.dat");
	std::ofstream ofr1("C:/Users/rossc/Work/Tempest_Spectral/test/r1.dat");
	std::ofstream ofda("C:/Users/rossc/Work/Tempest_Spectral/test/da.dat");
	std::ofstream ofdn("C:/Users/rossc/Work/Tempest_Spectral/test/dn.dat");
	std::ofstream ofpcd("C:/Users/rossc/Work/Tempest_Spectral/test/pcd.dat");

	R.S.simple_output(ofr);
	R1.S.simple_output(ofr1);
	DN.S.simple_output(ofdn);
	DA.S.simple_output(ofda);
	PCD.S.simple_output(ofpcd);
	return;
}

void write_responses(const TDEmResponse<double>& R, const TDEmResponse<double>& R1, const TDEmResponse<double>& DA, const TDEmResponse<double>& DN, const TDEmResponse<double>& PCD) {
	//using RT = TDEmResponse<cdouble>;
	std::ofstream ofr("C:/Users/rossc/Work/Tempest_Spectral/test/r.dat");
	std::ofstream ofr1("C:/Users/rossc/Work/Tempest_Spectral/test/r1.dat");
	std::ofstream ofda("C:/Users/rossc/Work/Tempest_Spectral/test/da.dat");
	std::ofstream ofdn("C:/Users/rossc/Work/Tempest_Spectral/test/dn.dat");
	std::ofstream ofpcd("C:/Users/rossc/Work/Tempest_Spectral/test/pcd.dat");

	R.S.simple_output(ofr);
	R1.S.simple_output(ofr1);
	DN.S.simple_output(ofdn);
	DA.S.simple_output(ofda);
	PCD.S.simple_output(ofpcd);
	return;
}

static void test_derivatives() {
	//using RT = TDEmResponse<double>;
	//fs::path stmpath = "../../examples/SkyTEM-BHMAR-2009/stmfiles/Skytem-LM.stm";
	//TDEmSystem T(stmpath);
	//TDEmGeometry G;
	//G.tx_height() = 30;
	//G.tx_roll() = 3;       G.tx_pitch() = 10; G.tx_yaw() = -5;
	//G.txrx_dx() = -12.62;  G.txrx_dy() = 12; G.txrx_dz() = +2.16;
	//G.rx_roll() = 4;       G.rx_pitch() = -3; G.rx_yaw() = 2;
	// The G.txrx_dy = 12 is so that we get some Y-component response
	
	//using RT = TDEmResponse<double>;
	//fs::path stmpath = "../../examples/Tempest-AusAEM-2020/stmfiles/Tempest-25.0Hz.stm";
	//TDEmSystem T(stmpath);

	using RT = TDEmResponse<cdouble>;
	fs::path stmpath = "C:/Users/rossc/Work/Tempest_Spectral/stmfiles/Tempest-Spectral.stm";
	SpectralAEMSystem T(stmpath);

	//Reference Geometry model
	TDEmGeometry G0;
	G0.tx_height() = 120;
	G0.tx_roll() = 3;    G0.tx_pitch() = 10; G0.tx_yaw() = -5;
	//G0.tx_roll() = 0;    G0.tx_pitch() = 0; G0.tx_yaw() = 0;
	G0.txrx_dx() = -40;  G0.txrx_dy() = 5; G0.txrx_dz() = -110;
	G0.rx_roll() = 4;    G0.rx_pitch() = -8; G0.rx_yaw() = 2;
	//G0.rx_roll() = 0;    G0.rx_pitch() = 0; G0.rx_yaw() = 0;
	const TDEmGeometry& G = G0;
	
	//Reference Earth model
	std::vector<double> c = { 0.01, 0.1, 0.01, 0.001};
	std::vector<double> t = { 10, 30, 20 };

	
	const Earth1D E(c, t);


	RT R;  // Forward response
	RT DA; // Analytic derivative
	RT DN; // Numerical derivative
	RT R1; // Numerical perturbation forward response
	RT PCD; // Percentage difference
	const double tol = 1e-7;
	CalculationType ctype;

	std::cout << "Conductivity derivatives" << std::endl;
	R = T.forward_model(E, G);
	for (size_t li = 0; li < c.size(); li++) {
		std::cout << "   Layer " << li << std::endl;
		ctype = CalculationType(CalculationType::Mode::DC, li);
		DA = T.derivative(G,ctype);

		std::vector<double> c1 = c;
		double delta = c[li] * 0.001;
		c1[li] = c[li] + delta;
		Earth1D E1(c1, t);
		R1 = T.forward_model(E1, G);
		DN = (R1 - R) / delta;
		PCD = percent_difference(DN, DA, tol);
		RT::display_max_abs_percent_difference(PCD);
		write_responses(R, R1, DA, DN, PCD);		
	};

	std::cout << "Thickness derivatives" << std::endl;
	for (size_t li = 0; li < c.size() - 1; li++) {
		std::cout << "   Layer " << li << std::endl;
		ctype = CalculationType(CalculationType::Mode::DT,li);
		DA = T.derivative(G,ctype);
		
		// Numerical derivative DN
		std::vector<double> t1 = t;
		double delta = t[li] * 0.01;
		delta = std::min(delta, 0.01);
		t1[li] = t[li] + delta;
		Earth1D E1(c, t1);
		R1 = T.forward_model(E1, G);
		DN = (R1 - R) / delta;
		PCD = percent_difference(DN, DA, tol);
		RT::display_max_abs_percent_difference(PCD);
		write_responses(R, R1, DA, DN, PCD);
	};

	TDEmGeometry G1; // Geometry perturbation
	double delta = 0.001; // 0.001 m
	
	std::cout << "DX derivative" << std::endl;
	ctype = CalculationType(CalculationType::Mode::DX);
	DA = T.derivative(G, ctype);
	G1 = G; G1.txrx_dx() += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA, tol);
	RT::display_max_abs_percent_difference(PCD);
	write_responses(R, R1, DA, DN, PCD);		

	std::cout << "DY derivative" << std::endl;
	ctype = CalculationType(CalculationType::Mode::DY);
	DA = T.derivative(G, ctype);
	G1 = G; G1.txrx_dy() += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA, tol); RT::display_max_abs_percent_difference(PCD);
	write_responses(R, R1, DA, DN, PCD);		

	std::cout << "DZ derivative" << std::endl;
	ctype = CalculationType(CalculationType::Mode::DZ);
	DA = T.derivative(G, ctype);
	G1 = G; G1.txrx_dz() += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA, tol); RT::display_max_abs_percent_difference(PCD);
	write_responses(R, R1, DA, DN, PCD);		

	std::cout << "DH derivative" << std::endl;
	ctype = CalculationType(CalculationType::Mode::DH);
	DA = T.derivative(G, ctype);
	G1 = G; // Tx moves independent of Rx
	G1.tx_height() += delta; G1.txrx_dz() -= delta;
	R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; 
	PCD = percent_difference(DN, DA, tol);
	RT::display_max_abs_percent_difference(PCD);
	write_responses(R, R1, DA, DN, PCD);		

	std::cout << "Tx height derivative" << std::endl;
	ctype = CalculationType(CalculationType::Mode::DTX_HEIGHT);
	DA = T.derivative(G, ctype);
	G1 = G; // Rx moves with Tx
	G1.tx_height() += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA, tol); RT::display_max_abs_percent_difference(PCD);
	write_responses(R, R1, DA, DN, PCD);		
	
	//Tx rotations
	delta = 0.001; //degree

	std::cout << "Tx roll derivative" << std::endl;
	ctype = CalculationType(CalculationType::Mode::DTX_ROLL);
	DA = T.derivative(G, ctype);
	G1 = G; G1.tx_roll() += delta; 
	R1 = T.forward_model(E, G1); 
	DN = (R1 - R) / delta;
	PCD = percent_difference(DN, DA, tol);
	//std::cout << R; std::cout << R1; std::cout << DN; std::cout << DA; std::cout << PCD;
	write_responses(R, R1, DA, DN, PCD); RT::display_max_abs_percent_difference(PCD);

	std::cout << "Tx pitch derivative" << std::endl;
	ctype = CalculationType(CalculationType::Mode::DTX_PITCH);
	DA = T.derivative(G, ctype);
	G1 = G; G1.tx_pitch() += delta;
	R1 = T.forward_model(E, G1);
	DN = (R1 - R) / delta;
	PCD = percent_difference(DN, DA, tol);
	//std::cout << R; std::cout << R1; std::cout << DN; std::cout << DA; std::cout << PCD;
	write_responses(R, R1, DA, DN, PCD); RT::display_max_abs_percent_difference(PCD);

	std::cout << "Tx yaw derivative" << std::endl;
	ctype = CalculationType(CalculationType::Mode::DTX_YAW);
	DA = T.derivative(G, ctype);
	G1 = G; G1.tx_yaw() += delta;
	R1 = T.forward_model(E, G1);
	DN = (R1 - R) / delta;
	PCD = percent_difference(DN, DA, tol);
	//std::cout << R; std::cout << R1; std::cout << DN; std::cout << DA; std::cout << PCD;
	write_responses(R, R1, DA, DN, PCD); RT::display_max_abs_percent_difference(PCD);

	//Rx rotations
	delta = 0.001; //degree

	std::cout << "Rx roll derivative" << std::endl;
	R = T.forward_model(E, G);
	ctype = CalculationType(CalculationType::Mode::DRX_ROLL);
	DA = T.derivative(G, ctype);
	G1 = G; G1.rx_roll() += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA, tol); RT::display_max_abs_percent_difference(PCD);
	write_responses(R, R1, DA, DN, PCD);

	std::cout << "Rx pitch derivative" << std::endl;
	R = T.forward_model(E, G);
	ctype = CalculationType(CalculationType::Mode::DRX_PITCH);
	DA = T.derivative(G, ctype);
	G1 = G; G1.rx_pitch() += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA, tol); RT::display_max_abs_percent_difference(PCD);
	write_responses(R, R1, DA, DN, PCD);

	std::cout << "Rx yaw derivative" << std::endl;
	R = T.forward_model(E, G);
	ctype = CalculationType(CalculationType::Mode::DRX_YAW);
	DA = T.derivative(G, ctype);
	G1 = G; G1.rx_yaw() += delta; R1 = T.forward_model(E, G1); DN = (R1 - R) / delta; PCD = percent_difference(DN, DA, tol); RT::display_max_abs_percent_difference(PCD);
	write_responses(R, R1, DA, DN, PCD);

 };

static void test_spectral() {
	fs::path stmpath = "C:/Users/rossc/Work/Tempest_Spectral/stmfiles/Tempest-Spectral.stm";
	SpectralAEMSystem S(stmpath);
	TDEmGeometry G;
	
	G.tx_height() = 120;
	G.tx_roll() = 0;	G.tx_pitch() = 0;	G.tx_yaw()  = 0;
	G.txrx_dx() = -110;	G.txrx_dy()  = 0;	G.txrx_dz() = -40;
	//G.rx_roll() = -10;	G.rx_pitch() = 10;	G.rx_yaw()  = 12;
	G.rx_roll() = 0;	G.rx_pitch() = 0;	G.rx_yaw() = 0;

	//std::vector<double> c = { 0.01, 0.1, 0.001 };
	//std::vector<double> c = { 0.01, 0.01, 0.01 };
	//std::vector<double> t = { 20, 30 };

	std::vector<double> c = { 0.01 };
	std::vector<double> t = {  };

	Earth1D E(c, t);
	TDEmResponse R = S.forward_model(E, G);
	TDEmVectorResponse<cdouble> TF = R.totalfield();
	std::ofstream ofst("C:/Users/rossc/Work/Tempest_Spectral/test/total.dat");
	std::ofstream ofsp("C:/Users/rossc/Work/Tempest_Spectral/test/primary.dat");
	std::ofstream ofss("C:/Users/rossc/Work/Tempest_Spectral/test/secondary.dat");

	TF.simple_output(std::cout);
	TF.simple_output(ofst);

	R.P.simple_output(ofsp);
	R.S.simple_output(ofss); 
	
};

static void test_simple() {		
	TDEmGeometry G;
	G.tx_height() = 30;
	G.tx_roll() = 0;	G.tx_pitch() = 0;	G.tx_yaw() = 0;
	G.txrx_dx() = -12;	G.txrx_dy() = -12;	G.txrx_dz() = 0;
	G.rx_roll() = 0;	G.rx_pitch() = 0;	G.rx_yaw() = 0;

	std::vector<double> c = { 0.1 };
	std::vector<double> t = {  };
	Earth1D E(c, t);

	std::vector<double> frequencies;
	for (double k = 0; k <= 6; k = k + 0.5) {
		frequencies.push_back(std::pow(10.0, k));
	}

	size_t numabscissa = 181;
	double modelling_loop_radius = 0.0;
	Vec3d tx_orientation(0,0,1);

	std::filesystem::path outdir = "C:/Users/rossc/Work/code/repos/deal_ii_tests/aem_solver/analytic";
	std::filesystem::create_directory(outdir);
	std::filesystem::path outfile = outdir / "1d_results.txt";
	std::ofstream ofs(outfile);
	const size_t np = 31;
	std::vector<Vec3d> rx_points(np);
	for (size_t pi = 0; pi < np; pi++) {
		Vec3d& p = rx_points[pi];
		p[0] = -150.0 + 10.0 * pi;
		p[1] = -12;
		p[2] = 30;
	}

	AEM::LEModeller S;
	S.initialise(frequencies, numabscissa, modelling_loop_radius);
	const double rxh = G.tx_height() + G.txrx_dz();
	S.set_earth(E);
	for (size_t pi = 0; pi < np; pi++) {
		const Vec3d rxpos = rx_points[pi];
		G.txrx_dx() = rxpos[0];
		G.txrx_dy() = rxpos[1];
		G.txrx_dz() = 0;
		//const double rxh = G.tx_height() + G.txrx_dz();

		S.set_geometry(tx_orientation, G.tx_height(), G.txrx_dx(), G.txrx_dy(), rxh);
		S.setup_computations();
		S.set_calculationtype(CMode::FM);
		double muzero = MUZERO<double>;
		for (size_t fi = 0; fi < frequencies.size(); fi++) {
			Vec3d  pf = 1e15 * muzero * S.primaryfield_inertial(tx_orientation);
			Vec3cd sf = 1e15 * muzero * S.secondaryfield_inertial(fi, tx_orientation);
			ofs << pi << ","
				<< rxpos[0] << ","
				<< rxpos[1] << ","
				<< rxpos[2] << ","
				<< frequencies[fi] << ","
				<< sf[0].real() << "," << sf[0].imag() << ","
				<< sf[1].real() << "," << sf[1].imag() << ","
				<< sf[2].real() << "," << sf[2].imag() << std::endl;
				//<< pf[0] << "," << pf[1] << "," << pf[2] << ","
		}
	}
	//prompttocontinue();
};

static int generate_synthetic_data() {
	using namespace IOManager;
	fs::path stmpath = "../stmfiles/Helitem-21m-25Hz-LM_25_w.stm";
	fs::path outfilepath = "../data/Helitem-21m-25Hz-LM_25_w.dat";
	makedirectory_for(outfilepath);

	TDEmSystem S(stmpath);

	double s = 0;
	TDEmGeometry G;
	double roll  = 5;
	double pitch = -10;
	double yaw   = 3;
	G.tx_height() = 30;
	G.tx_roll() = roll;   G.tx_pitch() = pitch;   G.tx_yaw()  = yaw;
	G.txrx_dx() = 0;      G.txrx_dy()  = 0;       G.txrx_dz() = 0;
	G.rx_roll() = roll;   G.rx_pitch() = pitch;   G.rx_yaw()  = yaw;
	std::cout << G.string();

	Earth1D E(3);
	E.conductivity[0] = 0.010;
	E.conductivity[1] = 0.100;
	E.conductivity[2] = 0.001;
	E.thickness[0] = 20;
	E.thickness[1] = 40;
	
	
	const size_t nw = S.nWindows();
	const size_t nl = E.nlayers();
	cASCIIOutputManager AM(outfilepath);
	cOutputField of_flight = cASCIIOutputManager::output_field("Flight", "Flight number", "", 1, cAsciiColumnFormat('I', 4, 0));
	cOutputField of_line = cASCIIOutputManager::output_field("Line", "Line number", "", 1, cAsciiColumnFormat('I', 7, 0));
	cOutputField of_fid  = cASCIIOutputManager::output_field("Fiducial", "Fiducial number", "", 1, cAsciiColumnFormat('I', 6, 0));
	cOutputField of_e = cASCIIOutputManager::output_field("Tx_Easting", "Transmitter Easting", "m", 1, cAsciiColumnFormat('F', 12, 2));
	cOutputField of_n = cASCIIOutputManager::output_field("Tx_Northing", "Transmitter Northing", "m", 1, cAsciiColumnFormat('F', 12, 2));
	cOutputField of_dtm = cASCIIOutputManager::output_field("DTM", "Ground elevation digital terrain model", "m", 1, cAsciiColumnFormat('F', 12, 2));

	std::vector<cOutputField> of_g;
	for (size_t gi = 0; gi < G.nelem(); gi++) {
		of_g.push_back(cASCIIOutputManager::output_field(G.element_name(gi), G.description(gi), G.units(gi), 1, cAsciiColumnFormat('F', 8, 3)));
	};

	cOutputField of_c = cASCIIOutputManager::output_field("Conductivity", "Conductivity", "S/m", nl, cAsciiColumnFormat('F', 10, 6));
	cOutputField of_t = cASCIIOutputManager::output_field("Thickness", "Thickness", "m", nl-1, cAsciiColumnFormat('F', 8, 2));

	cOutputField of_emx = cASCIIOutputManager::output_field("EMX", "EM X-component secondary field", "", nw, cAsciiColumnFormat('E', 14, 6));
	cOutputField of_emy = cASCIIOutputManager::output_field("EMY", "EM Y-component secondary field", "", nw, cAsciiColumnFormat('E', 14, 6));
	cOutputField of_emz = cASCIIOutputManager::output_field("EMZ", "EM Z-component secondary field", "", nw, cAsciiColumnFormat('E', 14, 6));
		
	AEM::TDEmResponse<double> R;
	bool status = AM.opendatafile();
	
	size_t np = 10;
	for (unsigned int pi = 0; pi < 10; pi++) {

		E.thickness[0] = 10 + 200 * (double)pi / (double)np;
		std::cout << E.thickness[0] << std::endl;
		R = S.forward_model(E, G);
		
		AM.begin_point_output();
		AM.writefield(pi, 100, of_flight);
		AM.writefield(pi, 10010, of_line);
		AM.writefield(pi, pi, of_fid);
		AM.writefield(pi, 200000+10*pi, of_e);
		AM.writefield(pi, 5000000.0, of_n);
		AM.writefield(pi, 0.0, of_dtm);

		for (size_t gi = 0; gi < G.nelem(); gi++) {
			AM.writefield(pi, G[gi], of_g[gi]);
		};

		AM.writefield(pi, E.conductivity, of_c);
		AM.writefield(pi, E.thickness, of_t);

		AM.writefield(pi, R.secondary(XCOMP), of_emx);
		AM.writefield(pi, R.secondary(YCOMP), of_emy);
		AM.writefield(pi, R.secondary(ZCOMP), of_emz);
		AM.end_point_output();

		if(pi==0) AM.end_first_record();  //Writes headers etc
	}

	for (size_t wi = 0; wi < nw; wi++) {
		std::cout << ixd(3) << wi;
		for (size_t ci = 0; ci < NCOMP; ci++) std::cout << exd(14, 6) << R.secondary(ci, wi);
		std::cout << std::endl;
	};

	return 0;
};

int main(int argc, char* argv[]) {
	try {
		skytem_example();
		//skytem_example_ip();
		//skytem_computation_time();
		//test_derivatives();
		//test_spectral();
		//test_simple();
		//generate_synthetic_data();
	}
	catch (std::exception& e) {
		std::cout << e.what();
	}
	return 0;
};

