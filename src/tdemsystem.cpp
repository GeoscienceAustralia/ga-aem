/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <complex>

#include "general_utils.h"
#include "file_utils.h"
#include "lem.h"
#include "tdemsystem.h"
#include "vector_utils.h"

using namespace std;

cTDEmSystem::cTDEmSystem()
{
	initialise();
};

cTDEmSystem::cTDEmSystem(std::string systemdescriptorfile)
{
	initialise();
	readsystemdescriptorfile(systemdescriptorfile);
};

cTDEmSystem::~cTDEmSystem()
{
	if (fftwplan_backward){
		fftw_destroy_plan(fftwplan_backward);
	}
};

void cTDEmSystem::initialise()
{
	xaxis = cVec(1.0, 0.0, 0.0);
	yaxis = cVec(0.0, 1.0, 0.0);
	zaxis = cVec(0.0, 0.0, 1.0);
	LEM.calculation_type = CT_FORWARDMODEL;
	LEM.rzerotype = RZM_PROPOGATIONMATRIX;
	fftwplan_backward = 0;
}

void cTDEmSystem::createwaveform()
{
	NumberOfFFTFrequencies = SamplesPerWaveform / 2 + 1;
	fft_frequency.resize(NumberOfFFTFrequencies);

	size_t N = SamplesPerWaveform;
	size_t NC = N;
	size_t NR = 2 * (N / 2 + 1);

	//Forward transform
	F_Waveform.resize(NC);
	double* in = (double*)&(T_Waveform[0]);
	fftw_complex* out = (fftw_complex*)&(F_Waveform[0]);
	fftw_plan fftwplan_forward = fftw_plan_dft_r2c_1d((int)N, in, out, FFTW_ESTIMATE);
	for (size_t k = 0; k < SamplesPerWaveform; k++){
		T_Waveform[k] /= (double)SamplesPerWaveform;
	}
	fftw_execute(fftwplan_forward);
	fftw_destroy_plan(fftwplan_forward);

	for (size_t k = 0; k < NumberOfFFTFrequencies; k++){
		fft_frequency[k] = calculate_fft_frequency(k);
	}

	bool convert_B_2_dBdT = false;
	bool convert_dBdT_2_B = false;
	if (WaveformType == WT_TX){
		if (OutputType==OT_DBDT){
			convert_B_2_dBdT = true;
		}
	}

	if (WaveformType == WT_RX){
		if (OutputType==OT_BFIELD){
			convert_dBdT_2_B = true;
		}
	}

	for (size_t k = 0; k < NumberOfFFTFrequencies; k += 2){
		F_Waveform[k] = 0.0;
	}

	Transfer.resize(NumberOfFFTFrequencies);
	Transfer = F_Waveform;

	for (size_t k = 0; k < NumberOfFFTFrequencies; k++){
		if (convert_B_2_dBdT == true){
			Transfer[k] *= cdouble(0.0, -TWOPI*fft_frequency[k]);
		}
		if (convert_dBdT_2_B == true){
			Transfer[k] *= cdouble(0.0, -1.0 / (TWOPI*fft_frequency[k]));
		}
		
		for (size_t fi = 0; fi < LowPassFilterCutoffFrequency.size(); fi++){			
			cdouble a = 1.0 / cdouble(1.0, fft_frequency[k] / LowPassFilterCutoffFrequency[fi]);
			Transfer[k] *= std::pow(a, LowPassFilterOrder[fi]);
		}
	}

	//Frequencies to be splined
	NumberOfSplinedFrequencies = NumberOfFFTFrequencies / 2;
	SplinedFrequencieslog10.resize(NumberOfSplinedFrequencies);
	for (size_t k = 0; k < NumberOfSplinedFrequencies; k++){
		SplinedFrequencieslog10[k] = log10(fabs(fft_frequency[k * 2 + 1]));
	}

	//Setup inverse transform work array	
	FFTWork.resize(NR);

#if defined MULTITHREADED
	//FFTW_MEASURE does not seem to be thread safe
	unsigned int FLAGS = FFTW_ESTIMATE;		
#else
	unsigned int FLAGS = FFTW_MEASURE;
#endif

	fftw_complex* invin = (fftw_complex*)&(FFTWork[0]);
	double* invout = (double*)&(FFTWork[0]);
	fftwplan_backward = fftw_plan_dft_c2r_1d((int)N, invin, invout, FLAGS);
}

double cTDEmSystem::calculate_fft_frequency(size_t index)
{
	double deltaF = 1.0 / ((double)SamplesPerWaveform*SampleInterval);

	double s = (double)index;
	if (index > (SamplesPerWaveform / 2)){
		s = (double)index - (double)SamplesPerWaveform;
	}
	return s*deltaF;
}

void cTDEmSystem::setup_splines()
{
	HxR = std::vector<double>(NumberOfDiscreteFrequencies);
	HxI = std::vector<double>(NumberOfDiscreteFrequencies);
	HyR = std::vector<double>(NumberOfDiscreteFrequencies);
	HyI = std::vector<double>(NumberOfDiscreteFrequencies);
	HzR = std::vector<double>(NumberOfDiscreteFrequencies);
	HzI = std::vector<double>(NumberOfDiscreteFrequencies);

	HxR_spline = std::vector<double>(NumberOfDiscreteFrequencies);
	HxI_spline = std::vector<double>(NumberOfDiscreteFrequencies);
	HyR_spline = std::vector<double>(NumberOfDiscreteFrequencies);
	HyI_spline = std::vector<double>(NumberOfDiscreteFrequencies);
	HzR_spline = std::vector<double>(NumberOfDiscreteFrequencies);
	HzI_spline = std::vector<double>(NumberOfDiscreteFrequencies);

	size_t ns = NumberOfSplinedFrequencies;
	h2_spline = std::vector<double>(ns);
	a_spline = std::vector<double>(ns);
	b_spline = std::vector<double>(ns);
	klo_spline = std::vector<size_t>(ns);
	khi_spline = std::vector<size_t>(ns);
	a3ma_spline = std::vector<double>(ns);
	b3mb_spline = std::vector<double>(ns);

	X_splined = std::vector<cdouble>(ns);
	Y_splined = std::vector<cdouble>(ns);
	Z_splined = std::vector<cdouble>(ns);

	setup_splineinterp(DiscreteFrequenciesLog10, SplinedFrequencieslog10);

	LEM.setfrequencies(DiscreteFrequencies);
}
void cTDEmSystem::spline(const std::vector<double>& x, const std::vector<double>& y, double yp1, double ypn, std::vector<double>& y2)
{
	double p, qn, sig, un;

	size_t n = y2.size();
	std::vector<double> u(n - 1);
	if (yp1 > 0.99e30)
		y2[0] = u[0] = 0.0;
	else {
		y2[0] = -0.5;
		u[0] = (3.0 / (x[1] - x[0]))*((y[1] - y[0]) / (x[1] - x[0]) - yp1);
	}
	for (size_t i = 1; i < n - 1; i++) {
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig*y2[i - 1] + 2.0;
		y2[i] = (sig - 1.0) / p;
		u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
		u[i] = (6.0*u[i] / (x[i + 1] - x[i - 1]) - sig*u[i - 1]) / p;
	}
	if (ypn > 0.99e30)
		qn = un = 0.0;
	else {
		qn = 0.5;
		un = (3.0 / (x[n - 1] - x[n - 2]))*(ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
	}
	y2[n - 1] = (un - qn*u[n - 2]) / (qn*y2[n - 2] + 1.0);

	for (size_t k = n - 1; k-- > 0;){
		//Note the unusual syntax for decrement of unsigned variable
		y2[k] = y2[k] * y2[k + 1] + u[k];
	}
	
	//for (size_t k = n - 2; k >= 0; k--){
	//	y2[k] = y2[k] * y2[k + 1] + u[k];
	//	if (k == 0)break;
	//}
}
void cTDEmSystem::setup_splineinterp(const std::vector<double>& xn, const std::vector<double>& xi)
{
	size_t n = xn.size();
	for (size_t fi = 0; fi < NumberOfSplinedFrequencies; fi++){
		size_t klo = 0;
		size_t khi = n - 1;

		while (khi - klo > 1) {
			size_t k = (khi + klo) >> 1;
			if (xn[k] > xi[fi]){
				khi = k;
			}
			else klo = k;
		}
		double h = xn[khi] - xn[klo];
		if (h == 0.0){
			errormessage("setup_splineinterp(): Bad xa input to routine splintsetup\n");
		}
		double a = (xn[khi] - xi[fi]) / h;
		double b = (xi[fi] - xn[klo]) / h;

		a_spline[fi] = a;
		b_spline[fi] = b;
		h2_spline[fi] = h*h;
		klo_spline[fi] = klo;
		khi_spline[fi] = khi;
		a3ma_spline[fi] = a*a*a - a;
		b3mb_spline[fi] = b*b*b - b;
	}

}
void cTDEmSystem::spline_interp()
{
	//y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;			
	for (size_t k = 0; k < NumberOfSplinedFrequencies; k++){
		const double& a = a_spline[k];
		const double& b = b_spline[k];
		const double& a3 = a3ma_spline[k];
		const double& b3 = b3mb_spline[k];
		const size_t& klo = klo_spline[k];
		const size_t& khi = khi_spline[k];
		const double& scale = (h2_spline[k]) / 6.0;

		if (XScale != 0.0){
			double rv, iv;
			rv = (a*HxR[klo] + b*HxR[khi] + scale * (a3*HxR_spline[klo] + b3*HxR_spline[khi])),
			iv = (a*HxI[klo] + b*HxI[khi] + scale * (a3*HxI_spline[klo] + b3*HxI_spline[khi]));
			X_splined[k] = cdouble(rv, iv);
		}

		if (YScale != 0.0){
			double rv, iv;
			rv = (a*HyR[klo] + b*HyR[khi] + scale * (a3*HyR_spline[klo] + b3*HyR_spline[khi]));
			iv = (a*HyI[klo] + b*HyI[khi] + scale * (a3*HyI_spline[klo] + b3*HyI_spline[khi]));
			Y_splined[k] = cdouble(rv, iv);
		}

		if (ZScale != 0.0){
			double rv, iv;
			rv = (a*HzR[klo] + b*HzR[khi] + scale * (a3*HzR_spline[klo] + b3*HzR_spline[khi]));
			iv = (a*HzI[klo] + b*HzI[khi] + scale * (a3*HzI_spline[klo] + b3*HzI_spline[khi]));
			Z_splined[k] = cdouble(rv, iv);
		}
	}
}

void cTDEmSystem::setearthproperties(const cEarth1D& E)
{
	LEM.setproperties(E);
}

void cTDEmSystem::setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness)
{
	LEM.setconductivitythickness(conductivity, thickness);
}

void cTDEmSystem::setgeometry(const cTDEmGeometry& g)
{
	//X = +ve in flight direction
	//Y = +ve on left wing
	//Z = +ve vertical up
	//ie different to Fugro convention

	//Steer left is positive yaw     X->Y axis
	//Left wing up is positive roll  Y->Z axis
	//Nose down is positive pitch	 Z->X axis

	TX_height = g.tx_height;
	TX_pitch = g.tx_pitch;
	TX_roll = g.tx_roll;
	TX_orientation = cVec(0.0, 0.0, 1.0);
	if (TX_pitch != 0.0) TX_orientation = TX_orientation.rotate(TX_pitch, yaxis);
	if (TX_roll != 0.0) TX_orientation = TX_orientation.rotate(TX_roll, xaxis);

	TX_RX_separation = cVec(g.txrx_dx, g.txrx_dy, g.txrx_dz);

	LEM.setgeometry(LEM.pitchrolldipole(TX_pitch, TX_roll), TX_height, TX_RX_separation.x, TX_RX_separation.y, TX_height + TX_RX_separation.z);

	RX_height = TX_height + g.txrx_dz;
	RX_pitch = g.rx_pitch;
	RX_roll = g.rx_roll;
};
void cTDEmSystem::setgeometry(const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
{	
	cTDEmGeometry g(tx_height,tx_roll,tx_pitch,tx_yaw,txrx_dx,txrx_dy,txrx_dz,rx_roll,rx_pitch,rx_yaw);
	setgeometry(g);
}
void cTDEmSystem::setupcomputations()
{
	for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++){
		LEM.setfrequencyabscissalayers(fi);
	}
}
void cTDEmSystem::setprimaryfields()
{
	LEM.setprimaryfields();
	PrimaryX = LEM.Fields.t.p.x;
	PrimaryY = LEM.Fields.t.p.y;
	PrimaryZ = LEM.Fields.t.p.z;

	if (LEM.calculation_type == CT_HDERIVATIVE){
		//This is because when H changes Z also changes
		//and DZ = DH
		//but they should be all zero anyway
		PrimaryX *= 2.0;
		PrimaryY *= 2.0;
		PrimaryZ *= 2.0;
	}

	if (Normalisation == NT_PPM_PEAKTOPEAK){
		PrimaryX *= 2.0;
		PrimaryY *= 2.0;
		PrimaryZ *= 2.0;	
	}

	if (OutputType == OT_DBDT){
		//Must convert to dB/dt. This happens implicitly for the secondary via the waveform.
		PrimaryX *= TX_PeakdIdT;
		PrimaryY *= TX_PeakdIdT;
		PrimaryZ *= TX_PeakdIdT;
	}

	if (RX_pitch != 0.0 || RX_roll != 0.0){
		cVec field = cVec(PrimaryX, PrimaryY, PrimaryZ);
		cVec rotatedfield = rotatetoreceiverorientation(field);
		PrimaryX = rotatedfield.x;
		PrimaryY = rotatedfield.y;
		PrimaryZ = rotatedfield.z;
	}
			
	PrimaryX *= XScale;
	PrimaryY *= YScale;
	PrimaryZ *= ZScale;

}
cVec cTDEmSystem::rotatetoreceiverorientation(cVec v)
{
	//Rotating in opposite sense because we are rotating the axes
	if (RX_roll != 0.0) v = v.rotate(-RX_roll, xaxis);
	if (RX_pitch != 0.0) v = v.rotate(-RX_pitch, yaxis);
	return v;
}
void cTDEmSystem::setsecondaryfields()
{
	//Computation for discrete frequencies 	
	for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++){
		LEM.dointegrals(fi);
		LEM.setsecondaryfields(fi);
		const cdouble& x = LEM.Fields.t.s.x;
		const cdouble& y = LEM.Fields.t.s.y;
		const cdouble& z = LEM.Fields.t.s.z;
		cVec vr = cVec(x.real(), y.real(), z.real());
		cVec vi = cVec(x.imag(), y.imag(), z.imag());

		vr = rotatetoreceiverorientation(vr);
		vi = rotatetoreceiverorientation(vi);


		if (LEM.calculation_type == CT_HDERIVATIVE){
			//This is because when H changes Z also changes
			//and DZ = DH
			vr *= 2.0;
			vi *= 2.0;
		}

		HxR[fi] = vr.x;
		HxI[fi] = vi.x;
		HyR[fi] = vr.y;
		HyI[fi] = vi.y;
		HzR[fi] = vr.z;
		HzI[fi] = vi.z;
	}

	//Spline discreet frequencies		
	if (XScale != 0.0){
		spline(DiscreteFrequenciesLog10, HxR, 1e-30, 1e-30, HxR_spline);
		spline(DiscreteFrequenciesLog10, HxI, 1e-30, 1e-30, HxI_spline);
	}
	if (YScale != 0.0){
		spline(DiscreteFrequenciesLog10, HyR, 1e-30, 1e-30, HyR_spline);
		spline(DiscreteFrequenciesLog10, HyI, 1e-30, 1e-30, HyI_spline);
	}
	if (ZScale != 0.0){
		spline(DiscreteFrequenciesLog10, HzR, 1e-30, 1e-30, HzR_spline);
		spline(DiscreteFrequenciesLog10, HzI, 1e-30, 1e-30, HzI_spline);
	}

	//Interpolate 	
	spline_interp();
	
	if (SaveDiagnosticFiles){
		write_discretefrequencies("diag_discretefrequencies.txt");
		write_splinedfrequencies("diag_splinedfrequencies.txt");
		write_frequencydomainwaveform("diag_frequencydomainwaveform.txt");
	}

	//Filter - splining only every second value (even index) of the Waveform filter is always zero	
	if (XScale != 0.0){
		size_t n = 0;
		FFTWork = Transfer;
		for (size_t k = 1; k < NumberOfFFTFrequencies; k += 2){
			FFTWork[k] *= X_splined[n];			
			n++;
		}
		//Inverse FFT		
		fftw_execute(fftwplan_backward);		
		computewindow((double*)FFTWork.data(), X);
		if (SaveDiagnosticFiles){
			write_timesseries("diag_xtimeseries.txt");
		}
		X *= XScale;
	}

	if (YScale != 0.0){
		size_t n = 0;
		FFTWork = Transfer;
		for (size_t k = 1; k < NumberOfFFTFrequencies; k += 2){
			FFTWork[k] = Y_splined[n];
			n++;
		}
		//Inverse FFT		
		fftw_execute(fftwplan_backward);		
		computewindow((double*)FFTWork.data(), Y);
		if (SaveDiagnosticFiles){
			write_timesseries("diag_ytimeseries.txt");
		}
		Y *= YScale;
	}

	if (ZScale != 0.0){
		size_t n = 0;
		FFTWork = Transfer;
		for (size_t k = 1; k < NumberOfFFTFrequencies; k += 2){
			FFTWork[k] *= Z_splined[n];
			n++;
		}		
		//Inverse FFT				
		fftw_execute(fftwplan_backward);		
		computewindow((double*)FFTWork.data(), Z);
		if (SaveDiagnosticFiles){
			write_timesseries("diag_ztimeseries.txt");
		}
		Z *= ZScale;
	}

	if (SaveDiagnosticFiles){
		write_windows("diag_windows.txt");
	}

}

void cTDEmSystem::initialise_windows()
{
	NumberOfWindows = (size_t)STM.getintvalue("Receiver.NumberOfWindows");

	WinSpec.resize(NumberOfWindows);
	X.resize(NumberOfWindows);
	Y.resize(NumberOfWindows);
	Z.resize(NumberOfWindows);

	//Read window times
	dmatrix wt = STM.getdoublematrix("Receiver.WindowTimes");
	size_t nw = wt.size();
	if (nw != NumberOfWindows){
		errormessage("cTDEmSystem::initialise_windows: the number of WindowTimes does not match the NumberOfWindows\n");
	}
	for (size_t i = 0; i < NumberOfWindows; i++){
		if (wt[i].size() != 2){
			errormessage("cTDEmSystem::initialise_windows: the number of WindowTimes must have exactly 2 columns (error in window %lu)\n", i + 1);
		}
		WinSpec[i].TimeLow = wt[i][0];
		WinSpec[i].TimeHigh = wt[i][1];
	}

	WindowWeightingScheme = STM.getstringvalue("Receiver.WindowWeightingScheme");
	if (strcasecmp(WindowWeightingScheme, "AreaUnderCurve") == 0){
		initialise_windows_area();
	}
	else if (strcasecmp(WindowWeightingScheme, "Boxcar") == 0){
		initialise_windows_boxcar();
	}
	else if (strcasecmp(WindowWeightingScheme, "LinearTaper") == 0){
		initialise_windows_lineartaper();
	}
	else{
		errormessage("cTDEmSystem::readsystemdescriptorfile: WindowWeightingScheme %s unknown (must be \"AreaUnderCurve\" or  \"Boxcar\" or \"LinearTaper\")\n", WindowWeightingScheme.c_str());
	}
}

void cTDEmSystem::initialise_windows_area()
{
	double tlow, thigh, t, tp, tn, tleft, tright;

	double dwt = WaveformTime[1] - WaveformTime[0];
	double eps = 1.0e-7;

	for (size_t w = 0; w < NumberOfWindows; w++){
		WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
		for (size_t s = 0; s < SamplesPerWaveform; s++){
			t = WaveformTime[s];
			if (t + eps >= WinSpec[w].TimeLow){
				WinSpec[w].SampleLow = s;
				break;
			}
		}
		
		for (size_t s = SamplesPerWaveform; s-- > 0;){
			//Note the unusual syntax for decrement of unsigned variable
			t = WaveformTime[s];
			if (t - eps <= WinSpec[w].TimeHigh){
				WinSpec[w].SampleHigh = s;
				break;
			}			
		}

		WinSpec[w].NumberOfSamples = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
		WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
		WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);
		for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++){
			WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
			WinSpec[w].Weight[k] = 0.0;
		}

		tlow = WinSpec[w].TimeLow;
		thigh = WinSpec[w].TimeHigh;
		double wsum = 0;
		for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++){
			size_t s = WinSpec[w].Sample[k];
			t = WaveformTime[s];
			tp = WaveformTime[s] - dwt;
			tn = WaveformTime[s] + dwt;
			tleft = max(tp, tlow);
			tright = min(tn, thigh);
			WinSpec[w].Weight[k] = 0.5*(t - tleft) + 0.5*(tright - t);
			WinSpec[w].Weight[k] /= WinSpec[w].TimeWidth;
			wsum += WinSpec[w].Weight[k];
		}
		for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++){
			WinSpec[w].Weight[k] /= wsum;
		}
	}
}

void cTDEmSystem::initialise_windows_boxcar()
{
	double eps = 1.0e-7;		
	for (size_t w = 0; w < NumberOfWindows; w++){
		WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
		for (size_t s = 0; s < SamplesPerWaveform; s++){
			double t = WaveformTime[s];
			if (t + eps >= WinSpec[w].TimeLow){
				WinSpec[w].SampleLow = s;
				break;
			}
		}
						
		for (size_t s = SamplesPerWaveform; s-- > 0;){			
			//Note the unusual syntax for decrement of unsigned variable
			double t = WaveformTime[s];
			if (t - eps <= WinSpec[w].TimeHigh){
				WinSpec[w].SampleHigh = s;
				break;
			}				
		}
		

		WinSpec[w].NumberOfSamples = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
		WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
		WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);
		double weightsum = 0.0;
		for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++){
			WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
			WinSpec[w].Weight[k] = 1.0;
			weightsum += WinSpec[w].Weight[k];
		}
		for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++){
			WinSpec[w].Weight[k] /= weightsum;
		}	
		//printf("%lu %lu\n", w, WinSpec[w].NumberOfSamples);
	}
}

void cTDEmSystem::initialise_windows_lineartaper()
{
	double eps = 1.0e-7;
	for (size_t w = 0; w < NumberOfWindows; w++){
		WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
		for (size_t s = 0; s < SamplesPerWaveform; s++){
			double t = WaveformTime[s];
			if (t + eps >= WinSpec[w].TimeLow){
				WinSpec[w].SampleLow = s;
				break;
			}
		}
		
		for (size_t s = SamplesPerWaveform; s-- > 0;){
			//Note the unusual syntax for decrement of unsigned variable
			double t = WaveformTime[s];
			if (t - eps <= WinSpec[w].TimeHigh){
				WinSpec[w].SampleHigh = s;
				break;
			}
		}

		size_t ns = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
		WinSpec[w].NumberOfSamples = ns * 3;
		WinSpec[w].SampleLow -= ns;
		WinSpec[w].SampleHigh += ns;
		WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
		WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);

		double weightsum = 0.0;
		for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++){
			WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
			if (k < ns){
				WinSpec[w].Weight[k] = (double)(k + 1) / (double)(ns + 1);
			}
			else if (k >= 2 * ns){
				WinSpec[w].Weight[k] = 1.0 - (double)((k + 1) - 2 * ns) / (double)(ns + 1);
			}
			else{
				WinSpec[w].Weight[k] = 1.0;
			}
			weightsum += WinSpec[w].Weight[k];
		}
		for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++){
			WinSpec[w].Weight[k] /= weightsum;
		}
	}
}

void cTDEmSystem::computewindow(const double* timeseries, std::vector<double>& W)
{
	for (size_t w = 0; w < NumberOfWindows; w++){
		W[w] = 0.0;
		for (size_t k = 0; k < WinSpec[w].Sample.size(); k++){
			W[w] += timeseries[WinSpec[w].Sample[k]] * WinSpec[w].Weight[k];
		}
	}
}

void cTDEmSystem::printwindows()
{
	printf("Primary   %15.8lf%15.8lf%15.8lf\n\n", PrimaryX, PrimaryY, PrimaryZ);
	printf("Window#             X               Y               Z\n");
	for (size_t w = 0; w < NumberOfWindows; w++){
		printf("%2lu        %15.8lf%15.8lf%15.8lf\n", w + 1, X[w], Y[w], Z[w]);
	}
}

void cTDEmSystem::write_windows(const std::string& path)
{
	FILE* fp = fileopen(path, "w");
	//printf("Primary   %15.8lf%15.8lf%15.8lf\n\n",PrimaryX,PrimaryY,PrimaryZ);
	//printf("Window#             X               Y               Z\n");
	for (size_t w = 0; w < NumberOfWindows; w++){
		fprintf(fp, "%2lu\t%20e\t%20e\t%15e%15e%15e\n", w + 1, WinSpec[w].TimeLow, WinSpec[w].TimeHigh, X[w], Y[w], Z[w]);
	}
	fclose(fp);
}

void cTDEmSystem::write_timedomainwaveform(const std::string& path)
{
	FILE* fp = fileopen(path, "w");
	for (size_t i = 0; i < SamplesPerWaveform; i++){
		fprintf(fp, "%20le\t%20le\n", WaveformTime[i], T_Waveform[i]);
	}
	fclose(fp);
}

void cTDEmSystem::write_frequencydomainwaveform(const std::string& path)
{
	FILE* fp = fileopen(path, "w");
	for (size_t i = 0; i < NumberOfFFTFrequencies; i++){
		fprintf(fp, "%15le\t%15le\t%15le\t%15le\t%15le\n", fft_frequency[i], F_Waveform[i].real(), F_Waveform[i].imag(), Transfer[i].real(), Transfer[i].imag());
	}
	fclose(fp);
}

void cTDEmSystem::write_discretefrequencies(const std::string& path)
{
	FILE* fp = fileopen(path, "w");
	for (size_t i = 0; i < NumberOfDiscreteFrequencies; i++){
		fprintf(fp, "%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n", DiscreteFrequencies[i], HxR[i], HxI[i], HyR[i], HyI[i], HzR[i], HzI[i]);
	}
	fclose(fp);
}

void cTDEmSystem::write_splinedfrequencies(const std::string& path)
{
	FILE* fp = fileopen(path, "w");
	for (size_t i = 0; i < NumberOfSplinedFrequencies; i++){
		double f = pow10(SplinedFrequencieslog10[i]);
		fprintf(fp, "%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n", f, X_splined[i].real(), X_splined[i].imag(), Y_splined[i].real(), Y_splined[i].imag(), Z_splined[i].real(), Z_splined[i].imag());
	}
	fclose(fp);
}

void cTDEmSystem::write_frequencyseries(const std::string& path)
{
	FILE* fp = fileopen(path, "w");
	for (size_t i = 0; i < NumberOfFFTFrequencies; i++){
		fprintf(fp, "%15le\t%15le\t%15le\n", fft_frequency[i], FFTWork[i].real(), FFTWork[i].imag());
	}
	fclose(fp);
}

void cTDEmSystem::write_timesseries(const std::string& path)
{
	FILE* fp = fileopen(path, "w");
	double* ts = (double*)&(FFTWork[0]);
	for (size_t i = 0; i < SamplesPerWaveform; i++){
		fprintf(fp, "%20.10le\t%20.10le\n", WaveformTime[i], ts[i]);
	}
	fclose(fp);
}

void cTDEmSystem::forwardmodel(const std::vector<double>& conductivity, const std::vector<double>& thickness, const cTDEmGeometry& geometry)
{
	setconductivitythickness(conductivity, thickness);
	setgeometry(geometry);
	LEM.calculation_type = CT_FORWARDMODEL;
	LEM.derivative_layer = INT_MAX;

	setupcomputations();
	setprimaryfields();
	setsecondaryfields();
}

void cTDEmSystem::drx_pitch(double xb, double zb, double p, double& dxbdp, double& dzbdp)
{
	//xi = (  xb*cosp  + zb*sinp);Inertial
	//zi = ( -xb*sinp  + zb*cosp);
	//xb = (  xi*cosp  - zi*sinp);As bird sees it
	//zb = (  xi*sinp  + zi*cosp);						

	double cosp = cos(D2R*p);
	double sinp = sin(D2R*p);

	double xi = (xb*cosp + zb*sinp);//convert back to real coordinate system
	double zi = (-xb*sinp + zb*cosp);

	dxbdp = D2R*(-xi*sinp - zi*cosp);
	dzbdp = D2R*(+xi*cosp - zi*sinp);

}
void cTDEmSystem::drx_roll(double yb, double zb, double r, double& dybdr, double& dzbdr)
{
	//yi = (  yb*cosr  - zb*sinr);Inertial
	//zi = (  yb*sinr  + zb*cosr);
	//yb = (  yi*cosr  + zi*sinr);As bird sees it
	//zb = ( -yi*sinr  + zi*cosr);						

	double cosr = cos(D2R*r);
	double sinr = sin(D2R*r);

	double yi = (yb*cosr - zb*sinr);//convert back to real coordinate system
	double zi = (yb*sinr + zb*cosr);

	dybdr = D2R*(-yi*sinr + zi*cosr);
	dzbdr = D2R*(-yi*cosr - zi*sinr);

}

void cTDEmSystem::drx_pitch(const std::vector<double>& xb, const std::vector<double>& zb, double p, std::vector<double>& dxbdp, std::vector<double>& dzbdp)
{
	//xi = (  xb*cosp  + zb*sinp);Inertial
	//zi = ( -xb*sinp  + zb*cosp);
	//xb = (  xi*cosp  - zi*sinp);As bird sees it
	//zb = (  xi*sinp  + zi*cosp);						

	double cosp = cos(D2R*p);
	double sinp = sin(D2R*p);

	//convert back to real coordinate system
	std::vector<double> xi = (xb *  cosp + zb * sinp);
	std::vector<double> zi = (xb * -sinp + zb * cosp);

	dxbdp = (xi * -sinp - zi * cosp)*D2R;
	dzbdp = (xi * cosp - zi * sinp)*D2R;

}
void  cTDEmSystem::drx_roll(const std::vector<double>& yb, const std::vector<double>& zb, double r, std::vector<double>& dybdr, std::vector<double>& dzbdr)
{
	//yi = (  yb*cosr  - zb*sinr);Inertial
	//zi = (  yb*sinr  + zb*cosr);
	//yb = (  yi*cosr  + zi*sinr);As bird sees it
	//zb = ( -yi*sinr  + zi*cosr);						

	double cosr = cos(D2R*r);
	double sinr = sin(D2R*r);

	//convert back to real coordinate system
	std::vector<double> yi = (yb*cosr - zb*sinr);
	std::vector<double> zi = (yb*sinr + zb*cosr);

	dybdr = (yi * -sinr + zi * cosr)*D2R;
	dzbdr = (yi * -cosr - zi * sinr)*D2R;

}

void cTDEmSystem::readsystemdescriptorfile(std::string systemdescriptorfile)
{
	STM = cBlock(systemdescriptorfile);
	SystemName = STM.getstringvalue("Name");
	SystemType = STM.getstringvalue("Type");
	
	if (strcasecmp(SystemType, "Time Domain") != 0){
		errormessage("cTDEmSystem::readsystemdescriptorfile(): System Type is not Time Domain\n");
	}
	BaseFrequency = STM.getdoublevalue("Transmitter.BaseFrequency");
	BasePeriod = 1.0 / BaseFrequency;
	SampleFrequency = STM.getdoublevalue("Transmitter.WaveformDigitisingFrequency");

	TX_NumberOfTurns = STM.getdoublevalue("Transmitter.NumberOfTurns");
	TX_PeakCurrent = STM.getdoublevalue("Transmitter.PeakCurrent");
	TX_LoopArea = STM.getdoublevalue("Transmitter.LoopArea");
	
	bool wavformdefined = false;

	if (wavformdefined == false){
		std::string path = STM.getstringvalue("Transmitter.WaveformReceived.File");
		if (isundefined(path) == false){
			sFilePathParts fpp = getfilepathparts(systemdescriptorfile);
			dmatrix wp = readwaveformfile(fpp.directory + path);
			if (wp.size() > 0){
				digitisewaveform(wp, WaveformTime, WaveformReceived);
				WaveformType = WT_RX;
				T_Waveform = WaveformReceived;
				wavformdefined = true;
			}
		}
	}

	if (wavformdefined == false){
		std::string path = STM.getstringvalue("Transmitter.WaveformCurrent.File");
		if (isundefined(path) == false){
			sFilePathParts fpp = getfilepathparts(systemdescriptorfile);
			dmatrix wp = readwaveformfile(fpp.directory + path);
			if (wp.size() > 0){
				digitisewaveform(wp, WaveformTime, WaveformCurrent);
				WaveformType = WT_TX;
				T_Waveform = WaveformCurrent;
				wavformdefined = true;				
			}
		}
	}

	if (wavformdefined == false){
		dmatrix wp = STM.getdoublematrix("Transmitter.WaveformCurrent");
		if (wp.size() > 0){
			digitisewaveform(wp, WaveformTime, WaveformCurrent);
			WaveformType = WT_TX;
			T_Waveform = WaveformCurrent;
			wavformdefined = true;			
		}
	}

	if (wavformdefined == false){
		dmatrix wp = STM.getdoublematrix("Transmitter.WaveformReceived");
		if (wp.size() > 0){
			digitisewaveform(wp, WaveformTime, WaveformReceived);
			WaveformType = WT_RX;
			T_Waveform = WaveformReceived;
			wavformdefined = true;
		}
	}



	if (wavformdefined == false){
		errormessage("cTDEmSystem::readsystemdescriptorfile(): The waveform is not defined\n");
	}

	initialise_windows();
	
	LEM.ModellingLoopRadius = STM.getdoublevalue("ForwardModelling.ModellingLoopRadius");
	if (isundefined(LEM.ModellingLoopRadius)){
		LEM.ModellingLoopRadius = 0.0;
	}

	std::string ot = STM.getstringvalue("ForwardModelling.OutputType");
	if (strcasecmp(ot, "B") == 0){
		OutputType = OT_BFIELD;
	}
	else if (strcasecmp(ot, "dB/dt") == 0){
		OutputType = OT_DBDT;
	}
	else{		
		errormessage("cTDEmSystem::readsystemdescriptorfile(): OutputType %s unknown (must be one of \"B\" or \"dB/dt\")\n", ot.c_str());
	}

	FrequenciesPerDecade = (size_t)STM.getintvalue("ForwardModelling.FrequenciesPerDecade");
	if (FrequenciesPerDecade < 5){
		warningmessage("cTDEmSystem::readsystemdescriptorfile: It is wise to use at least 5 frequencies per decade\n");
	}

	LEM.NumAbscissa = (size_t)STM.getintvalue("ForwardModelling.NumberOfAbsiccaInHankelTransformEvaluation");

	std::string n = STM.getstringvalue("ForwardModelling.SecondaryFieldNormalisation");
	if (strcasecmp(n, "None")==0){
		Normalisation = NT_NONE;
	}
	else if (strcasecmp(n, "PPM") == 0){
		Normalisation = NT_PPM;
	}
	else if (strcasecmp(n, "PPMPEAKTOPEAK") == 0){
		Normalisation = NT_PPM_PEAKTOPEAK;
	}
	else{	
		errormessage("cTDEmSystem::readsystemdescriptorfile(): Normalisation %s unknown (must be one of \"None,PPM,PPMPEAKTOPEAK\")\n", n.c_str());
	}
	
	SaveDiagnosticFiles = STM.getboolvalue("ForwardModelling.SaveDiagnosticFiles");
	
	if (WaveformTime.size() <= 2 || WaveformTime.size() != T_Waveform.size()){
		errormessage("cTDEmSystem::readsystemdescriptorfile(): The number of WaveformTime values must match number of WaveformCurrent/WaveformReceived values and also be more than two\n");
	}

	if (LEM.NumAbscissa < 17){
		warningmessage("cTDEmSystem::readsystemdescriptorfile(): It is wise to use at least 17 Absicca for integrating the Hankel Transforms");
	}

	LowPassFilterCutoffFrequency = STM.getdoublevector("Receiver.LowPassFilter.CutOffFrequency");
	LowPassFilterOrder = STM.getdoublevector("Receiver.LowPassFilter.Order");
	systeminitialise();
};

void cTDEmSystem::systeminitialise()
{
	createwaveform();
	setupdiscretefrequencies();
	setup_splines();
	setup_scaling();	
}

double cTDEmSystem::compute_peak_didt()
{
	double maxdidt = 0.0;
	for (size_t i = 1; i < WaveformCurrent.size(); i++){
		double dt = WaveformTime[i] - WaveformTime[i - 1];
		double dc = WaveformCurrent[i] - WaveformCurrent[i - 1];
		double didt = std::fabs(dc / dt);
		//double didt = (dc / dt);
		if (didt > maxdidt) maxdidt = didt;
	}
	return maxdidt;
}

void cTDEmSystem::setup_scaling(){

	TX_PeakdIdT = compute_peak_didt();
	double txscale = MUZERO*TX_LoopArea*TX_NumberOfTurns*TX_PeakCurrent;
	double xos = STM.getdoublevalue("ForwardModelling.XOutputScaling");
	double yos = STM.getdoublevalue("ForwardModelling.YOutputScaling");
	double zos = STM.getdoublevalue("ForwardModelling.ZOutputScaling");
	
	XScale  = txscale*xos;
	YScale  = txscale*yos;
	ZScale  = txscale*zos;
	
	if (Normalisation == NT_PPM || Normalisation == NT_PPM_PEAKTOPEAK){
		cBlock b = STM.findblock("ReferenceGeometry");
		if (b.Entries.size() == 0){
			errormessage("cTDEmSystem::setup_ppm_normalisation(): Must define a ReferenceGeometry for PPM or PPMPEAKTOPEAK normalisation\n");
		}
		NormalizationGeometry = cTDEmGeometry(b);
		setgeometry(NormalizationGeometry);
		setprimaryfields();

		double s = 1.0;
		if (Normalisation == NT_PPM){
			s *= 1.0e6;
		}
		else if (Normalisation == NT_PPM_PEAKTOPEAK){
			//PrimaryX *= 2.0;
			//PrimaryY *= 2.0;
			//PrimaryZ *= 2.0;
			s *= 1.0e6;
		}

		if (PrimaryX == 0.0) XScale = 0.0;
		else XScale  *= (s/PrimaryX);

		if (PrimaryY == 0.0) YScale = 0.0;
		else YScale  *= (s/PrimaryY);

		if (PrimaryZ == 0.0) ZScale = 0.0;					
		else ZScale  *=  (s/PrimaryZ);		
	}
		
}

void cTDEmSystem::setupdiscretefrequencies()
{
	double lf1 = log10(BaseFrequency);
	double lf2 = log10(SampleFrequency / 2);
	double dlf = 1.0 / FrequenciesPerDecade;

	lf1 = lf1 - 2.0*dlf;
	lf2 = lf2 + 2.0*dlf;

	size_t nf = (size_t)ceil((lf2 - lf1)*FrequenciesPerDecade);
	dlf = (lf2 - lf1) / double(nf - 1);

	NumberOfDiscreteFrequencies = nf;
	FrequencyLog10Spacing = dlf;
	DiscreteFrequencyLow = pow(10.0, lf1);
	DiscreteFrequencyHigh = pow(10.0, lf2);

	DiscreteFrequenciesLog10 = std::vector<double>(NumberOfDiscreteFrequencies);
	DiscreteFrequencies = std::vector<double>(NumberOfDiscreteFrequencies);
	for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++){
		DiscreteFrequenciesLog10[fi] = log10(DiscreteFrequencyLow) + FrequencyLog10Spacing*(double)fi;
		DiscreteFrequencies[fi] = pow(10.0, DiscreteFrequenciesLog10[fi]);
	}
}

void cTDEmSystem::digitisewaveform(const dmatrix& wp, std::vector<double>& t, std::vector<double>& v)
{
	double hp = 0.5 / BaseFrequency;
	SampleInterval = 1.0 / SampleFrequency;
	SamplesPerWaveform = (size_t)(SampleFrequency / BaseFrequency);

	t.resize(SamplesPerWaveform);
	v.resize(SamplesPerWaveform);

	size_t np = wp.size();

	if (wp[np - 1][0] - wp[0][0] < hp){
		errormessage("cTDEmSystem::digitisecurrentwaveform: One complete halfcycle of the waveform has not been specified\n");
	}

	std::vector<double> x(np);
	std::vector<double> y(np);
	for (size_t i = 0; i < np; i++){
		x[i] = wp[i][0];
		y[i] = wp[i][1];
	}


	for (size_t i = 0; i < SamplesPerWaveform / 2; i++){
		bool set = false;

		double time = (double)i*SampleInterval;
		t[i] = time;

		if (time >= x[0] && time <= x[np - 1]){
			v[i] = linearinterp(x, y, time);
			set = true;
		}
		else if ((time - hp) >= x[0] && (time - hp) <= x[np - 1]){
			v[i] = -linearinterp(x, y, time - hp);
			set = true;
		}
		else if ((time - 2 * hp) >= x[0] && (time - 2 * hp) <= x[np - 1]){
			v[i] = linearinterp(x, y, time - 2.0*hp);
			set = true;
		}
		else if ((time + hp) >= x[0] && (time + hp) <= x[np - 1]){
			v[i] = -linearinterp(x, y, time + hp);
			set = true;
		}
		else if ((time + 2 * hp) >= x[0] && (time + 2 * hp) <= x[np - 1]){
			v[i] = linearinterp(x, y, time + 2.0*hp);
			set = true;
		}

		if (set == false){
			errormessage("cTDEmSystem::Error in waveform - not all defined\n");
		}
	}

	for (size_t i = SamplesPerWaveform / 2; i < SamplesPerWaveform; i++){
		t[i] = hp + t[i - SamplesPerWaveform / 2];
		v[i] =     -v[i - SamplesPerWaveform / 2];
	}
}

dmatrix cTDEmSystem::readwaveformfile(const std::string& filename)
{
	FILE* fp = fileopen(filename, "r");
	if (fp == NULL){
		errormessage("cTDEmSystem::readwaveformfile: Unable to open waveformfile %s\n", filename.c_str());
		throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
	}

	int num;
	dmatrix w;
	double time, value;
	while ((num = fscanf(fp, "%lf %lf\n", &time, &value)) == 2){
		std::vector<double> v(2);
		v[0] = time;
		v[1] = value;
		w.push_back(v);
	}
	fclose(fp);
	return w;

}

void cTDEmSystem::forwardmodel(const cTDEmGeometry& G, const cEarth1D& E, cTDEmResponse& R)
{	
	setgeometry(G);
	setearthproperties(E);	
	setupcomputations();
	setprimaryfields();
	setsecondaryfields();

	R.PX = PrimaryX;
	R.PY = PrimaryY;
	R.PZ = PrimaryZ;
	R.SX = X;
	R.SY = Y;
	R.SZ = Z;	
}