/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/


#ifndef _fdemsystemclass_H
#define _fdemsystemclass_H

#include <cstring>
#include <vector>

#include "general_utils.h"
#include "blocklanguage.h"
#include "rollpitchyaw.h"
#include "earth1d.h"
#include "layeredearthmodeller.h"


//struct sFDEmInversionOptions {
//	bool solve_conductivity;
//	bool solve_thickness;
//	bool solve_birdheight;
//	double lambda;
//	double alpha;
//	double minimumphid;
//	double minimumimprovement;
//	int maxiterations;
//};
//
//struct sFDEmOutputOptions {
//	bool predictedbirdheight;
//
//	bool conductivity;
//	bool thickness;
//	bool depth;
//
//	bool observeddata;
//	bool predicteddata;
//};
//
//struct sFDEmResponse {
//	int     NC; //number of coilsets
//	double* SIP;//secondary inphase field
//	double* SQ;	//secondary quadrature field
//
//	double* EIP;//std dev inphase    noise estimates	
//	double* EQ; //std dev quadrature noise estimates	
//};

class cCoil {

public:
	Vec3d npos;//nominal position
	Vec3d pos;//position
	Vec3d naxis;//nominal axis
	Vec3d axis;//axis		

	void setrotation(const Mat3d& m) {				
		axis = m * naxis;
		pos  = m * npos;							
	}
};

class cFDEmCoilSet {

public:
	enum class Type { HCP, VCX, VCP, PER };

	double Frequency  = 0.0;
	double Separation = 0.0;
	Type Orientation  = Type::HCP;
	cCoil Tx;
	cCoil Rx;	

	cLayeredEarthModeller L;

	cFDEmCoilSet(const double& frequency, const cFDEmCoilSet::Type& orientation, const double& separation) {
		Frequency = frequency;
		L.setfrequency(frequency);

		Separation = separation;
		Tx.npos = Vec3d(Separation / 2.0, 0, 0);
		Rx.npos = Vec3d(-Separation / 2.0, 0, 0);

		Orientation = orientation;
		if (orientation == cFDEmCoilSet::Type::HCP) {
			Tx.naxis = Vec3d(0, 0, 1);
			Rx.naxis = Vec3d(0, 0, 1);
		}
		else if (orientation == cFDEmCoilSet::Type::VCX) {
			Tx.naxis = Vec3d(1, 0, 0);
			Rx.naxis = Vec3d(1, 0, 0);
		}
		else if (orientation == cFDEmCoilSet::Type::PER) {
			Tx.naxis = Vec3d(1, 0, 0);
			Rx.naxis = Vec3d(0, 1, 0);
		}
		else if (orientation == cFDEmCoilSet::Type::VCP) {
			Tx.naxis = Vec3d(0, 1, 0);
			Rx.naxis = Vec3d(0, 1, 0);
		}
		else {
			glog.errormsg(_SRC_, "Unknown coilset orientation\n");
		}
	}

	void setheightrotation(const double& height, const Mat3d& m) {
		Tx.setrotation(m);
		Rx.setrotation(m);
		Vec3d v = Tx.pos - Rx.pos;
		L.setxyzh(-v.x(), -v.y(), height + Tx.pos.z(), height + Rx.pos.z());		
	}
};

class cFDEmGeometry {

public:
	double Height=0.0;//Bird centre height above ground metres
	double Roll=0.0;//Bird roll degrees
	double Pitch=0.0;//Bird pitch degrees
	double Yaw=0.0;//Bird yaw degrees

	cFDEmGeometry() {};

	cFDEmGeometry(const double& height, const double& roll, const double& pitch, const double& yaw)
	{
		Height = height;
		Roll   = roll;
		Pitch  = pitch;
		Yaw    = yaw;
	};

	bool operator==(const cFDEmGeometry& rhs)
	{
		if (Height == rhs.Height &&
			Roll == rhs.Roll &&
			Pitch == rhs.Pitch && 
			Yaw == rhs.Yaw) return true;
		else return false;
	};

};

class cFDEmData {

public:
	std::vector<double> inphase;
	std::vector<double> quadrature;
};

class cFDEmNoiseModel {

public:
	enum class NoiseType { MAXPERCENTFLOOR, COMBINEPERCENTFLOOR };
	NoiseType noisetype;
	cvector   floor;
	cvector   percentage;
};

class cFDEmSystem {

private:
	
	cFDEmGeometry G;

	void setheightrollpitchyaw()
	{
		Eigen::Matrix<double, 3, 3> M = rollpitchyaw_matrix(D2R * G.Roll, D2R * G.Pitch, D2R * G.Yaw);
		for (size_t i = 0; i < NumberOfCoilSets(); i++) {
			CoilSets[i].setheightrotation(G.Height,M);
		}
	}

	
public:

	std::string SystemName;
	cBlock STM;
	
	std::vector<cFDEmCoilSet> CoilSets;

	cFDEmSystem() {};

	~cFDEmSystem() {};

	size_t NumberOfCoilSets() const {
		return CoilSets.size();
	}

	void readsystemdescriptorfile(std::string systemdescriptorfile)
	{
		STM.loadfromfile(systemdescriptorfile);
		SystemName = STM.getstringvalue("System.Name");
		std::vector<double> frequencies = STM.getdoublevector("System.Frequencies");
		std::vector<double> separations = STM.getdoublevector("System.Separations");
		std::vector<std::string> orientations = STM.getstringvector("System.Orientations");

		//NominalHeight = control.getdoublevalue("System.NominalHeight");
		std::vector<cFDEmCoilSet::Type> orientation;

		//Check system details	
		size_t nc = frequencies.size();
		if (separations.size() != nc || orientations.size() != nc) {
			glog.errormsg(_SRC_,"Error in system descriptor: Number of frequencies does not match number of seperations and/or orientations\n");
		}

		for (size_t i = 0; i < nc; i++) {
			if (strcasecmp("HCP", orientations[i]) == 0)orientation.push_back(cFDEmCoilSet::Type::HCP);
			else if (strcasecmp("VCX", orientations[i]) == 0)orientation.push_back(cFDEmCoilSet::Type::VCX);
			else if (strcasecmp("VCP", orientations[i]) == 0)orientation.push_back(cFDEmCoilSet::Type::VCP);
			else if (strcasecmp("PER", orientations[i]) == 0)orientation.push_back(cFDEmCoilSet::Type::PER);
			else {				
				glog.errormsg(_SRC_,"Error in system descriptor: %s is an unknown orientation\n", orientations[i].c_str());				
			}
		}

		//Add the coil sets				
		for (size_t i = 0; i < nc; i++) {			
			CoilSets.push_back(cFDEmCoilSet(frequencies[i], orientation[i], separations[i]));
		}
	};
	
	void setgeometry(const cFDEmGeometry& g)
	{		
		if (G == g) return;
		G = g;
		setheightrollpitchyaw();				
	}
	
	void setearth(const cEarth1D& e)
	{
		if (e.conductivity != CoilSets[0].L.getConductivity() || e.thickness != CoilSets[0].L.getThickness()) {
			for (size_t i = 0; i < NumberOfCoilSets(); i++) {
				CoilSets[i].L.setconductivitythickness(e.conductivity, e.thickness);
			}
		}
	}
	
	void setupcomputations()
	{
		for (size_t i = 0; i < NumberOfCoilSets(); i++) {
			CoilSets[i].L.setupcomputations();
		}
	}

	std::vector<double> p()
	{
		std::vector<double> v(NumberOfCoilSets());
		for (size_t i = 0; i < NumberOfCoilSets(); i++) {
			cFDEmCoilSet& C = CoilSets[i];
			v[i] = C.L.primary(C.Tx.axis, C.Rx.axis);
			if (C.Orientation == cFDEmCoilSet::Type::VCX) v[i] *= -1.0;
		}
		return v;
	}

	cvector s()
	{
		cvector v(NumberOfCoilSets());
		for (size_t i = 0; i < NumberOfCoilSets(); i++) {
			cFDEmCoilSet& C = CoilSets[i];
			v[i] = C.L.secondary(C.Tx.axis, C.Rx.axis);
			if (C.Orientation == cFDEmCoilSet::Type::VCX) v[i] *= -1.0;
		}
		return v;
	}

	cvector ppms()
	{
		cvector v(NumberOfCoilSets());
		for (size_t i = 0; i < NumberOfCoilSets(); i++) {
			cFDEmCoilSet& C = CoilSets[i];
			v[i] = C.L.ppm(C.Tx.axis, C.Rx.axis);
			if (C.Orientation == cFDEmCoilSet::Type::VCX) v[i] *= -1.0;
		}
		return v;
	}

	cvector dppms(const cLayeredEarthModeller::CalculationType& calculationtype, const size_t& derivativelayer)
	{
		cvector v(NumberOfCoilSets());
		for (size_t i = 0; i < NumberOfCoilSets(); i++) {
			cFDEmCoilSet& C = CoilSets[i];
			if (calculationtype == cLayeredEarthModeller::CalculationType::DB) {
				v[i]  = 2.0 * C.L.dppm(cLayeredEarthModeller::CalculationType::DH, derivativelayer, C.Tx.axis, C.Rx.axis);
			}
			else v[i] = C.L.dppm(calculationtype, derivativelayer, C.Tx.axis, C.Rx.axis);

			if (C.Orientation == cFDEmCoilSet::Type::VCX) v[i] *= -1.0;
		}
		return v;
	}
	
	cvector noiseestimates(const cvector& response, const cFDEmNoiseModel& noisemodel)
	{
		double anr, mnr, ani, mni;
		cvector v(NumberOfCoilSets());
		for (size_t i = 0; i < NumberOfCoilSets(); i++) {
			anr = noisemodel.floor[i].real();
			mnr = 0.01 * response[i].real() * noisemodel.percentage[i].real();
			ani = noisemodel.floor[i].imag();
			mni = 0.01 * response[i].imag() * noisemodel.percentage[i].imag();
			if (noisemodel.noisetype == cFDEmNoiseModel::NoiseType::MAXPERCENTFLOOR) {
				v[i] = cdouble(std::max(anr, mnr), std::max(ani, mni));
			}
			else if (noisemodel.noisetype == cFDEmNoiseModel::NoiseType::COMBINEPERCENTFLOOR) {
				v[i] = cdouble(std::sqrt(anr * anr + mnr * mnr), sqrt(ani * ani + mni * mni));
			}
			else {
				glog.errormsg(_SRC_,"unknown noise type\n");
			}
		}
		return v;
	}

	static std::vector<double> cv2dv(const cvector& cv)
	{
		std::vector<double> dv(cv.size() * 2);
		for (size_t i = 0; i < cv.size(); i++) {
			dv[i * 2] = cv[i].real();
			dv[i * 2 + 1] = cv[i].imag();
		}
		return dv;
	}

	static cvector dv2cv(const std::vector<double>& dv)
	{
		cvector cv(dv.size()/2);
		for (size_t i = 0; i < cv.size(); i++) {
			cv[i] = cdouble(dv[i * 2], dv[i * 2 + 1]);
		}
		return cv;
	}
};

#endif



