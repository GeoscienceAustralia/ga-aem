/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <mpi.h>


#include "radius_searcher.h"

#include "file_utils.h"
#include "file_formats.h"
#include "lem.h"
#include "tdemsystem.h"
#include "polygon.h"
#include "mpi_wrapper.h"
#include "petsc_wrapper.h"

#include "stopwatch.h"
#include "conductivity_logs.h"
#include "radius_searcher.h"
#include "inversion_line_searcher.h"
#include "stacktrace.h"

#define VERSION "1.0"

FILE* mylogfile = (FILE*)NULL;

class cSystemInfo;

class cComponentInfo;

enum eSmoothnessMethod { SM_1ST_DERIVATIVE, SM_2ND_DERIVATIVE };

class cLineTiles{

	std::vector<std::vector<double>> L;
	
public:

	cLineTiles(const double linelen, const double seglen, const double overlap){
		calculate_tiles(linelen, seglen, overlap);
	};

	void calculate_tiles(const double linelen, const double seglen, const double overlap){
		int n = std::round((linelen - overlap) / (seglen - overlap));		
		double newseglen = (linelen + 2.0*overlap*(n-1)) / (double)n;

		L.resize(n);				
		for (size_t i = 0; i < n; i++){
			L[i].resize(4);
			if (i == 0){		
				L[i][0] = 0.0;				
			}
			else{
				L[i][0] = L[i-1][0] + newseglen - 2*overlap;
			}
			
			L[i][1] = L[i][0] + overlap;
			L[i][2] = L[i][0] + newseglen - overlap;
			L[i][3] = L[i][0] + newseglen;

			if (i == 0){
				L[i][1] = L[i][0];
			}

			if (i == n-1){				
				L[i][2] = L[i][3];
			}			
		}
		std::cout << *this;
	};
	
	friend std::ostream& operator<<(std::ostream& os, const cLineTiles& t);
};

std::ostream& operator<<(std::ostream& os, const cLineTiles& t){
	for (size_t i = 0; i < t.L.size(); i++){
		os << t.L[i][0] << "\t";
		os << t.L[i][1] << "\t";
		os << t.L[i][2] << "\t";
		os << t.L[i][3] << std::endl;
	}
	return os;
};

class cField{

	std::vector<double> _data;
	size_t _ns;
	size_t _nb;

	bool defined;
	size_t nvals;
	std::string id;
	std::string def;
	size_t cn;
	bool cndef;
	bool ncndef;
	std::vector<double> defaultvec;


	double& el(const size_t& i, const size_t& j){
		return _data[i*_nb + j];
	};
	size_t cind(){ return cn - 1; }

public:

	cField(){
		defined = false;
		cndef = false;
		ncndef = false;
		nvals = 0;
	};

	cField(const cBlock& b, const std::string _id, size_t _size = 1){
		defined = false;
		cndef = false;
		ncndef = false;
		nvals = _size;

		id = _id;
		def = b.getstringvalue(_id);
		if (!isdefined(def)){
			nvals = 0;
			defined = false;
			return;
		}

		std::vector<std::string> t = tokenize(def);
		if (strcasecmp(t[0], "column") == 0){
			cndef = true;
			cn = (size_t)atoi(t[1].c_str());
		}
		else if (strcasecmp(t[0], "-column") == 0){
			ncndef = true;
			cn = (size_t)atoi(t[1].c_str());
		}
		else{
			for (size_t i = 0; i < t.size(); i++){
				defaultvec.push_back(atof(t[i].c_str()));
			}
			if (defaultvec.size() == 1){
				defaultvec.resize(nvals);
				for (size_t i = 1; i < defaultvec.size(); i++){
					defaultvec[i] = defaultvec[0];
				}
			}
		}
		defined = true;
	}

	void resize(const size_t ns, size_t nb = 1){
		_ns = ns;
		_nb = nb;
		_data.resize(ns*nb);
	};

	bool parse(const std::vector<std::string>& tokens, const size_t si){
		if (nvals == 0){
			return false;
		}
		else if (nvals == 1){
			el(si, 0) = get(tokens);
		}
		else{
			std::vector<double> v = getmulti(tokens);
			for (size_t bi = 0; bi < nvals; bi++){
				el(si, bi) = v[bi];
			}
		}
		return true;
	}

	double get(const std::vector<std::string>& tokens){
		double v;
		if (cndef) v = atof(tokens[cind()].c_str());
		else if (ncndef) v = -atof(tokens[cind()].c_str());
		else v = defaultvec[0];
		return v;
	}

	std::vector<double> getmulti(const std::vector<std::string>& tokens){
		std::vector<double> v(nvals);
		if (cndef){
			for (size_t k = 0; k < nvals; k++){
				v[k] = atof(tokens[cind() + k].c_str());
			}
		}
		else if (ncndef){
			for (size_t k = 0; k < nvals; k++){
				v[k] = -atof(tokens[cind() + k].c_str());
			}
		}
		else{
			for (size_t k = 0; k < nvals; k++){
				v[k] = defaultvec[k];
			}
		}
		return v;
	}

	double& operator()(const size_t& i, const size_t& j = 0){
		return el(i, j);
	};
};

class cComponentInfo {

public:
	
	bool Use = false;
	bool InvertTotalField = false;
	bool EstimateNoiseFromModel = false;
	size_t nw;
	size_t basedindex = 0;//sample data index for window 0
	cField fdp;
	cField fds;
	cField fdn;
	cField fdmn;
	cField fdan;	

	cComponentInfo() {};

	cComponentInfo(const cBlock& b, size_t nwindows, bool inverttotalfield)		
	{
		InvertTotalField = inverttotalfield;
		nw = nwindows;
		if (b.Entries.size() == 0){
			Use = false;
			return;
		}

		Use = b.getboolvalue("Use");
		if (Use == false)return;

		fdp = cField(b, "Primary");
		fds = cField(b, "Secondary", nwindows);

		EstimateNoiseFromModel = b.getboolvalue("EstimateNoiseFromModel");

		if (EstimateNoiseFromModel){
			fdmn = cField(b, "MultiplicativeNoise", nwindows);
			fdan = cField(b, "AdditiveNoise", nwindows);
		}
		else{
			fdn = cField(b, "Noise", nwindows);
		}

	}

	size_t ndata(){
		if (Use)return nw;
		else return 0;
	};

	bool allocate_data_arrays(const size_t nlocalsamples){
		if (Use == false)return true;
		fds.resize(nlocalsamples, nw);
		fdp.resize(nlocalsamples, 1);
		fdn.resize(nlocalsamples, nw);
		fdan.resize(nlocalsamples, nw);
		fdmn.resize(nlocalsamples, nw);
		return true;
	};

	bool parse(const std::vector<std::string> tokens, const size_t localsampleindex){
		if (Use == false)return true;
		fds.parse(tokens, localsampleindex);
		fdp.parse(tokens, localsampleindex);
		fdn.parse(tokens, localsampleindex);
		fdan.parse(tokens, localsampleindex);
		fdmn.parse(tokens, localsampleindex);
		return true;
	};

	std::vector<double> data(const size_t localsampleindex){
		std::vector<double> v;
		if (Use == false)return v;
		v.resize(nw);
		for (size_t wi = 0; wi < nw; wi++){
			v[wi] = fds(localsampleindex, wi);
			if (InvertTotalField){
				v[wi] += fdp(localsampleindex);
			}
		}
		return v;
	};

	std::vector<double> noise(const size_t localsampleindex){
		std::vector<double> v;
		if (Use == false)return v;

		v.resize(nw);
		if (EstimateNoiseFromModel){
			for (size_t wi = 0; wi < nw; wi++){
				double an = fdan(localsampleindex, wi);
				double pmn = fdmn(localsampleindex, wi);
				double s = fds(localsampleindex, wi);
				double mn = 0.01*pmn*s;
				v[wi] = sqrt(an*an + mn*mn);
			}
		}
		else{
			for (size_t wi = 0; wi < nw; wi++){
				v[wi] = fdn(localsampleindex, wi);
			}
		}
		return v;
	};

};

class cSystemInfo{

private:
	size_t nw;
	
public:
	cTDEmSystem T;
	bool InvertTotalField;	
	std::vector<cComponentInfo> Comp;
	cSystemInfo(){ 	};
	bool initialise(const cBlock& b){
		std::string stm = b.getstringvalue("SystemFile");
		T.readsystemdescriptorfile(stm);
		nw = T.NumberOfWindows;

		bool status;
		status = b.getvalue("InvertTotalField", InvertTotalField);
		if (status == false){
			InvertTotalField = false;
		}

		Comp.resize(3);
		Comp[0] = cComponentInfo(b.findblock("XComponent"), nw, InvertTotalField);
		Comp[1] = cComponentInfo(b.findblock("YComponent"), nw, InvertTotalField);
		Comp[2] = cComponentInfo(b.findblock("ZComponent"), nw, InvertTotalField);
		
		Comp[0].basedindex = 0;
		Comp[1].basedindex = 0;
		Comp[2].basedindex = 0;
		if (Comp[0].Use) Comp[1].basedindex += nw;		
		if (Comp[0].Use) Comp[2].basedindex += nw;
		if (Comp[1].Use) Comp[2].basedindex += nw;		
		return true;
	};

	size_t ndata(){
		return Comp[0].ndata() + Comp[1].ndata() + Comp[2].ndata();
	}

	inline size_t dindex(const size_t& component, const size_t& window){
		return Comp[component].basedindex + window;
	}

	bool allocate_data_arrays(const size_t nlocalsamples){
		for (size_t k = 0; k < Comp.size(); k++){
			Comp[k].allocate_data_arrays(nlocalsamples);
		}
		return true;
	}

	bool parse(const std::vector<std::string> tokens, const size_t localsampleindex){
		for (size_t ci = 0; ci < Comp.size(); ci++){
			Comp[ci].parse(tokens, localsampleindex);
		}
		return true;
	}

	std::vector<double> data(const size_t localsampleindex){
		std::vector<double> v;		
		for (size_t ci = 0; ci < Comp.size(); ci++){
			if (Comp[ci].Use){
				append(v, Comp[ci].data(localsampleindex));
			}
		}
		return v;
	}

	std::vector<double> noise(const size_t localsampleindex){
		std::vector<double> v;		
		for (size_t ci = 0; ci < Comp.size(); ci++){
			if (Comp[ci].Use){
				append(v, Comp[ci].noise(localsampleindex));
			}
		}
		return v;
	}

	bool forward_model(const std::vector<double>& conductivity, const std::vector<double>& thickness, const cTDEmGeometry& geometry){		
		T.setconductivitythickness(conductivity, thickness);
		T.setgeometry(geometry);
		T.LEM.calculation_type = CT_FORWARDMODEL;
		T.LEM.derivative_layer = INT_MAX;
		T.setupcomputations();
		T.setprimaryfields();
		T.setsecondaryfields();
		return true;
	}

	bool forward_model_and_derivatives(const std::vector<double>& conductivity, const std::vector<double>& thickness, const cTDEmGeometry& geometry, std::vector<double>& predicted, std::vector<std::vector<double>>& derivatives, const bool computederivatives, const std::vector<size_t> UGI){
		size_t nlayers = conductivity.size();
		T.setconductivitythickness(conductivity, thickness);
		T.setgeometry(geometry);
		T.LEM.calculation_type = CT_FORWARDMODEL;
		T.LEM.derivative_layer = INT_MAX;
		T.setupcomputations();
		T.setprimaryfields();
		T.setsecondaryfields();

		
		//Save for later derivative calculations
		std::vector<double> X = T.X;
		std::vector<double> Y = T.Y;
		std::vector<double> Z = T.Z;
		if (InvertTotalField){
			X += T.PrimaryX;
			Y += T.PrimaryY;
			Z += T.PrimaryZ;
		}
		
		predicted.resize(ndata());
		for (size_t ci = 0; ci < Comp.size(); ci++){
			if (Comp[ci].Use == false)continue;
			for (size_t wi = 0; wi < T.NumberOfWindows; wi++){
				predicted[dindex(ci,wi)] = T.secondary(ci, wi);
				if (Comp[ci].InvertTotalField){
					predicted[dindex(ci, wi)] += T.primary(ci);
				}				
			}
		}

		if (computederivatives == true){
			
			derivatives.resize(ndata());
			for (size_t di = 0; di < ndata(); di++){
				derivatives[di].resize(nlayers + UGI.size());
			}

			for (size_t li = 0; li < nlayers; li++){
				T.LEM.calculation_type = CT_CONDUCTIVITYDERIVATIVE;
				T.LEM.derivative_layer = li;
				T.setupcomputations();
				T.setprimaryfields();
				T.setsecondaryfields();
				
				for (size_t ci = 0; ci < Comp.size(); ci++){
					if (Comp[ci].Use == false)continue;
					for (size_t wi = 0; wi < T.NumberOfWindows; wi++){
						derivatives[dindex(ci, wi)][li] = T.secondary(ci, wi);
						if (Comp[ci].InvertTotalField){
							derivatives[dindex(ci, wi)][li] += T.primary(ci);
						}						
					}
				}
			}			

			for (size_t gi = 0; gi < UGI.size(); gi++){								
				if (cTDEmGeometry::elementtype(UGI[gi]) == GE_RX_PITCH){
					std::vector<double> dxbdp;
					std::vector<double> dzbdp;
					T.drx_pitch(X, Z, geometry.rx_pitch, dxbdp, dzbdp);					
					for (size_t ci = 0; ci < Comp.size(); ci++){
						if (Comp[ci].Use == false)continue;
						for (size_t wi = 0; wi < T.NumberOfWindows; wi++){
							if (ci == 0)      derivatives[dindex(ci, wi)][gi + nlayers] = dxbdp[wi];
							else if (ci == 1) derivatives[dindex(ci, wi)][gi + nlayers] = 0.0;
							else              derivatives[dindex(ci, wi)][gi + nlayers] = dzbdp[wi];							
						}
					}
				}
				else{
					T.LEM.calculation_type = cTDEmGeometry::derivativetype(UGI[gi]);
					T.LEM.derivative_layer = INT_MAX;
					T.setupcomputations();
					T.setprimaryfields();
					T.setsecondaryfields();						
					for (size_t ci = 0; ci < Comp.size(); ci++){
						if (Comp[ci].Use == false) continue;
						for (size_t wi = 0; wi < T.NumberOfWindows; wi++){
							derivatives[dindex(ci, wi)][gi + nlayers] = T.secondary(ci, wi);
							if (Comp[ci].InvertTotalField){
								derivatives[dindex(ci, wi)][gi + nlayers] += T.primary(ci);
							}														
						}
					}
				}				
			}
		}
		return true;
	}
};

class cEarthInfo{

private:
	size_t nlayers;

public:
	cField cref;
	cField tref;
	cField cstd;
	cField tstd;

	cEarthInfo(){};
	cEarthInfo(const cBlock& b){
		nlayers = b.getsizetvalue("NumberOfLayers");
		cref = cField(b, "ReferenceModel.Conductivity", nlayers);
		tref = cField(b, "ReferenceModel.Thickness", nlayers - 1);

		cstd = cField(b, "StdDevReferenceModel.Conductivity", nlayers);
		tstd = cField(b, "StdDevReferenceModel.Thickness", nlayers - 1);
	}

	size_t numlayers(){ return nlayers; }

	bool allocate_data_arrays(const size_t nlocalsamples){
		cref.resize(nlocalsamples, nlayers);
		tref.resize(nlocalsamples, nlayers - 1);
		cstd.resize(nlocalsamples, nlayers);
		tstd.resize(nlocalsamples, nlayers - 1);
		return true;
	}

	bool parse(const std::vector<std::string> tokens, const size_t localsampleindex){
		cref.parse(tokens, localsampleindex);
		tref.parse(tokens, localsampleindex);
		cstd.parse(tokens, localsampleindex);
		tstd.parse(tokens, localsampleindex);
		return true;
	}
};

class cInversionOptions{

public:
	double CorrelationRadius;
	double InverseDistancePower;
	double MinimumPhiD;
	double MinimumPercentageImprovement;
	size_t MaximumIterations;
	double AlphaV = 1;
	double AlphaH = 1;
	double AlphaB = 1;
	double AlphaR = 1;
	eSmoothnessMethod VerticalSmoothnessMethod = SM_2ND_DERIVATIVE;

	cInversionOptions(){};

	cInversionOptions(const cBlock& b){
		CorrelationRadius = b.getdoublevalue("CorrelationRadius");
		InverseDistancePower = b.getdoublevalue("InverseDistancePower");

		MinimumPhiD = b.getdoublevalue("MinimumPhiD");
		MinimumPercentageImprovement = b.getdoublevalue("MinimumPercentageImprovement");
		MaximumIterations = b.getsizetvalue("MaximumIterations");
		AlphaV = b.getdoublevalue("AlphaVertical");
		AlphaH = b.getdoublevalue("AlphaHorizontal");
		AlphaB = b.getdoublevalue("AlphaConductivityLogs");
		AlphaR = b.getdoublevalue("AlphaReferenceModel");

		std::string sm = b.getstringvalue("VerticalSmoothnessMethod");
		if (!isdefined(sm)){
			VerticalSmoothnessMethod = SM_2ND_DERIVATIVE;
		}
		else if (strcasecmp(sm, "Minimise1stDerivatives") == 0){
			VerticalSmoothnessMethod = SM_1ST_DERIVATIVE;
		}
		else if (strcasecmp(sm, "Minimize1stDerivatives") == 0){
			VerticalSmoothnessMethod = SM_1ST_DERIVATIVE;
		}
		else if (strcasecmp(sm, "Minimise2ndDerivatives") == 0){
			VerticalSmoothnessMethod = SM_2ND_DERIVATIVE;
		}
		else if (strcasecmp(sm, "Minimize2ndDerivatives") == 0){
			VerticalSmoothnessMethod = SM_2ND_DERIVATIVE;
		}
		else{
			rootmessage(mylogfile, "Unknown SmoothnessMethod %s\n", sm.c_str());
			std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw e;
		}
	}
};

class cInputOptions{

public:
	std::string DataFile;
	size_t HeaderLines = 0;
	size_t SubSample = 1;

	std::vector<cPolygon> IncludePolygons;
	std::vector<int> IncludeLines;
	std::vector<std::pair<int, int>> IncludeLineRanges;

	cInputOptions(){};

	cInputOptions(const cBlock& b){

		bool status;
		
		status = b.getvalue("HeaderLines", HeaderLines);
		if (status == false){
			HeaderLines = 0;
		}

		status = b.getvalue("Subsample", SubSample);
		if (status == false){
			SubSample = 1;
		}

		if (b.getvalue("DataFile", DataFile) == false){
			rootmessage(mylogfile, "Input DataFile was not specified\n");
			std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw e;
		}
		else if (exists(DataFile) == false){
			rootmessage(mylogfile, "Input DataFile %s not found\n", DataFile.c_str());
			std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw e;
		}

		std::string pfile;
		if (b.getvalue("IncludePolygon", pfile)){
			IncludePolygons.push_back(cPolygon(pfile));
		}

		std::string s;
		if (b.getvalue("IncludeLines", s)){
			parse_include_lines(s);
		}

	}

	void parse_include_lines(const std::string& s){
		std::pair<int, int> r;
		std::vector<std::string> t = tokenize(s);
		for (size_t i = 0; i < t.size(); i++){

			if (t[i][0] == ':' || t[i][t[i].size() - 1] == ':'){
				rootmessage(mylogfile, "Error bad token (%s) when parsing Input.IncludeLines\n", t[i].c_str());
				std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
				throw e;
			}
			else if (std::sscanf(t[i].c_str(), "%d:%d", &r.first, &r.second) == 2){
				IncludeLineRanges.push_back(r);
			}
			else if (std::sscanf(t[i].c_str(), "%d", &r.first) == 1){
				IncludeLines.push_back(r.first);
			}
			else{
				rootmessage(mylogfile, "Error bad token (%s) when parsing Input.IncludeLines\n", t[i].c_str());
				std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
				throw e;
			}
		}

	}

	bool is_line_included(const int& line){

		if (IncludeLines.size() == 0  && IncludeLineRanges.size() == 0) return true;

		auto rit = std::find_if(
			IncludeLineRanges.begin(),
			IncludeLineRanges.end(),
			[&line](const std::pair<int, int>& r)
		{
			return (line >= r.first && line <= r.second);
		}
		);
		if (rit != IncludeLineRanges.end())return true;
		
		auto lit = std::find(IncludeLines.begin(), IncludeLines.end(), line);
		if (lit != IncludeLines.end())return true;

		return false;
	}

	bool is_point_included(const cPoint& p)
	{
		bool insidestatus = true;
		for (size_t pi = 0; pi < IncludePolygons.size(); pi++){
			insidestatus = IncludePolygons[pi].isinside(p);
			if (insidestatus) break;
		}
		return insidestatus;
	};

};

class cOutputOptions{

public:
	std::string LogFile;
	std::string DataFile;
	bool PredictedData = true;
	bool ObservedData  = true;
	bool Noise         = true;
	bool PositiveLayerBottomDepths = false;
	bool NegativeLayerBottomDepths = false;
	bool InterfaceElevations = false;

	cOutputOptions(){};
	cOutputOptions(const cBlock& b){

		if (b.getvalue("LogFile", LogFile) == false){
			rootmessage(mylogfile, "Output LogFile was not specified\n");
			std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw e;
		}

		if (b.getvalue("DataFile", DataFile) == false){
			rootmessage(mylogfile, "Output DataFile was not specified\n");
			std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw e;
		}

		b.getvalue("PredictedData", PredictedData);
		b.getvalue("ObservedData", ObservedData);
		b.getvalue("Noise", Noise);
		b.getvalue("PositiveLayerBottomDepths", PositiveLayerBottomDepths);
		b.getvalue("NegativeLayerBottomDepths", NegativeLayerBottomDepths);
		b.getvalue("InterfaceElevations", InterfaceElevations);
	}
};

class cGeometryInfo {

public:

	bool   solve = false;
	cField ref;
	cField std;
	
	cGeometryInfo(){ };

	cGeometryInfo(const cBlock& b){
		ref  = cField(b, "Reference");
		b.getvalue("Solve", solve);
		if (solve){
			std = cField(b, "Uncertainty");
		}
	}

};

class cAllAtOnceInverter{

private:
	cBlock Control;
	
	cMpiEnv   mpienv;
	cMpiComm  mpicomm;
	cRadiusSearcher RS;

	cField fdsurvey;
	cField fddate;
	cField fdflight;
	cField fdline;
	cField fdfiducial;
	cField fdx;
	cField fdy;
	cField fdelevation;	
	std::vector<cGeometryInfo> G;
	std::vector<size_t> UGI;

	int mpisize;
	int mpirank;
	std::string mpipname;
	std::vector<cSystemInfo> T;
	cInputOptions InputOp;
	cOutputOptions OutputOp;
	cInversionOptions InversionOp;
	cEarthInfo E;

	std::vector<size_t> filerecordindex;
	size_t nlayers;
	size_t nchan;
	size_t nsamples;
	size_t ndata;

	size_t nparampersample;
	size_t nparam;
		
	size_t nlocaldata;
	size_t nlocalparam;
	cOwnership sown;

	cPetscDistVector dobs;
	cPetscDistVector dstd;
	cPetscDistVector mref;
	cPetscDistVector mstd;

	double ConductivityLogMaximumDistance = 500;
	double ConductivityLogPercentError = 5;
	std::vector<cConductivityLog> ConductivityLogs;
	cPetscDistVector clogref;//Conductivity logs ref
	cPetscDistVector clogstd;//Conductivity logs std

	double mLastPhiD = { 0.0 };
	double mLastLambda = { 0.0 };
	size_t mLastIteration = { 0 };

	cPetscDistMatrix J;//jacobian matrix	
	cPetscDistMatrix Wd;//data weights matrix			
	cPetscDistMatrix V;//vertical model regularization matrix
	cPetscDistMatrix H;//horizontal model regularization matrix
	cPetscDistMatrix B;//borehole log constraints
	cPetscDistMatrix Wb;
	cPetscDistMatrix P;//preconditioner matrix
	cPetscDistMatrix Wr;
	double lambda; //regularization parameter

public:

	cAllAtOnceInverter(int argc, char** argv)
	{
		if(argc < 2){
			rootmessage("%s\n", commandlinestring(argc, argv).c_str());
			rootmessage("%s\n", versionstring(VERSION, __TIME__, __DATE__).c_str());
			rootmessage("Usage: %s control_file_name\n", argv[0]);
			rootmessage("Too few command line arguments\n");
			exit(1);
		}
		else if(argc > 2){
			rootmessage("%s\n", commandlinestring(argc, argv).c_str());
			rootmessage("%s\n", versionstring(VERSION, __TIME__, __DATE__).c_str());
			rootmessage("Usage: %s control_file_name\n", argv[0]);
			rootmessage("Too many command line arguments\n");
			exit(1);
		}

		cStopWatch stopwatch;
		mpipname = mpienv.processor_name();
		mpicomm.set(MPI_COMM_WORLD);
		mpisize = mpicomm.size();
		mpirank = mpicomm.rank();

		std::string ControlFile = std::string(argv[1]);
		if(exists(ControlFile) == false){
			rootmessage("%s\n", commandlinestring(argc, argv).c_str());
			rootmessage("%s\n", versionstring(VERSION, __TIME__, __DATE__).c_str());
			rootmessage("Controlfile %s was not found\n", ControlFile.c_str());
			std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw(e);
		}

		Control = cBlock(ControlFile);
		OutputOp = cOutputOptions(Control.findblock("Output"));
		std::string s = strprint(".%04d", mpirank);
		
		//OutputOp.LogFile = insert_after_filename(OutputOp.LogFile, s);
		if (mpirank == 0){
			mylogfile = fileopen(OutputOp.LogFile, "w");
		}
		else{
			mylogfile = (FILE*)NULL;
		}

		rootmessage("Opening log file %s\n", OutputOp.LogFile.c_str());
		rootmessage(mylogfile, "Logfile opened on %s\n", timestamp().c_str());
		rootmessage(mylogfile, "Control file %s\n", Control.Filename.c_str());
		rootmessage(mylogfile, "Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		rootmessage(mylogfile, "Working directory %s\n", getcurrentdirectory().c_str());
		rootmessage(mylogfile, "Processes=%lu\tRank=%lu\n", mpisize, mpirank);
		rootmessage(mylogfile, "Processor name = %s\n", mpipname.c_str());
		if (mpirank == 0) Control.print();
		Control.write(mylogfile);		

		InputOp = cInputOptions(Control.findblock("Input"));
		InversionOp = cInversionOptions(Control.findblock("Options"));
		E = cEarthInfo(Control.findblock("Earth"));
		nlayers = E.numlayers();
		get_columns();
		get_geometry_columns();

		std::vector<cBlock> bv = Control.findblocks("EMSystem");
		T.resize(bv.size());
		for (size_t i = 0; i < bv.size(); i++){
			T[i].initialise(bv[i]);						
			std::string stmfile = bv[i].getstringvalue("SystemFile");
			std::string str = T[i].T.STM.get_as_string();
			rootmessage(mylogfile, "==============System file %s\n", stmfile.c_str());
			rootmessage(mylogfile, str.c_str());
		}
		rootmessage(mylogfile, "==========================================================================\n");
		nchan = calculate_nchan();
		rootmessage(mylogfile, "\nStarting setup\n");
		setup();
		rootmessage(mylogfile, "\nStarting iterations\n");
		iterate();
		rootmessage(mylogfile, "\nFinishing at at %s\n", timestamp().c_str());
		rootmessage(mylogfile, "Elapsed time = %.2lf\n", stopwatch.etimenow());
	};


	~cAllAtOnceInverter(){
		fclose(mylogfile);
	}

	bool get_columns(){

		cBlock b = Control.findblock("Input.Columns");
		fdsurvey = cField(b, "SurveyNumber");
		fddate = cField(b, "DateNumber");
		fdflight = cField(b, "FlightNumber");
		fdline = cField(b, "LineNumber");
		fdfiducial = cField(b, "FidNumber");
		fdx = cField(b, "Easting");
		fdy = cField(b, "Northing");
		fdelevation = cField(b, "GroundElevation");		
		return true;

	}

	bool get_geometry_columns(){

		G.resize(10);
		cBlock g = Control.findblock("Input.Geometry");
		for (size_t i = 0; i < G.size(); i++){
			std::string fname = cTDEmGeometry::fname(i);
			cBlock b = g.findblock(fname);	
			if (b.Name.size() == 0){				
				rootmessage(mylogfile,"Could not find block for geometry parameter %s\n", fname.c_str());
				std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);				
				throw(e);
			}
			G[i] = cGeometryInfo(b);
		}		
		UGI = unknown_geometry_indices();
		return true;
	}

	size_t count_closer_than(const std::vector<double> distance, const double& value){
		size_t nn = std::count_if(
			distance.begin(), distance.end(),
			[&value](const double& item) { return item <= value; }
		);
		return nn;
	}

	void read_conductivity_logs(){
		cBlock b = Control.findblock("ConductivityLogs");
		bool use = b.getboolvalue("Use");
		if (use){
			ConductivityLogMaximumDistance = b.getdoublevalue("MaximumDistance");
			ConductivityLogPercentError = b.getdoublevalue("PercentError");
			rootmessage(mylogfile, "Reading conductivity logs\n");
			std::string ldir = b.getstringvalue("Directory");
			auto flist = getfilelist(ldir, "con");
			for (size_t k = 0; k < flist.size(); k++){
				cConductivityLog clog(flist[k], true);
				cPoint p(clog.x, clog.y);
				if (InputOp.is_point_included(p)){
					std::vector<double> ndistances;
					std::vector<size_t> neighbours = RS.findneighbourstopoint(clog.x, clog.y, ndistances, ConductivityLogMaximumDistance);
					if (neighbours.size()>0){
						clog = cConductivityLog(flist[k], false);
						ConductivityLogs.push_back(clog);
					}
				}
			}
			rootmessage(mylogfile, "There are %lu conductivity logs available\n", flist.size());
			rootmessage(mylogfile, "There are %lu conductivity logs close enough to be included in the inversion\n", ConductivityLogs.size());
			for (size_t i = 0; i < ConductivityLogs.size(); i++){
				rootmessage(mylogfile, "%s\n", ConductivityLogs[i].infostring().c_str());
			}
		}
	}
	
	bool count_samples(){

		if (mpirank == 0){
			FILE* fp = fileopen(InputOp.DataFile, "r");
			if (fp == NULL){
				rootmessage(mylogfile, "Unable to open input DataFile %s\n", InputOp.DataFile.c_str());
				std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
				throw(e);
			}

			std::string s;
			for (size_t i = 0; i < InputOp.HeaderLines; i++){
				filegetline(fp, s);
			}

			size_t k = 0;
			while (filegetline(fp, s)){
				if (k % InputOp.SubSample == 0){
					std::vector<std::string> tokens = tokenize(s);
					int line = (int)fdline.get(tokens);
					if (InputOp.is_line_included(line)){
						cPoint p(fdx.get(tokens), fdy.get(tokens));
						if (InputOp.is_point_included(p)){
							filerecordindex.push_back(k);
						}
					}
				}
				k++;
			}
			fclose(fp);
		}
		mpicomm.bcast(filerecordindex);
		nsamples = filerecordindex.size();

		if (nsamples == 0){
			rootmessage(mylogfile, "There were no samples in the included lines and/or line ranges and/or polygon\n");
			std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw e;
		}

		ndata = calculate_ndata();
		nparampersample = calculate_nparampersample();		
		nparam = calculate_nparam();

		sown.set_petsc_default(mpisize, mpirank, nsamples);

		//Make sure the ownserships are intergral number of samples
		cOwnership down(sown.start*nchan, sown.end*nchan);
		cOwnership pown(sown.start*nparampersample, sown.end*nparampersample);
		nlocaldata  = down.nlocal();
		nlocalparam = pown.nlocal();
		return true;
	}

	bool allocate_data_arrays(){

		size_t nl = sown.nlocal();
		fdsurvey.resize(nl);
		fddate.resize(nl);
		fdflight.resize(nl);
		fdline.resize(nl);
		fdfiducial.resize(nl);
		fdx.resize(nl);
		fdy.resize(nl);
		fdelevation.resize(nl);

		for (size_t gi = 0; gi < G.size(); gi++){
			G[gi].ref.resize(nl);
			if (G[gi].solve){
				G[gi].std.resize(nl);
			}
		}
		for (size_t k = 0; k < T.size(); k++){
			T[k].allocate_data_arrays(nl);
		}
		E.allocate_data_arrays(nl);
		return true;
	};

	bool parse(const std::vector<std::string> tokens, const size_t localindex){
		fdsurvey.parse(tokens, localindex);
		fddate.parse(tokens, localindex);
		fdflight.parse(tokens, localindex);
		fdline.parse(tokens, localindex);
		fdfiducial.parse(tokens, localindex);
		fdx.parse(tokens, localindex);
		fdy.parse(tokens, localindex);
		fdelevation.parse(tokens, localindex);

		for (size_t gi = 0; gi < G.size(); gi++){
			G[gi].ref.parse(tokens, localindex);
			if (G[gi].solve){
				G[gi].std.parse(tokens, localindex);
			}
		}

		for (size_t ti = 0; ti < T.size(); ti++){
			T[ti].parse(tokens, localindex);
		}

		E.parse(tokens, localindex);
		return true;
	}

	bool read_data(){
		RS.x.resize(nsamples);
		RS.y.resize(nsamples);
		RS.elevation.resize(nsamples);

		allocate_data_arrays();

		FILE* fp = fileopen(InputOp.DataFile, "r");
		if (fp == NULL){
			rootmessage(mylogfile, "Unable to open input DataFile %s\n", InputOp.DataFile.c_str());
			std::string e = strprint("Error: exception thrown from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw e;
		}

		std::string s;
		size_t rec = 0;
		size_t gsi = 0;
		while (filegetline(fp, s)){
			if (rec == filerecordindex[gsi]){
				std::vector<std::string> tokens = tokenize(s);
				RS.x[gsi] = fdx.get(tokens);
				RS.y[gsi] = fdy.get(tokens);
				RS.elevation[gsi] = fdelevation.get(tokens);
				if (sown.owns(gsi)){
					size_t lsi = sown.localind(gsi);
					parse(tokens, lsi);
				}
				gsi++;
				if (gsi == nsamples)break;
			}
			rec++;
		}
		fclose(fp);

		RS.initialise(InversionOp.CorrelationRadius);
		return true;

	}

	void test_write_neighbours(){
		size_t k = 15;
		std::vector<double> ndistances;
		std::vector<size_t> n = RS.findneighbourstopoint(RS.x[k], RS.y[k], ndistances);
		//std::vector<size_t> n = RS.findneighbourstopoint(630325.2,6405732.5,ndistances,500.0);
		std::vector<double> xn(n.size());
		std::vector<double> yn(n.size());
		for (size_t i = 0; i < n.size(); i++){
			xn[i] = RS.x[n[i]];
			yn[i] = RS.y[n[i]];
		}

		if (mpirank == 0){
			write_xy("samples.txt", RS.x, RS.y);
			write_xy("neighbours.txt", xn, yn);
		}

	}

	void write_xy(const std::string& filename, const std::vector<double>& x, const std::vector<double>& y)
	{
		FILE* fp = fileopen(filename, "w");
		for (size_t i = 0; i < x.size(); i++){
			fprintf(fp, "%lf,%lf\n", x[i], y[i]);
		}
		fclose(fp);
	}

	size_t calculate_nchan(){
		size_t n = 0;
		for (size_t i = 0; i < T.size(); i++){
			n += T[i].ndata();
		}
		return n;
	}

	inline size_t calculate_ndata(){
		return nsamples*nchan;
	}

	inline size_t calculate_nparampersample(){
		return nlayers + UGI.size();				
	}

	size_t calculate_nparam(){
		return nsamples*nparampersample;
	}

	size_t dindex(const size_t& iglobalsample, const size_t& ichan){
		return iglobalsample*nchan + ichan;
	}

	size_t gpindex_c(const size_t& iglobalsample, const size_t& ilayer){
		return iglobalsample * nparampersample + ilayer;
	}

	size_t gpindex_g(const size_t& iglobalsample, const size_t& igparam){
		return iglobalsample*nparampersample + nlayers + igparam;
	}

	size_t localsampleindex(const size_t& _globalsampleindex){
		return _globalsampleindex - (size_t)sown.start;
	}

	size_t globalsampleindex(const size_t& _localsampleindex){
		return _localsampleindex + (size_t)sown.start;
	}

	std::vector<double> local_data(){
		std::vector<double> v;
		v.reserve(nlocaldata);
		for (size_t si = 0; si < (size_t)sown.nlocal(); si++){
			for (size_t ti = 0; ti < T.size(); ti++){
				append(v,T[ti].data(si));				
			}
		}
		return v;
	}

	std::vector<double> local_noise(){
		std::vector<double> v;
		v.reserve(nlocaldata);		
		for (size_t si = 0; si < (size_t)sown.nlocal(); si++){
			for (size_t ti = 0; ti < T.size(); ti++){
				append(v,T[ti].noise(si));				
			}
		}
		return v;
	}

	std::vector<double> local_mref(){
		std::vector<double> v(nlocalparam);
		size_t i = 0;
		for (size_t si = 0; si < (size_t)sown.nlocal(); si++){
			for (size_t li = 0; li < nlayers; li++){
				v[i] = std::log10(E.cref(si, li));
				i++;
			}
			for (size_t gi = 0; gi < UGI.size(); gi++){				
				v[i] = G[UGI[gi]].ref(si);
				i++;
			}
		}
		return v;
	}

	std::vector<double> local_mstd(){
		std::vector<double> v(nlocalparam);
		size_t i = 0;
		for (size_t si = 0; si < (size_t)sown.nlocal(); si++){
			for (size_t li = 0; li < nlayers; li++){
				v[i] = E.cstd(si, li);
				i++;
			}
			for (size_t gi = 0; gi < UGI.size(); gi++){				
				v[i] = G[UGI[gi]].std(si);
				i++;			
			}
		}
		return v;
	}

	void J_create(){

		rootmessage(mylogfile, "Creating matrix J\n");
		J.create_sparse("J", mpicomm, nlocaldata, nlocalparam, ndata, nparam);
		std::vector<PetscInt> d_nnz(J.nlocalrows(), 0);
		std::vector<PetscInt> o_nnz(J.nlocalrows(), 0);
		for (size_t si = (size_t)sown.start; si < (size_t)sown.end; si++){
			for (size_t ci = 0; ci < nchan; ci++){
				PetscInt di = dindex(si, ci);
				PetscInt lri = J.lri(di);
				for (size_t li = 0; li < nlayers; li++){
					PetscInt pi = gpindex_c(si, li);
					if (J.inownerdiagonalblock(di, pi)){
						d_nnz[lri]++;
					}
					else o_nnz[lri]++;;
				}
				for (size_t gi = 0; gi < UGI.size(); gi++){					
					PetscInt pi = gpindex_g(si, gi);
					if (J.inownerdiagonalblock(di, pi)){
						d_nnz[lri]++;
					}
					else o_nnz[lri]++;				
				}

			}
		}
		J.preallocate(d_nnz, o_nnz);
	}

	void J_set_nonzero_pattern(){
		rootmessage(mylogfile, "Assembling matrix G\n");
		PetscErrorCode ierr;
		for (size_t si = (size_t)sown.start; si < (size_t)sown.end; si++){
			for (size_t ci = 0; ci < nchan; ci++){
				PetscInt di = dindex(si, ci);
				for (size_t li = 0; li < nlayers; li++){
					PetscInt pi = gpindex_c(si, li);
					ierr = MatSetValue(J.mat(), di, pi, 0.0, INSERT_VALUES); CHKERR(ierr);
				}
				for (size_t gi = 0; gi < UGI.size(); gi++){
					PetscInt pi = gpindex_g(si, gi);
					ierr = MatSetValue(J.mat(), di, pi, 0.0, INSERT_VALUES); CHKERR(ierr);
				}
			}
		}
		J.assemble();
		rootmessage(mylogfile, "Finished assembling matrix G\n");
		return;
	}

	void Wd_create(){
		rootmessage(mylogfile, "Creating matrix Wd\n");
		Wd.create_diagonal_to_power("Wd", dstd, -2.0);
		Wd *= (1.0 / Wd.nglobalrows());
		rootmessage(mylogfile, "Finished creating matrix Wd\n");
		return;
	};

	void Wr_create(){
		rootmessage(mylogfile, "Creating matrix Wr\n");
		Wr.create_diagonal_to_power("Wr", mstd, -2.0);
		Wr *= (InversionOp.AlphaR / Wr.nglobalrows());
		rootmessage(mylogfile, "Finished creating matrix Wr\n");
		return;
	};

	void V_create_1st_derivative(){

		rootmessage(mylogfile, "Creating matrix V\n");
		
		size_t nglobalconstraints = (nlayers - 1)*nsamples;
		size_t nlocalconstraints = (nlayers - 1)*sown.nlocal();

		V.create_sparse("V", mpicomm, nlocalconstraints, nlocalparam, nglobalconstraints, nparam);
		V.preallocate(2, 0);
		PetscInt gri = V.gri(0);
		for (size_t si = (size_t)sown.start; si < (size_t)sown.end; si++){
			for (size_t li = 0; li < nlayers - 1; li++){
				PetscInt pa = gpindex_c(si, li);
				PetscInt pb = gpindex_c(si, li + 1);
				V.set(gri, pa, 1.0);
				V.set(gri, pb, -1.0);
				gri++;
			}
		}
		V.assemble();
		V *= std::sqrt(InversionOp.AlphaV / (double)V.nglobalrows());
		rootmessage(mylogfile, "Finished creating matrix V\n");
	};

	void V_create_2nd_derivative(){

		rootmessage(mylogfile, "Creating matrix V\n");
		size_t nglobalconstraints = (nlayers - 2)*nsamples;
		size_t nlocalconstraints = (nlayers - 2)*sown.nlocal();

		V.create_sparse("V", mpicomm, nlocalconstraints, nlocalparam, nglobalconstraints, nparam);
		V.preallocate(3, 0);
		PetscInt gri = V.gri(0);
		for (size_t si = (size_t)sown.start; si < (size_t)sown.end; si++){
			for (size_t li = 1; li < nlayers - 1; li++){
				PetscInt pa = gpindex_c(si, li - 1);
				PetscInt pb = gpindex_c(si, li);
				PetscInt pc = gpindex_c(si, li + 1);
				V.set(gri, pa, 1.0);
				V.set(gri, pb, -2.0);
				V.set(gri, pc, 1.0);
				gri++;
			}
		}
		V.assemble();
		V *= std::sqrt(InversionOp.AlphaV / (double)V.nglobalrows());
		rootmessage(mylogfile, "Finished creating matrix V\n");
	};

	void H_create_elevation(){

		rootmessage(mylogfile, "Creating matrix H\n");
		std::vector<double> t = get_thicknesses_ref(0);
		std::vector<double> d = get_interface_depths(t);
		d.push_back(d.back() + t.back());

		H.create_sparse("H", mpicomm, nlocalparam, nlocalparam, nparam, nparam);
		std::vector<PetscInt> d_nnz(H.nlocalrows(), 0);
		std::vector<PetscInt> o_nnz(H.nlocalrows(), 0);
		size_t nsum = 0;
		for (size_t gsi = (size_t)sown.start; gsi < (size_t)sown.end; gsi++){
			double elev = RS.elevation[gsi];
			std::vector<double> ndistance;
			std::vector<size_t> neighbours = RS.findneighbours(gsi, ndistance);
			nsum += neighbours.size();
			for (size_t li = 0; li < nlayers; li++){
				PetscInt gri = gpindex_c(gsi, li);
				PetscInt lri = H.lri(gri);
				d_nnz[lri]++;//non-zero for this samples layer
				for (size_t ni = 0; ni < neighbours.size(); ni++){
					size_t ngsi = neighbours[ni];
					double nelev = RS.elevation[ngsi];
					std::vector<double> nd = d + (elev - nelev);

					std::vector<double> fo = fractionaloverlaps(d[li], d[li + 1], nd);
					for (size_t nli = 0; nli < nlayers; nli++){
						if (fo[nli]>0){
							PetscInt gci = gpindex_c(ngsi, nli);
							if (H.inownerdiagonalblock(gri, gci)) d_nnz[lri]++;
							else o_nnz[lri]++;
						}
					}
				}
			}
		}
		H.preallocate(d_nnz, o_nnz);
		nsum = mpicomm.sum(nsum);
		rootmessage(mylogfile, "Average number of neighbours per sample = %lu\n", nsum / nsamples);

		//Set entries
		rootmessage(mylogfile, "Assembling matrix H\n");
		//loop over each sample
		int count = 0;
		for (size_t gsi = (size_t)sown.start; gsi < (size_t)sown.end; gsi++){
			double elev = RS.elevation[gsi];

			//Get neighbours and set weights
			std::vector<double> ndistances;
			std::vector<size_t> neighbours = RS.findneighbours(gsi, ndistances);
			std::vector<double> nweights(neighbours.size(), 0.0);
			for (size_t ni = 0; ni < neighbours.size(); ni++){
				//In case there are coincident samples make minimum distance 10 m
				if (ndistances[ni] < 10.0)ndistances[ni] = 10.0;

				nweights[ni] = std::pow(ndistances[ni], -1.0*InversionOp.InverseDistancePower);
			}
			nweights /= sum(nweights);

			//loop over this sample's layers
			for (size_t li = 0; li < nlayers; li++){

				PetscInt gri = gpindex_c(gsi, li);

				//Loop over each neighbour/layer
				double wsum = 0.0;
				std::vector<std::vector<double>> values(neighbours.size());
				for (size_t ni = 0; ni < neighbours.size(); ni++){
					size_t ngsi = neighbours[ni];
					double nelev = RS.elevation[ngsi];
					std::vector<double> nd = d + (elev - nelev);
					std::vector<double> fo = fractionaloverlaps(d[li], d[li + 1], nd);
					double s = sum(fo);
					if (s != 0.0) fo /= s;
					values[ni] = nweights[ni] * fo;
					wsum += sum(values[ni]);
				}

				//Set values for current sample/layer
				PetscInt gci = gpindex_c(gsi, li);
				PetscErrorCode ierr = MatSetValue(H.mat(), gri, gci, -wsum, INSERT_VALUES); CHKERR(ierr);
				//Set values for neighbours
				for (size_t ni = 0; ni < neighbours.size(); ni++){
					size_t ngsi = neighbours[ni];
					for (size_t nli = 0; nli < nlayers; nli++){
						if (values[ni][nli] == 0.0) continue;
						gci = gpindex_c(ngsi, nli);
						PetscErrorCode ierr = MatSetValue(H.mat(), gri, gci, values[ni][nli], INSERT_VALUES); CHKERR(ierr);
					}
				}
			}//loop over this samples layers
			count++;
		}//loop over samples
		H.assemble();
		H *= std::sqrt(InversionOp.AlphaH / (double)H.nglobalrows());
		rootmessage(mylogfile, "Finished assembling matrix H\n");
		return;
	}

	void B_create_elevation_new(){

		double clogerr = 0.5*(std::log10(100.0 + ConductivityLogPercentError) - std::log10(100.0 - ConductivityLogPercentError));
		rootmessage(mylogfile, "Creating matrix B\n");
		std::vector<double> lthick = get_thicknesses_ref(0);
		std::vector<double> ldepth = get_interface_depths(lthick);
		ldepth.push_back(ldepth.back() + lthick.back());
		PetscInt nconstraints = 0;
		for (size_t bi = 0; bi < ConductivityLogs.size(); bi++){
			cConductivityLog& clog = ConductivityLogs[bi];

			std::vector<double> ndistances;
			std::vector<size_t> neighbours = RS.findneighbourstopoint(clog.x, clog.y, ndistances, ConductivityLogMaximumDistance);
			for (size_t ni = 0; ni < neighbours.size(); ni++){
				size_t ngsi = neighbours[ni];
				double nelev = RS.elevation[ngsi];
				std::vector<double> nd = ldepth + (clog.z - nelev);
				for (size_t li = 0; li < nlayers; li++){
					if (clog.interval_has_overlap(nd[li], nd[li + 1])){
						nconstraints += 1;
					}
				}
			}
		}
		rootmessage(mylogfile, "nconductivitylogconstraints=%d\n", nconstraints);
		B.create_sparse("B", mpicomm, PETSC_DECIDE, nlocalparam, nconstraints, nparam);

		PetscInt nlocalconstraints = B.nlocalrows();
		PetscInt nglobalconstraints = B.nglobalrows();

		clogref.create("clogref", mpicomm, nlocalconstraints, nglobalconstraints);
		clogstd.create("clogstd", mpicomm, nlocalconstraints, nglobalconstraints);

		//Pre-allocate non zeros
		std::vector<PetscInt> d_nnz(B.nlocalrows(), 0);
		std::vector<PetscInt> o_nnz(B.nlocalrows(), 0);
		PetscInt gri = 0;
		for (size_t bi = 0; bi < ConductivityLogs.size(); bi++){
			cConductivityLog& clog = ConductivityLogs[bi];

			//Get neighbours
			std::vector<double> ndistances;
			std::vector<size_t> neighbours = RS.findneighbourstopoint(clog.x, clog.y, ndistances, ConductivityLogMaximumDistance);
			for (size_t ni = 0; ni < neighbours.size(); ni++){
				size_t ngsi = neighbours[ni];
				double nelev = RS.elevation[ngsi];
				std::vector<double> nd = ldepth + (clog.z - nelev);
				for (size_t li = 0; li < nlayers; li++){
					if (clog.interval_has_overlap(nd[li], nd[li + 1])){
						if (B.ownsrow(gri)){
							PetscInt lri = B.lri(gri);
							PetscInt gci = gpindex_c(ngsi, li);
							if (B.inownerdiagonalblock(gri, gci)) d_nnz[lri]++;
							else o_nnz[lri]++;
						}
						gri++;
					}
				}
			}
		}
		B.preallocate(d_nnz, o_nnz);

		//Set values
		double* clogref_loc = clogref.getlocalarray();
		double* clogstd_loc = clogstd.getlocalarray();		
		gri = 0;
		for (size_t bi = 0; bi < ConductivityLogs.size(); bi++){
			cConductivityLog& clog = ConductivityLogs[bi];

			std::vector<double> ndistances;
			std::vector<size_t> neighbours = RS.findneighbourstopoint(clog.x, clog.y, ndistances, ConductivityLogMaximumDistance);
			
			std::vector<double> nweights(neighbours.size());
			for (size_t ni = 0; ni < neighbours.size(); ni++){
				size_t ngsi = neighbours[ni];
				double nelev = RS.elevation[ngsi];
				std::vector<double> nd = ldepth + (clog.z - nelev);
				
				//In case there are coincident positions make minimum distance 10 m
				if (ndistances[ni] < 10.0) ndistances[ni] = 10.0;
				nweights[ni] = std::pow(ndistances[ni], -1.0*InversionOp.InverseDistancePower);

				for (size_t li = 0; li < nlayers; li++){
					PetscInt gci = gpindex_c(ngsi, li);					
					if (clog.interval_has_overlap(nd[li], nd[li + 1])){
						if (B.ownsrow(gri)){
							PetscInt lri = B.lri(gri);														
							
							//Set rhs vector's value
							size_t npoints;
							double linear_mean, log10_mean;
							clog.interval_means(nd[li], nd[li + 1], npoints, linear_mean, log10_mean);
							clogref_loc[lri] = std::log10(log10_mean);
							clogstd_loc[lri] = clogerr;//Converted to log10 decades
							
							//Set matrix entry
							B.set(gri, gci, 1.0);
						}//is owner
						gri++;
					}//has overlap					
				}//layer loop
			}//neighbour loop
		}//bore loop
		clogref.restorelocalarray(clogref_loc);
		clogstd.restorelocalarray(clogstd_loc);
		B.assemble();
		Wb.create_diagonal_to_power("Wb", clogstd, -2.0);
		Wb *= (InversionOp.AlphaB / Wb.nglobalrows());
		rootmessage(mylogfile, "Finished creating matrix B\n");
	}

	void B_create_elevation(){

		double clogerr = 0.5*(std::log10(100.0 + ConductivityLogPercentError) - std::log10(100.0 - ConductivityLogPercentError));
		rootmessage(mylogfile, "Creating matrix B\n");
		std::vector<double> t = get_thicknesses_ref(0);
		std::vector<double> d = get_interface_depths(t);
		d.push_back(d.back() + t.back());
		PetscInt nconstraints = 0;
		for (size_t bi = 0; bi < ConductivityLogs.size(); bi++){
			cConductivityLog& clog = ConductivityLogs[bi];
			std::vector<double> ndistances;
			std::vector<size_t> neighbours = RS.findneighbourstopoint(clog.x, clog.y, ndistances, ConductivityLogMaximumDistance);
			if (neighbours.size() == 0)continue;
			for (size_t li = 0; li < nlayers; li++){
				bool hasoverlap = clog.interval_has_overlap(d[li], d[li + 1]);
				if (hasoverlap) nconstraints += 1;
			}
		}
		rootmessage(mylogfile, "nconductivitylogconstraints=%d\n", nconstraints);
		B.create_sparse("B", mpicomm, PETSC_DECIDE, nlocalparam, nconstraints, nparam);

		PetscInt nlocalconstraints = B.nlocalrows();
		PetscInt nglobalconstraints = B.nglobalrows();

		//Pre-allocate non zeros
		std::vector<PetscInt> d_nnz(B.nlocalrows(), 0);
		std::vector<PetscInt> o_nnz(B.nlocalrows(), 0);
		PetscInt gri = 0;
		for (size_t bi = 0; bi < ConductivityLogs.size(); bi++){
			cConductivityLog& clog = ConductivityLogs[bi];

			//Get neighbours
			std::vector<double> ndistances;
			std::vector<size_t> neighbours = RS.findneighbourstopoint(clog.x, clog.y, ndistances, ConductivityLogMaximumDistance);
			if (neighbours.size() == 0)continue;

			for (size_t li = 0; li < nlayers; li++){
				bool hasoverlap = clog.interval_has_overlap(d[li], d[li + 1]);
				if (hasoverlap){
					if (B.ownsrow(gri)){
						PetscInt lri = B.lri(gri);
						for (size_t ni = 0; ni < neighbours.size(); ni++){
							size_t ngsi = neighbours[ni];
							double nelev = RS.elevation[ngsi];
							std::vector<double> nd = d + (clog.z - nelev);
							std::vector<double> fo = fractionaloverlaps(d[li], d[li + 1], nd);
							for (size_t nli = 0; nli < fo.size(); nli++){
								if (fo[nli]>0){
									PetscInt gci = gpindex_c(ngsi, nli);
									if (B.inownerdiagonalblock(gri, gci)) d_nnz[lri]++;
									else o_nnz[lri]++;
								}
							}
						}
					}
					gri++;
				}
			}
		}
		B.preallocate(d_nnz, o_nnz);

		//Set values
		clogref.create("clogref", mpicomm, nlocalconstraints, nglobalconstraints);
		clogstd.create("clogstd", mpicomm, nlocalconstraints, nglobalconstraints);
		double* clogref_loc = clogref.getlocalarray();
		double* clogstd_loc = clogstd.getlocalarray();
		gri = -1;
		for (size_t bi = 0; bi < ConductivityLogs.size(); bi++){
			cConductivityLog& clog = ConductivityLogs[bi];
			std::vector<double> ndistances;
			std::vector<size_t> neighbours = RS.findneighbourstopoint(clog.x, clog.y, ndistances, ConductivityLogMaximumDistance);
			if (neighbours.size() == 0)continue;

			//Set neighbour weights
			std::vector<double> nweights(neighbours.size());
			for (size_t ni = 0; ni < neighbours.size(); ni++){
				//In case there are coincident positions make minimum distance 10 m
				if (ndistances[ni] < 10.0) ndistances[ni] = 10.0;

				nweights[ni] = std::pow(ndistances[ni], -1.0*InversionOp.InverseDistancePower);
			}
			nweights /= sum(nweights);

			for (size_t li = 0; li < nlayers; li++){
				size_t npoints;
				double linear_mean, log10_mean;
				bool hasoverlap = clog.interval_means(d[li], d[li + 1], npoints, linear_mean, log10_mean);
				if (hasoverlap == false) continue;
				gri++;
				if (B.ownsrow(gri) == false) continue;

				PetscInt lri = B.lri(gri);
				//Loop over each neighbour/layer
				double wsum = 0.0;
				std::vector<std::vector<double>> values(neighbours.size());
				for (size_t ni = 0; ni < neighbours.size(); ni++){
					size_t ngsi = neighbours[ni];
					double nelev = RS.elevation[ngsi];
					std::vector<double> nd = d + (clog.z - nelev);
					std::vector<double> fo = fractionaloverlaps(d[li], d[li + 1], nd);
					double s = sum(fo);
					if (s != 0.0) fo /= s;
					values[ni] = nweights[ni] * fo;
					wsum += sum(values[ni]);
				}

				//if wsum is tiny there is very little overlap with any neighbour
				if (wsum < 0.001){
					clogref_loc[lri] = 0.0;
					clogstd_loc[lri] = clogerr;//Converted to log10 decades				
					continue;
				}

				//Set rhs vector's value
				clogref_loc[lri] = std::log10(log10_mean);
				clogstd_loc[lri] = clogerr;//Converted to log10 decades				

				//Set values for neighbours
				for (size_t ni = 0; ni < neighbours.size(); ni++){
					size_t ngsi = neighbours[ni];
					for (size_t nli = 0; nli < nlayers; nli++){
						if (values[ni][nli] == 0.0) continue;
						values[ni][nli] /= wsum;
						PetscInt gci = gpindex_c(ngsi, nli);
						B.set(gri, gci, values[ni][nli]);
					}
				}
			}//loop over bore's layers'
		}//loop over bores	
		clogref.restorelocalarray(clogref_loc);
		clogstd.restorelocalarray(clogstd_loc);
		B.assemble();
		Wb.create_diagonal_to_power("Wb", clogstd, -2.0);
		if (Wb.nglobalrows() > 0){
			Wb *= (InversionOp.AlphaB / Wb.nglobalrows());
		}
		rootmessage(mylogfile, "Finished creating matrix B\n");
	}

	void report_matrix_memory_usage(){
		rootmessage(mylogfile, "J matrix global memory %.3lf MiB\n", J.globalmemory() / 1e6);
		rootmessage(mylogfile, "Wd matrix global memory %.3lf MiB\n", Wd.globalmemory() / 1e6);
		rootmessage(mylogfile, "V matrix global memory %.3lf MiB\n", V.globalmemory() / 1e6);
		rootmessage(mylogfile, "H matrix global memory %.3lf MiB\n", H.globalmemory() / 1e6);
		rootmessage(mylogfile, "B matrix global memory %.3lf MiB\n", B.globalmemory() / 1e6);
		rootmessage(mylogfile, "Wb matrix global memory %.3lf MiB\n", Wb.globalmemory() / 1e6);
		rootmessage(mylogfile, "Wr matrix global memory %.3lf MiB\n", Wr.globalmemory() / 1e6);
		rootmessage(mylogfile, "P matrix global memory %.3lf MiB\n", P.globalmemory() / 1e6);

		double total = J.globalmemory() + Wd.globalmemory() + V.globalmemory() + H.globalmemory() + B.globalmemory() + Wb.globalmemory() + Wr.globalmemory() + P.globalmemory();
		rootmessage(mylogfile, "Total matrix global memory %.3lf MiB\n", total / 1e6);
	}

	std::vector<double> get_conductivity_model(const size_t& localsampleindex, const double* mlocal){
		std::vector<double> c(nlayers);
		size_t pi = nparampersample * localsampleindex;
		for (size_t li = 0; li < nlayers; li++){
			c[li] = std::pow(10.0, mlocal[pi+li]);
		}
		return c;
	}

	std::vector<double> get_thicknesses_ref(const size_t localsampleindex){
		std::vector<double> thickness(nlayers - 1);
		for (size_t i = 0; i < nlayers - 1; i++){
			thickness[i] = E.tref(localsampleindex, i);
		}
		return thickness;
	}

	std::vector<double> get_interface_depths(const std::vector<double>& thickness){
		std::vector<double> d(thickness.size() + 1);
		d[0] = 0;
		for (size_t i = 0; i < thickness.size(); i++){
			d[i + 1] = d[i] + thickness[i];
		}
		return d;
	}

	cTDEmGeometry get_geometry_ref(const size_t localsampleindex){
		cTDEmGeometry geom;
		for (size_t gi = 0; gi < G.size(); gi++){
			geom[gi] = G[gi].ref(localsampleindex);
		}
		return geom;
	}

	cTDEmGeometry get_geometry_model(const size_t& localsampleindex, const double* mlocal)
	{
		cTDEmGeometry geom = get_geometry_ref(localsampleindex);						
		size_t pi = nparampersample * localsampleindex;
		for (size_t gi = 0; gi < UGI.size(); gi++){			
			geom[UGI[gi]] = mlocal[pi + nlayers + gi];			
		}
		return geom;	
	}

	std::vector<size_t> unknown_geometry_indices(){
		std::vector<size_t> indices;
		for (size_t gi = 0; gi < G.size(); gi++){
			if (G[gi].solve) indices.push_back(gi);
		}
		return indices;
	}

	bool forwardmodel_and_jacobian(const cPetscDistVector& m, cPetscDistVector& g, const bool computejacobian){
		
		static double natlog10 = std::log(10.0);		
		std::vector<double> predicted;
		std::vector<std::vector<double>> derivatives;

		const double* mlocal = m.getlocalreadonlyarray();
		//cOwnership  mdist = m.ownership();

		double* glocal = g.getlocalarray();
		cOwnership gdist = g.ownership();
		//rootmessage(mylogfile, "Starting forward modelling\n");
		for (size_t si = (size_t)sown.start; si < (size_t)sown.end; si++){
			size_t lsi = sown.localind(si);
			//size_t gpi = gpindex_c(si, 0);
			//size_t lpi = mdist.localind(gpi);
			
			std::vector<double> conductivity = get_conductivity_model(lsi,mlocal);
			std::vector<double> thickness    = get_thicknesses_ref(lsi);
			cTDEmGeometry       geometry     = get_geometry_model(lsi,mlocal);
			
			size_t gdi = dindex(si, 0);
			size_t ldi = gdist.localind(gdi);
			for (size_t ti = 0; ti < T.size(); ti++){

				T[ti].forward_model_and_derivatives(conductivity, thickness, geometry, predicted, derivatives, computejacobian, UGI);

				for (size_t k = 0; k < predicted.size(); k++){
					glocal[ldi + k] = predicted[k];
				}

				if (computejacobian){
					for (size_t k = 0; k < derivatives.size(); k++){
						for (size_t li = 0; li < nlayers; li++){
							//multiply by natural log(10) as parameters are in log base 10 units
							const double val = natlog10 * conductivity[li] * derivatives[k][li];
							J.set(gdi + k, gpindex_c(si,li), val);
						}
						for (size_t gi = 0; gi < UGI.size(); gi++){
							const double val = derivatives[k][nlayers+gi];
							J.set(gdi + k, gpindex_g(si,gi), val);
						}
					}
				}
				gdi += T[ti].ndata();
				ldi += T[ti].ndata();
			}
		}
		m.restorelocalreadonlyarray(mlocal);
		g.restorelocalarray(glocal);
		if (computejacobian) J.assemble();
		
		//rootmessage(mylogfile, "Finished forward modelling\n");
		return true;
	}

	bool setup(){

		count_samples();
		rootmessage(mylogfile, "nsamples=%lu\n", nsamples);
		rootmessage(mylogfile, "nchannels=%lu\n", nchan);
		rootmessage(mylogfile, "ndata=%lu\n", ndata);
		rootmessage(mylogfile, "nlayers=%lu\n", nlayers);
		rootmessage(mylogfile, "nparampersample=%lu\n", nparampersample);
		rootmessage(mylogfile, "nparam=%lu\n", nparam);

		read_data();
		read_conductivity_logs();

		dobs.create("dobs", mpicomm, nlocaldata, ndata);
		dobs.set_local(local_data());
		//dobs.writetextfile("dobs.vec");

		dstd.create("dstd", mpicomm, nlocaldata, ndata);
		dstd.set_local(local_noise());
		//dstd.writetextfile("dstd.vec");

		mref.create("mref", mpicomm, nlocalparam, nparam);
		mref.set_local(local_mref());
		//mref.writetextfile("mref.vec");

		mstd.create("mstd", mpicomm, nlocalparam, nparam);
		mstd.set_local(local_mstd());
		//mstd.writetextfile("mstd.vec");

		Wd_create();
		//Wd.writetextfile("Wd.smat");

		Wr_create();
		//Wr.writetextfile("Wr.smat");

		B_create_elevation();
		//B.writetextfile("B.smat");
		//Wb.writetextfile("Wb.smat");
		//clogref.writetextfile("clogref.vec");
		//clogstd.writetextfile("clogstd.vec");

		if (InversionOp.VerticalSmoothnessMethod == SM_1ST_DERIVATIVE){
			V_create_1st_derivative();
		}
		else if (InversionOp.VerticalSmoothnessMethod == SM_2ND_DERIVATIVE){
			V_create_2nd_derivative();
		}
		//V.writetextfile("V.smat");

		H_create_elevation();
		//H.writetextfile("H.smat");
		//H.printinfostring();
		//mref.writetextfile("mref.vec");


		J_create();
		J_set_nonzero_pattern();
		//J.writetextfile("J.smat");

		rootmessage(mylogfile, "Creating preconditioner\n");
		P.create_identity("P", mpicomm, nlocalparam, nparam);
		//P.writetextfile("P.smat");		
		rootmessage(mylogfile, "Finished creating preconditioner\n");

		report_matrix_memory_usage();
		
		return true;
	};

	static void shellmatrixmult(Mat shellmat, Vec x, Vec b)
	{
		void* context;
		PetscErrorCode ierr = MatShellGetContext(shellmat, &context); CHKERR(ierr);
		cAllAtOnceInverter* pA = (cAllAtOnceInverter*)context;
		pA->rawvecmult(x, b);
		return;
	}

	PetscErrorCode rawvecmult(Vec& xVec, Vec& bVec)
	{
		// Ax=b
		cPetscDistVector x(xVec);
		cPetscDistVector b(bVec);
		b = (J ^ (Wd*(J*x)));
		if (InversionOp.AlphaV > 0.0)b += ((V ^ (V*x)) *= lambda);
		if (InversionOp.AlphaH > 0.0)b += ((H ^ (H*x)) *= lambda);
		if (InversionOp.AlphaB > 0.0)b += ((B ^ (Wb*(B*x))) *= lambda);
		if (InversionOp.AlphaR > 0.0)b += ((Wr*x) *= lambda);
		return 0;
	}

	cPetscDistVector cg_solve(cPetscDistShellMatrix& A, const cPetscDistVector& m, const cPetscDistVector& g){

		cPetscDistVector b("b", mpicomm, nlocalparam, nparam);
		rootmessage(mylogfile, "Starting CG solve\n");
		cStopWatch sw;

		b = J ^ (Wd*(dobs - g + J*m));
		if (InversionOp.AlphaB > 0) b += lambda*((B ^ (Wb^clogref)));
		if (InversionOp.AlphaR > 0) b += lambda*((Wr^mref));

		cPetscDistVector mtrial = A.solve_CG(P, b, m);
		rootmessage(mylogfile, "Finished CG solve time=%lf\n", sw.etimenow());
		rootmessage(mylogfile, A.convergence_summary().c_str());
		return (mtrial - m);
	};

	PetscErrorCode rawvecmult_dm(Vec& xVec, Vec& bVec)
	{
		// Ax=b
		cPetscDistVector x(xVec);
		cPetscDistVector b(bVec);
		b = (J ^ (Wd*(J*x)));
		if (InversionOp.AlphaV > 0.0)b += ((V ^ (V*x)) *= lambda);
		if (InversionOp.AlphaH > 0.0)b += ((H ^ (H*x)) *= lambda);
		if (InversionOp.AlphaB > 0.0)b += ((B ^ (Wb*(B*x))) *= lambda);
		if (InversionOp.AlphaR > 0.0)b += ((Wr*x) *= lambda);
		return 0;
	}
	
	cPetscDistVector cg_solve_dm(cPetscDistShellMatrix& A, const cPetscDistVector& m, const cPetscDistVector& g){

		cPetscDistVector b("b", mpicomm, nlocalparam, nparam);
		rootmessage(mylogfile, "Starting CG solve\n");
		cStopWatch sw;

		b = J ^ (Wd*(dobs - g));
		if (InversionOp.AlphaV > 0) b -= lambda*((V ^ (V*m)));
		if (InversionOp.AlphaH > 0) b -= lambda*((H ^ (H*m)));
		if (InversionOp.AlphaB > 0) b -= lambda*((B ^ (B*m)));
		if (InversionOp.AlphaB > 0) b += lambda*((B ^ (Wb^clogref)));
		if (InversionOp.AlphaR > 0) b += lambda*((Wr^ (mref-m)));

		cPetscDistVector mtrial = A.solve_CG(P, b, m);
		rootmessage(mylogfile, "Finished CG solve time=%lf\n", sw.etimenow());
		rootmessage(mylogfile, A.convergence_summary().c_str());
		return (mtrial - m);
	};

	double PhiD(const cPetscDistVector& g){
		return Wd.vtAv(dobs - g);
	}

	double PhiV(const cPetscDistVector& m){
		return (V*m).l2_norm();
	}

	double PhiH(const cPetscDistVector& m){
		return (H*m).l2_norm();
	}

	double PhiB(const cPetscDistVector& m){
		return Wb.vtAv(B*m - clogref);
	}

	double PhiR(const cPetscDistVector& m){
		return Wr.vtAv(m - mref);
	}
	
	void find_stepfactor(const cPetscDistVector& m, const cPetscDistVector& dm, const double& currentphid, const double& targetphid, double& bestsf, double& bestphid, double& improvement){

		rootmessage(mylogfile, "Finding step factor\n");
		cStopWatch sw;		
		cPetscDistVector gtrial = dobs;
		cPetscDistVector mtrial = m;
		cInversionLineSearcher LS(currentphid, targetphid);

		double sf;
		while (LS.next(sf)){
			mtrial = m + sf*dm;
			forwardmodel_and_jacobian(mtrial, gtrial, false);
			double phid = PhiD(gtrial);
			LS.addtrial(sf, phid);
		}
		LS.nearestindex(bestsf, bestphid);
		rootmessage(mylogfile, "Find stepfactor time=%lf\n", sw.etimenow());
		improvement = 100.0*(currentphid - bestphid) / currentphid;
		rootmessage(mylogfile, "Step factor = %.5lf\n", bestsf);
		rootmessage(mylogfile, "Found PhiD  = %.5lf\n", bestphid);
		rootmessage(mylogfile, "Improvement = %.5lf%%\n", improvement);
		
		//if (mpirank == 0){
		//	std::string stepsfile = strprint("output//steps//steps_%02llu.txt", mLastIteration);
		//	LS.writetextfile(stepsfile);
		//}

		return;
	}

	void iterate(){
		
		cPetscDistVector m = mref;
		m.setname("m");

		cPetscDistShellMatrix A("A", mpicomm, nlocalparam, nlocalparam, nparam, nparam, (void*)this);
		A.set_multiply_function_vec((void*)shellmatrixmult);

		lambda = 1.0;

		bool keepgoing = true;
		size_t iteration = 1;
		while (keepgoing){			
			rootmessage(mylogfile, "\n\nItaration = %lu\n", iteration);
			cStopWatch sw;
			cPetscDistVector g("g", mpicomm, nlocaldata, ndata);
			forwardmodel_and_jacobian(m, g, true);
			rootmessage(mylogfile, "Forward modelling time=%lf\n", sw.etimenow());

			double phiv = PhiV(m); double phih = PhiH(m);
			double phib = PhiB(m); double phir = PhiR(m);
			double phid = PhiD(g);
			double phi = phid + phiv + phih + phib + phir;
			double targetphid = phid * 0.7;
			if (targetphid < InversionOp.MinimumPhiD) targetphid = InversionOp.MinimumPhiD;
			
			log_iteration_msg(lambda, phi, phiv, phih, phib, phir, phid, targetphid);
						
			cPetscDistVector dm = cg_solve(A, m, g);									
			double bestsf, bestphid, improvement;
			find_stepfactor(m, dm, phid, targetphid, bestsf, bestphid, improvement);
			
			iteration++;
			if (improvement > 0.0){
				m += bestsf*dm;
				forwardmodel_and_jacobian(m, g, false);

				mLastPhiD = bestphid;
				mLastLambda = lambda;
				mLastIteration = iteration;
				write_results(OutputOp.DataFile, m, g);				
			}

			if (improvement < InversionOp.MinimumPercentageImprovement){
				keepgoing = false;
				//lambda = lambda * 0.7;
			}
			if (iteration > InversionOp.MaximumIterations){
				keepgoing = false;
			}
			if (bestphid < InversionOp.MinimumPhiD){
				keepgoing = false;
			}
		}
	};

	void log_iteration_msg(const double& lam, const double& phi, const double& phiv, const double& phih, const double& phib,	const double& phir, const double& phid, const double& targetphid){
		rootmessage(mylogfile, "Current Lambda = %lf\n", lam);
		rootmessage(mylogfile, "Current Phi  = %lf\n", phi);
		rootmessage(mylogfile, "Current PhiV = %lf\n", phiv);
		rootmessage(mylogfile, "Current PhiH = %lf\n", phih);
		rootmessage(mylogfile, "Current PhiB = %lf\n", phib);
		rootmessage(mylogfile, "Current PhiR = %lf\n", phir);
		rootmessage(mylogfile, "Current PhiD = %lf\n", phid);
		rootmessage(mylogfile, "Target  PhiD = %lf\n", targetphid);
	}
	
	void write_results(const std::string& filename, const cPetscDistVector& m, const cPetscDistVector& g){

		for (int p = 0; p < mpisize; p++){
			if (p == mpirank){
				if (mpirank == 0){
					FILE* fp = fileopen(filename, "w");
					fclose(fp);
				}
				append_my_results(filename, m, g);
			}
			m.mpibarrier();
		}

	}

	void append_my_results(const std::string& filename, const cPetscDistVector& m, const cPetscDistVector& g){
		
		const double* mlocal = m.getlocalreadonlyarray();
		const double* glocal = g.getlocalreadonlyarray();
		//cOwnership mdist = m.ownership();
		cOwnership gdist = g.ownership();
		std::vector<double> samplephid = get_sample_phid(g);

		cOutputFileInfo OI;
		std::string buf;
		FILE* fp = fileopen(filename, "a");
		for (size_t lsi = 0; lsi < (size_t)sown.nlocal(); lsi++){
			size_t gsi = sown.globalind(lsi);
			//size_t lpi = mdist.localind(gpindex_c(gsi, 0));
			std::vector<double> conductivity = get_conductivity_model(lsi, mlocal);
			std::vector<double> thickness    = get_thicknesses_ref(lsi);
			thickness.push_back(thickness.back());
			cTDEmGeometry gref = get_geometry_ref(lsi);
			cTDEmGeometry ginv = get_geometry_model(lsi, mlocal);


			OI.addfield("survey", 'I', 12, 0);
			OI.setcomment("Survey number");
			buf += strprint("%12d", (int)fdsurvey(lsi));

			OI.addfield("date", 'I', 12, 0);
			OI.setcomment("Date number");
			buf += strprint("%12d", (int)fddate(lsi));

			OI.addfield("flight", 'I', 12, 0);
			OI.setcomment("Flight number, IntrepidFlightNumber");
			buf += strprint("%12d", (int)fdflight(lsi));

			OI.addfield("line", 'I', 12, 0);
			OI.setcomment("Line number, IntrepidLineNumber");
			buf += strprint("%12d", (int)fdline(lsi));

			OI.addfield("fiducial", 'F', 12, 2);
			OI.setcomment("Fiducial number, IntrepidFiducial");
			buf += strprint("%12.2lf", fdfiducial(lsi));

			//Location
			OI.addfield("easting", 'F', 10, 1);
			OI.setunits("m"); OI.setcomment("IntrepidX");
			buf += strprint("%10.1lf", fdx(lsi));

			OI.addfield("northing", 'F', 10, 1);
			OI.setunits("m"); OI.setcomment("IntrepidY");
			buf += strprint("%10.1lf", fdy(lsi));

			OI.addfield("elevation", 'F', 10, 2);
			OI.setunits("m"); OI.setcomment("Ground elevation relative to sea-level");
			buf += strprint("%10.2lf", fdelevation(lsi));
					
			for (size_t gi = 0; gi < G.size(); gi++){
				OI.addfield(gref.fname(gi), 'F', 9, 2);
				OI.setunits(gref.units(gi));
				OI.setcomment(gref.description(gi));
				buf += strprint("%9.2lf", gref[gi]);
			}
			
			for (size_t gi = 0; gi < UGI.size(); gi++){
				OI.addfield("inverted_"+ginv.fname(UGI[gi]), 'F', 9, 2);
				OI.setunits(ginv.units(UGI[gi]));
				OI.setcomment("Inverted " + ginv.description(UGI[gi]));
				buf += strprint("%9.2lf", ginv[UGI[gi]]);
			}

			OI.addfield("nlayers", 'I', 4, 0);
			OI.setcomment("Number of layers");
			buf += strprint("%4lu", nlayers);

			OI.addfield("conductivity", 'E', 15, 6, nlayers);
			OI.setunits("S/m"); OI.setcomment("Layer conductivity");
			for (size_t li = 0; li < nlayers; li++){
				buf += strprint("%15.6le", conductivity[li]);
			}

			OI.addfield("thickness", 'F', 9, 2, nlayers);
			OI.setunits("m"); OI.setcomment("Layer thickness");
			for (size_t li = 0; li < nlayers; li++){
				buf += strprint("%9.2lf", thickness[li]);
			}

			if (OutputOp.PositiveLayerBottomDepths){
				OI.addfield("depth_bottom", 'F', 9, 2, nlayers);
				OI.setunits("m"); OI.setcomment("Depth to bottom of layer");
				double tsum = 0.0;
				for (size_t i = 0; i < nlayers; i++){
					buf += strprint("%9.2lf", tsum);
					tsum += thickness[i];
				}
			}

			if (OutputOp.NegativeLayerBottomDepths){
				OI.addfield("depth_bottom_negative", 'F', 9, 2, nlayers);
				OI.setunits("m"); OI.setcomment("Negative of depth to bottom of layer");
				double tsum = 0.0;
				for (size_t i = 0; i < nlayers; i++){
					tsum += thickness[i];
					buf += strprint("%9.2lf", -tsum);
				}
			}

			if (OutputOp.InterfaceElevations){
				OI.addfield("elevation_interfaces", 'F', 9, 2, nlayers);
				OI.setunits("m"); OI.setcomment("Elevation of interfaces");
				double etop = fdelevation(lsi);
				for (size_t i = 0; i < nlayers; i++){
					buf += strprint("%9.2lf", etop);
					etop -= thickness[i];
				}
			}

			char cid[3] = { 'X', 'Y', 'Z' };
			if (OutputOp.ObservedData){								
				for (size_t si = 0; si < T.size(); si++){
					cSystemInfo& S = T[si];					
					std::string sys = strprint("EMSystem_%lu_", si + 1);
					for (size_t ci = 0; ci < S.Comp.size(); ci++){
						cComponentInfo& C = S.Comp[ci];
						if (C.Use == false)continue;
						if (S.InvertTotalField){
							OI.addfield("observed_" + sys + cid[ci] + "P", 'E', 15, 6);
							OI.setcomment("Observed " + sys + cid[ci] + "-component primary field");
							buf += strprint("%15.6le", C.fdp(lsi, 0));
						}							
						OI.addfield("observed_" + sys + cid[ci] + "S", 'E', 15, 6, C.nw);
						OI.setcomment("Observed " + sys + cid[ci] + "-component secondary field windows");
						for (size_t w = 0; w < C.nw; w++){
							buf += strprint("%15.6le", C.fds(lsi,w));
						}						
					}
				}
			}

			if (OutputOp.Noise){				
				for (size_t si = 0; si < T.size(); si++){
					cSystemInfo& S = T[si];
					std::string sys = strprint("EMSystem_%lu_", si + 1);
					for (size_t ci = 0; ci < S.Comp.size(); ci++){
						cComponentInfo& C = S.Comp[ci];
						if (C.Use == false)continue;						
						OI.addfield("noise_" + sys + cid[ci] + "S", 'E', 15, 6, C.nw);
						OI.setcomment("Estimated noise " + sys + cid[ci] + "-component secondary field windows");
						for (size_t w = 0; w < C.nw; w++){
							buf += strprint("%15.6le", C.fdn(lsi, w));
						}
					}
				}
			}

			if (OutputOp.PredictedData){				
				size_t ldi = gdist.localind(dindex(gsi, 0));
				for (size_t si = 0; si < T.size(); si++){
					cSystemInfo& S = T[si];
					S.forward_model(conductivity, thickness, ginv);
					std::string sys = strprint("EMSystem_%lu_", si + 1);
					for (size_t ci = 0; ci < S.Comp.size(); ci++){
						cComponentInfo& C = S.Comp[ci];
						if (C.Use == false)continue;
						if (S.InvertTotalField){
							OI.addfield("predicted_" +  sys + cid[ci] + "P", 'E', 15, 6);
							OI.setcomment("Predicted " + sys + cid[ci] + "-component primary field");
							if (ci == 0) buf += strprint("%15.6le", S.T.PrimaryX);
							else if (ci == 1) buf += strprint("%15.6le", S.T.PrimaryY);
							else              buf += strprint("%15.6le", S.T.PrimaryZ);

							OI.addfield("predicted_" + sys + cid[ci] + "S", 'E', 15, 6, C.nw);
							OI.setcomment("Predicted " + sys + cid[ci] + "-component secondary field windows");
							for (size_t w = 0; w < C.nw; w++){
								if      (ci == 0) buf += strprint("%15.6le", S.T.X[w]);
								else if (ci == 1) buf += strprint("%15.6le", S.T.Y[w]);
								else              buf += strprint("%15.6le", S.T.Z[w]);
							}
						}
						else{
							OI.addfield("predicted_" + sys + cid[ci] + "S", 'E', 15, 6, C.nw);
							OI.setcomment("Predicted " + sys + cid[ci] + "-component secondary field windows");
							for (size_t w = 0; w < C.nw; w++){
								buf += strprint("%15.6le", glocal[ldi]);
								ldi++;
							}
						}
					}
				}
			}
			

			//Inversion parameters
			OI.addfield("SamplePhiD", 'E', 15, 6);
			OI.setcomment("Normalised data misfit for this sample");
			buf += strprint("%15.6le", samplephid[lsi]);

			OI.addfield("PhiD", 'E', 15, 6);
			OI.setcomment("Normalised data misfit");
			buf += strprint("%15.6le", mLastPhiD);

			OI.addfield("Lambda", 'E', 15, 6);
			OI.setcomment("Lambda regularization parameter");
			buf += strprint("%15.6le", mLastLambda);

			OI.addfield("Iterations", 'I', 4, 0);
			OI.setcomment("Number of iterations");
			buf += strprint("%4lu", mLastIteration);

			//Carriage return
			buf += strprint("\n");
			if ((int)lsi == sown.nlocal() - 1 || buf.size() >= 2048){
				fprintf(fp, buf.c_str());
				fflush(fp);
				buf.resize(0);
			}

			OI.lockfields();
			if (lsi == 0 && mpirank == 0){
				sFilePathParts fpp = getfilepathparts(filename);
				std::string hdrfile = fpp.directory + fpp.prefix + ".hdr";
				OI.write_simple_header(hdrfile);

				std::string aseggdffile = fpp.directory + fpp.prefix + ".dfn";
				OI.write_aseggdf_header(aseggdffile);
			}
		}
		m.restorelocalreadonlyarray(mlocal);
		g.restorelocalreadonlyarray(glocal);
		fclose(fp);
	};

	std::vector<double> get_sample_phid(const cPetscDistVector& g){

		cPetscDistVector nr2 = (g - dobs) / dstd;
		nr2.pow(2.0);

		const double* nr2local = nr2.getlocalreadonlyarray();
		cOwnership    nr2dist = nr2.ownership();
		std::vector<double> samplephid((size_t)sown.nlocal());
		for (size_t lsi = 0; lsi < (size_t)sown.nlocal(); lsi++){
			size_t gsi = sown.globalind(lsi);
			size_t ldi = nr2dist.localind(dindex(gsi, 0));
			double sum = 0.0;
			for (size_t i = 0; i < nchan; i++){
				sum += nr2local[ldi + i];
			}
			samplephid[lsi] = sum / nchan;
		}
		nr2.restorelocalreadonlyarray(nr2local);
		return samplephid;
	}

};

int main(int argc, char** argv)
{
	mylogfile = (FILE*) NULL;
	PetscErrorCode ierr;
	try{
		ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

		if(argc < 2){
			rootmessage("%s\n", commandlinestring(argc, argv).c_str());
			rootmessage("%s\n", versionstring(VERSION, __TIME__, __DATE__).c_str());
			rootmessage("Usage: %s control_file_name\n", argv[0]);
			rootmessage("Too few command line arguments\n");
		}
		else if(argc > 2){
			rootmessage("%s\n", commandlinestring(argc, argv).c_str());
			rootmessage("%s\n", versionstring(VERSION, __TIME__, __DATE__).c_str());
			rootmessage("Usage: %s control_file_name\n", argv[0]);
			rootmessage("Too many command line arguments\n");
		}
		else{
			cAllAtOnceInverter(argc, argv);
		}
		ierr = PetscFinalize(); CHKERRQ(ierr);
	}
	catch (const std::string msg){
		rootmessage(mylogfile, "%s", msg.c_str());
		ierr = PetscFinalize(); CHKERRQ(ierr);
	}
	catch (const std::runtime_error e){
		rootmessage(mylogfile, "%s", e.what());
		ierr = PetscFinalize(); CHKERRQ(ierr);
	}
	catch (const std::exception e){
		rootmessage(mylogfile, "%s", e.what());
		ierr = PetscFinalize(); CHKERRQ(ierr);
	}

	fflush(stdout);
	return 0;
};

