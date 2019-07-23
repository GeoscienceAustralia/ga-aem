/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _rjmcmc1dtdeminverter_H
#define _rjmcmc1dtdeminverter_H

#include "gaaem_version.h"
#include "blocklanguage.h"
#include "fielddefinition.h"
#include "tdemsystem.h"
#include "file_formats.h"
#include "rjmcmc1d.h"

struct sTDEmSystemInfo{
	cTDEmSystem T;	
	size_t ncomps;
	size_t nwindows;
	size_t nchans;
	
	bool useX;
	bool useY;
	bool useZ;
	bool useTotal;
	bool reconstructPrimary;
	bool estimateNoiseFromModel;

	double x_multiplicativenoise;
	double y_multiplicativenoise;
	double z_multiplicativenoise;
	std::vector<double> x_additivenoise;
	std::vector<double> y_additivenoise;
	std::vector<double> z_additivenoise;
	double  oPX,oPY,oPZ;
	std::vector<double> oSX,oSY,oSZ;	
	std::vector<double> oEX,oEY,oEZ;	
	cFieldDefinition fd_oPX,fd_oPY,fd_oPZ;
	cFieldDefinition fd_oSX,fd_oSY,fd_oSZ;
	cFieldDefinition fd_oEX,fd_oEY,fd_oEZ;
};

class rjmcmc1dTDEmInverter : public rjMcMC1DSampler{
	
	public:

	rjmcmc1dTDEmInverter(const std::string& executable, const std::string& controlfile, size_t size, size_t rank);
	virtual ~rjmcmc1dTDEmInverter();

		
	size_t nsystems;
	std::vector<sTDEmSystemInfo> SV;
		
	std::vector<rjMcMCNuisance> ntemplate;
	std::vector<std::string> ninitial;
	cTDEmGeometry  IG;	
	cOutputFileInfo OI;
		
	void initialise(const std::string& executable, const std::string& controlfilename);
	void loadoptions();
	
	void initialise_systems();
	void getcolumnnumbers();
	void initialise_sampler();
	
	bool readnextrecord();
	bool readnextrecord_thisprocess();	
	void parsecurrentrecord();
	void set_data();
	void set_nuisance();
	
	void sample();
	
	std::string results_string();
	std::string prefixstring();
	void writemapstofile();
	void writechainstofile();

	cTDEmGeometry getgeometry(const rjMcMC1DModel& m);
	std::vector<double> collect(const sTDEmSystemInfo& S, const cTDEmSystem& T);	
	std::vector<double> forwardmodel(const rjMcMC1DModel& m);		
	double misfit(const std::vector<double>& g);
	double misfit(const rjMcMC1DModel& m);
		
	bool SaveMaps;
	int  SaveMapsRate;
	bool SaveChains;
	int  SaveChainsRate;
		
	cBlock Control;
	std::string LogFile;	
	double memoryusedatstart;

	std::string InputDataFile;	
	FILE*  fp_indata;		

	size_t Outputcolumn; //output column number
	std::string OutputDirectory;
	std::string OutputDataFile;	
	std::string MapsDirectory;
	std::string ChainsDirectory;
		
	size_t HeaderLines;		
	size_t SubSample;	
	size_t FirstRecord;
	size_t LastRecord;
	size_t CurrentRecord;	
	std::string CurrentRecordString;
	std::vector<std::string> CurrentRecordFields;	
								
	size_t surveynumber;
	size_t datenumber;
	size_t flightnumber;
	size_t linenumber;	
	double fidnumber;
	double timenumber;
	double xord;
	double yord;
	double elevation;
	double altimeter;

	cFieldDefinition fd_surveynumber;
	cFieldDefinition fd_datenumber;	
	cFieldDefinition fd_flightnumber;	
	cFieldDefinition fd_linenumber;
	cFieldDefinition fd_fidnumber;	
	cFieldDefinition fd_xord;
	cFieldDefinition fd_yord;
	cFieldDefinition fd_elevation;	
	cFieldDefinition fd_geometry[10];			
};

#endif
