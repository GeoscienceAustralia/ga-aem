/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstring>

#include "gaaem_version.h"
#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "earth1d.h"
#include "lem.h"
#include "tdemsystem.h"

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

int process(std::string controlfile);
int parseinputrecord(const char* record, cTDEmGeometry& G, cEarth1D& E);
int writeoutputrecord(const bool& csvoutput, FILE* fout, FILE* fhdr, size_t recnum, const cTDEmSystem& T, const cTDEmResponse& R);
int writehdr(FILE* fhdr, const size_t& nw);
int writehdrentry(FILE* fhdr, const char* s, size_t& colnum, size_t nbands);
int writecsvheader(FILE* fout, const size_t& nw);

int main(int argc, char* argv[])
{			
	if (argc < 2){
		printf("Usage: %s control_file_name\n", argv[0]);
		return EXIT_FAILURE;
	}
	else if (argc > 2){
		printf("**Error: Too many command line arguments\n");
		printf("Usage: %s control_file_name\n", argv[0]);
		return EXIT_FAILURE;
	}
	else{
		printf("Program 'gaforwardmodeltdem'\n");
		printf("Geoscience Australia's Airborne Electromagnetic Layered Earth Forward Modelling\n\n");
		printf("Working directory: %s\n", getcurrentdirectory().c_str());
		printf("%s\n", commandlinestring(argc, argv).c_str());
		printf("%s\n", versionstring(GAAEM_VERSION, __TIME__, __DATE__).c_str());		
		std::string controlfilename = argv[1];
		process(controlfilename);
		return EXIT_SUCCESS;
	}
}

int process(std::string controlfilename)
{
	cBlock C;
	fixseparator(controlfilename);
	printf("Loading control file %s\n", controlfilename.c_str());
	C.loadfromfile(controlfilename);

	std::string inputfile  = C.getstringvalue("Control.InputModelFile");
	std::string outputfile = C.getstringvalue("Control.OutputDataFile");
	std::string outputhdr  = C.getstringvalue("Control.OutputDataHeader");
	
	bool csvoutput = false;
	std::string ext = extractfileextension(outputfile);
	if (strcasecmp(ext, ".csv") == 0){
		csvoutput = true;
	}
	
	std::string ipmodel = C.getstringvalue("Control.IPModel");
	cLEM::IPType iptype = cLEM::IPType::NONE;
	if (strcasecmp(ipmodel, ud_string()) == 0){
		iptype = cLEM::IPType::NONE;
	}
	else if (strcasecmp(ipmodel, "none")==0){
		iptype = cLEM::IPType::NONE;
	}
	else if(strcasecmp(ipmodel, "colecole")==0){
		iptype = cLEM::IPType::COLECOLE;
	}
	else if (strcasecmp(ipmodel, "pelton")==0){
		iptype = cLEM::IPType::PELTON;
	}
	else{
		printf("Unknown IPModel %s: use none colecole or peltion\n", ipmodel.c_str());
		return EXIT_FAILURE;
	}

	std::string sysfile = C.getstringvalue("Control.SystemFile");
	printf("Opening AEM system file %s\n", sysfile.c_str());
	cTDEmSystem T(sysfile.c_str());
	T.LEM.iptype = (cLEM::IPType)iptype;

	printf("Opening input file %s\n", inputfile.c_str());
	FILE* fin = fileopen(inputfile, "r");
	printf("Opening output data file %s\n", outputfile.c_str());
	FILE* fout = fileopen(outputfile, "w");

	
	printf("Opening output header file %s\n", outputhdr.c_str());
	FILE* fhdr = fileopen(outputhdr, "w");
	writehdr(fhdr, T.NumberOfWindows);
	
	cTDEmResponse R;	
	size_t recnum = 1;
	char* CurrentRecordStr = new char[5001];
	while (fgets(CurrentRecordStr, 5000, fin) != NULL){
		printf("Processing record %zu\n", recnum);
		cTDEmGeometry G;
		cEarth1D E;
		parseinputrecord(CurrentRecordStr, G, E);
		printf("%s", CurrentRecordStr);						
		T.forwardmodel(G, E, R);
		if (recnum == 1){
			writecsvheader(fout, R.SX.size());
		}
		writeoutputrecord(csvoutput, fout, fhdr, recnum, T, R);
		recnum++;
	};
	printf("End of input\n");	
	delete[]CurrentRecordStr;
	fclose(fin);
	fclose(fout);
	fclose(fhdr);
	return 0;
}

int parseinputrecord(const char* record, cTDEmGeometry& G, cEarth1D& E)
{
	std::vector<double> v = getdoublevector(record," ,\t\r\n");
	G.tx_height = v[0];
	G.tx_roll = v[1];
	G.tx_pitch = v[2];
	G.tx_yaw = v[3];
	G.txrx_dx = v[4];
	G.txrx_dy = v[5];
	G.txrx_dz = v[6];
	G.rx_roll = v[7];
	G.rx_pitch = v[8];
	G.rx_yaw = v[9];


	size_t nlayers = (size_t)v[10];
	E.conductivity.resize(nlayers);
	E.thickness.resize(nlayers - 1);
	
	size_t k = 11;
	for (size_t i = 0; i < nlayers; i++){
		E.conductivity[i] = v[k];
		k++;
	}

	for (size_t i = 0; i < nlayers - 1; i++){
		E.thickness[i] = v[k];
		k++;
	}

	if (v.size() > k){
		E.chargeability.resize(nlayers);
		E.frequencydependence.resize(nlayers);
		E.timeconstant.resize(nlayers);

		for (size_t i = 0; i < nlayers; i++){
			E.chargeability[i] = v[k];
			k++;
		}		

		for (size_t i = 0; i < nlayers; i++){
			E.timeconstant[i] = v[k];
			k++;
		}

		for (size_t i = 0; i < nlayers; i++){
			E.frequencydependence[i] = v[k];
			k++;
		}
	}

	return 0;
}

int writehdr(FILE* fhdr, const size_t& nw)
{		
	size_t colnum = 1;
	writehdrentry(fhdr, "XP", colnum, 1);
	writehdrentry(fhdr, "YP", colnum, 1);
	writehdrentry(fhdr, "ZP", colnum, 1);
	writehdrentry(fhdr, "XS", colnum, nw);
	writehdrentry(fhdr, "YS", colnum, nw);
	writehdrentry(fhdr, "ZS", colnum, nw);
	return 0;
}

int writehdrentry(FILE* fhdr, const char* s, size_t& colnum, size_t nbands)
{
	if (nbands == 1){
		fprintf(fhdr, "%zu\t%s\n", colnum, s);
	}
	else{
		fprintf(fhdr, "%zu-%zu\t%s\n", colnum, colnum + nbands - 1, s);
	}
	colnum = colnum + nbands;
	return 0;
}

int writecsvheader(FILE* fout, const size_t& nw)
{	
	char delim = ',';	
	fprintf(fout, "XP%c",delim);
	fprintf(fout, "YP%c",delim);
	fprintf(fout, "ZP%c",delim);	
	for (size_t i = 0; i < nw; i++)fprintf(fout, "XS[%02zu]%c", i+1, delim);
	for (size_t i = 0; i < nw; i++)fprintf(fout, "YS[%02zu]%c", i+1, delim);
	for (size_t i = 0; i < nw; i++)fprintf(fout, "ZS[%02zu]%c", i+1, delim);
	fprintf(fout, "\n");
	return 0;
}

int writeoutputrecord(const bool& csvoutput, FILE* fout, FILE* fhdr, size_t recnum, const cTDEmSystem& T, const cTDEmResponse& R)
{
	char delim = ' ';
	if (csvoutput)delim = ',';	
	fprintf(fout, " %15g%c", R.PX,delim);	
	fprintf(fout, " %15g%c", R.PY, delim);	
	fprintf(fout, " %15g%c", R.PZ, delim);

	size_t nw = R.SX.size();	
	for (size_t i = 0; i < nw; i++)fprintf(fout, " %15g%c", R.SX[i], delim);	
	for (size_t i = 0; i < nw; i++)fprintf(fout, " %15g%c", R.SY[i], delim);	
	for (size_t i = 0; i < nw; i++)fprintf(fout, " %15g%c", R.SZ[i], delim);
	fprintf(fout, "\n");
	return 0;
}


