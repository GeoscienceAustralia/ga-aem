/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstring>

#include "gaaem_version.hpp"
#include "general_utils.hpp"
#include "file_utils.hpp"
#include "blocklanguage.hpp"
#include "earth1d.hpp"
#include "lem.hpp"
#include "tdemsystem.hpp"

class cLogger glog; //The global instance of the log file manager

using namespace AEM;

static int parseinputrecord(const char* record, TDEmGeometry& G, Earth1D& E)
{
	std::vector<double> v = getdoublevector(record, " ,\t\r\n");

	const size_t nf = v.size();
	if (nf < 11) {
		glog.errormsg(_SRC_, "There should be at least 11 columns per record\n");
	};

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


	size_t nlayers = (size_t) v[10];

	size_t n1 = 11 + 2 * nlayers - 1;
	size_t n2 = 11 + 5 * nlayers - 1;
	if (nf != n1 && nf != n2) {
		glog.errormsg(_SRC_, "For %d layers there should be either %zu (no IP) or %zu (for IP) columns per record\n",nlayers,n1,n2);
	};

	if (v.size() != 11 + 2 * nlayers - 1) {
	}

	E.conductivity.resize(nlayers);
	E.thickness.resize(nlayers - 1);

	size_t k = 11;
	for (size_t i = 0; i < nlayers; i++) {
		E.conductivity[i] = v[k];
		k++;
	}

	for (size_t i = 0; i < nlayers - 1; i++) {
		E.thickness[i] = v[k];
		k++;
	}

	if (nf > k) {
		E.chargeability.resize(nlayers);
		for (size_t i = 0; i < nlayers; i++) {
			E.chargeability[i] = v[k];
			k++;
		}

		E.timeconstant.resize(nlayers);
		for (size_t i = 0; i < nlayers; i++) {
			E.timeconstant[i] = v[k];
			k++;
		}

		E.frequencydependence.resize(nlayers);
		for (size_t i = 0; i < nlayers; i++) {
			E.frequencydependence[i] = v[k];
			k++;
		}
	}

	return 0;
}

static int writehdrentry(std::ofstream& ofs, const char* s, size_t& colnum, size_t nbands)
{
	if (nbands == 1) {
		ofs << strprint("%zu\t%s\n", colnum, s);
	}
	else {
		ofs << strprint("%zu-%zu\t%s\n", colnum, colnum + nbands - 1, s);
	}
	colnum = colnum + nbands;
	return 0;
}

static int writehdr(std::ofstream& ofs, const size_t& nw)
{
	size_t colnum = 1;
	writehdrentry(ofs, "XP", colnum, 1);
	writehdrentry(ofs, "YP", colnum, 1);
	writehdrentry(ofs, "ZP", colnum, 1);
	writehdrentry(ofs, "XS", colnum, nw);
	writehdrentry(ofs, "YS", colnum, nw);
	writehdrentry(ofs, "ZS", colnum, nw);
	return 0;
}

static int writecsvheader(std::ofstream& ofs, const size_t& nw)
{
	char delim = ',';
	ofs << strprint("XP%c", delim);
	ofs << strprint("YP%c", delim);
	ofs << strprint("ZP%c", delim);
	for (size_t i = 0; i < nw; i++) ofs << strprint("XS[%02zu]%c", i + 1, delim);
	for (size_t i = 0; i < nw; i++) ofs << strprint("YS[%02zu]%c", i + 1, delim);
	for (size_t i = 0; i < nw; i++) ofs << strprint("ZS[%02zu]%c", i + 1, delim);
	ofs << std::endl;
	return 0;
}

static int writeoutputrecord(const bool& csvoutput, std::ofstream& ofsout, size_t recnum, const TDEmSystem& T, const TDEmResponse& R) {
	char delim = ' ';
	if (csvoutput) delim = ',';
	ofsout << strprint(" %15g%c", R.primary(XCOMP), delim);
	ofsout << strprint(" %15g%c", R.primary(YCOMP), delim);
	ofsout << strprint(" %15g%c", R.primary(ZCOMP), delim);

	const size_t nw = R.size();
	for (size_t i = 0; i < nw; i++) ofsout << strprint(" %15g%c", R.secondary(XCOMP, i), delim);
	for (size_t i = 0; i < nw; i++) ofsout << strprint(" %15g%c", R.secondary(YCOMP, i), delim);
	for (size_t i = 0; i < nw; i++) ofsout << strprint(" %15g%c", R.secondary(ZCOMP, i), delim);
	ofsout << std::endl;
	return 0;
}

static int process(std::string controlfilename)
{
	cBlock C;
	fixseparator(controlfilename);
	glog.logmsg("Loading control file %s\n", controlfilename.c_str());
	C.loadfromfile(controlfilename);

	std::string inputfile = C.getstringvalue("Control.InputModelFile");
	std::string outputfile = C.getstringvalue("Control.OutputDataFile");
	std::string outputhdr = C.getstringvalue("Control.OutputDataHeader");

	bool csvoutput = false;
	std::string ext = extractfileextension(outputfile);
	if (strcasecmp(ext, ".csv") == 0) {
		csvoutput = true;
	}

	std::string ipmodel = C.getstringvalue("Control.IPModel");
	AEM::IPType iptype = AEM::IPType::NONE;
	if (strcasecmp(ipmodel, undefinedvalue<std::string>()) == 0) {
		iptype = AEM::IPType::NONE;
	}
	else if (strcasecmp(ipmodel, "none") == 0) {
		iptype = AEM::IPType::NONE;
	}
	else if (strcasecmp(ipmodel, "colecole") == 0) {
		iptype = AEM::IPType::COLECOLE;
	}
	else if (strcasecmp(ipmodel, "pelton") == 0) {
		iptype = AEM::IPType::PELTON;
	}
	else {
		glog.errormsg(_SRC_,"Unknown IPModel %s: use none, colecole, or peltion\n", ipmodel.c_str());
	}

	std::string sysfile = C.getstringvalue("Control.SystemFile");
	glog.logmsg("Opening AEM system file %s\n", sysfile.c_str());
	TDEmSystem T(sysfile.c_str());
	T.lem().set_iptype((AEM::IPType)iptype);

	glog.logmsg("Opening input file %s\n", inputfile.c_str());
	std::ifstream ofsin = ifstream_ex(inputfile);
	glog.logmsg("Opening output data file %s\n", outputfile.c_str());
	std::ofstream ofsout = ofstream_ex(outputfile);

	glog.logmsg("Opening output header file %s\n", outputhdr.c_str());
	std::ofstream ofshdr = ofstream_ex(outputhdr);
	writehdr(ofshdr, T.nWindows());

	size_t recnum = 1;
	std::string CurrentRecord;
	while (filegetline_ifs(ofsin, CurrentRecord)) {
		trim_inplace(CurrentRecord);
		if (CurrentRecord.size() == 0) {
			recnum++;
			continue;
		}
		glog.logmsg("Processing record %zu: ", recnum);
		TDEmGeometry G;
		Earth1D E;
		parseinputrecord(CurrentRecord.c_str(), G, E);
		glog.logmsg("%s\n", CurrentRecord.c_str());


		TDEmResponse R = T.forward_model(E, G);

		if (recnum == 1) writecsvheader(ofsout, R.size());
		writeoutputrecord(csvoutput, ofsout, recnum, T, R);
		recnum++;
	};
	glog.logmsg("End of input\n");
	return 0;
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		glog.logmsg("Usage: %s control_file_name\n", argv[0]);
		return EXIT_FAILURE;
	}
	else if (argc > 2) {
		glog.logmsg("**Error: Too many command line arguments\n");
		glog.logmsg("Usage: %s control_file_name\n", argv[0]);
		return EXIT_FAILURE;
	}
	else {
		try {
			glog.logmsg("Program 'gaforwardmodeltdem'\n");
			glog.logmsg("Geoscience Australia's Airborne Electromagnetic Layered Earth Forward Modelling\n\n");
			glog.logmsg("Working directory: %s\n", getcurrentdirectory().c_str());
			glog.logmsg("%s\n", commandlinestring(argc, argv).c_str());
			glog.logmsg("%s\n", versionstring(GAAEM_VERSION, __TIME__, __DATE__).c_str());
			std::string controlfilename = argv[1];
			process(controlfilename);
		}
		catch (std::exception& e) {
			std::cout << e.what();
			return EXIT_FAILURE;
		}
		return EXIT_SUCCESS;
	}
}
