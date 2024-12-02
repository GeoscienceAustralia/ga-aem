/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>

#include "gaaem_version.hpp"
#include "undefinedvalues.hpp"
#include "file_formats.hpp"
#include "general_types.hpp"
#include "general_utils.hpp"
#include "file_utils.hpp"
#include "blocklanguage.hpp"
#include "geometry3d.hpp"
#include "stopwatch.hpp"
#include "ctlinedata.hpp"
#include "filesplitter.hpp"

#include "ticpp.h"
using namespace ticpp;

class cLogger glog; //The global instance of the log file manager

class cSGridCreator {

private:
	const cCTLineData& D;

	cBlock mControl;
	bool   binary = true;
	std::string sgriddir;
	std::string sgridprefix;
	std::string sgridsuffix;


	//bool usedepthextent;
	//double depthextent;

	bool usecellalignment = false;
	double cellwidth = 0;

	bool isconstantthickness;
	std::vector<double> constantthickness;

	double NullInputConductivity;
	double NullOutputProperty;
	double NullBelowElevation;
	double NullBelowDepth;

public:

	cSGridCreator(const cBlock& control, const cCTLineData& _D) : D(_D)
	{
		mControl = control;
		getsgridoptions();
	};

	void process() {
		if (usecellalignment) {
			create_sgrid_prop_alignment_cells();
		}
		else {
			create_sgrid_prop_alignment_points();
		}
		savexml();
	}

	void getsgridoptions() {

		cBlock b = mControl.findblock("SGrid");
		binary = b.getboolvalue("Binary");
		sgriddir = b.getstringvalue("OutDir");
		fixseparator(sgriddir);
		add_trailing_separator(sgriddir);

		sgridprefix = "";
		sgridsuffix = "";
		b.getvalue("Prefix", sgridprefix);
		b.getvalue("Suffix", sgridsuffix);

		NullBelowElevation = b.getdoublevalue("NullBelowElevation");
		if (!isdefined(NullBelowElevation)) {
			NullBelowElevation = std::numeric_limits<double>::lowest();
		}

		NullBelowDepth = b.getdoublevalue("NullBelowDepth");
		if (!isdefined(NullBelowDepth)) {
			NullBelowDepth = std::numeric_limits<double>::max();;
		}

		NullOutputProperty = b.getdoublevalue("NullOutputProperty");
		if (!isdefined(NullOutputProperty)) {
			NullOutputProperty = -999;
		}

		usecellalignment = b.getboolvalue("UseCellAlignment");
		if (usecellalignment) {
			cellwidth = b.getdoublevalue("CellWidth");
		}
	}

	void create_sgrid_prop_alignment_cells() {

		std::string sgriddatafile = sgridname() + ".sg.data";
		std::string sgriddatapath = sgriddir + sgridname() + ".sg.data";
		std::string sgridhdrpath = sgriddir + sgridname() + ".sg";
		std::ofstream ofs = ofstream_ex(sgriddatapath);
		ofs << strprint("*\n");
		ofs << strprint("*   X   Y   Z  Conductivity I   J   K\n");
		ofs << strprint("*\n");

		for (int wi = 0; wi < 2; wi++) {
			for (int li = 0; li <= D.nlayers; li++) {
				for (int si = 0; si <= D.nsamples; si++) {
					double xc, yc, zc, ec, c0;
					cVec v;
					if (si > 0 && si < D.nsamples) {
						ec = (D.e[si - 1] + D.e[si]) / 2.0;
						xc = (D.x[si - 1] + D.x[si]) / 2.0;
						yc = (D.y[si - 1] + D.y[si]) / 2.0;
						zc = (D.z[si - 1][li] + D.z[si][li]) / 2.0;
						v = cVec(D.x[si] - D.x[si - 1], D.y[si] - D.y[si - 1], 0.0);
					}
					else if (si == 0) {//first cell
						ec = D.e[0] - (D.e[1] - D.e[0]) / 2.0;
						xc = D.x[0] - (D.x[1] - D.x[0]) / 2.0;
						yc = D.y[0] - (D.y[1] - D.y[0]) / 2.0;
						zc = D.z[0][li] - (D.z[1][li] - D.z[0][li]) / 2.0;
						v = cVec(D.x[1] - D.x[0], D.y[1] - D.y[0], 0.0);
					}
					else if (si == D.nsamples) {//last cell
						ec = D.e[D.nsamples - 1] + (D.e[D.nsamples - 1] - D.e[D.nsamples - 2]) / 2.0;
						xc = D.x[D.nsamples - 1] + (D.x[D.nsamples - 1] - D.x[D.nsamples - 2]) / 2.0;
						yc = D.y[D.nsamples - 1] + (D.y[D.nsamples - 1] - D.y[D.nsamples - 2]) / 2.0;
						zc = D.z[D.nsamples - 1][li] + (D.z[D.nsamples - 1][li] - D.z[D.nsamples - 2][li]) / 2.0;
						v = cVec(D.x[D.nsamples - 1] - D.x[D.nsamples - 2], D.y[D.nsamples - 1] - D.y[D.nsamples - 2], 0.0);
					}
					else {
						printf("Error\n");
					}

					double ang = 90;
					if (wi == 1)ang = -90;
					cVec vr = (cellwidth / 2.0) * v.rotate(ang, cVec(0.0, 0.0, 1.0)).unit();
					xc = xc + vr.x;
					yc = yc + vr.y;

					c0 = NullOutputProperty;
					if (wi == 0 && li < D.nlayers && si < D.nsamples) {
						c0 = D.c[si][li];
						if (c0 == NullInputConductivity) {
							c0 = NullOutputProperty;
						}
						else if (c0 <= 0.0) {
							c0 = NullOutputProperty;
						}
						else {
							c0 = D.c[si][li];
						}

						if (zc < NullBelowElevation) {
							c0 = NullOutputProperty;
						}

						double dc = ec - zc;
						if (dc > NullBelowDepth) {
							c0 = NullOutputProperty;
						}
					}

					ofs << strprint("%8.1f %9.1f %7.1f %10.6f %4d %4d %4d\n", xc, yc, zc, c0, si, li, wi);
				}
			}
		}

		std::ofstream ofs_hdr = ofstream_ex(sgridhdrpath);
		ofs_hdr << strprint("GOCAD SGrid 1\n");
		ofs_hdr << strprint("HEADER {\n");
		ofs_hdr << strprint("name:%s\n", sgridname().c_str());
		ofs_hdr << strprint("painted:true\n");
		ofs_hdr << strprint("*painted*variable:Conductivity\n");
		ofs_hdr << strprint("cage:false\n");
		ofs_hdr << strprint("volume:true\n");
		ofs_hdr << strprint("*volume*grid:false\n");
		ofs_hdr << strprint("*volume*transparency_allowed:false\n");
		ofs_hdr << strprint("*volume*points:false\n");
		ofs_hdr << strprint("shaded_painted:false\n");
		ofs_hdr << strprint("precise_painted:true\n");
		ofs_hdr << strprint("*psections*grid:false\n");
		ofs_hdr << strprint("*psections*solid:true\n");
		ofs_hdr << strprint("dead_cells_faces:false\n");
		ofs_hdr << strprint("}\n");

		ofs_hdr << strprint("\n");
		ofs_hdr << strprint("AXIS_N %zu %zu %d\n", D.nsamples + 1, D.nlayers + 1, 2);
		ofs_hdr << strprint("PROP_ALIGNMENT CELLS\n");
		ofs_hdr << strprint("ASCII_DATA_FILE %s\n", sgriddatafile.c_str());

		ofs_hdr << strprint("\n");
		ofs_hdr << strprint("PROPERTY 1 Conductivity\n");
		ofs_hdr << strprint("PROP_UNIT 1 S/m\n");
		ofs_hdr << strprint("PROP_NO_DATA_VALUE 1 -999\n");

		ofs_hdr << strprint("\n");
		ofs_hdr << strprint("END\n");
	}

	void create_sgrid_prop_alignment_points() {

		static bool flipendian = !isbigendian();
		//must flip binary bytes to bigendian if not natively bigendian (MSBFIRST)

		std::string sgridhdrpath = sgriddir + sgridname() + ".sg";

		std::ofstream ofs_asciidata;
		std::ofstream ofs_points;
		std::ofstream ofs_property;

		std::string asciidatapath = sgriddir + sgridname() + ".sg.data";
		std::string pointspath = sgriddir + sgridname() + "_points@@";
		std::string proppath = sgriddir + sgridname() + "_Conductivity@@";

		if (binary) {
			ofs_points = ofstream_ex(pointspath, std::ios_base::out | std::ios_base::binary);
			ofs_points = ofstream_ex(proppath, std::ios_base::out | std::ios_base::binary);
		}
		else {
			ofs_asciidata = ofstream_ex(asciidatapath);
			ofs_asciidata << strprint("*\n");
			ofs_asciidata << strprint("*   X   Y   Z  Conductivity  I   J   K\n");
			ofs_asciidata << strprint("*\n");
		};

		for (int li = 0; li < D.nlayers; li++) {
			for (int si = 0; si < D.nsamples; si++) {
				double xc, yc, zc, ec, c0;
				cVec v;
				if (si >= 0 && si <= D.nsamples) {
					ec = D.e[si];
					xc = D.x[si];
					yc = D.y[si];
					zc = (D.z[si][li] + D.z[si][li + 1]) / 2.0;
				}
				else {
					glog.errormsg(_SRC_,"Error: sample out of range\n");
				}

				c0 = NullOutputProperty;
				if (li < D.nlayers && si < D.nsamples) {
					c0 = D.c[si][li];
					if (c0 == NullInputConductivity) {
						c0 = NullOutputProperty;
					}
					else if (c0 <= 0.0) {
						c0 = 1e-6;
					}
					else {
						c0 = D.c[si][li];
					}

					if (zc < NullBelowElevation) {
						c0 = NullOutputProperty;
					}

					double dc = ec - zc;
					if (dc > NullBelowDepth) {
						c0 = NullOutputProperty;
					}
				}

				if (binary) {
					float fxc = (float)xc;
					float fyc = (float)yc;
					float fzc = (float)zc;
					float fc0 = (float)c0;
					if (flipendian) {
						swap_endian(&fxc, 1);
						swap_endian(&fyc, 1);
						swap_endian(&fzc, 1);
						swap_endian(&fc0, 1);
					}

					ofs_points.write(reinterpret_cast<char*>(&fxc), sizeof(float));
					ofs_points.write(reinterpret_cast<char*>(&fyc), sizeof(float));
					ofs_points.write(reinterpret_cast<char*>(&fzc), sizeof(float));
					ofs_points.write(reinterpret_cast<char*>(&fc0), sizeof(float));
				}
				else {
					ofs_asciidata << strprint("%8.1f %9.1f %7.1f %10.6f %4d %4d %4d\n", xc, yc, zc, c0, si, li, 0);
				}
			}
		}

		std::ofstream ofs_hdr = ofstream_ex(sgridhdrpath);
		ofs_hdr << strprint("GOCAD SGrid 1\n");
		ofs_hdr << strprint("HEADER {\n");
		ofs_hdr << strprint("name:%s\n", sgridname().c_str());
		ofs_hdr << strprint("painted:true\n");
		ofs_hdr << strprint("*painted*variable:Conductivity\n");
		ofs_hdr << strprint("ascii:on\n");
		ofs_hdr << strprint("double_precision_binary:off\n");
		ofs_hdr << strprint("cage:false\n");
		ofs_hdr << strprint("volume:true\n");
		ofs_hdr << strprint("*volume*grid:false\n");
		ofs_hdr << strprint("*volume*transparency_allowed:false\n");
		ofs_hdr << strprint("*volume*points:false\n");
		ofs_hdr << strprint("shaded_painted:false\n");
		ofs_hdr << strprint("precise_painted:true\n");
		ofs_hdr << strprint("*psections*grid:false\n");
		ofs_hdr << strprint("*psections*solid:true\n");
		ofs_hdr << strprint("dead_cells_faces:false\n");
		ofs_hdr << strprint("}\n");

		ofs_hdr << strprint("\n");
		ofs_hdr << strprint("AXIS_N %zu %zu %d\n", D.nsamples, D.nlayers, 1);
		ofs_hdr << strprint("PROP_ALIGNMENT POINTS\n");

		if (binary) {
			ofs_hdr << strprint("POINTS_FILE %s\n", extractfilename(pointspath).c_str());
		}
		else {
			ofs_hdr << strprint("ASCII_DATA_FILE %s\n", extractfilename(asciidatapath).c_str());
		}

		ofs_hdr << strprint("\n");
		ofs_hdr << strprint("PROPERTY 1 Conductivity\n");
		ofs_hdr << strprint("PROP_UNIT 1 S/m\n");
		ofs_hdr << strprint("PROP_NO_DATA_VALUE 1 -999\n");
		if (binary) {
			ofs_hdr << strprint("PROP_FILE 1 %s\n", extractfilename(proppath).c_str());
			ofs_hdr << strprint("PROP_ESIZE 1 4\n");
			ofs_hdr << strprint("PROP_ETYPE 1 IEEE\n");
			ofs_hdr << strprint("PROP_ALIGNMENT 1 POINTS\n");
			ofs_hdr << strprint("PROP_FORMAT 1 RAW\n");
			ofs_hdr << strprint("PROP_OFFSET 1 0\n");
		}


		ofs_hdr << strprint("\n");
		ofs_hdr << strprint("END\n");
	}

	std::string sgridname() {
		std::string name = sgridprefix + strprint("%d", D.linenumber) + sgridsuffix;
		return name;
	}

	std::string sgridhdrname() {
		std::string s = sgridname() + ".sg";
		return s;
	}

	std::string sgridhdrpath() {
		std::string s = sgriddir + sgridhdrname();
		return s;
	}

	std::string xmlname() {
		std::string s = sgridname() + ".xml";
		return s;
	}

	std::string xmlpath() {
		std::string s = sgriddir + xmlname();
		fixseparator(s);
		return s;
	}

	void savexml()
	{
		cBlock b = mControl.findblock("XML");
		if (b.Entries.size() == 0) return;
		std::string cp = b.getstringvalue("DataCachePrefix");
		if (cp[cp.size() - 1] != '/') cp += '/';
		std::string datacachename = cp + sgridhdrname();

		Document doc(xmlpath());
		std::string ver = "1.0";
		std::string enc = "UTF-8";
		std::string std = "yes";
		Declaration dec(ver, enc, std);
		doc.InsertEndChild(dec);

		Element l("Layer");
		l.SetAttribute("layerType", "VolumeLayer");
		l.SetAttribute("version", "1");
		l.InsertEndChild(Element("DisplayName", sgridname()));
		l.InsertEndChild(Element("URL", sgridhdrname()));
		l.InsertEndChild(Element("DataFormat", "GOCAD SGrid"));
		l.InsertEndChild(Element("DataCacheName", datacachename));
		l.InsertEndChild(Element("CoordinateSystem", D.inputdatumprojection));
		doc.InsertEndChild(l);
		doc.SaveFile();
	}

};

void save_dataset_xml(const std::string xmlpath, const std::string datasetname, const std::vector<std::string> names, const std::vector<std::string> urls)
{
	makedirectory_for(xmlpath);
	try
	{
		Element a, b;
		Document doc(xmlpath);

		Element dl("DatasetList");

		Element d("Dataset");
		d.SetAttribute("name", datasetname);

		for (size_t i = 0; i < names.size(); i++) {
			Element l("Layer");
			l.SetAttribute("name", names[i]);
			l.SetAttribute("url", urls[i]);
			d.InsertEndChild(l);
		}
		dl.InsertEndChild(d);
		doc.InsertEndChild(dl);
		doc.SaveFile();
	}
	catch (ticpp::Exception& ex)
	{
		std::cout << ex.what();
	}
}

int main(int argc, char** argv)
{
	try {
		if (argc >= 2) {
			glog.logmsg("Executing %s %s\n", argv[0], argv[1]);
			glog.logmsg("Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
			glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
		}
		else {
			glog.logmsg("Executing %s\n", argv[0]);
			glog.logmsg("Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
			glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
			glog.logmsg("Error: Not enough input arguments\n");
			glog.logmsg("Usage: %s controlfilename\n", argv[0]);
			return 0;
		}

		cBlock b(argv[1]);
		cBlock ib = b.findblock("Input");
		cBlock sb = b.findblock("Section");
		cBlock xb = b.findblock("XML");

		std::string headerfile;
		if (ib.getvalue("DfnFile", headerfile) == true) {
			glog.logmsg("Headerfile = %s\n", headerfile.c_str());
			glog.logmsg("Note: in future please use 'HeaderFile = ...' instead of 'DfnFile = ...'\n");
		}
		else if (ib.getvalue("HeaderFile", headerfile) == true) {
			glog.logmsg("Headerfile = %s\n", headerfile.c_str());
		}
		else {
			glog.logmsg("No Headerfile defined, columns numbers to be used'\n");
		}

		int linefieldindex = -1;
		std::string linefieldname;
		if (ib.getvalue("line", linefieldname) == false) {
			glog.logmsg("Error: you must define a line field name or column number for the line number field using 'Line = ...'\n");
			return 0;
		}

		cRange<int> r;
		std::vector<cAsciiColumnField> fields;
		bool status = cCTLineData::parse_column_range(linefieldname, r);
		if (status == true) {
			linefieldindex = r.from;
		}
		else {
			if (cHDRHeader::is_of_format(headerfile)) {
				cHDRHeader H(headerfile);
				fields = H.getfields();
				linefieldindex = H.column_range_by_name(linefieldname).from;
			}
			else if (cASEGGDF2Header::is_of_format(headerfile)) {
				cASEGGDF2Header A(headerfile);
				fields = A.getfields();
				linefieldindex = A.column_range_by_name(linefieldname).from;
			}
		}

		if (linefieldindex < 0) {
			glog.logmsg("Error: cannot find the line field %s\n", linefieldname.c_str());
			return 0;
		}

		std::string infiles = ib.getstringvalue("DataFiles");
		std::vector<std::string> filelist = DirectoryAccess::getfilelist_multi_pattern(infiles);
		if (filelist.size() == 0) {
			glog.logmsg("Error: no data files found matching %s\n", infiles.c_str());
			return 0;
		}

		std::vector<std::string> names;
		std::vector<std::string> urls;
		std::string xmldir;
		cStopWatch stopwatch;
		for (size_t i = 0; i < filelist.size(); i++) {
			glog.logmsg("Processing file %s %3zu of %3zu\n", filelist[i].c_str(), i + 1, filelist.size());

			std::string datafile = filelist[i];
			cFileSplitter FS(datafile, 0, linefieldindex);
			std::vector<std::string> L;
			while (FS.getnextgroup(L) > 0) {
				cCTLineData D(ib, fields);
				D.load(L);
				glog.logmsg("Line %d\n", D.linenumber);
				cSGridCreator S(b, D);
				S.process();
				if (xb.Entries.size() > 0) {
					names.push_back(S.sgridname());
					urls.push_back(S.xmlname());
					xmldir = extractfiledirectory(S.xmlpath());
				}
			}
		}
		if (xb.Entries.size() > 0) {
			std::string datasetname = xb.getstringvalue("DatasetName");
			std::string datasetxml = xmldir + datasetname + ".xml";
			save_dataset_xml(datasetxml, datasetname, names, urls);
		}
		printf("Done ... \nElapsed time = %.3lf seconds\n", stopwatch.etimenow());
	}
	catch (ticpp::Exception& e) {
		std::cout << e.what();
	}
	catch (std::runtime_error& e) {
		std::cout << e.what();
	}
	return 0;
}
