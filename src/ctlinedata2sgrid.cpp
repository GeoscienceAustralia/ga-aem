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

#include "gaaem_version.h"
#include "undefinedvalues.h"
#include "file_formats.h"
#include "general_types.h"
#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "geometry3d.h"
#include "stopwatch.h"
#include "ctlinedata.h"
#include "filesplitter.h"

#include "ticpp.h"
using namespace ticpp;

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

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

		addtrailingseparator(sgriddir);

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
		FILE* fp_data = fileopen(sgriddatapath, "w");
		fprintf(fp_data, "*\n");
		fprintf(fp_data, "*   X   Y   Z  Conductivity I   J   K\n");
		fprintf(fp_data, "*\n");

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

					fprintf(fp_data, "%8.1f %9.1f %7.1f %10.6f %4d %4d %4d\n", xc, yc, zc, c0, si, li, wi);
				}
			}
		}
		fclose(fp_data);

		FILE* fp_hdr = fileopen(sgridhdrpath, "w");

		fprintf(fp_hdr, "GOCAD SGrid 1\n");
		fprintf(fp_hdr, "HEADER {\n");
		fprintf(fp_hdr, "name:%s\n", sgridname().c_str());
		fprintf(fp_hdr, "painted:true\n");
		fprintf(fp_hdr, "*painted*variable:Conductivity\n");
		fprintf(fp_hdr, "cage:false\n");
		fprintf(fp_hdr, "volume:true\n");
		fprintf(fp_hdr, "*volume*grid:false\n");
		fprintf(fp_hdr, "*volume*transparency_allowed:false\n");
		fprintf(fp_hdr, "*volume*points:false\n");
		fprintf(fp_hdr, "shaded_painted:false\n");
		fprintf(fp_hdr, "precise_painted:true\n");
		fprintf(fp_hdr, "*psections*grid:false\n");
		fprintf(fp_hdr, "*psections*solid:true\n");
		fprintf(fp_hdr, "dead_cells_faces:false\n");
		fprintf(fp_hdr, "}\n");

		fprintf(fp_hdr, "\n");
		fprintf(fp_hdr, "AXIS_N %zu %zu %d\n", D.nsamples + 1, D.nlayers + 1, 2);
		fprintf(fp_hdr, "PROP_ALIGNMENT CELLS\n");
		fprintf(fp_hdr, "ASCII_DATA_FILE %s\n", sgriddatafile.c_str());

		fprintf(fp_hdr, "\n");
		fprintf(fp_hdr, "PROPERTY 1 Conductivity\n");
		fprintf(fp_hdr, "PROP_UNIT 1 S/m\n");
		fprintf(fp_hdr, "PROP_NO_DATA_VALUE 1 -999\n");

		fprintf(fp_hdr, "\n");
		fprintf(fp_hdr, "END\n");
		fclose(fp_hdr);

	}

	void create_sgrid_prop_alignment_points() {

		static bool flipendian = !isbigendian();
		//must flip binary bytes to bigendian if not natively bigendian (MSBFIRST)

		std::string sgridhdrpath = sgriddir + sgridname() + ".sg";

		FILE* fp_asciidata = (FILE*)NULL;
		FILE* fp_points = (FILE*)NULL;
		FILE* fp_property = (FILE*)NULL;

		std::string asciidatapath = sgriddir + sgridname() + ".sg.data";
		std::string pointspath = sgriddir + sgridname() + "_points@@";
		std::string proppath = sgriddir + sgridname() + "_Conductivity@@";

		if (binary) {
			fp_points = fileopen(pointspath, "w+b");
			fp_property = fileopen(proppath, "w+b");
		}
		else {
			fp_asciidata = fileopen(asciidatapath, "w");
			fprintf(fp_asciidata, "*\n");
			fprintf(fp_asciidata, "*   X   Y   Z  Conductivity  I   J   K\n");
			fprintf(fp_asciidata, "*\n");
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
					printf("Error\n");
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

					fwrite(&fxc, sizeof(float), 1, fp_points);
					fwrite(&fyc, sizeof(float), 1, fp_points);
					fwrite(&fzc, sizeof(float), 1, fp_points);
					fwrite(&fc0, sizeof(float), 1, fp_property);
				}
				else {
					fprintf(fp_asciidata, "%8.1f %9.1f %7.1f %10.6f %4d %4d %4d\n", xc, yc, zc, c0, si, li, 0);
				}
			}
		}
		if (binary) {
			fclose(fp_points);
			fclose(fp_property);
		}
		else {
			fclose(fp_asciidata);
		}


		FILE* fp_hdr = fileopen(sgridhdrpath, "w");
		fprintf(fp_hdr, "GOCAD SGrid 1\n");
		fprintf(fp_hdr, "HEADER {\n");
		fprintf(fp_hdr, "name:%s\n", sgridname().c_str());
		fprintf(fp_hdr, "painted:true\n");
		fprintf(fp_hdr, "*painted*variable:Conductivity\n");
		fprintf(fp_hdr, "ascii:on\n");
		fprintf(fp_hdr, "double_precision_binary:off\n");
		fprintf(fp_hdr, "cage:false\n");
		fprintf(fp_hdr, "volume:true\n");
		fprintf(fp_hdr, "*volume*grid:false\n");
		fprintf(fp_hdr, "*volume*transparency_allowed:false\n");
		fprintf(fp_hdr, "*volume*points:false\n");
		fprintf(fp_hdr, "shaded_painted:false\n");
		fprintf(fp_hdr, "precise_painted:true\n");
		fprintf(fp_hdr, "*psections*grid:false\n");
		fprintf(fp_hdr, "*psections*solid:true\n");
		fprintf(fp_hdr, "dead_cells_faces:false\n");
		fprintf(fp_hdr, "}\n");

		fprintf(fp_hdr, "\n");
		fprintf(fp_hdr, "AXIS_N %zu %zu %d\n", D.nsamples, D.nlayers, 1);
		fprintf(fp_hdr, "PROP_ALIGNMENT POINTS\n");

		if (binary) {
			fprintf(fp_hdr, "POINTS_FILE %s\n", extractfilename(pointspath).c_str());
		}
		else {
			fprintf(fp_hdr, "ASCII_DATA_FILE %s\n", extractfilename(asciidatapath).c_str());
		}

		fprintf(fp_hdr, "\n");
		fprintf(fp_hdr, "PROPERTY 1 Conductivity\n");
		fprintf(fp_hdr, "PROP_UNIT 1 S/m\n");
		fprintf(fp_hdr, "PROP_NO_DATA_VALUE 1 -999\n");
		if (binary) {
			fprintf(fp_hdr, "PROP_FILE 1 %s\n", extractfilename(proppath).c_str());
			fprintf(fp_hdr, "PROP_ESIZE 1 4\n");
			fprintf(fp_hdr, "PROP_ETYPE 1 IEEE\n");
			fprintf(fp_hdr, "PROP_ALIGNMENT 1 POINTS\n");
			fprintf(fp_hdr, "PROP_FORMAT 1 RAW\n");
			fprintf(fp_hdr, "PROP_OFFSET 1 0\n");
		}


		fprintf(fp_hdr, "\n");
		fprintf(fp_hdr, "END\n");
		fclose(fp_hdr);
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
	makedirectorydeep(extractfiledirectory(xmlpath));
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
		std::vector<std::string> filelist = cDirectoryAccess::getfilelist(infiles);
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
