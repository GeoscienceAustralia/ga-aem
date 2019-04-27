/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>

#if defined _OPENMP
#include <omp.h>
#endif

#pragma warning( push )
#pragma warning (disable: 4251)
#include <gdal_priv.h>
#pragma warning( pop )
#include <gdal_alg.h>
#include <ogr_spatialref.h>

#include "general_utils.h"
#include "file_utils.h"
#include "vector_utils.h"
#include "asciicolumnfile.h"
#include "blocklanguage.h"
#define VERSION "1.0"
class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

using namespace std;

class cGridOptions{

public:	
	double cellsize;
	double xmin;
	double xmax;
	double ymin;
	double ymax;	
	GUInt32 nx(){
		return (GUInt32)((xmax - xmin) / cellsize);
	};
	GUInt32 ny(){
		return (GUInt32)((ymax - ymin) / cellsize);
	};
	GUInt32 ncells(){
		return nx()*ny();
	};
		
	double searchradius;	
	bool   gridinlog10space;
	std::string datum;
	std::string projection;	

	static float nullcellvalue(){
		return -32768;
	}

	static bool isnullcellvalue(const float& v){
		if(v == nullcellvalue())return true;
		else return false;
	}

	void setextents(const std::vector<double>& x, const std::vector<double>& y){

		auto result = std::minmax_element(x.begin(), x.end());
		xmin = *result.first;
		xmax = *result.second;
		result = std::minmax_element(y.begin(), y.end());
		ymin = *result.first;
		ymax = *result.second;
		
		xmin = rounddownnearest(xmin, cellsize) - cellsize / 2.0;
		xmax = roundupnearest(xmax, cellsize) + cellsize / 2.0;
		ymin = rounddownnearest(ymin, cellsize) - cellsize / 2.0;
		ymax = roundupnearest(ymax, cellsize) + cellsize / 2.0;
		
	}

};

std::vector<double> depth_slices(const std::vector<double>& conductivity, const std::vector<double>& interfacedepth, const std::vector<double>& dslicetop, const std::vector<double>& dslicebottom)
{
	double nullvalue = cGridOptions::nullcellvalue();
	size_t nslices = dslicetop.size();
	size_t nlayers = conductivity.size();
	std::vector<double> dslice(nslices, nullvalue);
	for (size_t si = 0; si < nslices; si++){
		double top = dslicetop[si];
		double bot = dslicebottom[si];
		double minimumintersection = (bot - top) / 2.0;

		if (bot < 0.0){
			dslice[si] = nullvalue;
			continue;
		}

		if (top < 0.0)top = 0.0;


		size_t itop = nlayers - 1;
		for (size_t li = 0; li < nlayers - 1; li++){
			if (top >= interfacedepth[li] && top <= interfacedepth[li + 1]){
				itop = li;
				break;
			}
		}

		size_t ibot = nlayers - 1;
		for (size_t li = itop; li < nlayers - 1; li++){
			if (bot >  interfacedepth[li] && bot <= interfacedepth[li + 1]){
				ibot = li;
				break;
			}
		}

		double tsum = 0.0;
		double csum = 0.0;
		for (size_t li = itop; li <= ibot; li++){
			double t = overlap(top, bot, interfacedepth[li], interfacedepth[li + 1]);
			csum += t*conductivity[li];
			tsum += t;
		}
		dslice[si] = csum / tsum;
		if (tsum < minimumintersection) dslice[si] = nullvalue;
	}
	return dslice;
}

std::vector<double> elevation_slices(const std::vector<double>& conductivity, const std::vector<double>& interfacedepth, const double& elevation, const std::vector<double>& eslicetop, const std::vector<double>& eslicebottom)
{
	std::vector<double> dslicetop = elevation - eslicetop;
	std::vector<double> dslicebottom = elevation - eslicebottom;
	std::vector<double> eslice = depth_slices(conductivity, interfacedepth, dslicetop, dslicebottom);
	return eslice;
}

class cSliceGridder{

private:
	cBlock Control;
	std::string DataFile;
	std::string DfnFile;
	cGridOptions GO;
	bool AutoExtents;
	
	size_t HeaderLines;
	size_t SubSample;	
	size_t nlayers;	

	std::vector<double> xdata;
	std::vector<double> ydata;
	std::vector<double> edata;
	std::vector<std::vector<double>> cdata;
	std::vector<std::vector<double>> ddata;

public:

	cSliceGridder(const std::string& controlfile){
		Control.loadfromfile(controlfile);				
		getgriddingoptions();
		readdatafile();		
		gridlayers();
		griddepthslices();
		gridelevationslices();
	}

	void getgriddingoptions(){

		cBlock b = Control.findblock("Options");
		GO.cellsize = b.getdoublevalue("CellSize");
		AutoExtents = b.getboolvalue("AutoExtents");
		if (AutoExtents == false){
			GO.xmin = b.getdoublevalue("XMin");
			GO.xmax = b.getdoublevalue("XMax");
			GO.ymin = b.getdoublevalue("YMin");
			GO.ymax = b.getdoublevalue("YMax");
		}
		GO.searchradius = b.getdoublevalue("SearchRadius");		
		GO.gridinlog10space = b.getboolvalue("GridInLog10Space");
		GO.datum = b.getstringvalue("Datum");
		GO.projection = b.getstringvalue("Projection");
	}

	std::vector<std::vector<double>> interpretstrides(const std::string& str){
		double x1, dx, x2;
		sscanf(str.c_str(), "%lf:%lf:%lf", &x1, &dx, &x2);
		std::vector<std::vector<double>> m;

		if (x1 <= x2){
			dx = abs(dx);
			for (double x = x1; x < x2; x += dx){
				std::vector<double> v(2);
				v[0] = x;
				v[1] = x + dx;
				m.push_back(v);
			}
		}
		else if (x1 > x2){
			dx = -abs(dx);
			for (double x = x1; x > x2; x += dx){
				std::vector<double> v(2);
				v[0] = x;
				v[1] = x + dx;
				m.push_back(v);
			}
		}
		return m;
	}

	std::vector<std::vector<double>> getslices(const cBlock& b){
		std::vector<std::vector<double>> slices;
		std::string ss = b.getstringvalue("Slices");
		if (isdefined(ss)){
			slices = interpretstrides(ss);
		}
		else{
			slices = b.getdoublematrix("Slices");
		}
		return slices;
	}

	bool readdatafile()	{		
		cBlock b = Control.findblock("Input");
		DataFile = b.getstringvalue("DataFile");
		DfnFile = b.getstringvalue("DfnFile");
		SubSample = b.getsizetvalue("Subsample");
		HeaderLines = b.getsizetvalue("HeaderLines");
		
		if (!isdefined(SubSample))SubSample = 1;
		if (!isdefined(HeaderLines))HeaderLines = 0;

		double t1 = gettime();
		printf("Reading input data file %s\n",DataFile.c_str());
		cAsciiColumnFile A;		
		A.parse_aseggdf2_header(DfnFile);

		std::vector<std::vector<double>> exclude = b.getdoublematrix("ExcludeLineRanges");
		std::string lfield = b.getstringvalue("LineNumber");
		size_t findex_l = A.fieldindexbyname(lfield);		
		if (exclude.size() > 0){
			if (findex_l == A.nullfieldindex()){
				printf("Error could not find %s\n", lfield.c_str());
				exit(1);
			}
		}

		std::string xfield = b.getstringvalue("Easting");				
		size_t findex_x = A.fieldindexbyname(xfield);
		if (findex_x == A.nullfieldindex()){
			printf("Error could not find %s\n", xfield.c_str());
			exit(1);
		}

		std::string  yfield = b.getstringvalue("Northing");		
		size_t findex_y = A.fieldindexbyname(yfield);
		if (findex_y == A.nullfieldindex()){
			printf("Error could not find %s\n", yfield.c_str());
			exit(1);
		}

		std::string efield = b.getstringvalue("Elevation");
		size_t findex_e = A.fieldindexbyname(efield);
		if (findex_e == A.nullfieldindex()){
			printf("Error could not find %s\n", efield.c_str());
			exit(1);
		}

		bool haveconductivity = false;
		std::string cfield = b.getstringvalue("Conductivity");		
		std::string rfield = b.getstringvalue("Resistivity");
		size_t findex_c;
		if (isdefined(cfield)){
			haveconductivity = true;
			findex_c = A.fieldindexbyname(cfield);
			if (findex_c == A.nullfieldindex()){
				printf("Error could not find Conductivity field: %s\n", cfield.c_str());
				exit(1);
			}
		}
		else if (isdefined(rfield)){
			findex_c = A.fieldindexbyname(rfield);
			if (findex_c == A.nullfieldindex()){
				printf("Error could not find Resistivity field: %s\n", rfield.c_str());
				exit(1);
			}
		}
		else{
			printf("Error could not find either a Conductivity or Resistivity field\n");
			exit(1);
		}
		nlayers = A.fields[findex_c].nbands;

		bool havedepth = false;
		bool havethickness = false;
		bool haveelevationtop = false;
		std::string dfield = b.getstringvalue("DepthTop");		
		std::string tfield = b.getstringvalue("Thickness");
		std::string etopfield = b.getstringvalue("ElevationTop");
		size_t findex_d;		
		if (isdefined(dfield)){
			havedepth = true;
			findex_d = A.fieldindexbyname(dfield);
			if (findex_d == A.nullfieldindex()){
				printf("Error could not find DepthTop field: %s\n", dfield.c_str());
				exit(1);
			}
		}
		else if (isdefined(tfield)){
			havethickness = true;
			findex_d = A.fieldindexbyname(tfield);
			if (findex_d == A.nullfieldindex()){
				printf("Error could not find Thickness field: %s\n", tfield.c_str());
				exit(1);
			}
		}
		else if (isdefined(etopfield)){
			haveelevationtop = true;
			findex_d = A.fieldindexbyname(etopfield);
			if (findex_d == A.nullfieldindex()){
				printf("Error could not find ElevationTop field: %s\n", etopfield.c_str());
				exit(1);
			}
		}
		else{
			printf("Error could not find either a DepthTop or Thickness field\n");
			exit(1);
		}
		
		double cscale = 1.0;
		double rscale = 1.0;
		if (haveconductivity){
			cscale = b.getdoublevalue("ConductivityScaling");
			if (!isdefined(cscale)){
				cscale = 1.0;
			}
		}
		else{
			rscale = b.getdoublevalue("ResistivityScaling");
			if (!isdefined(rscale)){
				rscale = 1.0;
			}
		}		
		
		double ta, tb;		
		ta = gettime();
		size_t nrecords = countlines(DataFile);		
		tb = gettime();
		printf("Record count time %lf records=%zu\n", tb - ta, nrecords);
		
		size_t nreserve = nrecords / SubSample;
		A.openfile(DataFile);	

		//Skip header lines
		for (size_t k = 0; k < HeaderLines; k++){
			A.readnextrecord();
		}

		xdata.reserve(nreserve);
		ydata.reserve(nreserve);
		edata.reserve(nreserve);

		cdata.resize(nlayers);		
		ddata.resize(nlayers);
		for (size_t i = 0; i < nlayers; i++){
			cdata[i].reserve(nreserve);
			ddata[i].reserve(nreserve);
		}
				
		while (A.readnextrecord()){
			A.parserecord();

			bool excludedline = false;
			if (exclude.size() > 0){
				double line;
				A.getfield(findex_l, line);
				for (size_t k = 0; k < exclude.size(); k++){
					if (line >= exclude[k][0] && line <= exclude[k][1]){
						excludedline = true;
						break;
					}
				}
			}
			if (excludedline)continue;

			double x, y, e;
			A.getfield(findex_x, x);
			xdata.push_back(x);

			A.getfield(findex_y, y);
			ydata.push_back(y);

			A.getfield(findex_e, e);
			edata.push_back(e);

			std::vector<double> v;
			A.getfield(findex_c, v);
			if (haveconductivity){
				for (int i = 0; i < (int) nlayers; i++){
					cdata[i].push_back(cscale*v[i]);
				}
			}
			else{
				for (int i = 0; i < (int) nlayers; i++){
					cdata[i].push_back(1.0/(rscale*v[i]));
				}
			}
			
			std::vector<double> w;
			A.getfield(findex_d, w);
			if (havedepth){
				for (int i = 0; i < (int) nlayers; i++){
					ddata[i].push_back(w[i]);
				}
			}
			else if(havethickness){
				double sum = 0.0;
				ddata[0].push_back(sum);
				for (int i = 0; i < (int) nlayers-1; i++){
					sum += w[i];
					ddata[i+1].push_back(sum);
				}
			}
			else if (haveelevationtop){
				double sum = 0.0;
				ddata[0].push_back(sum);				
				for (int i = 0; i < (int) nlayers - 1; i++){
					sum += w[i]-w[i+1];
					ddata[i + 1].push_back(sum);
				}
			}
			else{
				printf("Error do not have a depth-top/thickness/elevation-top field\n");
				exit(1);
			}

			//Skip due to subsampling
			for (size_t k = 0; k < SubSample-1; k++){
				A.readnextrecord();
			}
		}
		A.closefile();
		double t2 = gettime();
		printf("Reading input data elapsed time = %.2f\n", t2 - t1);

		if (AutoExtents){
			GO.setextents(xdata, ydata);
		}
		return true;
	}

	void gridlayers(){

		cBlock b = Control.findblock("Layers");
		bool create = b.getboolvalue("Create");
		if (create == false)return;
		printf("Gridding layers\n");

		std::string outdir = b.getstringvalue("OutputDir");
		outdir = outdir + pathseparatorstring();
		fixseparator(outdir);
		if (!exists(outdir)) makedirectorydeep(outdir);
		
		for (size_t i = 0; i < nlayers; i++){		
			double t1 = gettime();			
			double dtop = ddata[i][0];
			std::string gridfilename;
			if (i < nlayers - 1){				
				double dbot = ddata[i+1][0];
				gridfilename = outdir + strprint("layer_%02d_%05.1lf_%05.1lf_m", i + 1, dtop, dbot);
			}
			else{
				gridfilename = outdir + strprint("layer_%02d_%05.1lf_Inf_m", i + 1, dtop);
			}			
			gdalgrid(cdata[i], gridfilename);
			double t2 = gettime();
			printf("Gridding %s elapsed time = %.2f\n", gridfilename.c_str(), t2 - t1);
		}		
	}

	void griddepthslices(){

		cBlock b = Control.findblock("DepthSlices");
		bool create = b.getboolvalue("Create");
		if (create == false)return;
		printf("Gridding depth slices\n");

		std::string outdir = b.getstringvalue("OutputDir");
		outdir = outdir + pathseparatorstring();
		fixseparator(outdir);
		if (!exists(outdir)) makedirectorydeep(outdir);
		
		std::vector<std::vector<double>> slices = getslices(b);
		size_t nslices = slices.size();
		std::vector<double> dslicetop(nslices);
		std::vector<double> dslicebot(nslices);
		for (size_t i = 0; i < nslices; i++){
			if (slices[i][0] > slices[i][1]){
				std::swap(slices[i][0], slices[i][1]);
			}
			dslicetop[i] = slices[i][0];
			dslicebot[i] = slices[i][1];
		}

		std::vector<std::vector<double>> dsdata(nslices);
		for (size_t si = 0; si < nslices; si++){
			dsdata[si].resize(xdata.size());
		}

		std::vector<double> c(nlayers);
		std::vector<double> d(nlayers);
		for (size_t k = 0; k < xdata.size(); k++){			

			for (size_t li = 0; li < nlayers; li++){
				c[li] = cdata[li][k];
				d[li] = ddata[li][k];
			}
						
			std::vector<double> ds = depth_slices(c, d, dslicetop, dslicebot);
			for (size_t si = 0; si < nslices; si++){
				dsdata[si][k] = ds[si];
			}
		}
		
		for (size_t si = 0; si < nslices; si++){
			double t1 = gettime();
			std::string gridfilename;						
			gridfilename = outdir + strprint("dslice_%02d_%05.1lf_%05.1lf_m", si + 1, slices[si][0], slices[si][1]);
			gdalgrid(dsdata[si], gridfilename);
			double t2 = gettime();
			printf("Gridding %s elapsed time = %.2f\n", gridfilename.c_str(), t2 - t1);
		}
	}
	
	void gridelevationslices(){

		cBlock b = Control.findblock("ElevationSlices");
		bool create = b.getboolvalue("Create");
		if (create == false)return;
		printf("Gridding elevation slices\n");

		std::string outdir = b.getstringvalue("OutputDir");
		outdir = outdir + pathseparatorstring();
		fixseparator(outdir);
		if (!exists(outdir)) makedirectorydeep(outdir);
		
		std::vector<std::vector<double>> slices = getslices(b);
		size_t nslices = slices.size();
		std::vector<double> eslicetop(nslices);
		std::vector<double> eslicebot(nslices);
		for (size_t i = 0; i < nslices; i++){
			if (slices[i][0] < slices[i][1]){
				std::swap(slices[i][0], slices[i][1]);
			}
			eslicetop[i] = slices[i][0];
			eslicebot[i] = slices[i][1];
		}

		std::vector<std::vector<double>> esdata(nslices);
		for (size_t si = 0; si < nslices; si++){
			esdata[si].resize(xdata.size());
		}

		std::vector<double> c(nlayers);
		std::vector<double> d(nlayers);
		for (size_t k = 0; k < xdata.size(); k++){

			for (size_t li = 0; li < nlayers; li++){
				c[li] = cdata[li][k];
				d[li] = ddata[li][k];
			}			
			
			std::vector<double> es = elevation_slices(c, d, edata[k], eslicetop, eslicebot);
			for (size_t si = 0; si < nslices; si++){
				esdata[si][k] = es[si];
			}
		}

		for (size_t si = 0; si < nslices; si++){
			double t1 = gettime();
			std::string gridfilename;
			gridfilename = outdir + strprint("eslice_%02d_%05.1lf_%05.1lf_m", si + 1, slices[si][0], slices[si][1]);
			gdalgrid(esdata[si], gridfilename);
			double t2 = gettime();
			printf("Gridding %s elapsed time = %.2f\n", gridfilename.c_str(), t2 - t1);
		}
	}

	void gdalgrid(const std::vector<double>& zdata, const std::string& gridfilename){	
				
		GUInt32  npoints = (GUInt32)xdata.size();
		std::vector<double> x, y, z;
		x.reserve(npoints);
		y.reserve(npoints);
		z.reserve(npoints);		
		if (GO.gridinlog10space){			
			for (size_t k = 0; k < npoints; k++){
				if (zdata[k]>0 && zdata[k] != GO.nullcellvalue()){
					x.push_back(xdata[k]);
					y.push_back(ydata[k]);
					z.push_back(log10(zdata[k]));
				}
			}			
		}
		else{
			for (size_t k = 0; k < npoints; k++){
				if (zdata[k]>0 && zdata[k] != GO.nullcellvalue()){
					x.push_back(xdata[k]);
					y.push_back(ydata[k]);
					z.push_back(zdata[k]);
				}
			}
		}		

		GDALGridInverseDistanceToAPowerOptions *poOptions = new GDALGridInverseDistanceToAPowerOptions();
		poOptions->dfPower = 2;
		poOptions->dfRadius1 = GO.searchradius;
		poOptions->dfRadius2 = GO.searchradius;
		poOptions->dfNoDataValue = GO.nullcellvalue();
				
		std::vector<float> griddata(GO.ncells());
		GDALGridCreate(GGA_InverseDistanceToAPower, poOptions,
			npoints, x.data(), y.data(), z.data(),
			GO.xmin, GO.xmax,
			GO.ymin, GO.ymax,
			GO.nx(), GO.ny(),
			GDT_Float32, griddata.data(), (GDALProgressFunc)NULL, NULL);
		delete poOptions;
		
		if (GO.gridinlog10space){			
			for(size_t k = 0; k < GO.ncells(); k++){
				if (griddata[k] != GO.nullcellvalue()){
					griddata[k] = powf(10.0, griddata[k]);
				}
			}			
		}

		GDALDriver*  pdriver = GetGDALDriverManager()->GetDriverByName("ERS");
		GDALDataset* pdataset = pdriver->Create(gridfilename.c_str(), GO.nx(), GO.ny(), 1, GDT_Float32, NULL);
		double adfGeoTransform[6] = { GO.xmin, GO.cellsize, 0, GO.ymax, 0, -GO.cellsize };
		pdataset->SetGeoTransform(adfGeoTransform);
		
		OGRSpatialReference SRS;
		std::string units = "Meters";
		SRS.importFromERM(GO.projection.c_str(),GO.datum.c_str(),units.c_str());
		char* wkt = new char[2048];
		SRS.exportToWkt(&wkt);
		pdataset->SetProjection(wkt);
		CPLFree(wkt);

		
		//Write griddata to file
		//Because the buffer is south-up rather than north-down order
		//point rasterio to last line with negative line increment		
		int nLineSpace = -(int)GO.nx()*sizeof(float);
		float* plastline = &(griddata[(GO.ny() - 1)*GO.nx()]);
		pdataset->RasterIO(GF_Write, 0, 0, GO.nx(), GO.ny(), plastline, GO.nx(), GO.ny(), GDT_Float32, 1, 0, 0, nLineSpace, 0);
		GDALRasterBand* pband = pdataset->GetRasterBand(1);
		pband->SetNoDataValue(GO.nullcellvalue());
		GDALClose(pdataset);				
	}
};

int main(int argc, char** argv)
{
	if (argc >= 2){
		glog.logmsg("Executing %s %s\n", argv[0], argv[1]);
		glog.logmsg("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
		#if defined _OPENMP
		glog.logmsg("Maximum number of OpenMP Threads is %d\n", omp_get_max_threads());
		#endif
	}
	else{
		glog.logmsg("Executing %s %s\n", argv[0], argv[1]);
		glog.logmsg("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
		#if defined _OPENMP
		glog.logmsg("Maximum number of OpenMP Threads is %d\n", omp_get_max_threads());
		#endif
		glog.logmsg("Error: Not enough input arguments\n");
		glog.logmsg("Usage: %s controlfilename\n", argv[0]);
		return 0;;
	}

	GDALAllRegister();
	double t1 = gettime();
	cSliceGridder S(argv[1]);	
	double t2 = gettime();
	printf("Done ... Elapsed time = %.2lf seconds\n", t2 - t1);
	return 0;
}



 