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
#include <valarray>
#include <cstring>
#include <iostream>
#include <iomanip>

#include "gaaem_version.h"
#include "general_types.h"
#include "general_utils.h"
#include "general_types.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "geometry3d.h"
#include "crs.h"
#include "gdal_utils.h"
#include "stretch.h"
#include "colormap.h"
#include "gdiplus_utils.h"
#include "stopwatch.h"
#include "filesplitter.h"
#include "asciicolumnfile.h"
#include "ctlinedata.h"

#include "ticpp.h"
using namespace ticpp;

#include "RamerDouglasPeucker.h"
using namespace RDP;

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

template<typename T> 
class cBIL {
	
private:

	std::vector<T> data;
	const int nhpixels;
	const int nvpixels;

public:

	cBIL(const int& _nhpixels, const int& _nvpixels) : 
		nhpixels(_nhpixels), nvpixels(_nvpixels)
	{		
		data.resize(nhpixels * nvpixels, (T)0.0);
	}

	const T nodata_value() {
		return (T) -32767;
	}
			
	std::string type_string(const int& v) { return "int"; }
	std::string type_string(const float& v) { return "float"; }
	std::string type_string(const double& v) { return "double"; }

	std::string type_string() {
		return type_string(data[0]);
	}

	void SetPixel(const size_t& i, const size_t& j, const T& val)
	{
		data[(i - 1)*nhpixels + j];
	}

	void save(const std::string& filename)
	{
		sFilePathParts fpp(filename);
		std::string bilpath = fpp.directory + pathseparatorstring() + fpp.prefix + ".bil";
		std::string hdrpath = fpp.directory + pathseparatorstring() + fpp.prefix + ".hdr";
		savebil(bilpath);
		savehdr(hdrpath);
	}

	void savebil(const std::string& bilfilename)
	{				
		std::ofstream ofs(bilfilename, std::ios::out | std::ios::binary);
		ofs.write((char*)(&data[0]), sizeof(T)*data.size());
	}

	void savehdr(const std::string& hdrfilename)
	{
		std::ofstream ofs(hdrfilename);
		ofs << "ncols " << nhpixels << std::endl;
		ofs << "nrows " << nvpixels << std::endl;
		ofs << "cellsize " << 1 << std::endl;
		ofs << "xllcorner " << 0 << std::endl;
		ofs << "yllcorner " << 0 << std::endl;
		ofs << "nodata_value " << nodata_value() << std::endl;
		ofs << "nbits " << 8*sizeof(T) << std::endl;
		ofs << "pixeltype " << type_string() << std::endl;
		if (isbigendian()) ofs << "byteorder " << "msb" << std::endl;		
		else ofs << "byteorder " << "lsb" << std::endl;				
	}
};

class cCurtainImageSection {

private:

	const cCTLineData& D;	
	size_t sequence_number = 0;
	
	int nhpixels = 0;
	int nvpixels = 0;	
	double h0 = 0.0;//Centre of left pixel column
	double h1 = 0.0;//Centre of right pixel column
	double dh = 0.0;//Horizontal pixel size
	double v0 = 0.0;//Centre of bottom pixel row
	double v1 = 0.0;//Centre of top pixel row
	double dv = 0.0;//Vertical pixel size
	std::vector<double> imagefid;
	std::vector<double> imagex;
	std::vector<double> imagey;	
	std::vector<double> imageelevation;

	double gratdiv = 0.0;

	std::string outdir;
	std::string prefix;
	std::string suffix;

	cColorMap cmap;
	cStretch stretch;
	std::vector<double> cbarticks;

	Color BkgColor;
	Color AirColor;
	Color NullsColor;
	
	bool autozsectiontop=true;
	bool autozsectionbot=true;

	bool   spreadfade = 0.0;
	double spreadfadelowclip = 0.0;
	double spreadfadehighclip = 0.0;
	double lowspreadfade = 0.0;
	double highspreadfade = 0.0;

	double geometrytolerance = 0.0;
	int    tilesize = 0;
	std::string datasetname;
	std::string datacachename;

	//std::string imageformat = "png";
	//std::string worldfileextension = "pngw";

	std::string imageformat = "jpg";
	std::string worldfileextension = "jgw";
public:
	
	cCurtainImageSection(const cCTLineData& _D) : D(_D)
	{		
		BkgColor = Color(0, 255, 255, 255);
		AirColor = Color(0, 255, 255, 255);
		NullsColor = Color(255, 128, 128, 128);
	}

	std::string image_format() {
		return imageformat;
	}

	std::string worldfile_extension() {
		return worldfileextension;
	}

	void set_sequence_number(const size_t& seqn){
		sequence_number = seqn;
	}

	void getoptions(const cBlock& b){
		
		spreadfade = b.getboolvalue("SpreadFade");
		if(spreadfade){
			spreadfadelowclip   = b.getdoublevalue("Log10SpreadLowClip");
			spreadfadehighclip  = b.getdoublevalue("Log10SpreadHighClip");
			lowspreadfade       = b.getdoublevalue("LowSpreadFade");
			highspreadfade      = b.getdoublevalue("HighSpreadFade");	
		}

		v1 = b.getdoublevalue("ElevationTop");
		if(isdefined(v1)) autozsectiontop = false;					
		else autozsectiontop = true;
		
		v0 = b.getdoublevalue("ElevationBottom");
		if (isdefined(v0)) autozsectionbot = false;			
		else autozsectionbot = true;

		outdir = b.getstringvalue("OutDir");
		addtrailingseparator(outdir);

		prefix = b.getstringvalue("Prefix");
		suffix = b.getstringvalue("Suffix");

		if (!b.getvalue("HorizontalResolution", dh)){
			std::string msg("HorizontalResolution must be defined\n");
			throw(std::runtime_error(msg));
		}

		if (!b.getvalue("VerticalResolution", dv)){		
			std::string msg("VerticalResolution must be defined\n");
			throw(std::runtime_error(msg));
		}

		cBlock cs = b.findblock("ColourStretch");
		cmap      = cColorMap(cs);
		stretch   = cStretch(cs);
		cbarticks = cs.getdoublevector("ColourBarTicks");

		std::vector<int> bkgclr = cs.getintvector("BackgroundColour");
		std::vector<int> airclr = cs.getintvector("AirColour");
		std::vector<int> nullsclr = cs.getintvector("NullsColour");

		if (bkgclr.size() == 4) BkgColor = Color(bkgclr[0], bkgclr[1], bkgclr[2], bkgclr[3]);
		if (airclr.size() == 4) AirColor = Color(airclr[0], airclr[1], airclr[2], airclr[3]);
		if (nullsclr.size() == 4) NullsColor = Color(nullsclr[0], nullsclr[1], nullsclr[2], nullsclr[3]);
				
		gratdiv = b.getdoublevalue("ElevationGridDivision");

		geometrytolerance = b.getdoublevalue("GeometryTolerance");
		tilesize = b.getintvalue("TileSize");		
		datasetname   = b.getstringvalue("DatasetName");
		datacachename = b.getstringvalue("DataCacheName");
	}
	
	void process(){					
		calculateextents();

		Bitmap bm(nhpixels, nvpixels, PixelFormat32bppARGB);		
		createbitmap(bm);
		generatesectiondata(bm);
		//drawgraticule();		
		//fixtransparenttext();
		saveimage(bm);
		save_world_file();
		savegeometry();		
		saveribbontilerbatchcommand();		

		//cBIL<float> bil(nhpixels, nvpixels);
		//generatesectiondata(bil);
		//saveimage(bil);
		//savegeometry();
		//saveribbontilerbatchcommand();
	}

	void calculateextents()
	{												
		double hlength = D.linedistance.back();		
		hlength = roundupnearest(hlength,dh);			
		h0 = 0.0;
		h1 = hlength;
		nhpixels  = 1+(int)(hlength/dh);
		
		double zmin =  DBL_MAX;
		double zmax = -DBL_MAX;
		for(int si=0; si<D.nsamples; si++){
			if(D.z[si][D.nlayers] < zmin)zmin = D.z[si][D.nlayers];
			if (D.z[si][0]       > zmax)zmax = D.z[si][0];
		}

		if(autozsectionbot==true) v0 = zmin;
		if(autozsectiontop==true) v1 = zmax;		
		double vlength = v1-v0;
		//Adjust extents to nearest pixel
		vlength = roundupnearest(vlength,dv);		
		v0 = v1-vlength;
		nvpixels = 1 + (int)(vlength/dv);
	}

	int wh2ix(double wh){
		int ix = (int)roundnearest((wh-h0)/dh,1.0);
		return ix;
	}

	int wv2iy(double wv){		
		int iy = (int)roundnearest((v1-wv)/dv,1.0);
		return iy;
	}
		
	void createbitmap(Bitmap& bm) {
		for (int i = 0; i < nhpixels; i++) {
			for (int j = 0; j < nvpixels; j++) {
				bm.SetPixel(i, j, Color(255, 255, 255, 255));
			}
		}
	}

	double fadevalue(const double& spread){
		if(spread <= spreadfadelowclip)return  lowspreadfade;
		if(spread >= spreadfadehighclip)return highspreadfade;
		return lowspreadfade + (spread - spreadfadelowclip)/(spreadfadehighclip - spreadfadelowclip) * (highspreadfade - lowspreadfade);		
	}

	void fadecolor(Color& clr, const double& fade){
		//fade=0 gives original color
		//fade=1 gives white
		double alpha = 1.0-fade;
		unsigned char w = (unsigned char)(255.0*fade);		
		unsigned char r = w + (unsigned char)(clr.GetRed()   * alpha);
		unsigned char g = w + (unsigned char)(clr.GetGreen() * alpha);
		unsigned char b = w + (unsigned char)(clr.GetBlue()  * alpha);		
		clr = Color(255,r,g,b);				
	}

	void generatesectiondata(Bitmap& bm){

		int vp0 = wv2iy(v0);
		int vp1 = wv2iy(v1);

		imagefid.resize(nhpixels);
		imagex.resize(nhpixels);
		imagey.resize(nhpixels);
		imageelevation.resize(nhpixels);
				
		int k = 0;
		for (int i = 0; i < nhpixels; i++){
			double hp = h0 + (double)i*dh;
			while (k < D.linedistance.size() - 1 && std::abs(D.linedistance[k] - hp) > std::abs(D.linedistance[k + 1] - hp)){
				k++;
			}

			bool is_in_gap = false;
			if (std::abs(D.linedistance[k] - hp) > (2.0 * dh)) {
				is_in_gap = true;
			}
			
			int a = k;
			int b = k + 1;
			if (k == D.linedistance.size()-1){
				a = k-1; b = k;
			}
			imagefid[i]       = linearinterp(D.linedistance[a], D.fid[a], D.linedistance[b], D.fid[b], hp);
			imagex[i]         = linearinterp(D.linedistance[a], D.x[a],   D.linedistance[b], D.x[b], hp);
			imagey[i]         = linearinterp(D.linedistance[a], D.y[a],   D.linedistance[b], D.y[b], hp);;
			imageelevation[i] = linearinterp(D.linedistance[a], D.e[a],   D.linedistance[b], D.e[b], hp);
			

			/*if (i > 10270) {
				std::cout << i << " ";
				std::cout << k << " ";
				std::cout << std::setprecision(10) << D.linedistance[k] << " ";
				std::cout << std::setprecision(10) << D.x[k] << " ";
				std::cout << std::setprecision(10) << D.y[k] << " ";
				std::cout << std::setprecision(10) << D.y[k+1] << " ";
				std::cout << std::endl;
			}*/

			for(int j=vp1; j<=vp0; j++){				
				bm.SetPixel(i,j,BkgColor);

				if (is_in_gap) {
					bm.SetPixel(i, j, AirColor);
					continue;
				}

				double zp = v1 - (double)j*dv;
				if(zp>D.e[k]){
					bm.SetPixel(i,j,AirColor);					
				}
				else{					
					for (size_t li = 0; li<D.nlayers; li++){
						if (zp < D.z[k][li] && zp >= D.z[k][li + 1]){
							double conductivity = D.c[k][li];
							if(conductivity < 0.0){
								bm.SetPixel(i,j,NullsColor);
							}
							else{							
								int ind = stretch.index(conductivity);
								Color clr(255,cmap.r[ind],cmap.g[ind],cmap.b[ind]);

								if(spreadfade){
									double conductivity_p10 = D.cp10[k][li];
									double conductivity_p90 = D.cp90[k][li];
									double spread = log10(conductivity_p90) - log10(conductivity_p10);
									double fade   = fadevalue(spread);
									fadecolor(clr,fade);
								}								
								bm.SetPixel(i,j,clr);
							}
							break;
						}						
					}										
				}//layer loop
			}//v pixel loop
		}//h pixel loop		
	}

	template<typename T>
	void generatesectiondata(cBIL<T>& b) {

		T BkgValue(0.0);
		T AirValue(0.0);
		T NullsValue(0.0);

		int vp0 = wv2iy(v0);
		int vp1 = wv2iy(v1);

		int mini = 0;
		for (int i = 0; i <= nhpixels; i++) {
			double hp = h0 + (double)i*dh;

			while (mini < D.linedistance.size() - 1 && 
				std::abs(D.linedistance[mini] - hp) > 
				std::abs(D.linedistance[mini + 1] - hp)) {
				mini++;
			}
			

			for (int j = vp1; j <= vp0; j++) {
				b.SetPixel(i, j, BkgValue);				

				double zp = v1 - (double)j*dv;
				if (zp > D.e[mini]) {
					b.SetPixel(i, j, AirValue);
				}
				else {
					for (int li = 0; li < D.nlayers; li++) {
						if (zp < D.z[mini][li] && zp >= D.z[mini][li + 1]) {
							double conductivity = D.c[mini][li];
							if (conductivity < 0.0){
								b.SetPixel(i, j, NullsValue);
							}
							else {
								//int ind = stretch.index(conductivity);
								//Color clr(255, cmap.r[ind], cmap.g[ind], cmap.b[ind]);
								//if (spreadfade) {
								//	double conductivity_p10 = D.cp10[mini][li];
								//	double conductivity_p90 = D.cp90[mini][li];
								//	double spread = log10(conductivity_p90) - log10(conductivity_p10);
								//	double fade = fadevalue(spread);
								//	fadecolor(clr, fade);
								//}
								b.SetPixel(i, j, (T)conductivity);
							}
							break;
						}
					}
				}//layer loop
			}//v pixel loop
		}//h pixel loop		
	}

	void drawgraticule(Bitmap& bm){

		if (gratdiv == undefinedvalue<double>()) return;
		if (gratdiv <= 0.0) return;

		Graphics gr(&bm);
		gr.SetCompositingMode(CompositingModeSourceOver);						
		gr.SetTextRenderingHint(TextRenderingHintAntiAlias);		

		int zg1=(int)roundnearest(v0,(int)gratdiv);
		int zg2=(int)roundnearest(v1,(int)gratdiv);

		Pen blackpen(Color::Black,0);		
		SolidBrush blackbrush(Color::Black);
		
		//Font font(L"Arial", getfontsize(), FontStyleRegular, UnitPoint);				
		Font font(L"Arial", 12, FontStyleBold, UnitPoint);

		StringFormat textformat;
		textformat.SetAlignment(StringAlignmentFar);
		textformat.SetLineAlignment(StringAlignmentCenter);
		SizeF layoutsize(32767,32767);
		SizeF textsize;				
		gr.MeasureString(L"-0000 m",-1,&font,layoutsize,&textformat,&textsize);

		for(int zg=zg1; zg<=zg2; zg+=(int)gratdiv){
			gr.DrawLine(&blackpen,wh2ix(h0),wv2iy(zg),wh2ix(h1),wv2iy(zg));

			wchar_t s[20];				
			swprintf(s,20,L"%5d m",zg);						
			gr.MeasureString(s,-1,&font,layoutsize,&textformat,&textsize);

			int tx = wh2ix(h1);
			int ty = wv2iy(zg);
			if(ty > nvpixels)continue;
			if(ty-textsize.Height < 0)continue;	
			PointF txpos((Gdiplus::REAL)tx, (Gdiplus::REAL)ty);
			gr.DrawString(s, -1, &font, txpos, &textformat, &blackbrush);
		}		    				
	}

	void createcolorbar(){
		makedirectorydeep(extractfiledirectory(colorbarfile()));
		std::string title = "Conductivity (S/m)";
		Bitmap* bm = cGDIplusHelper::colorbar(cmap, stretch, title, cbarticks);
		cGDIplusHelper::saveimage(bm, colorbarfile());
		delete bm;
	}

	void fixtransparenttext(Bitmap& bm){
		int hp0 = wh2ix(h0);
		int hp1 = wh2ix(h1);
		int vp1 = wv2iy(v0);
		int vp2 = wv2iy(v1);
		for(int i=hp0; i<=hp1; i++){		
			for(int j=vp2; j<=vp1; j++){
				Color c;
				bm.GetPixel(i,j,&c);				
				bm.SetPixel(i, j, Color(255, c.GetR(), c.GetG(), c.GetB()));
			}
		}		
	}

	std::string getdatasetname(){
		std::string s = datasetname;
		return s;
	}

	std::string basename(){
		std::string bn = prefix + strprint("%d", D.linenumber) + suffix;
		return bn;
	}

	std::string tilesetdir_nod(){
		std::string s = datasetname + "\\";;
		return s;
	}

	std::string tilesetdir(){
		std::string s = outdir + tilesetdir_nod();
		return s;
	}

	std::string imagefile_nod(){
		std::string s = image_format() + "\\" + basename() + "." + image_format();
		return s;
	}
	
	std::string imagefile(){
		std::string s = outdir + imagefile_nod();
		return s;
	}
		
	std::string worldfile_nod() {
		std::string s = image_format() + "\\" + basename() + "." + worldfile_extension();
		return s;
	}

	std::string worldfile() {
		std::string s = outdir + worldfile_nod();
		return s;
	}

	std::string xyfile() {
		std::string s = outdir + "xy\\" + basename() + ".xy";
		return s;
	}

	std::string bilfile_nod() {
		std::string s = "bil\\" + basename() + ".bil";
		return s;
	}

	std::string bilfile() {
		std::string s = outdir + bilfile_nod();
		return s;
	}

	std::string xmlname(){
		std::string s = basename() + ".xml";
		return s;
	}

	std::string xmlpath(){
		std::string s = tilesetdir() + xmlname();
		return s;
	}

	std::string datasetxmlpath(){
		std::string s = tilesetdir() + datasetname + ".xml";
		return s;
	}

	std::string ribbontilerbatchfilepath(){
		std::string s = outdir + "\\run_ribbon_tiler.bat";
		return s;
	}

	std::string colorbarfile(){
		std::string s = tilesetdir() + "colourbar." + image_format();
		return s;
	}

	std::string colorbarfile_nod(){
		std::string s = "colourbar." + image_format();
		return s;
	}
	
	void saveimage(Bitmap& bm){		
		cGDIplusHelper::saveimage(&bm, imagefile());
	}

	template<typename T>
	void saveimage(cBIL<T>& bil) {
		bil.savebil(bilfile());
		bil.savehdr(bilfile());
	}
	
	void savegeometry(){
		
		std::vector<RDP::Point> pl(D.x.size());
		std::vector<RDP::Point> plout;
		for (size_t i = 0; i < D.x.size(); i++){
			pl[i] = RDP::Point(D.x[i], D.y[i]);
		}
		
		RDP::RamerDouglasPeucker(pl, geometrytolerance, plout);
		
		std::string ppath = outdir + basename() + ".points.dat";
		std::string pfpath = outdir + basename() + ".points_filtered.dat";
		//savepoints(pl, ppath);
		//savepoints(plout, pfpath);

		saveimage_geometry();
				
		std::vector<double> x(plout.size());
		std::vector<double> y(plout.size());		
		for (size_t i = 0; i < plout.size(); i++){
			x[i] = plout[i].first;
			y[i] = plout[i].second;
		}

		std::vector<double> longitude(plout.size());
		std::vector<double> latitude(plout.size());		
		en2lonlat(x, y, longitude, latitude);
		std::string pfpathll = outdir + basename() + ".points_filtered_geodetic.dat";		
		savexml(longitude, latitude);
	}

	void en2lonlat(const std::vector<double>& e_in, const std::vector<double>& n_in, std::vector<double>& lon_out, std::vector<double>& lat_out)
	{
		int inepsgcode = cCRS::epsgcode(D.inputdatumprojection);
		if (inepsgcode < 0) {
			std::string msg = strprint("Invalid DatumProjection %s was specified\n", D.inputdatumprojection.c_str()) + _SRC_;
			throw(std::runtime_error(msg));
		}		
		int outepsgcode = cCRS::epsgcode("WGS84|GEODETIC");
		//transform(inepsgcode, e_in, n_in, outepsgcode, lon_out, lat_out);
		transform(inepsgcode, e_in, n_in, outepsgcode, lat_out, lon_out);
	}

	void savepoints(const std::vector<RDP::Point>& p, const std::string& filename){
		FILE* fp = fileopen(filename, "w");
		for (size_t i = 0; i < p.size(); i++){
			fprintf(fp, "%10lf %10lf\n", p[i].first,p[i].second);
		}				
		fclose(fp);
	};

	void savepoints(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename)
	{
		FILE* fp = fileopen(filename, "w");
		for (size_t i = 0; i < x.size(); i++){
			fprintf(fp, "%10.2lf %10.2lf\n", x[i], y[i]);
		}
		fclose(fp);
	};

	void saveimage_geometry()
	{
		std::vector<double> imagelon(nhpixels);
		std::vector<double> imagelat(nhpixels);		
		en2lonlat(imagex, imagey, imagelon, imagelat);
		std::string xypath = outdir + "geometry\\" + basename() + ".path.txt";
		FILE* fp = fileopen(xypath,"w");
		for (size_t i = 0; i < nhpixels; i++) {
			double hp = dh/2.0 + dh * i;
			fprintf(fp,"%10d %10zu %10.2lf %10.2lf %10.2lf %10.2lf %12.6lf %12.6lf %10.2lf\n",
				D.linenumber, i + 1, hp, imagefid[i], imagex[i], imagey[i], imagelon[i], imagelat[i], imageelevation[i]);
		}
		fclose(fp);

		double ulx = h0 - dh/2.0 + dh/2.0;
		double uly = v1 + dv/2.0; 
		double lrx = h1 + dh/2.0 + dh/2.0;
		double lry = v0 - dv/2.0;
		std::string xyext = outdir + "geometry\\" + basename() + ".extent.txt";
		fp = fileopen(xyext, "w");
		fprintf(fp, "%10d %10d %10d %10d %10d", D.linenumber, 0, 0, nhpixels, -nvpixels);
		fprintf(fp, " %10.2lf %10.2lf %10.2lf %10.2lf\n", ulx, uly, lrx, lry);
		fclose(fp);
	};

	void save_world_file()
	{
		//Line 1 : A: pixel size in the x - direction in map units / pixel
		//Line 2 : D : rotation about y - axis
		//Line 3 : B : rotation about x - axis
		//Line 4 : E : pixel size in the y - direction in map units, almost always negative[3]
		//Line 5 : C : x - coordinate of the center of the upper left pixel
		//Line 6 : F : y - coordinate of the center of the upper left pixel		
		std::ofstream ofs(worldfile());
		ofs << dh  << std::endl;
		ofs << 0.0 << std::endl;
		ofs << 0.0 << std::endl;
		ofs << -dv << std::endl;
		ofs << h0+dh/2.0  << std::endl;
		ofs << v1  << std::endl;
	}
	
	void savexml(const std::vector<double> longitude, const std::vector<double> latitude)
	{				
		makedirectorydeep(extractfiledirectory(xmlpath()));
		try
		{
			//Levels
			int nlevels = nlevels = levelCount(nhpixels, nvpixels, tilesize);
			while(nlevels == 0){
				//only to get around bug in the ribbon tiler java code
				tilesize /= 2;
				nlevels = levelCount(nhpixels, nvpixels, tilesize);				
			} 

			Element a, b;
			Document doc(xmlpath());
			std::string ver = "1.0";
			std::string enc = "UTF-8";
			std::string std = "yes";
			Declaration dec(ver, enc, std);
			doc.InsertEndChild(dec);			

			Element l("Layer");
			l.SetAttribute("version", "1");
			l.SetAttribute("layerType", "CurtainImageLayer");
			l.InsertEndChild(Element("DisplayName", basename()));
			l.InsertEndChild(Element("Legend", colorbarfile_nod()));
			
			bool local = true;
			std::string url = "http://www.ga.gov.au/apps/world-wind/tiles.jsp";
			if (local)  url   = "./";

			a = Element("Service");
			a.SetAttribute("serviceName", "DelegatorTileService");
			a.InsertEndChild(Element("URL", url));
			l.InsertEndChild(a);

			a = Element("Delegates");
			if(local) a.InsertEndChild(Element("Delegate", "LocalRequester"));
			a.InsertEndChild(Element("Delegate", "TransparentColorTransformer(255,255,255,0.2)"));
			a.InsertEndChild(Element("Delegate", strprint("ResizeTransformer(%d,%d)",tilesize,tilesize)));
			l.InsertEndChild(a);
						
			//Expects timestamps in the format “dd MM yyyy HH:mm:ss Z”			
			std::string tf = "%d %m %Y %H:%M:%S +11:00";
			std::string timestampstr = timestring(tf);
			l.InsertEndChild(Element("LastUpdate", timestampstr));
			
			std::string dn = basename();
			std::string dc = datacachename + "/" + basename();
			l.InsertEndChild(Element("DatasetName", dn));
			l.InsertEndChild(Element("DataCacheName", dc));
			
			//Image formats
			l.InsertEndChild(Element("ImageFormat", "image/" + image_format()));
			l.InsertEndChild(Element("FormatSuffix", "." + image_format()));
			
			a = Element("AvailableImageFormats");
			a.InsertEndChild(Element("ImageFormat", "image/" + image_format()));
			l.InsertEndChild(a);

			a = Element("NumLevels");
			a.SetAttribute("count", nlevels);
			a.SetAttribute("numEmpty", "0");
			l.InsertEndChild(a);
			
			//Tilesize			
			a = Element("TileSize");
			b = Element("Dimension");
			b.SetAttribute("width", strprint("%d",tilesize));
			b.SetAttribute("height", strprint("%d", tilesize));
			a.InsertEndChild(b);
			l.InsertEndChild(a);

			//Fullsize
			a = Element("FullSize");
			b = Element("Dimension");
			b.SetAttribute("width", nhpixels);
			b.SetAttribute("height", nvpixels);
			a.InsertEndChild(b);
			l.InsertEndChild(a);

			l.InsertEndChild(getxmlpathelement(longitude, latitude));

			l.InsertEndChild(Element("CurtainTop", v1));
			l.InsertEndChild(Element("CurtainBottom", v0));
			l.InsertEndChild(Element("FollowTerrain", "false"));

			l.InsertEndChild(Element("Subsegments",2));
			l.InsertEndChild(Element("UseTransparentTextures","true"));
			l.InsertEndChild(Element("ForceLevelZeroLoads","true"));
			l.InsertEndChild(Element("RetainLevelZeroTiles","image/dds"));
			l.InsertEndChild(Element("UseMipMaps","true"));
			l.InsertEndChild(Element("DetailHint",0.5));
			
			doc.InsertEndChild(l);
			doc.SaveFile();
		}
		catch (ticpp::Exception& ex)
		{
			std::cout << ex.what();
		}
	}		

	int levelCount(const int& width, const int& height, const int& tilesize)
	{
		//This is lifted direct from 
		//https://github.com/GeoscienceAustralia/ga-worldwind-suite/blob/master/Tiler/src/main/java/au/gov/ga/worldwind/tiler/ribbon/RibbonTiler.java
		double xCount = width / (double)tilesize;
		double yCount = height / (double)tilesize;
		int levels = 0;	
		while (4.0 * xCount * yCount >= 1)
		{
			levels++;
			xCount /= 2.0;
			yCount /= 2.0;
		}
		return levels;
	}
	
	Element getxmlpathelement(const std::vector<double>& longitude, const std::vector<double>& latitude)
	{
		Element p("Path");
		for (size_t i = 0; i < longitude.size(); i++){
			Element a("LatLon");
			a.SetAttribute("units", "degrees");
			a.SetAttribute("latitude", strprint("%10.6lf", latitude[i]));
			a.SetAttribute("longitude", strprint("%10.6lf", longitude[i]));
			p.InsertEndChild(a);
		}
		return p;
	}
	
	void saveribbontilerbatchcommand()
	{						
		std::ios_base::openmode mode = std::ios::app;
		if (sequence_number == 0) mode = std::ios::trunc;
		std::ofstream file(ribbontilerbatchfilepath().c_str(), mode);
		
		std::string s;
		if (sequence_number == 0) s += strprint("@echo off\n\n");
		
		s += strprint("call ribbon.bat");		
		s += strprint(" -tilesize %d",tilesize);
		s += strprint(" -noLayerDef");
		s += strprint(" -format %s", image_format().c_str());
		s += strprint(" -source %s", imagefile_nod().c_str());
		s += strprint(" -output %s", tilesetdir_nod().c_str());
		s += strprint("\n");				
		file << s.c_str();		
	}		
	
	static void appendpause(const std::string filename)
	{
		std::ios_base::openmode mode = std::ios::app;		
		std::ofstream file(filename, mode);		
		file << "\npause\n";
	}
};


void save_dataset_xml(const std::string xmlpath, 
	const std::string datasetname,
	const std::vector<std::string> names,
	const std::vector<std::string> urls
	)
{
	makedirectorydeep(extractfiledirectory(xmlpath));
	try
	{		
		Element a, b;
		Document doc(xmlpath);
		
		//std::string ver = "1.0";
		//std::string enc = "UTF-8";
		//std::string std = "yes";
		//Declaration dec(ver, enc, std);
		//doc.InsertEndChild(dec);
		
		Element dl("DatasetList");

		Element d("Dataset");
		d.SetAttribute("name", datasetname);
		//d.SetAttribute("info", "http://www.ga.gov.au/eftf");

		for (size_t i = 0; i < names.size(); i++){
			Element l("Layer");
			l.SetAttribute("name", names[i]);
			l.SetAttribute("url", urls[i]);
			//l.SetAttribute("icon", "icon.png");
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
	ULONG_PTR token;
	try{	
		token = cGDIplusHelper::start();

		if (argc >= 2){
			glog.logmsg("Executing %s %s\n", argv[0], argv[1]);
			glog.logmsg("Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
			glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
		}
		else{
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
		std::string dfnfile = ib.getstringvalue("DfnFile");
		std::string infiles = ib.getstringvalue("DataFiles");
		std::vector<std::string> filelist = cDirectoryAccess::getfilelist(infiles);
		cStopWatch stopwatch;
		size_t sequence_number = 0;
		std::string tilerbatchfile;
		std::vector<std::string> names;
		std::vector<std::string> urls;
		std::string datasetxml;
		std::string datasetname;

		for (size_t i = 0; i < filelist.size(); i++){
			std::printf("Processing file %s %3zu of %3zu\n", filelist[i].c_str(), i + 1, filelist.size());

			std::string datafile = filelist[i];			

			cCTLineData dummy(ib, dfnfile);
			cRange<size_t> lcol = dummy.getcolumns("line");

			cFileSplitter FS(datafile, 0, lcol.from);
			std::vector<std::string> L;
			while (FS.getnextgroup(L) > 0){				
				cCTLineData D(ib, dfnfile);
				D.load(L);
				std::printf("Line %d\n", D.linenumber);

				cCurtainImageSection C(D);
				C.set_sequence_number(sequence_number);				
				C.getoptions(sb);
				C.process();				
				if (sequence_number == 0) C.createcolorbar();
				sequence_number++;
				tilerbatchfile = C.ribbontilerbatchfilepath();
				
				names.push_back(C.basename());
				urls.push_back(C.xmlname());
				datasetxml  = C.datasetxmlpath();
				datasetname = C.getdatasetname();
			}
		}
				
		save_dataset_xml(datasetxml,datasetname,names,urls);
		cCurtainImageSection::appendpause(tilerbatchfile);		
		printf("Done ... \nElapsed time = %.3lf seconds\n", stopwatch.etimenow());
		cGDIplusHelper::stop(token);		
	}
	catch (ticpp::Exception& e){
		std::cout << e.what();
		cGDIplusHelper::stop(token);
	}
	catch (std::runtime_error& e){
		std::cout << e.what();
		cGDIplusHelper::stop(token);
	}	
	return 0;
}



