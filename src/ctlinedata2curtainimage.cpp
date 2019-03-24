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
#include <cstring>

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

#define VERSION "1.0"

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

class cCurtainImageSection{

private:

	const cCTLineData& D;	
	size_t sequence_number = 0;
	Bitmap* pBitmap;
	int nhpixels;
	int nvpixels;
	double hlength;
	double vlength;
	double h0, h1, dh;
	double v0, v1, dv;

	double gratdiv;

	std::string outdir;
	std::string prefix;
	std::string suffix;

	cColorMap cmap;
	cStretch stretch;
	std::vector<double> cbarticks;

	bool autozsectiontop;
	bool autozsectionbot;

	bool   spreadfade;
	double spreadfadelowclip;
	double spreadfadehighclip;
	double lowspreadfade;
	double highspreadfade;

	double geometrytolerance;
	int    tilesize;
	std::string datasetname;
	std::string datacachename;

public:
	
	

	cCurtainImageSection(const cCTLineData& _D) : D(_D)
	{

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
		
		gratdiv = b.getdoublevalue("ElevationGridDivision");

		geometrytolerance = b.getdoublevalue("GeometryTolerance");
		tilesize = b.getintvalue("TileSize");		
		datasetname   = b.getstringvalue("DatasetName");
		datacachename = b.getstringvalue("DataCacheName");
	}
	
	void process(){					
		calculateextents();
		createbitmap();
		generatesectiondata();
		//drawgraticule();		
		//fixtransparenttext();
		saveimage();
		savegeometry();		
		saveribbontilerbatchcommand();
		deletebitmap();
	}

	void calculateextents()
	{												
		hlength = D.linedistance.back();
		hlength = roundupnearest(hlength,dh);	
		h0 = 0.0;
		h1 = hlength;
		nhpixels  = 1 + (int)(hlength/dh);
		
		double zmin =  DBL_MAX;
		double zmax = -DBL_MAX;
		for(int si=0; si<D.nsamples; si++){
			if(D.z[si][D.nlayers] < zmin)zmin = D.z[si][D.nlayers];
			if (D.z[si][0]       > zmax)zmax = D.z[si][0];
		}

		if(autozsectionbot==true) v0 = zmin;
		if(autozsectiontop==true) v1 = zmax;		
		vlength = v1-v0;
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
	
	void createbitmap(){
		pBitmap = new Bitmap(nhpixels, nvpixels, PixelFormat32bppARGB);		
		for(int i=0; i<nhpixels; i++){
			for(int j=0; j<nvpixels; j++){				
				pBitmap->SetPixel(i,j,Color(255,255,255,255));				
			}
		}		
	}

	void deletebitmap(){
		delete pBitmap;
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

	void generatesectiondata(){

		Color BkgColor(255,128,128,128);
		Color AirColor(0,255,255,255);
		Color NullsColor(255,128,128,128);
		
		int vp0 = wv2iy(v0);
		int vp1 = wv2iy(v1);

		int mini = 0;
		for (int i = 0; i <= nhpixels; i++){
			double hp = h0 + (double)i*dh;

			while (mini < D.linedistance.size() - 1 && std::abs(D.linedistance[mini] - hp) > std::abs(D.linedistance[mini + 1] - hp)){
				mini++;
			}
			
			for(int j=vp1; j<=vp0; j++){				
				pBitmap->SetPixel(i,j,BkgColor);

				double zp = v1 - (double)j*dv;
				if(zp>D.e[mini]){
					pBitmap->SetPixel(i,j,AirColor);					
				}
				else{					
					for (int li = 0; li<D.nlayers; li++){
						if (zp < D.z[mini][li] && zp >= D.z[mini][li + 1]){
							double conductivity = D.c[mini][li];
							if(conductivity < 0.0){
								pBitmap->SetPixel(i,j,NullsColor);
							}
							else{							
								int ind = stretch.index(conductivity);
								Color clr(255,cmap.r[ind],cmap.g[ind],cmap.b[ind]);

								if(spreadfade){
									double conductivity_p10 = D.cp10[mini][li];
									double conductivity_p90 = D.cp90[mini][li];
									double spread = log10(conductivity_p90) - log10(conductivity_p10);
									double fade   = fadevalue(spread);
									fadecolor(clr,fade);
								}								
								pBitmap->SetPixel(i,j,clr);
							}
							break;
						}						
					}										
				}//layer loop
			}//v pixel loop
		}//h pixel loop		
	}

	void drawgraticule(){

		if (gratdiv == ud_double()) return;
		if (gratdiv <= 0.0) return;

		Graphics gr(pBitmap);
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

	void fixtransparenttext(){	
		int hp0 = wh2ix(h0);
		int hp1 = wh2ix(h1);
		int vp1 = wv2iy(v0);
		int vp2 = wv2iy(v1);
		for(int i=hp0; i<=hp1; i++){		
			for(int j=vp2; j<=vp1; j++){
				Color c;
				pBitmap->GetPixel(i,j,&c);				
				pBitmap->SetPixel(i, j, Color(255, c.GetR(), c.GetG(), c.GetB()));
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

	std::string jpegfile_nod(){
		std::string s = "jpeg\\" + basename() + ".jpg";
		return s;
	}

	std::string jpegfile(){
		std::string s = outdir + jpegfile_nod();
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
		std::string s = tilesetdir() + "colourbar.jpg";
		return s;
	}

	std::string colorbarfile_nod(){
		std::string s = "colourbar.jpg";
		return s;
	}
	
	void saveimage(){
		cGDIplusHelper::saveimage(pBitmap, jpegfile());
	}
	
	void savegeometry(){
		
		//std::string txtpath = outdir + basename() + ".txt";
		//FILE* fp = fileopen(txtpath, "w");
		//fprintf(fp, "top %lf\n", zsectiontop);
		//fprintf(fp, "bottom %lf\n", zsectionbot);
		//fprintf(fp, "hlength %lf\n", hlength);
		//fprintf(fp, "vlength %lf\n", vlength);
		//fprintf(fp, "nhpixels %d\n", nhpixels);
		//fprintf(fp, "nvpixels %d\n", nvpixels);				
		//fclose(fp);

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
		
		std::vector<double> x(plout.size());
		std::vector<double> y(plout.size());
		std::vector<double> longitude(plout.size());
		std::vector<double> latitude(plout.size());
		for (size_t i = 0; i < plout.size(); i++){
			x[i] = plout[i].first;
			y[i] = plout[i].second;
		}
		
		int inepsgcode = cCRS::epsgcode(D.inputdatumprojection);
		if (inepsgcode < 0){
			std::string msg = strprint("Invalid DatumProjection %s was specified\n", D.inputdatumprojection.c_str()) + _SRC_;
			throw(std::runtime_error(msg));
		}

		int outepsgcode = cCRS::epsgcode("WGS84|GEODETIC");
		transform(inepsgcode, x, y, outepsgcode, longitude, latitude);
		std::string pfpathll = outdir + basename() + ".points_filtered_geodetic.dat";
		//savepoints(longitude, latitude, pfpathll);
		savexml(longitude, latitude);
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
			fprintf(fp, "%10lf %10lf\n", x[i], y[i]);
		}
		fclose(fp);
	};
	
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
			l.InsertEndChild(Element("ImageFormat", "image/jpg"));			
			l.InsertEndChild(Element("FormatSuffix", ".jpg"));
			
			a = Element("AvailableImageFormats");
			a.InsertEndChild(Element("ImageFormat", "image/jpg"));			
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
		//s += strprint(" -copySource");
		s += strprint(" -noLayerDef");
		s += strprint(" -source %s", jpegfile_nod().c_str());
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
			glog.logmsg("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
			glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
		}
		else{
			glog.logmsg("Executing %s\n", argv[0]);
			glog.logmsg("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
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
			//std::string dfnfile = extractfiledirectory(datafile);
			//dfnfile += "inversion.output.dfn";

			cCTLineData dummy(ib, dfnfile);
			cRange<int> lcol = dummy.getcolumns("line");

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



