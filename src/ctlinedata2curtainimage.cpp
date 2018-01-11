/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <windows.h>
#include <gdiplus.h>
#include <GdiPlusImageCodec.h>
using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")

#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cstring>
#include "general_types.h"
#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "geometry3d.h"
#include "colormap.h"
#include "RamerDouglasPeucker.h"
#include "gdal_utils.h"

#include "ticpp.h"
using namespace ticpp;

using namespace RDP;


#define VERSION "1.0"

int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes

	ImageCodecInfo* pImageCodecInfo = NULL;

	GetImageEncodersSize(&num, &size);
	if (size == 0)
		return -1;  // Failure

	pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
	if (pImageCodecInfo == NULL)
		return -1;  // Failure

	GetImageEncoders(num, size, pImageCodecInfo);

	for (UINT j = 0; j < num; ++j)
	{
		if (wcscmp(pImageCodecInfo[j].MimeType, format) == 0)
		{
			*pClsid = pImageCodecInfo[j].Clsid;
			free(pImageCodecInfo);
			return j;  // Success
		}
	}

	free(pImageCodecInfo);
	return -1;  // Failure
}

class cCurtainImageSection{

private:

	Bitmap* pBitmap;
	cColorMap cmap;
	bool drawcbar;
	std::vector<double> cbarticks;
	int nhpixels;
	int nvpixels;
	std::vector<double> linedistance;
	double hlength;
	double vlength;
	double h0,h1,dh;
	double v0,v1,dv;
	//double x0,x1,x2,dx;
	//double y0,y1,y2,dy;
	double gratdiv;
	double geometrytolerance;

	std::string outdir;
	std::string prefix;
	std::string suffix;

	int nlayers;
	int nsamples;
	int linenumber;


	double LowClip;
	double HighClip;
	bool   Log10Stretch;

	bool SaveJPG;
	bool SavePNG;
	bool SaveEMF;

	double elevation_median;
	bool autozsectiontop;
	bool autozsectionbot;
	double zsectiontop;
	double zsectionbot;

public:
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> e;
	std::vector<std::vector<double>> c;
	std::vector<std::vector<double>> z;
	std::vector<std::vector<double>> cp10;
	std::vector<std::vector<double>> cp90;
	bool spreadfade;	
	double spreadfadelowclip;
	double spreadfadehighclip;
	double lowspreadfade;
	double highspreadfade;
	
	void getoptions(cBlock b){
		
		spreadfade = b.getboolvalue("SpreadFade");
		if(spreadfade){
			spreadfadelowclip   = b.getdoublevalue("Log10SpreadLowClip");
			spreadfadehighclip  = b.getdoublevalue("Log10SpreadHighClip");
			lowspreadfade       = b.getdoublevalue("LowSpreadFade");
			highspreadfade      = b.getdoublevalue("HighSpreadFade");	
		}


		double v = b.getdoublevalue("ElevationTop");
		if(isundefined(v)){
			autozsectiontop = true;
			zsectiontop = v;
		}
		else{
			autozsectiontop = false;
			zsectiontop = v;
		}

		v = b.getdoublevalue("ElevationBottom");
		if (isundefined(v)){
			autozsectionbot = true;
			zsectionbot = v;
		}
		else{
			autozsectionbot = false;
			zsectionbot = v;
		}

		outdir = b.getstringvalue("OutDir");
		prefix = b.getstringvalue("Prefix");
		suffix = b.getstringvalue("Suffix");

		dh = b.getdoublevalue("HorizontalResolution");
		double vxg   = b.getdoublevalue("VerticalExaggeration");
		dv = dh/vxg;

		std::string cmapname = b.getstringvalue("ColourMap");
		cmap.set(cmapname);

		drawcbar  = b.getboolvalue("ColourBar");
		if(drawcbar){
			cbarticks = b.getdoublevector("ColourBarTicks");
		}

		LowClip  = b.getdoublevalue("LowClip");
		HighClip = b.getdoublevalue("HighClip");
		Log10Stretch  = b.getboolvalue("Log10Stretch");

		gratdiv = b.getdoublevalue("ElevationGridDivision");

		geometrytolerance = b.getdoublevalue("GeometryTolerance");

		SaveJPG = b.getboolvalue("SaveJPG");
		SavePNG = b.getboolvalue("SavePNG");
		SaveEMF = b.getboolvalue("SaveEMF");
	}

	void readdatafile(const cBlock& input, const std::string filename)
	{
		int subsample = input.getintvalue("Subsample");
		if (isundefined(subsample))subsample = 1;

		std::string lstr = input.getstringvalue("Line");
		std::string xstr = input.getstringvalue("Easting");
		std::string ystr = input.getstringvalue("Northing");
		std::string estr = input.getstringvalue("Elevation");
		
		bool isresistivity = false;
		std::string crstr = input.getstringvalue("Conductivity");		
		if (isundefined(crstr)){
			crstr = input.getstringvalue("Resistivity");
			isresistivity = true;
		}

		std::string cp10str = input.getstringvalue("Conductivity_p10");
		std::string cp90str = input.getstringvalue("Conductivity_p90");								

		double cscale=1.0;
		std::string cunits = input.getstringvalue("InputConductivityUnits");
		if(isundefined(cunits)){
			cscale = 1.0;
		}
		else if(strcasecmp(cunits,"S/m") == 0){
			cscale = 1.0;
		}
		else if(strcasecmp(cunits,"mS/m") == 0){
			cscale = 0.001;
		}		
		else{
			message("Unknown InputConductivityUnits %s\n",cunits.c_str());			
		}

		int lcol, xcol, ycol, ecol;
		int crcol1, crcol2, tcol1, tcol2;
		int cp10col1, cp10col2;
		int cp90col1, cp90col2;
		sscanf(lstr.c_str(), "Column %d", &lcol); lcol--;
		sscanf(xstr.c_str(), "Column %d", &xcol); xcol--;
		sscanf(ystr.c_str(), "Column %d", &ycol); ycol--;
		sscanf(estr.c_str(), "Column %d", &ecol); ecol--;

		sscanf(crstr.c_str(), "Column %d-%d", &crcol1, &crcol2); crcol1--; crcol2--;
		nlayers = crcol2 - crcol1 + 1;

		if(spreadfade){
			sscanf(cp10str.c_str(), "Column %d-%d", &cp10col1, &cp10col2); cp10col1--; cp10col2--;
			sscanf(cp90str.c_str(), "Column %d-%d", &cp90col1, &cp90col2); cp90col1--; cp90col2--;
		}

		bool isconstantthickness = false;
		std::vector<double> constantthickness;
		std::string tstr = input.getstringvalue("Thickness");
		if(sscanf(tstr.c_str(), "Column %d-%d", &tcol1, &tcol2) == 2){
			tcol1--; tcol2--;	
			isconstantthickness = false;
		}
		else{			
			constantthickness = input.getdoublevector("Thickness");
			tcol1 = 0; tcol2 = 0;
			isconstantthickness = true;
			if(constantthickness.size() == 0){
				errormessage("Thickness not set\n");
			}
			else if(constantthickness.size() > 1 && constantthickness.size() < nlayers - 1){
				errormessage("Thickness not set correctly\n");
			}
			else if(constantthickness.size() == 1){
				constantthickness = std::vector<double>(nlayers-1, constantthickness[0]);
			}
			else{
				//all good
			}

		}

		FILE* fp = fileopen(filename,"r");
		std::string str;
		std::vector<std::vector<double>> M;

		int k=0;
		while (filegetline(fp, str)){			
			if(k%subsample == 0){
				M.push_back(getdoublevector(str.c_str(), " "));
			}	
			k++;
		}
		fclose(fp);


		nsamples   = (int)M.size();
		linenumber = (int)M[0][lcol];		

		x.resize(nsamples);
		y.resize(nsamples);
		e.resize(nsamples);
		z.resize(nsamples);
		c.resize(nsamples);
		if(spreadfade){
			cp10.resize(nsamples);
			cp90.resize(nsamples);		
		}
		for (int si = 0; si < nsamples; si++){
			x[si] = M[si][xcol];
			y[si] = M[si][ycol];
			e[si] = M[si][ecol];

			c[si].resize(nlayers);
			for (int li = 0; li < nlayers; li++){
				c[si][li] = M[si][crcol1 + li];				
				if(c[si][li] > 0.0){
					if(isresistivity)c[si][li] = 1.0 / c[si][li];
					c[si][li] *= cscale;
				}
			}

			if(spreadfade){
				cp10[si].resize(nlayers);
				for (int li = 0; li < nlayers; li++){
					cp10[si][li] = M[si][cp10col1 + li];
					if(cp10[si][li] > 0.0){
						cp10[si][li] *= cscale;
					}
				}

				cp90[si].resize(nlayers);
				for (int li = 0; li < nlayers; li++){
					cp90[si][li] = M[si][cp90col1 + li];
					if(cp90[si][li] > 0.0){
						cp90[si][li] *= cscale;
					}
				}
			}

			z[si].resize(nlayers+1);			
			z[si][0] = e[si];
			for (int li = 0; li < nlayers; li++){				

				double t;
				if (li < nlayers - 1){
					if (isconstantthickness == true){
						t = constantthickness[li];
					}
					else{
						t = M[si][tcol1 + li];
					}
				}		
				else{
					if (isconstantthickness == true){
						t = constantthickness[li-1];
					}
					else{
						t = M[si][tcol1 + li - 1];
					}					
				}	

				z[si][li + 1] = z[si][li] - t;

			}
		}		
	}

	void process(){		
		//savegeometry();
		//return;

		calculateextents();				
		createbitmap();
		generatesectiondata();
		//drawgraticule();
		//drawcolorbar();		
		//fixtransparenttext();
		saveimage();
		savegeometry();
		deletebitmap();
	}

	void calculateextents()
	{				
		//bestfitlineendpoints(x,y,x1,y1,x2,y2);
		//double angle = atan2(y2-y1,x2-x1);
		//if(angle > PIONTWO || angle <= -PIONTWO){
		//	std::swap(x1,x2);
		//	std::swap(y1,y2);			
		//}

		
		linedistance.resize(x.size());
		linedistance[0] = 0.0;
		for (size_t i = 1; i < x.size(); i++){
			linedistance[i] = linedistance[i - 1] + distance(x[i - 1], y[i - 1], x[i], y[i]);
		}
		double rawlength = linedistance.back();
		hlength = roundupnearest(rawlength,dh);

		//x2 = x1 + (x2-x1)*hlength/rawlength;
		//y2 = y1 + (y2-y1)*hlength/rawlength;

		int hmarginpixels = getmarginpixels();
		h0 = 0.0;		
		h1 = hlength;

		nhpixels  = 1 + (int)((h1-h0)/dh);
		//dx = dh*(x2-x1)/(h2-h1);
		//dy = dh*(y2-y1)/(h2-h1);		
		//x0 = x1 - (x2-x1)*h1/hlength;
		//y0 = y1 - (y2-y1)*h1/hlength;

		double zmin =  DBL_MAX;
		double zmax = -DBL_MAX;
		for(int si=0; si<nsamples; si++){
			if(z[si][nlayers] < zmin)zmin = z[si][nlayers];
			if(z[si][0]       > zmax)zmax = z[si][0];
		}

		if(autozsectionbot==false){
			zmin = zsectionbot;
		}
		else{
			zsectionbot = zmin;;
		}

		if(autozsectiontop==false){
			zmax = zsectiontop;
		}
		else{
			zsectiontop = zmax;
		}

		elevation_median = median(&e[0],e.size());

		double rawheight = zmax-zmin;
		vlength = roundupnearest(rawheight,dv);		
		v0 = zmin;		
		v1 = vlength;
		nvpixels = 1 + (int)(vlength/dv);		
	}

	int wh2lx(double wh){
		int lx = (int)roundnearest((wh-h0)/dh,1.0);
		return lx;
	}

	int wv2ly(double wv){		
		int ly = (int)roundnearest((v1-wv)/dv,1.0);
		return ly;
	}

	int getmarginpixels(){
		return 0;

		//int graticulemargin = 65;
		//int colourbarmargin = 95;

		//int margin = graticulemargin;
		//if(drawcbar) margin += colourbarmargin;
		//return margin;
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

		//int hp1 = wh2lx(h1);
		//int hp2 = wh2lx(h2);
		int vp1 = wv2ly(v0);
		int vp2 = wv2ly(v1);

		int mini = 0;
		for (int i = 0; i <= nhpixels; i++){
			double hp = h0 + (double)i*dh;

			while(mini < linedistance.size()-1 && std::abs(linedistance[mini] - hp) > std::abs(linedistance[mini + 1] - hp)){
				mini++;
			}
			
			
			/*
			double mind=DBL_MAX;
			int mini;
			for(int si=0; si<nsamples; si++){
				double d = distance(0.0,0.0,xp-x[si],yp-y[si]);
				if(d<mind){
					mini = si;
					mind = d;
				}
			}
			*/


			for(int j=vp2; j<=vp1; j++){				
				pBitmap->SetPixel(i,j,BkgColor);

				double zp = v1 - (double)j*dv;
				if(zp>e[mini]){
					pBitmap->SetPixel(i,j,AirColor);					
				}
				else{					
					for(int li=0; li<nlayers; li++){
						if(zp < z[mini][li] && zp >= z[mini][li+1]){
							double conductivity = c[mini][li];							
							if(conductivity < 0.0){
								pBitmap->SetPixel(i,j,NullsColor);
							}
							else{							
								int ind;
								if(Log10Stretch) ind = log10stretch(conductivity,LowClip,HighClip);							
								else ind = linearstretch(conductivity,LowClip,HighClip);

								Color clr(255,cmap.r[ind],cmap.g[ind],cmap.b[ind]);

								if(spreadfade){
									double conductivity_p10 = cp10[mini][li];
									double conductivity_p90 = cp90[mini][li];
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

	REAL getfontsize(){
		//REAL fontsize = (REAL)abs((wv2ly(0.33*gratdiv) - wv2ly(0)));
		REAL fontsize = 8.5;
		return fontsize;
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
		Font font(L"Arial", getfontsize(), FontStyleRegular, UnitPoint);				

		StringFormat textformat;
		textformat.SetAlignment(StringAlignmentFar);
		textformat.SetLineAlignment(StringAlignmentCenter);
		SizeF layoutsize(32767,32767);
		SizeF textsize;				
		gr.MeasureString(L"-0000 m",-1,&font,layoutsize,&textformat,&textsize);

		for(int zg=zg1; zg<=zg2; zg+=(int)gratdiv){
			gr.DrawLine(&blackpen,wh2lx(h0),wv2ly(zg),wh2lx(h1),wv2ly(zg));

			wchar_t s[20];				
			swprintf(s,20,L"%5d m",zg);						
			gr.MeasureString(s,-1,&font,layoutsize,&textformat,&textsize);

			int tx = wh2lx(h1);
			int ty = wv2ly(zg);
			if(ty > nvpixels)continue;
			if(ty-textsize.Height < 0)continue;	
			PointF txpos((REAL)tx,(REAL)ty);
			gr.DrawString(s, -1, &font, txpos, &textformat, &blackbrush);
		}		    				
	}

	void drawcolorbar(){

		if(drawcbar==false)return;

		Pen blackpen(Color::Black, 0);
		SolidBrush blackbrush(Color::Black);

		Graphics gr(pBitmap);

		TextRenderingHint hint = gr.GetTextRenderingHint();		
		gr.SetTextRenderingHint(TextRenderingHintAntiAlias);

		int ph1 = wv2ly(v0);
		int ph2 = wv2ly(v1);
		ph1 = 25; ph2 = 40;

		int pvtop = wv2ly(v0+(v1-v0)*0.95);
		int pvbot = wv2ly(v0+(v1-v0)*0.05);
		for (int j = pvtop; j <= pvbot; j++){			
			int ind = 255*(j - pvbot) / (pvtop - pvbot);
			for (int i = ph1; i <= ph2; i++){
				pBitmap->SetPixel(i, j, Color(255, cmap.r[ind], cmap.g[ind], cmap.b[ind]));
			}
		}
		gr.DrawRectangle(&blackpen, ph1, pvtop, ph2-ph1, pvbot-pvtop);

		Font font(L"Arial", getfontsize(), FontStyleRegular, UnitPoint);
		StringFormat textformat;		
		textformat.SetAlignment(StringAlignmentNear);
		textformat.SetLineAlignment(StringAlignmentCenter);
		SizeF layoutsize(32767, 32767);
		SizeF textsize;
		gr.MeasureString(L"0.0001 m", -1, &font, layoutsize, &textformat, &textsize);


		for (int i = 0; i<cbarticks.size(); i++){			
			if(cbarticks[i] < LowClip)continue;
			if(cbarticks[i] > HighClip)continue;

			int tickv;
			if (Log10Stretch){
				int ind = log10stretch(cbarticks[i], LowClip, HighClip);
				tickv = pvbot + ind*(pvtop - pvbot) / 255;
			}
			else{
				int ind = linearstretch(cbarticks[i], LowClip, HighClip);
				tickv = pvbot + ind*(pvtop - pvbot) / 255;
			}

			wchar_t s[20];
			swprintf(s, 20,L"%5.3lf", cbarticks[i]);

			PointF p;
			p.X = (REAL)ph2+3;
			p.Y = (REAL)tickv;
			gr.DrawString(s, -1, &font, p, &textformat, &blackbrush);
			gr.DrawLine(&blackpen,ph2-2,tickv,ph2+2,tickv);
		}		

		PointF p;
		p.X = (REAL)ph1-2;
		p.Y = (REAL)(pvtop + pvbot) / 2;

		textformat.SetAlignment(StringAlignmentCenter);
		textformat.SetLineAlignment(StringAlignmentFar);
		gr.TranslateTransform(p.X,p.Y);
		gr.RotateTransform(-90);		
		gr.DrawString(L"Conductivity (S/m)", -1, &font, PointF(0,0), &textformat, &blackbrush);

	}

	void fixtransparenttext(){	
		int hp0 = wh2lx(h0);
		int hp1 = wh2lx(h1);
		int vp1 = wv2ly(v0);
		int vp2 = wv2ly(v1);
		for(int i=hp0; i<=hp1; i++){		
			for(int j=vp2; j<=vp1; j++){
				Color c;
				pBitmap->GetPixel(i,j,&c);				
				pBitmap->SetPixel(i, j, Color(255, c.GetR(), c.GetG(), c.GetB()));
			}
		}		
	}

	std::string basename(){
		std::string bn = prefix + strprint("%d", linenumber) + suffix;
		return bn;
	}

	void saveimage(){
		std::string bname = basename();
		makedirectorydeep(outdir);

		if(SavePNG){		
			std::string imagepath = outdir + bname + ".png";						
			wchar_t wcimagepath[200];
			size_t len;			
			mbstowcs_s(&len, wcimagepath, 200, imagepath.c_str(), _TRUNCATE);

			CLSID pngClsid;
			int ret = GetEncoderClsid(L"image/png", &pngClsid);			
			Status result = pBitmap->Save(wcimagepath, &pngClsid, NULL);			
		}

		if(SaveJPG){			
			std::string imagepath = outdir + bname + ".jpg";
			wchar_t wcimagepath[500];
			size_t len;
			mbstowcs_s(&len, wcimagepath, 500, imagepath.c_str(), _TRUNCATE);

			CLSID jpgClsid;
			GetEncoderClsid(L"image/jpeg", &jpgClsid);
			Status result = pBitmap->Save(wcimagepath, &jpgClsid, NULL);			
		}

	}

	void savegeometry(){
		
		std::string txtpath = outdir + basename() + ".txt";
		FILE* fp = fileopen(txtpath, "w");
		fprintf(fp, "top %lf\n", zsectiontop);
		fprintf(fp, "bottom %lf\n", zsectionbot);
		fprintf(fp, "hlength %lf\n", hlength);
		fprintf(fp, "vlength %lf\n", vlength);
		fprintf(fp, "nhpixels %d\n", nhpixels);
		fprintf(fp, "nvpixels %d\n", nvpixels);				
		fclose(fp);


		std::vector<RDP::Point> pl(x.size());
		std::vector<RDP::Point> plout;
		for (size_t i = 0; i < x.size(); i++){
			pl[i]=RDP::Point(x[i],y[i]);
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

		int inepsgcode = 28353;//GDA94,MGA53		
		int outepsgcode = 4283;//GDA94,Geodetic
		
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

	void savepoints(
		const std::vector<double>& x,
		const std::vector<double>& y,
		const std::string& filename)
	{
		FILE* fp = fileopen(filename, "w");
		for (size_t i = 0; i < x.size(); i++){
			fprintf(fp, "%10lf %10lf\n", x[i], y[i]);
		}
		fclose(fp);
	};
	
	void savexml(const std::vector<double> longitude, const std::vector<double> latitude)
	{
		std::string xmlpath = outdir + basename() + ".xml";
		try
		{
			Element a, b;
			Document doc(xmlpath);
			std::string ver = "1.0";
			std::string enc = "UTF-8";
			std::string std = "yes";
			Declaration dec(ver, enc, std);
			doc.InsertEndChild(dec);
			doc.InsertEndChild(Element("CreatedBy", "ctlinedata2curtainimage"));

			Element l("Layer");
			l.SetAttribute("version", "1");
			l.SetAttribute("layerType", "CurtainImageLayer");
			l.InsertEndChild(Element("DisplayName", basename()));
			
			bool local = true;
			std::string url = "http://www.ga.gov.au/apps/world-wind/tiles.jsp";
			if (local)  url   = ".\\";

			a = Element("Service");
			a.SetAttribute("serviceName", "DelegatorTileService");
			a.InsertEndChild(Element("URL", url));
			l.InsertEndChild(a);

			a = Element("Delegates");
			if(local) a.InsertEndChild(Element("Delegate", "LocalRequester"));
			a.InsertEndChild(Element("Delegate", "TransparentColorTransformer(255, 255, 255, 0.2)"));
			a.InsertEndChild(Element("Delegate", "ResizeTransformer(512, 512)"));
			l.InsertEndChild(a);
			l.InsertEndChild(Element("LastUpdate", timestamp()));
			
			std::string datasetname  = "TestAEM/AEM_Lines/"+basename();
			std::string datasetcache = "cache/TestAEM/AEM_Lines/"+basename();
			l.InsertEndChild(Element("DatasetName", datasetname));
			l.InsertEndChild(Element("DataCacheName", datasetcache));
			
			a = Element("AvailableImageFormats");
			a.InsertEndChild(Element("ImageFormat", "image/jpg"));
			l.InsertEndChild(a);
			l.InsertEndChild(Element("FormatSuffix", ".jpg"));

			a = Element("NumLevels");
			a.SetAttribute("count", "5");
			a.SetAttribute("numEmpty", "0");
			l.InsertEndChild(a);
			
			//Tilesize
			a = Element("TileSize");
			b = Element("Dimension");
			b.SetAttribute("width", "512");
			b.SetAttribute("height", "512");
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

			l.InsertEndChild(Element("CurtainTop", zsectiontop));
			l.InsertEndChild(Element("CurtainBottom", zsectionbot));
			l.InsertEndChild(Element("FollowTerrain", "false"));

			l.InsertEndChild(Element("Subsegments",2));
			l.InsertEndChild(Element("UseTransparentTextures",true));
			l.InsertEndChild(Element("ForceLevelZeroLoads",true));
			l.InsertEndChild(Element("RetainLevelZeroTiles","image/dds"));
			l.InsertEndChild(Element("UseMipMaps",true));
			l.InsertEndChild(Element("DetailHint",0.5));
			
			doc.InsertEndChild(l);
			doc.SaveFile();
		}
		catch (ticpp::Exception& ex)
		{
			std::cout << ex.what();
		}
	}

	Element getxmlpathelement(const std::vector<double>& longitude, const std::vector<double>& latitude)
	{
		Element p("Path");
		for (size_t i = 0; i < longitude.size(); i++){
			Element a("LatLon");
			a.SetAttribute("units", "degrees");
			a.SetAttribute("latitude", strprint("%10.6lf", longitude[i]));
			a.SetAttribute("longitude", strprint("%10.6lf", latitude[i]));
			p.InsertEndChild(a);
		}
		return p;
	}
};

int main(int argc, char** argv)
{	
	GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

	if (argc >= 2){
		message("Executing %s %s\n", argv[0], argv[1]);
		message("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		message("Working directory %s\n", getcurrentdirectory().c_str());
	}
	else{
		message("Executing %s\n", argv[0]);
		message("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		message("Working directory %s\n", getcurrentdirectory().c_str());
		message("Error: Not enough input arguments\n");
		message("Usage: %s controlfilename\n",argv[0]);		
		return 0;
	}

	cBlock b(argv[1]);
	cBlock inputblock   = b.findblock("Input");	
	cBlock sectionblock = b.findblock("Section");	
	std::string infiles = inputblock.getstringvalue("DataFiles");

	std::vector<std::string> filelist =  cDirectoryAccess::getfilelist(infiles);
	double t1 = gettime();	
	for (size_t i = 0; i < filelist.size(); i++){	
		std::printf("Processing file %s %3zu of %3zu\n", filelist[i].c_str(),i+1,filelist.size());
		cCurtainImageSection S;
		S.getoptions(sectionblock);		
		S.readdatafile(inputblock,filelist[i].c_str());		
		S.process();		
	}
	double t2 = gettime();
	printf("Done ... Elapsed time = %.2lf seconds\n", t2 - t1);
	
	GdiplusShutdown(gdiplusToken);
	return 0;
}



