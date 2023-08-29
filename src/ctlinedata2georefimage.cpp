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

#include "gaaem_version.h"
#include "undefinedvalues.h"
#include "general_types.h"
#include "general_utils.h"
#include "colormap.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "geometry3d.h"
#include "stretch.h"
#include "file_formats.h"
#include "filesplitter.h"
#include "ctlinedata.h"
#include "gdiplus_utils.h"
#include "stopwatch.h"


class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

class cGeorefSection {

private:

	const cCTLineData& D;

	Bitmap* pBitmap;
	cColorMap m_cmap;
	bool m_drawcbar;
	std::vector<double> m_cbarticks;
	int nhpixels;
	int nvpixels;
	double hlength;
	double vlength;
	double h0, h1, h2, dh;
	double v0, v1, v2, dv;
	double x0, x1, x2, dx;
	double y0, y1, y2, dy;
	double gratdiv;

	std::string outdir;
	std::string prefix;
	std::string suffix;

	double LowClip;
	double HighClip;
	bool   Log10Stretch;

	bool SaveJPG;
	bool SavePNG;

	Color BkgColor;
	Color AirColor;
	Color NullsColor;

	double elevation_median;
	bool autozsectiontop;
	bool autozsectionbot;
	double zsectiontop;
	double zsectionbot;
	
	//std::string imageformat = "jpg";
	//std::string worldfileextension = "jgw";

public:

	cGeorefSection(const cCTLineData& _D) : D(_D)
	{
		BkgColor = Color(0, 255, 255, 255);
		AirColor = Color(0, 255, 255, 255);
		NullsColor = Color(255, 128, 128, 128);
	};	
	bool spreadfade;
	double spreadfadelowclip;
	double spreadfadehighclip;
	double lowspreadfade;
	double highspreadfade;

	void getoptions(cBlock b) {

		spreadfade = b.getboolvalue("SpreadFade");
		if (spreadfade) {
			spreadfadelowclip = b.getdoublevalue("Log10SpreadLowClip");
			spreadfadehighclip = b.getdoublevalue("Log10SpreadHighClip");
			lowspreadfade = b.getdoublevalue("LowSpreadFade");
			highspreadfade = b.getdoublevalue("HighSpreadFade");
		}


		double v = b.getdoublevalue("ElevationTop");
		if (!isdefined(v)) {
			autozsectiontop = true;
			zsectiontop = v;
		}
		else {
			autozsectiontop = false;
			zsectiontop = v;
		}

		v = b.getdoublevalue("ElevationBottom");
		if (!isdefined(v)) {
			autozsectionbot = true;
			zsectionbot = v;
		}
		else {
			autozsectionbot = false;
			zsectionbot = v;
		}

		outdir = b.getstringvalue("OutDir");
		prefix = b.getstringvalue("Prefix");
		suffix = b.getstringvalue("Suffix");

		dh = b.getdoublevalue("HorizontalResolution");
		double vxg = b.getdoublevalue("VerticalExaggeration");
		dv = dh / vxg;

		std::string cmapname = b.getstringvalue("ColourMap");
		m_cmap = cColorMap(cmapname);

		m_drawcbar = b.getboolvalue("ColourBar");
		if (m_drawcbar) {
			m_cbarticks = b.getdoublevector("ColourBarTicks");
		}

		LowClip = b.getdoublevalue("LowClip");
		HighClip = b.getdoublevalue("HighClip");
		Log10Stretch = b.getboolvalue("Log10Stretch");

		gratdiv = b.getdoublevalue("ElevationGridDivision");

		SaveJPG = b.getboolvalue("SaveJPG");
		SavePNG = b.getboolvalue("SavePNG");

		BkgColor = Color(0, 255, 255, 255);
		AirColor = Color(0, 255, 255, 255);
		NullsColor = Color(255, 128, 128, 128);

		std::vector<int> bkgclr = b.getintvector("BackgroundColour");
		std::vector<int> airclr = b.getintvector("AirColour");
		std::vector<int> nullsclr = b.getintvector("NullsColour");

		if (bkgclr.size() == 4) BkgColor = Color(bkgclr[0], bkgclr[1], bkgclr[2], bkgclr[3]);
		if (airclr.size() == 4) AirColor = Color(airclr[0], airclr[1], airclr[2], airclr[3]);
		if (nullsclr.size() == 4) NullsColor = Color(nullsclr[0], nullsclr[1], nullsclr[2], nullsclr[3]);

	}

	/*
	void readdatafile(const cBlock& input, const std::string filename)
	{
		int subsample = input.getintvalue("Subsample");
		if (!isdefined(subsample))subsample = 1;

		std::string lstr = input.getstringvalue("Line");
		std::string xstr = input.getstringvalue("Easting");
		std::string ystr = input.getstringvalue("Northing");
		std::string estr = input.getstringvalue("Elevation");

		bool isresistivity = false;
		std::string crstr = input.getstringvalue("Conductivity");
		if (!isdefined(crstr)) {
			crstr = input.getstringvalue("Resistivity");
			isresistivity = true;
		}

		std::string cp10str = input.getstringvalue("Conductivity_p10");
		std::string cp90str = input.getstringvalue("Conductivity_p90");

		double cscale = 1.0;
		std::string cunits = input.getstringvalue("InputConductivityUnits");
		if (!isdefined(cunits)) {
			cscale = 1.0;
		}
		else if (strcasecmp(cunits, "S/m") == 0) {
			cscale = 1.0;
		}
		else if (strcasecmp(cunits, "mS/m") == 0) {
			cscale = 0.001;
		}
		else {
			glog.logmsg("Unknown InputConductivityUnits %s\n", cunits.c_str());
		}

		int lcol, xcol, ycol, ecol;
		int crcol1, crcol2, tcol1, tcol2;
		int cp10col1, cp10col2;
		int cp90col1, cp90col2;
		std::sscanf(lstr.c_str(), "Column %d", &lcol); lcol--;
		std::sscanf(xstr.c_str(), "Column %d", &xcol); xcol--;
		std::sscanf(ystr.c_str(), "Column %d", &ycol); ycol--;
		std::sscanf(estr.c_str(), "Column %d", &ecol); ecol--;

		std::sscanf(crstr.c_str(), "Column %d-%d", &crcol1, &crcol2); crcol1--; crcol2--;
		nlayers = crcol2 - crcol1 + 1;

		if (spreadfade) {
			std::sscanf(cp10str.c_str(), "Column %d-%d", &cp10col1, &cp10col2); cp10col1--; cp10col2--;
			std::sscanf(cp90str.c_str(), "Column %d-%d", &cp90col1, &cp90col2); cp90col1--; cp90col2--;
		}

		bool isconstantthickness = false;
		std::vector<double> constantthickness;
		std::string tstr = input.getstringvalue("Thickness");
		if (std::sscanf(tstr.c_str(), "Column %d-%d", &tcol1, &tcol2) == 2) {
			tcol1--; tcol2--;
			isconstantthickness = false;
		}
		else {
			constantthickness = input.getdoublevector("Thickness");
			tcol1 = 0; tcol2 = 0;
			isconstantthickness = true;
			if (constantthickness.size() == 0) {
				glog.errormsg(_SRC_, "Thickness not set\n");
			}
			else if (constantthickness.size() > 1 && constantthickness.size() < nlayers - 1) {
				glog.errormsg(_SRC_, "Thickness not set correctly\n");
			}
			else if (constantthickness.size() == 1) {
				constantthickness = std::vector<double>(nlayers - 1, constantthickness[0]);
			}
			else {
				//all good
			}

		}

		FILE* fp = fileopen(filename, "r");
		std::string str;
		std::vector<std::vector<double>> M;

		int k = 0;
		while (filegetline(fp, str)) {
			if (k % subsample == 0) {
				M.push_back(getdoublevector(str.c_str(), " "));
			}
			k++;
		}
		fclose(fp);


		nsamples = (int)M.size();
		linenumber = (int)M[0][lcol];

		x.resize(nsamples);
		y.resize(nsamples);
		e.resize(nsamples);
		z.resize(nsamples);
		c.resize(nsamples);
		if (spreadfade) {
			cp10.resize(nsamples);
			cp90.resize(nsamples);
		}
		for (int si = 0; si < nsamples; si++) {
			x[si] = M[si][xcol];
			y[si] = M[si][ycol];
			e[si] = M[si][ecol];

			c[si].resize(nlayers);
			for (int li = 0; li < nlayers; li++) {
				c[si][li] = M[si][crcol1 + li];
				if (c[si][li] > 0.0) {
					if (isresistivity)c[si][li] = 1.0 / c[si][li];
					c[si][li] *= cscale;
				}
			}

			if (spreadfade) {
				cp10[si].resize(nlayers);
				for (int li = 0; li < nlayers; li++) {
					cp10[si][li] = M[si][cp10col1 + li];
					if (cp10[si][li] > 0.0) {
						cp10[si][li] *= cscale;
					}
				}

				cp90[si].resize(nlayers);
				for (int li = 0; li < nlayers; li++) {
					cp90[si][li] = M[si][cp90col1 + li];
					if (cp90[si][li] > 0.0) {
						cp90[si][li] *= cscale;
					}
				}
			}

			z[si].resize(nlayers + 1);
			z[si][0] = e[si];
			for (int li = 0; li < nlayers; li++) {

				double t;
				if (li < nlayers - 1) {
					if (isconstantthickness == true) {
						t = constantthickness[li];
					}
					else {
						t = M[si][tcol1 + li];
					}
				}
				else {
					if (isconstantthickness == true) {
						t = constantthickness[li - 1];
					}
					else {
						t = M[si][tcol1 + li - 1];
					}
				}

				z[si][li + 1] = z[si][li] - t;

			}
		}
	}
	*/

	void process() {
		calculateextents();
		createbitmap();
		generatesectiondata();
		drawgraticule();
		drawcolorbar();
		fixtransparenttext();
		save();
		deletebitmap();
	}

	void calculateextents()
	{
		bestfitlineendpoints(D.x, D.y, x1, y1, x2, y2);
		if (std::abs(x2 - x1) >= std::abs(y2 - y1)) {
			//if predominantly EW
			if (x1 > x2) {
				std::swap(x1, x2);
				std::swap(y1, y2);
			}
		}
		else {
			//if predominantly NS
			if (y1 > y2) {
				std::swap(x1, x2);
				std::swap(y1, y2);
			}
		}

		double rawlength = distance(x1, y1, x2, y2);
		double hlength = roundupnearest(rawlength, dh);

		x2 = x1 + (x2 - x1) * hlength / rawlength;
		y2 = y1 + (y2 - y1) * hlength / rawlength;

		int hmarginpixels = getmarginpixels();
		h0 = 0.0;
		h1 = hmarginpixels * dh;
		h2 = h1 + hlength;

		nhpixels = 1 + (int)((h2 - h0) / dh);
		dx = dh * (x2 - x1) / (h2 - h1);
		dy = dh * (y2 - y1) / (h2 - h1);
		x0 = x1 - (x2 - x1) * h1 / hlength;
		y0 = y1 - (y2 - y1) * h1 / hlength;

		double zmin = (std::numeric_limits<double>::max)();
		double zmax = (std::numeric_limits<double>::lowest)();
		for (int si = 0; si < D.nsamples; si++) {
			if (D.z[si][D.nlayers] < zmin)zmin = D.z[si][D.nlayers];
			if (D.z[si][0] > zmax)zmax = D.z[si][0];
		}

		if (autozsectionbot == false) {
			zmin = zsectionbot;
		}

		if (autozsectiontop == false) {
			zmax = zsectiontop;
		}

		elevation_median = median(D.e.data(), D.e.size());

		double rawheight = zmax - zmin;
		vlength = roundupnearest(rawheight, dv);
		v0 = zmin;
		v1 = zmin;
		v2 = v1 + vlength;
		nvpixels = 1 + (int)(vlength / dv);
	}

	int wh2lx(double wh) {
		int lx = (int)roundnearest((wh - h0) / dh, 1.0);
		return lx;
	}

	int wv2ly(double wv) {
		int ly = (int)roundnearest((v2 - wv) / dv, 1.0);
		return ly;
	}

	int getmarginpixels() {
		int graticulemargin = 65;
		int colourbarmargin = 95;

		int margin = graticulemargin;
		if (m_drawcbar) margin += colourbarmargin;
		return margin;
	}

	void createbitmap() {
		pBitmap = new Bitmap(nhpixels, nvpixels, PixelFormat32bppARGB);
		for (int i = 0; i < nhpixels; i++) {
			for (int j = 0; j < nvpixels; j++) {
				pBitmap->SetPixel(i, j, Color(255, 255, 255, 255));
			}
		}
	}

	void deletebitmap() {
		delete pBitmap;
	}

	double fadevalue(const double& spread) {
		if (spread <= spreadfadelowclip)return  lowspreadfade;
		if (spread >= spreadfadehighclip)return highspreadfade;
		return lowspreadfade + (spread - spreadfadelowclip) / (spreadfadehighclip - spreadfadelowclip) * (highspreadfade - lowspreadfade);
	}

	void fadecolor(Color& clr, const double& fade) {
		//fade=0 gives original color
		//fade=1 gives white
		double alpha = 1.0 - fade;
		unsigned char w = (unsigned char)(255.0 * fade);
		unsigned char r = w + (unsigned char)(clr.GetRed() * alpha);
		unsigned char g = w + (unsigned char)(clr.GetGreen() * alpha);
		unsigned char b = w + (unsigned char)(clr.GetBlue() * alpha);
		clr = Color(255, r, g, b);
	}

	void generatesectiondata() {

		//Unit vec in section direction
		cVec u(dx, dy, 0.0);
		u.unitise();

		int hp1 = wh2lx(h1);
		int hp2 = wh2lx(h2);
		int vp1 = wv2ly(v1);
		int vp2 = wv2ly(v2);
		for (int i = hp1; i <= hp2; i++) {
			double xp = x0 + (double)i * dx;
			double yp = y0 + (double)i * dy;

			double mind = (std::numeric_limits<double>::max)();
			int mini = -1;
			for (int si = 0; si < D.nsamples; si++) {
				//double d = distance(0.0,0.0,xp-x[si],yp-y[si]);
				//Projected distance along section to sample
				cVec v(xp - D.x[si], yp - D.y[si], 0.0);
				double d = std::abs(u.dot(v));
				if (d < mind) {
					mini = si;
					mind = d;
				}
			}

			for (int j = vp2; j <= vp1; j++) {
				pBitmap->SetPixel(i, j, BkgColor);

				double zp = v2 - (double)j * dv;
				if (zp > D.e[mini]) {
					pBitmap->SetPixel(i, j, AirColor);
				}
				else {
					if (mind > (dh * 2)) {
						//if more than 2 pixels away then leave it null
						continue;
					}
					for (int li = 0; li < D.nlayers; li++) {
						if (zp < D.z[mini][li] && zp >= D.z[mini][li + 1]) {
							double conductivity = D.c[mini][li];
							if (conductivity < 0.0) {
								pBitmap->SetPixel(i, j, NullsColor);
							}
							else {
								int ind;
								if (Log10Stretch) ind = cStretch::log10stretch(conductivity, LowClip, HighClip);
								else ind = cStretch::linearstretch(conductivity, LowClip, HighClip);

								Color clr(255, m_cmap.r[ind], m_cmap.g[ind], m_cmap.b[ind]);

								if (spreadfade) {
									double conductivity_p10 = D.cp10[mini][li];
									double conductivity_p90 = D.cp90[mini][li];
									double spread = log10(conductivity_p90) - log10(conductivity_p10);
									double fade = fadevalue(spread);
									fadecolor(clr, fade);
								}
								pBitmap->SetPixel(i, j, clr);
							}
							break;
						}
					}
				}//layer loop
			}//v pixel loop
		}//h pixel loop		
	}

	REAL getfontsize() {
		//REAL fontsize = (REAL)abs((wv2ly(0.33*gratdiv) - wv2ly(0)));
		REAL fontsize = 8.5;
		return fontsize;
	}

	void drawgraticule() {

		Graphics gr(pBitmap);
		gr.SetCompositingMode(CompositingModeSourceOver);
		gr.SetTextRenderingHint(TextRenderingHintAntiAlias);

		int zg1 = (int)roundnearest(v1, (int)gratdiv);
		int zg2 = (int)roundnearest(v2, (int)gratdiv);

		Pen blackpen(Color::Black, 0);
		SolidBrush blackbrush(Color::Black);
		Font font(L"Arial", getfontsize(), FontStyleRegular, UnitPoint);

		StringFormat textformat;
		textformat.SetAlignment(StringAlignmentFar);
		textformat.SetLineAlignment(StringAlignmentCenter);
		SizeF layoutsize(32767, 32767);
		SizeF textsize;
		gr.MeasureString(L"-0000 m", -1, &font, layoutsize, &textformat, &textsize);

		for (int zg = zg1; zg <= zg2; zg += (int)gratdiv) {
			gr.DrawLine(&blackpen, wh2lx(h1), wv2ly(zg), wh2lx(h2), wv2ly(zg));

			wchar_t s[20];
			swprintf(s, 20, L"%5d m", zg);
			gr.MeasureString(s, -1, &font, layoutsize, &textformat, &textsize);

			int tx = wh2lx(h1);
			int ty = wv2ly(zg);
			if (ty > nvpixels)continue;
			if (ty - textsize.Height < 0)continue;
			PointF txpos((REAL)tx, (REAL)ty);
			gr.DrawString(s, -1, &font, txpos, &textformat, &blackbrush);
		}
	}

	void drawcolorbar() {

		if (m_drawcbar == false)return;

		Pen blackpen(Color::Black, 0);
		SolidBrush blackbrush(Color::Black);

		Graphics gr(pBitmap);

		TextRenderingHint hint = gr.GetTextRenderingHint();
		gr.SetTextRenderingHint(TextRenderingHintAntiAlias);

		int ph1 = wv2ly(v1);
		int ph2 = wv2ly(v2);
		ph1 = 25; ph2 = 40;

		int pvtop = wv2ly(v1 + (v2 - v1) * 0.95);
		int pvbot = wv2ly(v1 + (v2 - v1) * 0.05);
		for (int j = pvtop; j <= pvbot; j++) {
			int ind = 255 * (j - pvbot) / (pvtop - pvbot);
			for (int i = ph1; i <= ph2; i++) {
				pBitmap->SetPixel(i, j, Color(255, m_cmap.r[ind], m_cmap.g[ind], m_cmap.b[ind]));
			}
		}
		gr.DrawRectangle(&blackpen, ph1, pvtop, ph2 - ph1, pvbot - pvtop);

		Font font(L"Arial", getfontsize(), FontStyleRegular, UnitPoint);
		StringFormat textformat;
		textformat.SetAlignment(StringAlignmentNear);
		textformat.SetLineAlignment(StringAlignmentCenter);
		SizeF layoutsize(32767, 32767);
		SizeF textsize;
		gr.MeasureString(L"0.0001 m", -1, &font, layoutsize, &textformat, &textsize);


		for (int i = 0; i < m_cbarticks.size(); i++) {
			if (m_cbarticks[i] < LowClip)continue;
			if (m_cbarticks[i] > HighClip)continue;

			int tickv;
			if (Log10Stretch) {
				int ind = cStretch::log10stretch(m_cbarticks[i], LowClip, HighClip);
				tickv = pvbot + ind * (pvtop - pvbot) / 255;
			}
			else {
				int ind = cStretch::linearstretch(m_cbarticks[i], LowClip, HighClip);
				tickv = pvbot + ind * (pvtop - pvbot) / 255;
			}

			wchar_t s[20];
			swprintf(s, 20, L"%5.3lf", m_cbarticks[i]);

			PointF p;
			p.X = (REAL)ph2 + 3;
			p.Y = (REAL)tickv;
			gr.DrawString(s, -1, &font, p, &textformat, &blackbrush);
			gr.DrawLine(&blackpen, ph2 - 2, tickv, ph2 + 2, tickv);
		}

		PointF p;
		p.X = (REAL)ph1 - 2;
		p.Y = (REAL)(pvtop + pvbot) / 2;

		textformat.SetAlignment(StringAlignmentCenter);
		textformat.SetLineAlignment(StringAlignmentFar);
		gr.TranslateTransform(p.X, p.Y);
		gr.RotateTransform(-90);
		gr.DrawString(L"Conductivity (S/m)", -1, &font, PointF(0, 0), &textformat, &blackbrush);

	}

	void fixtransparenttext() {
		int hp0 = wh2lx(h0);
		int hp1 = wh2lx(h1);
		int vp1 = wv2ly(v1);
		int vp2 = wv2ly(v2);
		for (int i = hp0; i <= hp1; i++) {
			for (int j = vp2; j <= vp1; j++) {
				Color c;
				pBitmap->GetPixel(i, j, &c);
				pBitmap->SetPixel(i, j, Color(255, c.GetR(), c.GetG(), c.GetB()));
			}
		}
	}

	std::string basename() const {
		std::string s = prefix + strprint("%d", D.linenumber) + suffix;
		return s;
	}

	std::string imagefile_nod(const std::string& imageformat) const {
		std::string s = basename() + "." + imageformat;
		return s;
	}

	std::string imagefile(const std::string& imageformat) const {
		std::string s = outdir + imagefile_nod(imageformat);
		return s;
	}

	std::string worldfile_nod(const std::string& worldfileextension) const {
		std::string s = basename() + "." + worldfileextension;
		return s;
	}

	std::string worldfile(const std::string& worldfileextension) const {
		std::string s = outdir + worldfile_nod(worldfileextension);
		return s;
	}


	void saveimage(Bitmap& bm, const std::string imageformat) {
		cGDIplusHelper::saveimage(&bm, imagefile(imageformat));
	}

	void save() {
		makedirectorydeep(outdir);
		if (SavePNG) {
			saveimage(*pBitmap,"png");
			saveworldfile(worldfile("pngw"));
		}

		if (SaveJPG) {
			saveimage(*pBitmap, "jpg");
			saveworldfile(worldfile("jgw"));
		}
	}

	void saveworldfile(std::string& worldfilepath) {

		//Line 1: A: x component of the pixel width (x-scale)
		//Line 2: D: y component of the pixel width (y-skew)
		//Line 3: B: x component of the pixel height (x-skew)
		//Line 4: E: y component of the pixel height (y-scale), almost always negative
		//Line 5: C: x-coordinate of the center of the upper left pixel
		//Line 6: F: y-coordinate of the center of the upper left pixel

		double angle = atan2(y2 - y1, x2 - x1);
		double vshift = (v2 - elevation_median);
		double ix0 = x0 - (dh / dv) * vshift * sin(angle);
		double iy0 = y0 + (dh / dv) * vshift * cos(angle);

		FILE* fp = fileopen(worldfilepath, "w");
		fprintf(fp, "%lf\n", dh * cos(angle));
		fprintf(fp, "%lf\n", dh * sin(angle));
		fprintf(fp, "%lf\n", dh * sin(angle));
		fprintf(fp, "%lf\n", -dh * cos(angle));
		fprintf(fp, "%lf\n", ix0);
		fprintf(fp, "%lf\n", iy0);
		fclose(fp);
	}

};

int main(int argc, char** argv)
{
	unsigned long long gdikey;
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
		std::string headerfile;
		if (ib.getvalue("DfnFile", headerfile) == true) {
			glog.logmsg("Headerfile = %s'\n", headerfile.c_str());
			glog.logmsg("Note: in future please use 'HeaderFile = ...' instead of 'DfnFile = ...'\n");
		}
		else if (ib.getvalue("HeaderFile", headerfile) == true) {
			glog.logmsg("Headerfile = %s'\n", headerfile.c_str());
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

		if (linefieldindex <= 0) {
			glog.logmsg("Error: cannot find the line field %s\n", linefieldname.c_str());
			return 0;
		}

		std::string infiles = ib.getstringvalue("DataFiles");
		std::vector<std::string> filelist = cDirectoryAccess::getfilelist(infiles);
		if (filelist.size() == 0) {
			glog.logmsg("Error: no data files found matching %s\n", infiles.c_str());
			return 0;
		}

		gdikey = cGDIplusHelper::start();
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

				cGeorefSection S(D);
				S.getoptions(sb);
				//S.readdatafile(ib, filelist[i].c_str());
				S.process();
			}
		}
		printf("Done ... \nElapsed time = %.3lf seconds\n", stopwatch.etimenow());
		cGDIplusHelper::stop(gdikey);
	}
	catch (std::runtime_error& e) {
		std::cout << e.what();
		cGDIplusHelper::stop(gdikey);
	}
	return 0;
}

