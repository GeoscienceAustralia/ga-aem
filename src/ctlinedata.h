/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _ctlinedata_H
#define _ctlinedata_H

#include <vector>

#include "blocklanguage.h"
#include "file_formats.h"

class cCTLineData {

private:
	cBlock B;
	std::vector<cAsciiColumnField> fields;

public:
	int linenumber;
	size_t nlayers;
	size_t nsamples;
	bool spreadfade = false;
	std::string inputdatumprojection;
	std::vector<double> fid;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> e;	
	std::vector<double> linedistance;
	std::vector<std::vector<double>> c;
	std::vector<std::vector<double>> z;
	std::vector<std::vector<double>> cp10;
	std::vector<std::vector<double>> cp90;

	cCTLineData(const cBlock& b, const std::vector<cAsciiColumnField>& _fields) {
		B = b;
		fields = _fields;
	}

	cCTLineData(const cBlock& b, const std::string& headerfile){
		B = b;
		cASEGGDF2Header A(headerfile);
		if (A.isvalid()) {
			fields = A.getfields();
		}
		else {
			cHDRHeader H(headerfile);
			if (H.isvalid()) {
				fields = H.getfields();
			}
		}
	}

	void load(const std::vector<std::string>& L)
	{

		inputdatumprojection = B.getstringvalue("DatumProjection");

		int subsample = B.getintvalue("Subsample");
		if (isdefined(subsample) == false)subsample = 1;

		cRange<int> lcol = getcolumns("Line");
		cRange<int> fcol = getcolumns("Fiducial");
		cRange<int> xcol = getcolumns("Easting");
		cRange<int> ycol = getcolumns("Northing");
		cRange<int> ecol = getcolumns("Elevation");
		
		bool isresistivity = false;
		cRange<int> crcol = getcolumns("Conductivity");
		if (crcol.valid() == false){
			crcol = getcolumns("Resistivity");
			isresistivity = true;
		}

		if (crcol.valid() == false){
			std::string msg = strprint("Either a conductivity or resistivity field must be specified\n") + _SRC_;
			throw(std::runtime_error(msg));
		}
		nlayers = crcol.to - crcol.from + 1;

		double cscale = 1.0;
		std::string cunits = B.getstringvalue("InputConductivityUnits");
		if (isdefined(cunits) == false){
			cscale = 1.0;
		}
		else if (strcasecmp(cunits, "S/m") == 0){
			cscale = 1.0;
		}
		else if (strcasecmp(cunits, "mS/m") == 0){
			cscale = 0.001;
		}
		else{
			glog.logmsg("Unknown InputConductivityUnits %s\n", cunits.c_str());
		}

		cRange<int> cp10col, cp90col;
		if (spreadfade){
			cp10col = getcolumns("Conductivity_p10");
			cp10col = getcolumns("Conductivity_p90");
		}

		//Thickness
		bool isconstantthickness = false;
		std::vector<double> constantthickness;						
		cRange<int> tcol = getcolumns("Thickness");
		if (tcol.valid() == true){
			isconstantthickness = false;
		}
		else{
			isconstantthickness = true;
			constantthickness = B.getdoublevector("Thickness");
			if (constantthickness.size() == 0){
				std::string msg("Thickness not set\n");
				throw(std::runtime_error(msg));
			}
			else if (constantthickness.size() == 1) {
				constantthickness = std::vector<double>(nlayers - 1, constantthickness[0]);
			}
			else if (constantthickness.size() > 1 && constantthickness.size() < nlayers - 1){
				std::string msg("Thickness not set correctly\n");
				throw(std::runtime_error(msg));
			}			
			else{
				//all good
			}
		}

		std::vector<std::vector<double>> M;
		for (size_t i = 0; i<L.size(); i++){
			if (i % subsample == 0){
				M.push_back(getdoublevector(L[i].c_str(), " "));
			}
		}

		nsamples = M.size();
		linenumber = (int)M[0][lcol.from];

		fid.resize(nsamples);
		x.resize(nsamples);
		y.resize(nsamples);
		e.resize(nsamples);
		z.resize(nsamples);
		c.resize(nsamples);
		if (spreadfade){
			cp10.resize(nsamples);
			cp90.resize(nsamples);
		}

		for (int si = 0; si < nsamples; si++){
			if(fcol.from >=0 ) fid[si] = M[si][fcol.from];
			x[si] = M[si][xcol.from];
			y[si] = M[si][ycol.from];
			e[si] = M[si][ecol.from];
			
			c[si].resize(nlayers);
			for (int li = 0; li < nlayers; li++){
				c[si][li] = M[si][crcol.from + li];
				if (c[si][li] > 0.0){
					if (isresistivity)c[si][li] = 1.0 / c[si][li];
					c[si][li] *= cscale;
				}
			}

			if (spreadfade){
				cp10[si].resize(nlayers);
				for (int li = 0; li < nlayers; li++){
					cp10[si][li] = M[si][cp10col.from + li];
					if (cp10[si][li] > 0.0){
						cp10[si][li] *= cscale;
					}
				}

				cp90[si].resize(nlayers);
				for (int li = 0; li < nlayers; li++){
					cp90[si][li] = M[si][cp90col.from + li];
					if (cp90[si][li] > 0.0){
						cp90[si][li] *= cscale;
					}
				}
			}

			z[si].resize(nlayers + 1);
			z[si][0] = e[si];
			for (int li = 0; li < nlayers; li++){

				double t;
				if (li < nlayers - 1){
					if (isconstantthickness == true){
						t = constantthickness[li];
					}
					else{
						t = M[si][tcol.from + li];
					}
				}
				else{
					if (isconstantthickness == true){
						t = constantthickness[li - 1];
					}
					else{
						t = M[si][tcol.from + li - 1];
					}
				}
				z[si][li + 1] = z[si][li] - t;
			}
		}

		linedistance.resize(x.size());
		linedistance[0] = 0.0;
		for (size_t i = 1; i < x.size(); i++){
			linedistance[i] = linedistance[i - 1] + distance(x[i - 1], y[i - 1], x[i], y[i]);
		}

	}

	int field_index_by_name(const std::string& fieldname) const {
		int index = field_index_by_name_impl(fields, fieldname);
		return index;
	};

	static bool parse_column_range(const std::string& token, cRange<int>& r) {
		int status = sscanf(token.c_str(), "Column %d-%d", &r.from, &r.to);
		if (status == 1) {
			r.to = r.from;
			r.from--; r.to--;
			return true;
		}
		else if (status == 2) {
			r.from--; r.to--;
			return true;
		}
		else return false;
	}

	cRange<int> getcolumns(const std::string& token) {
		return getcolumns_block_fields(B, token);
	};

	cRange<int> getcolumns_block_fields(const cBlock& b, const std::string &fieldkey){

		std::string fieldvalue = b.getstringvalue(fieldkey);
		cRange<int> r(-1,-1);

		if (fields.size() == 0) {
			bool status = parse_column_range(fieldvalue, r);
			return r;
		}
		else {
			int findex = field_index_by_name(fieldvalue);
			if (findex >= 0 && findex < fields.size()) {
				r.from = fields[findex].startcol();
				r.to = fields[findex].endcol();
				return r;
			}
		}
		return r;//invalid;
	}
	

};

#endif
