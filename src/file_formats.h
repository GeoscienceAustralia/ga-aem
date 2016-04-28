/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _file_formats_H
#define _file_formats_H

#include <cstring>
#include <vector>
#include "general_utils.h"

class cAsciiColumnField{

public:
	size_t order;
	size_t startcolumn;
	std::string name;
	std::string expandedname;

	size_t nbands;

	char   fmttype;
	size_t fmtwidth;
	size_t fmtdecimals;
	
	std::string units;

	bool hasnullvalue;
	std::string nullvaluestr;
	double nullvalue;
	std::string comment;

	cAsciiColumnField(){
		initialise();
	};

	cAsciiColumnField(const size_t _order, const size_t _startcolumn, const std::string _name, const char _fmttype, const size_t _fmtwidth, const size_t _fmtdecimals, const size_t _nbands = 1){
		initialise();
		order = _order;
		startcolumn = _startcolumn;
		name = _name;
		fmttype = _fmttype;
		fmtwidth = _fmtwidth;
		fmtdecimals = _fmtdecimals;
		nbands = _nbands;
	};

	void initialise(){
		nbands = 1;
		fmtdecimals = 0;
		hasnullvalue = false;
		nullvalue = -32767;
	};

	std::string simple_header_record(){

		std::string s;
		if (nbands == 1){
			s = strprint("%lu\t%s\n", startcolumn + 1, name.c_str());
		}
		else{
			s = strprint("%lu-%lu\t%s\n", startcolumn + 1, startcolumn + nbands, name.c_str());
		}
		return s;
	}

	std::string i3_header_record(){

		std::string channelname = name;
		if (nbands > 1){
			channelname = strprint("%s{%lu}", name.c_str(), nbands);
		}

		std::string typestr;
		if (fmttype == 'I'){
			typestr = "INTEGER";
		}
		else if (fmttype == 'F'){
			typestr = "DOUBLE";
		}

		std::string s;
		s += strprint("DATA\t%lu, %lu, NORMAL, , , ,\n", startcolumn, fmtwidth);
		s += strprint("CHAN\t%s, %s, NORMAL, %lu, %lu, LABEL = \"%s\"\n", channelname.c_str(), typestr.c_str(), fmtwidth, fmtdecimals, name.c_str());
		return s;
	}

	std::string aseggdf_header_record(){
		std::string fmtstr;
		if (nbands > 1) fmtstr += strprint("%lu", nbands);
		fmtstr += strprint("%c", fmttype);
		fmtstr += strprint("%lu", fmtwidth);
		if (fmttype != 'I')fmtstr += strprint(".%lu", fmtdecimals);

		std::string s;
		s += strprint("DEFN %lu ST=RECD,RT=; %s : %s", order + 1, name.c_str(), fmtstr.c_str());

		char comma = ',';
		char colon = ':';
		std::string o = " :";
		if (units.size() > 0){
			if (o[o.size() - 1] != colon && o[o.size() - 1] != comma) o += comma;
			o += strprint(" UNITS = %s ", units.c_str());
		}

		if (nullvaluestr.size() > 0){
			if (o[o.size() - 1] != colon && o[o.size() - 1] != comma) o += comma;
			o += strprint(" NULL = %s ", nullvaluestr.c_str());
		}

		if (comment.size() > 0){
			if (o[o.size() - 1] != colon && s[s.size() - 1] != comma) o += comma;
			o += strprint(" %s", comment.c_str());
		}

		if (o.size() > 2) s += o;
		s += strprint("\n");
		return s;
	}

	void print(){		
		printf("\n");
		printf(" name=%s", name.c_str());
		printf(" order=%lu", order);
		printf(" startcolumn=%lu", startcolumn);		
		printf(" bands=%lu", nbands);
		printf(" type=%c", fmttype);
		printf(" width=%lu", fmtwidth);
		printf(" decimals=%lu", fmtdecimals);
		printf(" units=%s", units.c_str());
		printf(" nullvalue=%s", nullvaluestr.c_str());
		printf(" nullvalue=%lf", nullvalue);
		printf(" expandedname=%s", expandedname.c_str());
		printf(" comment=%s", comment.c_str());
		printf("\n");
	}

	std::string datatype(){
		std::string dtype;
		if (fmttype == 'I' || fmttype == 'i'){
			dtype = "Integer";
		}
		else if (fmttype == 'F' || fmttype == 'f'){
			dtype = "Real";
		}
		else if (fmttype == 'E' || fmttype == 'e'){
			dtype = "Real";
		}
		else{
			dtype = "Unknown";
			warningmessage("Unknown data type <%c>\n", fmttype);
		}
		return dtype;
	};

	bool isinteger(){
		if (strcasecmp(datatype(), "Integer") == 0)return true;
		return false;
	}

	bool isreal(){
		if (strcasecmp(datatype(), "Real") == 0)return true;
		return false;
	}

	bool isnull(const double v){
		if (hasnullvalue){
			if (v == nullvalue)return true;
		}
		return false;
	}
};

class cOutputFileInfo{

	size_t lastfield;
	size_t lastcolumn;
	bool   allowmorefields;

	public:

	std::vector<cAsciiColumnField> fields;

	cOutputFileInfo(){
		lastfield  = 0;
		lastcolumn = 0;
		allowmorefields = true;
	}

	void lockfields(){		
		allowmorefields = false;
	}

	void addfield(const std::string _name, const char _form, const size_t _width, const size_t _decimals, const size_t _nbands = 1){
		if (allowmorefields){
			cAsciiColumnField cf(lastfield,lastcolumn, _name, _form, _width, _decimals, _nbands);
			fields.push_back(cf);
			lastfield++;
			lastcolumn += _nbands;
		}
	}

	void setunits(const std::string _units){
		if (allowmorefields){			
			fields[lastfield-1].units = _units;			
		}
	}

	void setnullvalue(const std::string _nullvaluestr){
		if (allowmorefields){
			fields[lastfield - 1].nullvaluestr = _nullvaluestr;
		}
	}

	void setcomment(const std::string _comment){
		if (allowmorefields){
			fields[lastfield-1].comment = _comment;
		}
	}

	void write_simple_header(const std::string pathname){
		FILE* fp = fileopen(pathname.c_str(), "w");		
		for (size_t i = 0; i < fields.size(); i++){
			std::string s = fields[i].simple_header_record();
			fprintf(fp, s.c_str());
		}
		fclose(fp);
	};

	void write_PAi3_header(const std::string pathname){
		FILE* fp = fileopen(pathname.c_str(), "w");
		fprintf(fp, "[IMPORT ARCHIVE]\n");
		fprintf(fp, "FILEHEADER\t1\n");
		fprintf(fp, "RECORDFORM\tFIXED\n");
		fprintf(fp, "SKIPSTRING\t\"/\"\n");
		for (size_t i = 0; i < fields.size(); i++){
			std::string s = fields[i].i3_header_record();
			fprintf(fp, s.c_str());
		}	
		fclose(fp);
	};

	void write_aseggdf_header(const std::string pathname){
		FILE* fp = fileopen(pathname.c_str(), "w");		
		fprintf(fp,"DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A76\n");
		for (size_t i = 0; i < fields.size(); i++){
			std::string s = fields[i].aseggdf_header_record();
			fprintf(fp, s.c_str());
		}
		fprintf(fp, "DEFN %lu ST=RECD,RT=;END DEFN\n", fields.size()+1);
		fclose(fp);
	};

};

#endif

 