/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _inputmanager_H_
#define _inputmanager_H_

#include "samplebunch.h"
#include "asciicolumnfile.h"
#include "fielddefinition.h"
#if defined HAVE_NETCDF
	#include "geophysics_netcdf.h"
#endif

class cInputManager {

public:
	enum class IOType { ASCII, NETCDF, NONE };

protected:		
	std::string DataFileName;
	IOType iotype = IOType::NONE;
	size_t Subsample = 1;
	size_t Record = 0;

public:
	cInputManager() {};

	virtual ~cInputManager() {};

	void initialise(const cBlock& b) {
		DataFileName = b.getstringvalue("DataFile");
		fixseparator(DataFileName);
	}

	void set_subsample_rate(const size_t& subsamplerate) {
		Subsample = subsamplerate;		
	}

	static bool isnetcdf(const cBlock& b) {
		std::string fname = b.getstringvalue("DataFile");
		std::string ext = extractfileextension(fname);
		if (strcasecmp(ext, ".nc") == 0) {
			return true;
		}
		return false;
	}
	
	size_t subsamplerate() { return Subsample; }

	virtual bool is_record_valid() { return true; }
	
	virtual bool get_bunch(cSampleBunch& bunch, const cFieldDefinition& fd, const int& pointindex, const int& bunchsize, const int& bunchsubsample)
	{
		glog.errormsg(_SRC_+"\nfunction not yet implemented\n");
		return false;
	};

	virtual bool load_record(const size_t& record) = 0;

	virtual bool parse_record() { return true; }
	
	virtual bool get_acsiicolumnfield(const std::string& fname, cAsciiColumnField& c) const {
		glog.errormsg("get_acsiicolumnfield() not yet implemented\n");
		return true;
	}

	virtual bool set_variant_type(const std::string& fname, cVrnt& vnt) const {
		glog.errormsg("set_variant_type() not yet implemented\n");
		return true;
	}
		
	const std::string& datafilename() { return DataFileName; }

	const size_t& record() const { return Record; }
			
	bool readvnt(cFDVar& fdv, const size_t n=1)
	{
		cFieldDefinition& fd  = fdv.first;
		cVrnt& vnt = fdv.second;

		auto ReadVisitor = [&](auto& t) {
			bool s = read(fd, t, n);
		};		
		std::visit(ReadVisitor, vnt);		
		return true;
	}

	bool readfdvnt(cFdVrnt& fdv, const size_t n = 1)
	{		
		auto ReadVisitor = [&](auto& t) {
			bool s = read(fdv.fd, t, n);
		};
		std::visit(ReadVisitor, fdv.vnt);
		return true;
	}


	template<typename T>
	bool read(const cFieldDefinition& fd, T& v, const size_t n=1)
	{
		if (fd.definitiontype() == cFieldDefinition::TYPE::UNAVAILABLE) {
			v = undefinedvalue(v);
			return false;
		}
		else if (fd.definitiontype() == cFieldDefinition::TYPE::NUMERIC) {
			v = (T)fd.numericvalue[0];
			return true;
		}

		std::vector<T> vec;
		file_read(fd,vec,1);
		v = vec[0];

		if (iotype != IOType::ASCII) {//Don't flip if ASCII reader as it already does this - to be fixed
			if (fd.flip) { v = -1 * v; } {
				fd.applyoperator(v);
			}
		}
		return true;
	}

	template<typename T>
	bool read(const cFieldDefinition& fd, std::vector<T>& vec, const size_t n)
	{
		vec.resize(n);
		if (fd.definitiontype() == cFieldDefinition::TYPE::NUMERIC) {
			size_t deflen = fd.numericvalue.size();
			if (deflen != 1 && deflen != n) {
				std::ostringstream oss;
				oss << "Mismatch in field sizes for '<" << fd.keyname << ">' : expected " << n << " but got " << deflen << std::endl;
				glog.errormsg(oss.str());
			}
			

			if (deflen == 1) {
				for (size_t i = 0; i < n; i++) {
					vec[i] = (T)fd.numericvalue[0];					
				}
			}
			else {				
				for (size_t i = 0; i < n; i++) {					
					vec[i] = (T) fd.numericvalue[i];
				}
			}
			return true;
		}

		bool status = file_read(fd, vec, n);

		if (iotype != IOType::ASCII) {//Don't flip if ASCII reader as it already does this - to be fixed
			for (size_t i = 0; i < vec.size(); i++) {
				if (fd.flip) { vec[i] = -vec[i]; }
				fd.applyoperator(vec[i]);
			}
		}
		return true;
	}
	
	//virtual template classes are not allowed - therefore repeated
	virtual bool file_read(const cFieldDefinition& fd, std::vector<char>& vec, const size_t n) = 0;
	virtual bool file_read(const cFieldDefinition& fd, std::vector<int>& vec, const size_t n) = 0;
	virtual bool file_read(const cFieldDefinition& fd, std::vector<float>& vec, const size_t n) = 0;
	virtual bool file_read(const cFieldDefinition& fd, std::vector<double>& vec, const size_t n) = 0;	
};

class cASCIIInputManager : public cInputManager {

private:
	cAsciiColumnFile AF;	

public:

	std::string HeaderFileName;
	std::size_t HeaderLines;
		
	cASCIIInputManager(const cBlock& b) {
		cInputManager::initialise(b);
		initialise(b);
	}

	~cASCIIInputManager() {	};

	void initialise(const cBlock& b)
	{		
		iotype = IOType::ASCII;
		HeaderFileName = b.getstringvalue("DfnFile");
		if(!isdefined(HeaderFileName)){
			HeaderFileName = b.getstringvalue("HeaderFile");
		}
		
		if(isdefined(HeaderFileName)) {
			fixseparator(HeaderFileName);
			if (!exists(HeaderFileName)) {
				std::string msg = _SRC_;
				msg += strprint("\n\tD'oh! the specified header file (%s) does not exist\n", HeaderFileName.c_str());
				throw(std::runtime_error(msg));
			}
		}

		fixseparator(DataFileName);
		if (!exists(DataFileName)) {
			std::string msg = _SRC_;
			msg += strprint("\n\tD'Oh! the specified data file (%s) does not exist\n", DataFileName.c_str());
			throw(std::runtime_error(msg));
		}
		
		AF.openfile(DataFileName);

		if (isdefined(HeaderFileName)) {			
			glog.logmsg(0, "Parsing input HeaderFile %s\n", HeaderFileName.c_str());
			std::string ext = extractfileextension(HeaderFileName);
			if (strcasecmp(ext,".dfn")==0){
				AF.read_dfn(HeaderFileName);
				AF.headertype = cAsciiColumnFile::HeaderType::DFN;
				AF.parsetype  = cAsciiColumnFile::ParseType::FIXEDWIDTH;
			}
			else if (strcasecmp(ext,".csv")==0) {
				AF.parse_csv_header(HeaderFileName);
				AF.headertype = cAsciiColumnFile::HeaderType::CSV;
				AF.parsetype = cAsciiColumnFile::ParseType::FIXEDWIDTH;
			}
			else if (strcasecmp(ext,".csvh")==0) {
				AF.parse_csv_header(HeaderFileName);
				AF.headertype = cAsciiColumnFile::HeaderType::CSV;
				AF.parsetype = cAsciiColumnFile::ParseType::FIXEDWIDTH;
			}
			else {
				std::string msg = _SRC_;
				msg += strprint("\n\tD'oh! the specified header file (%s) is not .dfn or .csv or .csvh\n", HeaderFileName.c_str());
				throw(std::runtime_error(msg));
			}
		}
		else{
			AF.headertype = cAsciiColumnFile::HeaderType::NONE;
			AF.parsetype  = cAsciiColumnFile::ParseType::DELIMITED;
		}

		size_t headerlines = b.getsizetvalue("Headerlines");
		if (!isdefined(headerlines)) { headerlines = 0; }
		HeaderLines = headerlines;				
	}

	bool is_record_valid() {
		bool status =  AF.is_record_valid();
		if (status == false) {
			std::string msg;
			msg = strprint("Skipping non-valid record at line %zu of Input DataFile %s\n", record(), datafilename().c_str());
			glog.logmsg(msg);
			std::cerr << msg;

			msg = strprint("%s\n", recordstring().c_str());
			glog.logmsg(msg);
			std::cerr << msg;
			return false;
		}
		return true;
	}
	
	bool load_record(const size_t& n)
	{
		Record = n;
		return AF.load_record(n+HeaderLines);
	}

	bool parse_record() {
		size_t n = AF.parse_record();
		if (n <= 1) return false;
		return true;		
	}
	
	template<typename T>
	bool get_one(const cFieldDefinition& fd, const size_t& pointindex, T& val) {
		if (load_record(pointindex)) {
			if (parse_record()) {
				if (read(fd, val)) {
					return true;
				};
			}
		}
		return false;
	}
	
	bool get_bunch(cSampleBunch& bunch, const cFieldDefinition& fd, const int& pointindex, const int& bunchsize, const int& bunchsubsample)
	{							
		int line;
		bool status = get_one(fd, pointindex, line);
		if (status == false) return false;

		int pn = pointindex;
		int ln = line;

		int pa = pointindex - bunchsubsample*((bunchsize-1)/2);
		while (pa < 0) pa += bunchsubsample;
		pn = pointindex - bunchsubsample;
		if (pn < pa) pn = pa;

		while(pn>pa){
			bool status = get_one(fd, pn, ln);
			if (status && ln == line) {				
				pn -= bunchsubsample;
			}
			else{
				pa = pn + bunchsubsample;
				break;					
			}													
		}

		int pb = pa + bunchsubsample * (bunchsize - 1);
		pn = pointindex;
		ln = line;
		while (pn <= pb) {
			bool status = get_one(fd, pn, ln);
			if (status && ln == line) {
				pn += bunchsubsample;
			}
			else {
				pb = pn - bunchsubsample;
				pa = pb - bunchsubsample * (bunchsize - 1);
				break;
			}
		}

		std::vector<std::size_t> indices = increment((size_t)bunchsize, (size_t)pa,(size_t)bunchsubsample);
		bunch=cSampleBunch(indices,pointindex);
		return true;
	};

	const std::string& recordstring() const { return AF.currentrecord_string(); }

	const std::vector<std::string>& fields() const { return AF.currentrecord_columns(); }
	
	bool file_read(const cFieldDefinition& fd, std::vector<char>& vec, const size_t n) { return file_read_impl(fd, vec, n); }
	bool file_read(const cFieldDefinition& fd, std::vector<int>& vec, const size_t n) { return file_read_impl(fd, vec, n); }
	bool file_read(const cFieldDefinition& fd, std::vector<float>& vec, const size_t n) { return file_read_impl(fd, vec, n); }
	bool file_read(const cFieldDefinition& fd, std::vector<double>& vec, const size_t n) { return file_read_impl(fd, vec, n); }

	template<typename T>
	bool file_read_impl(const cFieldDefinition& fd, std::vector<T>& vec, const size_t n)
	{									
		bool status = AF.getvec_fielddefinition(fd, vec, n);
		return status;	
	}	
	
	bool get_acsiicolumnfield(const std::string& fname, cAsciiColumnField& c) const {		
		int findex = AF.fieldindexbyname(fname);
		if (findex >= 0) {
			c = AF.fields[findex];
			return true;
		}
		return false;		
	}

	bool set_variant_type(const std::string& fname, cVrnt& vnt) const {
		cAsciiColumnField c;
		bool status = get_acsiicolumnfield(fname, c);
		if (status) {			
			c.set_variant_type(vnt);
			return true;
		}
		else {
			glog.errormsg("Could not find field %s",fname.c_str());
			return false;
		}		
	}

};

#if defined HAVE_NETCDF
class cNetCDFInputManager : public cInputManager {

private:	
	cGeophysicsNcFile NC;

public:

	cNetCDFInputManager(const cBlock& b) {
		cInputManager::initialise(b);
		initialise(b);
	}

	~cNetCDFInputManager() {	};

	void initialise(const cBlock& b)
	{						
		glog.logmsg(0, "Opening Input DataFile %s\n", DataFileName.c_str());		
		iotype = IOType::NETCDF;
		NC.open(DataFileName, netCDF::NcFile::FileMode::read);		
	}

	bool load_record(const size_t& n)
	{
		Record = n;
		#if defined HAVE_NETCDF
			if (Record > NC.ntotalsamples()) return false;
		#endif
		return true;
	}

	bool file_read(const cFieldDefinition& fd, std::vector<char>& vec, const size_t n) { return file_read_impl(fd, vec, n); }
	bool file_read(const cFieldDefinition& fd, std::vector<int>& vec, const size_t n) { return file_read_impl(fd, vec, n); }
	bool file_read(const cFieldDefinition& fd, std::vector<float>& vec, const size_t n) { return file_read_impl(fd, vec, n); }
	bool file_read(const cFieldDefinition& fd, std::vector<double>& vec, const size_t n) { return file_read_impl(fd, vec, n); }

	template<typename T>
	bool file_read_impl(const cFieldDefinition& fd, std::vector<T>& v, const size_t n)
	{
		NC.getDataByPointIndex(fd.varname, Record, v);
		return true;
	}	
};
#endif

#endif

 