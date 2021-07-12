/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _inputmanager_H_
#define _inputmanager_H_

#include "asciicolumnfile.h"
#include "fielddefinition.h"
#if defined HAVE_NETCDF
	#include "geophysics_netcdf.h"
#endif

class cASCIIInputManager;
class cNetCDFInputManager;

class cInputManager {

public:
	enum class IOType { ASCII, NETCDF, NONE };

protected:	
	bool AtStart = true;
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

		Subsample = b.getsizetvalue("Subsample");
		if (!isdefined(Subsample)) { Subsample = 1; }
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
		
	virtual bool readnextrecord() = 0;

	virtual bool parserecord() { return true; }
	

	const std::string& datafilename() { return DataFileName; }

	const size_t& record() const { return Record; }
	
	template<typename T>
	bool read(const cFieldDefinition& cd, T& v)
	{
		if (cd.definitiontype() == UNAVAILABLE) {
			v = undefinedvalue(v);
			return false;
		}
		else if (cd.definitiontype() == NUMERIC) {
			v = (T)cd.numericvalue[0];
			return true;
		}

		std::vector<T> vec;
		bool status = file_read(cd,vec,1);
		v = vec[0];

		if (iotype != IOType::ASCII) {//Don't flip if ASCII reader as it already does this - to be fixed
			if (cd.flip) { v = -1 * v; } {
				cd.applyoperator(v);
			}
		}
		return true;
	}

	template<typename T>
	bool read(const cFieldDefinition& cd, std::vector<T>& vec, const size_t n)
	{
		vec.resize(n);
		if (cd.definitiontype() == NUMERIC) {
			size_t deflen = cd.numericvalue.size();
			for (size_t i = 0; i < n; i++) {
				if (deflen == 1)vec[i] = (T)cd.numericvalue[0];
				else            vec[i] = (T)cd.numericvalue[i];
			}
			return true;
		}

		bool status = file_read(cd, vec, n);

		if (iotype != IOType::ASCII) {//Don't flip if ASCII reader as it already does this - to be fixed
			for (size_t i = 0; i < vec.size(); i++) {
				if (cd.flip) { vec[i] = -vec[i]; }
				cd.applyoperator(vec[i]);
			}
		}
		return true;
	}

	template<typename T>
	bool file_read(const cFieldDefinition& cd, std::vector<T>& vec, const size_t n)
	{		
		if (iotype == IOType::ASCII) {
			return ((cASCIIInputManager*)this)->file_read(cd, vec, n);
		}
		return ((cNetCDFInputManager*)this)->file_read(cd, vec, n);
	}	
};

class cASCIIInputManager : public cInputManager {

private:
	cAsciiColumnFile AF;
	static bool contains_non_numeric_characters(const std::string& str, const size_t& startpos=0)
	{				
		const static std::string validchars = "0123456789.+-eE ,\t\r\n";
		size_t pos = str.find_first_not_of(validchars,startpos);
		if (pos == std::string::npos) return false;
		else return true;
	}

public:

	std::string HeaderFileName;
		
	cASCIIInputManager(const cBlock& b) {
		cInputManager::initialise(b);
		initialise(b);
	}

	~cASCIIInputManager() {	};

	void initialise(const cBlock& b)
	{		
		HeaderFileName = b.getstringvalue("DfnFile");
		fixseparator(HeaderFileName);		
		
		iotype = IOType::ASCII;
		AF.openfile(DataFileName);
		if (isdefined(HeaderFileName)) {
			glog.logmsg(0, "Parsing Input DfnFile %s\n", HeaderFileName.c_str());
			AF.parse_aseggdf2_header(HeaderFileName);
		}

		size_t headerlines = b.getsizetvalue("Headerlines");
		if (!isdefined(headerlines)) { headerlines = 0; }
		for (size_t k = 0; k < headerlines; k++) {
			AF.readnextrecord();
		}				
	}

	bool is_record_valid() {
		size_t startpos = AF.RT_string.size();
		bool nonnumeric = contains_non_numeric_characters(recordstring(),startpos);
		if (nonnumeric) {
			glog.logmsg("Skipping non-numeric record at line %zu of Input DataFile %s\n", record(), datafilename().c_str());
			glog.logmsg("\n%s\n\n", recordstring().c_str());
			return false;
		}
		return true;
	}
	
	bool readnextrecord()
	{
		bool status = true;
		if (AtStart == true) {
			status = AF.readnextrecord();
			if (status == false)return false;
			AtStart = false;
			Record = 0;
		}
		else {
			AF.skiprecords(Subsample - 1);
			status = AF.readnextrecord();
			if (status == false) return false;
			Record += Subsample;
		}		
		return true;
	}

	bool parserecord() {
		size_t n = AF.parserecord();
		if (n <= 1) return false;
		return true;		
	}
	
	const std::string& recordstring() const { return AF.currentrecord_string(); }

	const std::vector<std::string>& fields() const { return AF.currentrecord_columns(); }

	template<typename T>
	bool file_read(const cFieldDefinition& cd, std::vector<T>& vec, const size_t n)
	{		
		bool status = cd.getvalue(AF, vec, n);		
		return status;	
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

	bool readnextrecord()
	{
		bool status = true;		
		if (AtStart == true) {
			AtStart = false;
			Record = 0;
		}
		else {
			Record += Subsample;
			#if defined HAVE_NETCDF
				if (Record > NC.ntotalsamples()) return false;
			#endif
		}				
		return true;
	}

	template<typename T>
	bool file_read(const cFieldDefinition& cd, std::vector<T>& v, const size_t n)
	{
		NC.getDataByPointIndex(cd.varname, Record, v);
		return true;
	}

};
#endif

#endif

 