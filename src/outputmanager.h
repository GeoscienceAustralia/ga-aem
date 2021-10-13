/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _outputmanager_H_
#define _outputmanager_H_

#include <list>
#include "asciicolumnfile.h"
#include "fielddefinition.h"


#if defined _MPI_ENABLED
	#include "mpi_wrapper.h"
#endif

#if defined HAVE_NETCDF
	#include "geophysics_netcdf.h"
#endif

class cASCIIOutputManager;
class cNetCDFOutputManager;
class cOutputManager;
class cOutputField {
	
	public:				
		std::shared_ptr<cGeophysicsVar> var;
		std::string name;//name
		std::string description;//description
		std::string units;//units	
		size_t bands = 0;//number of bands
		nc_type ncstoragetype ;//binary storage tyrpe
		std::string ncdimname;//dimension names		
		char notation = 0;//ascii notation form I, F, E
		size_t width = 0 ;//ascii width
		size_t precision = 0;//ascii number of decimals places
				
		cOutputField(			
			const std::string& _name,//name
			const std::string& _description,//description
			const std::string& _units,//units	
			const size_t& _bands,//number of bands
			const nc_type& _ncstoragetype,//binary storage tyrpe
			const std::string& _ncdimname,//dimension names		
			const char& _notation,//ascii notation I, F, E
			const size_t& _width,//ascii width
			const size_t& _precision//ascii number of decimals places			
			)
		{					
			name = _name;
			description = _description;
			units = _units;
			bands = _bands;
			ncstoragetype = _ncstoragetype;
			ncdimname = _ncdimname;
			notation = _notation;
			width = _width;
			precision = _precision;			
		}	

		//~cOutputField() {};
};

class cOutputManager {
		
public:
	int Size = 1;
	int Rank = 0;
	enum class IOType { ASCII, NETCDF, NONE };

protected:		
	std::string DataFileName;	
	IOType iotype = IOType::NONE;	
	std::list<cOutputField> flist;

public:
	
	cOutputManager() {};

	virtual ~cOutputManager() {};

	void initialise(const cBlock& b, const int& size, const int& rank) {
		Size = size;
		Rank = rank;
		DataFileName = b.getstringvalue("DataFile");
		fixseparator(DataFileName);
		std::string suffix = stringvalue(Rank, ".%04d");
		DataFileName = insert_after_filename(DataFileName, suffix);		
	}

	static bool isnetcdf(const cBlock& b) {
		std::string fname = b.getstringvalue("DataFile");
		std::string ext = extractfileextension(fname);
		if (strcasecmp(ext, ".nc") == 0) {
			return true;
		}
		return false;
	}
	
	const std::string& datafilename() { return DataFileName; }

	virtual bool opendatafile(const std::string& srcfile, const size_t& subsample) = 0;
	
	std::shared_ptr<cOutputField> getfield(const std::string& name) {
		std::shared_ptr<cOutputField> f;
		for(auto it = flist.begin(); it != flist.end(); it++) {			
			if (strcasecmp(it->name, name) == 0) {
				f = std::make_shared<cOutputField>(*it);
				return f;
			}
		}
		return f;
	}

	virtual std::shared_ptr<cOutputField> addfield(
		const std::string& _name,//name
		const std::string& _description,//description
		const std::string& _units,//units	
		const size_t& _bands,//number of bands
		const nc_type& _ncstoragetype,//binary storage tyrpe
		const std::string& _ncdimname,//dimension names		
		const char& _fmtchar,//ascii form I, F, E
		const size_t& _width,//ascii width
		const size_t& _precision//ascii number of decimals places			
	) = 0;

	template <typename T> 
	bool writefield(
		const int& pointindex,//point index of sample in the file
		const T& vals,//values to be written
		const std::string& _name,//name
		const std::string& _description,//description
		const std::string& _units,//units	
		const size_t& _bands,//number of bands
		const nc_type& _ncstoragetype,//binary storage tyrpe
		const std::string& _ncdimname,//dimension names		
		const char& _fmtchar,//ascii form I, F, E, A
		const size_t& _width,//ascii width
		const size_t& _precision//ascii number of decimals places					
	) 
	{
		std::shared_ptr<cOutputField> f = getfield(_name);
		if (!f) {
			f = addfield(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _precision);
		}
		write(vals, f, pointindex);
		return true;
	}

	virtual void begin_point_output() = 0;
	virtual void end_point_output() = 0;
	virtual bool end_first_record() {
		return true;
	};

	virtual bool addvar(cOutputField& of) { return true; }

	//Cannot templatize these virtual functions
	virtual bool write(const int& val, std::shared_ptr<cOutputField> of, const int& pointindex) = 0;
	virtual bool write(const size_t& val, std::shared_ptr<cOutputField> of, const int& pointindex) = 0;
	virtual bool write(const float& val, std::shared_ptr<cOutputField> of, const int& pointindex) = 0;
	virtual bool write(const double& val, std::shared_ptr<cOutputField> of, const int& pointindex) = 0;

	virtual bool write(const std::vector<int>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) = 0;	
	virtual bool write(const std::vector<size_t>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) = 0;
	virtual bool write(const std::vector<float>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) = 0;
	virtual bool write(const std::vector<double>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) = 0;	
};

class cASCIIOutputManager : public cOutputManager {

private:
	std::ofstream fstrm;		
	std::ostringstream buf;
	cOutputFileInfo OI;	
	void set_formatflags(std::shared_ptr<cOutputField> of){
		if (of->notation == 'e' || of->notation == 'E') buf << std::scientific;
		else buf << std::fixed;
	}

	bool SaveDFNHeader = false;
	bool SaveCSVHeader = false;
	bool SaveHDRHeader = false;

public:
	
	cASCIIOutputManager(const cBlock& b, const int& size, const int& rank) {
		cOutputManager::initialise(b,size,rank);
		initialise(b);
	}

	~cASCIIOutputManager() {	};

	void initialise(const cBlock& b)
	{						
		SaveDFNHeader = b.getboolvalue("SaveDFNHeader");
		SaveCSVHeader = b.getboolvalue("SaveCSVHeader");
		SaveHDRHeader = b.getboolvalue("SaveHDRHeader");
	}
	
	bool opendatafile(const std::string& srcfile, const size_t& subsample) {	
		glog.logmsg(0,"Opening Output ASCII DataFile %s\n", DataFileName.c_str());				
		fstrm.open(DataFileName, std::ofstream::out);				
		return fstrm.is_open();
	};
	
	std::shared_ptr<cOutputField> addfield(
		const std::string& _name,//name
		const std::string& _description,//description
		const std::string& _units,//units	
		const size_t& _bands,//number of bands
		const nc_type& _ncstoragetype,//binary storage tyrpe
		const std::string& _ncdimname,//dimension names		
		const char& _fmtchar,//ascii form I, F, E
		const size_t& _width,//ascii width
		const size_t& _precision//ascii number of decimals places			
	) {
		std::shared_ptr<cOutputField> f = std::make_shared<cOutputField>(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _precision);
		OI.addfield(_name, _fmtchar, _width, _precision, _bands);
		OI.setunits(_units);
		OI.setdescription(_description);
		return f;
	}

	void begin_point_output() {		
		buf.str(std::string());	//empty buffer
	};

	void end_point_output() {

		buf << std::endl; //Carriage return		
		fstrm << buf.str() << std::flush; // Write to file
	};

	bool end_first_record(){
		if(Rank == 0) {
			write_headers();
		}
		return true;
	}

	void write_headers(){		
		sFilePathParts fpp = getfilepathparts(datafilename());

		if (SaveDFNHeader) {
			std::string aseggdffile = fpp.directory + fpp.prefix + ".dfn";
			OI.write_aseggdf_header(aseggdffile);
		}

		if (SaveCSVHeader) {
			std::string csvfile = fpp.directory + fpp.prefix + ".csvh";
			OI.write_csv_header(csvfile);
		}

		if (SaveHDRHeader) {
			std::string hdrfile = fpp.directory + fpp.prefix + ".hdr";
			OI.write_simple_header(hdrfile);
		}		
	};
		
	bool write(const int& val, const std::shared_ptr<cOutputField> of, const int& pointindex) {		
		return write_scalar(val, of, pointindex);
	}
	bool write(const size_t& val, const std::shared_ptr<cOutputField> of, const int& pointindex) {
		return write_scalar(val, of, pointindex);
	}
	bool write(const float& val, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return write_scalar(val, of, pointindex);
	}
	bool write(const double& val, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return write_scalar(val, of, pointindex);
	}

	bool write(const std::vector<int>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {		
		return write_vector(vals, of, pointindex);
	}
	bool write(const std::vector<size_t>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return write_vector(vals, of, pointindex);
	}
	bool write(const std::vector<float>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return write_vector(vals, of, pointindex);
	}
	bool write(const std::vector<double>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return write_vector(vals, of, pointindex);
	}

	template <typename T> 
	bool write_scalar(const T& val, std::shared_ptr<cOutputField> of, const int& pointindex) {		
		set_formatflags(of);
		buf.width(of->width);
		buf.precision(of->precision);
		buf << val;		
		return true;
	}

	template <typename T>
	bool write_vector(const std::vector<T>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {		
		set_formatflags(of);
		for (size_t i = 0; i < vals.size(); i++) {
			buf.width(of->width);
			buf.precision(of->precision);
			buf << vals[i];
		}
		return true;
	}
	
};

#if defined HAVE_NETCDF
class cNetCDFOutputManager : public cOutputManager {

private:	
	cGeophysicsNcFile NC;

public:

	cNetCDFOutputManager(const cBlock& b, const int& size, const int& rank) {
		cOutputManager::initialise(b, size, rank);
		initialise(b);
	}

	~cNetCDFOutputManager() {	};

	void initialise(const cBlock& b)
	{						
		glog.logmsg(0, "Opening Output DataFile %s\n", DataFileName.c_str());		
		iotype = IOType::NETCDF;		
	}

	bool opendatafile(const std::string& srcfile, const size_t& subsample) {		
		//std::vector<std::string> include_varnames = { "proj_client", "flight", "flight_index", "longitude", "latitude", "emz_nonhprg" };
		//std::vector<std::string> exclude_varnames = { "easting", "northing" };
		std::vector<std::string> include_varnames;
		std::vector<std::string> exclude_varnames;

		
		int rank = 0;
#if defined _MPI_ENABLED
		rank = cMpiEnv::world_rank();
#endif
		if (rank == 0){
			glog.logmsg(0, "Creating Output NetCDF DataFile %s\n", DataFileName.c_str());
			cGeophysicsNcFile inncfile(srcfile, NcFile::FileMode::read);			
			cGeophysicsNcFile outncfile(datafilename(),NcFile::FileMode::replace);
			outncfile.subsample(srcfile, subsample, include_varnames, exclude_varnames);
		}	

#if defined _MPI_ENABLED
		cMpiEnv::world_barrier();
#endif

		glog.logmsg(0, "Opening Output NetCDF DataFile %s\n", DataFileName.c_str());
		NC.open(datafilename(), NcFile::FileMode::write);
		return true;
	};

	std::shared_ptr<cOutputField> addfield(
		const std::string& _name,//name
		const std::string& _description,//description
		const std::string& _units,//units	
		const size_t& _bands,//number of bands
		const nc_type& _ncstoragetype,//binary storage tyrpe
		const std::string& _ncdimname,//dimension names		
		const char& _fmtchar,//ascii form I, F, E
		const size_t& _width,//ascii width
		const size_t& _precision//ascii number of decimals places			
	) {	
		std::shared_ptr<cOutputField> of = getfield(_name);
		if (!of) {
			flist.push_back(cOutputField(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _precision));
			addvar(flist.back());
			of = std::make_shared<cOutputField>(flist.back());			
		}		
		return of;
	}

	bool writefield(
		const std::string& _name,//name
		const std::string& _description,//description
		const std::string& _units,//units	
		const size_t& _bands,//number of bands
		const nc_type& _ncstoragetype,//binary storage tyrpe
		const std::string& _ncdimname,//dimension names		
		const char& _fmtchar,//ascii form I, F, E
		const size_t& _width,//ascii width
		const size_t& _precision,//ascii number of decimals places			
		const int& pointindex,
		const size_t& vals
	) {
		std::shared_ptr<cOutputField> of = getfield(_name);		
		if (!of) {
			of = addfield(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _precision);
		}
		write(of, pointindex, vals);
		return true;
	}
	
	void set_manager(cOutputManager* m) {

	};

	void begin_point_output() { };

	void end_point_output() { };

	bool addvar(cOutputField& of) {
		if (NC.hasVar(of.name) == false) {
			NcDim dim;
			if (of.ncdimname.size()) {
				dim = NC.addDim(of.ncdimname, of.bands);
			}

			of.var = std::make_shared<cSampleVar>(NC.addSampleVar(of.name, of.ncstoragetype, dim));
			of.var->add_description(of.description);
			of.var->add_units(of.units);
		}
		else {
			if (NC.isLineVar(NC.getVar(of.name))) {
				cLineVar v  = NC.getLineVar(of.name);
				of.var = std::make_shared<cLineVar>(v);				
			}
			else if (NC.isSampleVar(NC.getVar(of.name))) {
				cSampleVar v = NC.getSampleVar(of.name);
				of.var = std::make_shared<cSampleVar>(v);
			}
				
		}
		return true;
	}
	
	template <typename T>
	bool write(std::shared_ptr<cOutputField> of, const int& pointindex, const T& val) {
		return of->var->putRecord(pointindex, val);
	}

	virtual bool write(const int& val, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return of->var->putRecord(pointindex, val);		
	}
	virtual bool write(const size_t& val, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return of->var->putRecord(pointindex, val);
	}
	virtual bool write(const float& val, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return of->var->putRecord(pointindex, val);
	}
	virtual bool write(const double& val, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return of->var->putRecord(pointindex, val);
	}

	virtual bool write(const std::vector<int>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return of->var->putRecord(pointindex, vals);		
	}
	virtual bool write(const std::vector<size_t>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return of->var->putRecord(pointindex, vals);
	}
	virtual bool write(const std::vector<float>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return of->var->putRecord(pointindex, vals);
	}
	virtual bool write(const std::vector<double>& vals, std::shared_ptr<cOutputField> of, const int& pointindex) {
		return of->var->putRecord(pointindex, vals);		
	}
		
};
#endif

#endif

 