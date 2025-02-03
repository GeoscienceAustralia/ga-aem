/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <list>
#include <iterator>
#include <optional>
#include "asciicolumnfile.hpp"
#include "fielddefinition.hpp"


#ifdef ENABLE_MPI
#include "mpi_wrapper.hpp"
#endif

#ifdef HAVE_NETCDF
#include "geophysics_netcdf.hpp"
using namespace GeophysicsNetCDF;
#endif

//template<class T>
//std::ostream& operator<<(std::ostream& stream, const std::vector<T>& values) {
//	std::copy(begin(values), end(values), std::ostream_iterator<T>(stream, ""));
//	//c++20 std::ranges::copy(values, std::ostream_iterator<T>(stream, ""));	
//	return stream;
//}

namespace IOManager {

	enum class BinaryStorageType { FLOAT, DOUBLE, INT, UINT };//binary storagetype
	auto constexpr ST_INT = BinaryStorageType::INT;
	auto constexpr ST_UINT = BinaryStorageType::UINT;
	auto constexpr ST_FLOAT = BinaryStorageType::FLOAT;
	auto constexpr ST_DOUBLE = BinaryStorageType::DOUBLE;

	class cOutputField {

	public:
		//Base
		std::string name;//name		
		size_t nbands = 0;//number of bands

		//Atts
		cKeyVecCiStr atts;

		//Netcdf
#ifdef HAVE_NETCDF
		std::shared_ptr<GVar> var;
		nc_type nctype() const {
			if (btype == BinaryStorageType::FLOAT) return NC_FLOAT;
			else if (btype == BinaryStorageType::DOUBLE) return NC_DOUBLE;
			else if (btype == BinaryStorageType::INT) return NC_INT;
			else if (btype == BinaryStorageType::UINT) return NC_UINT;
			else {
				//Won't get here but silence the compiler warning
				glog.errormsg(_SRC_, "Unknown binary storage tyep\n");
			}
			return NC_FLOAT;
		}
#endif

		// Binary
		BinaryStorageType btype = BinaryStorageType::DOUBLE;// NC_DOUBLE;//binary storage tyrpe		
		std::string ncdimname;//dimension names		

		//Ascii
		cAsciiColumnField asciicolumn;

	public:

		cOutputField() {};

		cOutputField(const std::string& _name, const std::string& _description, const std::string& _units, const size_t& _bands, const cOutputField::BinaryStorageType& _ncstoragetype, const std::string& _ncdimname, const char& _fmtchar, const size_t& _width, const size_t& _decimals) {
			initialise(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _decimals);
		};

		cOutputField(const std::string& _name, const std::string& _description, const std::string& _units, const size_t& _bands, const cOutputField::BinaryStorageType& _ncstoragetype, const std::string& _ncdimname, const cAsciiColumnFormat& fmt) {
			initialise(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, fmt.fmtchar, fmt.width, fmt.decimals);
		};

		cOutputField(const cAsciiColumnField& c) {
			initialise(c);
		};

		void initialise(
			const std::string& _name,//name
			const std::string& _description,//description
			const std::string& _units,//units	
			const size_t& _nbands,//number of bands
			const BinaryStorageType& _ncstoragetype,//binary storage tyrpe
			const std::string& _ncdimname,//dimension names		
			const char& _fmtchar,//ascii notation I, F, E
			const size_t& _width,//ascii width
			const size_t& _decimals//ascii number of decimals places			
		)
		{
			name = _name;
			nbands = _nbands;
			atts.add(cAsciiColumnField::DESC, _description);
			atts.add(cAsciiColumnField::UNITS, _units);
			btype = _ncstoragetype;
			ncdimname = _ncdimname;

			asciicolumn.name = name;
			asciicolumn.nbands = nbands;
			asciicolumn.atts = atts;
			asciicolumn.fmtchar = _fmtchar;
			asciicolumn.width = _width;
			asciicolumn.decimals = _decimals;
		}

		void initialise(const cAsciiColumnField& c) {
			name = c.name;
			atts = c.atts;
			nbands = c.nbands;

			asciicolumn.name = name;
			asciicolumn.nbands = nbands;
			asciicolumn.atts = atts;
			asciicolumn.fmtchar = c.fmtchar;
			asciicolumn.width = c.width;
			asciicolumn.decimals = c.decimals;

			btype = IOManager::BinaryStorageType::DOUBLE;
			ncdimname = DN_NONE;
		}

	};

	class cOutputManager {

	protected:
		bool firstpointwritten = false;

	public:
		int Size = 1;
		int Rank = 0;
		enum class IOType { ASCII, NETCDF, NONE };

	protected:
		using spcOutputField = std::shared_ptr<cOutputField>;

		std::string DataFileName;
		IOType iotype = IOType::NONE;
		std::list<spcOutputField> flist;

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

		spcOutputField getfield(const std::string& name) {
			spcOutputField f;
			for (auto it = flist.begin(); it != flist.end(); it++) {
				if (strcasecmp((*it)->name, name) == 0) {
					return *it;
				}
			}
			return f;
		}

		virtual spcOutputField addfield(const cOutputField& f_uninitialised) = 0;

		virtual spcOutputField addfield(
			const std::string& _name,//name
			const std::string& _description,//description
			const std::string& _units,//units	
			const size_t& _bands,//number of bands
			const BinaryStorageType& _ncstoragetype,//binary storage tyrpe
			const std::string& _ncdimname,//dimension names		
			const char& _fmtchar,//ascii form I, F, E
			const size_t& _width,//ascii width
			const size_t& _decimals//ascii number of decimals places			
		) = 0;

		virtual bool writevrnt(const int& pointindex, const cVrnt& vrnt, const cAsciiColumnField& c) {
			std::string msg = func_nyi_msg(__FUNCTION__);
			glog.errormsg(_SRC_, msg);
			return false;
		};

		template <typename T>
		bool writefield(const int& pointindex, const T& vals, const cOutputField& of) {
			spcOutputField f;
			if (firstpointwritten == false) {
				f = addfield(of);
			}
			else {
				f = getfield(of.name);
			}
			write(vals, f, pointindex);
			return true;
		}

		template <typename T>
		bool writefield(
			const int& pointindex,//point index of sample in the file
			const T& vals,//values to be written
			const std::string& _name,//name
			const std::string& _description,//description
			const std::string& _units,//units	
			const size_t& _bands,//number of bands
			const BinaryStorageType& _ncstoragetype,//binary storage tyrpe
			const std::string& _ncdimname,//dimension names		
			const char& _fmtchar,//ascii form I, F, E, A
			const size_t& _width,//ascii width
			const size_t& _decimals//ascii number of decimals places					
		)
		{
			spcOutputField f;
			if (firstpointwritten == false) {
				if (!f) {
					f = addfield(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _decimals);
				};
			}
			else {
				f = getfield(_name);
			}
			write(vals, f, pointindex);
			return true;
		}

		virtual void begin_point_output() = 0;
		virtual void end_point_output() = 0;
		virtual bool end_first_record() {
			return true;
		};

		//virtual bool addvar(cOutputField& of) { return true; }

		//Cannot templatize these virtual functions
		virtual bool write(const int& val, const spcOutputField& of, const int& pointindex) = 0;
		virtual bool write(const size_t& val, const spcOutputField& of, const int& pointindex) = 0;
		virtual bool write(const float& val, const spcOutputField& of, const int& pointindex) = 0;
		virtual bool write(const double& val, const spcOutputField& of, const int& pointindex) = 0;
		virtual bool write(const std::vector<int>& vals, const spcOutputField& of, const int& pointindex) = 0;
		virtual bool write(const std::vector<size_t>& vals, const spcOutputField& of, const int& pointindex) = 0;
		virtual bool write(const std::vector<float>& vals, const spcOutputField& of, const int& pointindex) = 0;
		virtual bool write(const std::vector<double>& vals, const spcOutputField& of, const int& pointindex) = 0;

		template<typename T>
		bool write(const T& val, const spcOutputField& of, const int& pointindex) {
			std::string msg = template_func_nyi_msg<T>(__FUNCTION__);
			glog.errormsg(_SRC_, msg);
			return false;
		}

		//bool write(const char& val, const spcOutputField& of, const int& pointindex) {
		//	std::string msg = template_func_nyi_msg<const char&>(__FUNCTION__);
		//	glog.errormsg(_SRC_, msg);
		//	return false;
		//}

		//bool write(const std::complex<double>& val, const spcOutputField& of, const int& pointindex) {
		//	std::string msg = template_func_nyi_msg<const std::complex<double>&>(__FUNCTION__);
		//	glog.errormsg(_SRC_, msg);
		//	return false;
		//}

		//bool write(const std::vector<std::complex<double>>& vals, const spcOutputField& of, const int& pointindex) {
		//	std::string msg = template_func_nyi_msg<const std::vector<std::complex<double>>&>(__FUNCTION__);
		//	glog.errormsg(_SRC_, msg);
		//	return false;
		//}

		//bool write(const std::vector<char>& vals, const spcOutputField& of, const int& pointindex) {
		//	std::string msg = template_func_nyi_msg<const std::vector<char>&>(__FUNCTION__);
		//	glog.errormsg(_SRC_, msg);
		//	return false;
		//}

		static std::vector<std::string> preferred_sort_order() {
			static const std::vector<std::string> porder = {
				cAsciiColumnField::NULLSTR,
				cAsciiColumnField::UNITS,
				cAsciiColumnField::DATUM,
				cAsciiColumnField::PROJECTION,
				cAsciiColumnField::LONGNAME,
				cAsciiColumnField::DESC
			};
			return porder;
		}

		void sort_field_atts() {
			const std::vector<std::string> porder = preferred_sort_order();
			for (const auto& f : flist) {
				//auto& c = f->atts;
				f->atts = f->atts.preferred_sort(porder);
			}
		};

	};

	class cASCIIOutputManager : public cOutputManager {

	private:
		std::ofstream filestream;
		std::ostringstream buffer;
		void set_formatflags(const spcOutputField& of) {
			if (of->asciicolumn.fmtchar == 'e' || of->asciicolumn.fmtchar == 'E') {
				buffer << std::scientific;
			}
			else {
				buffer << std::fixed;
			}
		}
		bool SaveDFNHeader = true;
		bool SaveCSVHeader = true;
		bool SaveHDRHeader = true;
		bool SaveI3Header = true;

	public:

		cASCIIOutputManager(const cBlock& b, const int& size, const int& rank) {
			cOutputManager::initialise(b, size, rank);
			initialise(b);
		}

		~cASCIIOutputManager() {};

		void initialise(const cBlock& b)
		{
			bool status;
			if (b.getvalue("SaveDFNHeader", status)) {
				SaveDFNHeader = status;
			};

			if (b.getvalue("SaveCSVHeader", status)) {
				SaveCSVHeader = status;
			};

			if (b.getvalue("SaveHDRHeader", status)) {
				SaveHDRHeader = status;
			};

			if (b.getvalue("SaveI3Header", status)) {
				SaveHDRHeader = status;
			};
		}

		bool opendatafile(const std::string& srcfile, const size_t& subsample) {
			glog.logmsg(0, "Opening Output ASCII DataFile %s\n", DataFileName.c_str());
			filestream.open(DataFileName, std::ofstream::out);
			return filestream.is_open();
		};

		spcOutputField add_smartptr(cOutputField of) {
			of.asciicolumn.name = of.name;
			of.asciicolumn.nbands = of.nbands;
			of.asciicolumn.atts = of.atts;

			of.asciicolumn.fileorder = flist.size();
			of.asciicolumn.startcolumn = 0;
			of.asciicolumn.startchar = 0;
			if (flist.size() > 0) {
				of.asciicolumn.startcolumn = 1 + flist.back()->asciicolumn.endcol();
				of.asciicolumn.startchar = 1 + flist.back()->asciicolumn.endchar();
			}
			flist.push_back(std::make_shared<cOutputField>(of));
			return flist.back();
		}

		spcOutputField addfield(const cOutputField& f_uninitialised) {
			spcOutputField sp = getfield(f_uninitialised.name);
			if (sp) {
				if (firstpointwritten == false) {//Already exists and on first record so must be duplicate
					glog.errormsg(_SRC_, "Conflict in output field names. Field %s has already been added to the output file\n", f_uninitialised.name.c_str());
				}
			}
			else {
				sp = add_smartptr(f_uninitialised);
			}
			return sp;
		};

		spcOutputField addfield(
			const std::string& _name,//name
			const std::string& _description,//description
			const std::string& _units,//units	
			const size_t& _bands,//number of bands
			const BinaryStorageType& _ncstoragetype,//binary storage tyrpe
			const std::string& _ncdimname,//dimension names		
			const char& _fmtchar,//ascii form I, F, E
			const size_t& _width,//ascii width
			const size_t& _decimals//ascii number of decimals places			
		) {
			//todo
			spcOutputField sp = getfield(_name);
			if (!sp) {
				cOutputField f(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _decimals);
				sp = add_smartptr(f);
			}
			else {
				//Already exists and on first record so must be duplicate
				if (firstpointwritten == false) {
					glog.errormsg(_SRC_, "Conflict in output field names. Field %s has already been added to the output file\n", _name.c_str());
				}
			}
			return sp;
		}

		virtual spcOutputField addfield(cAsciiColumnField c) {
			//todo
			spcOutputField sp = getfield(c.name);
			if (!sp) {
				cOutputField f(c);
				sp = add_smartptr(f);
			}
			else {
				//Already exists and on firts record so must be duplicate
				if (firstpointwritten == false) {
					glog.errormsg(_SRC_, "Conflict in output field names. Field %s has already been added to the output file\n", c.name.c_str());
				}
			}
			return sp;
		}

		void begin_point_output() {
			buffer.str(std::string());	//empty buffer
		};

		void end_point_output() {
			buffer << std::endl; //Carriage return		
			filestream << buffer.str() << std::flush; // Write to file
		};

		bool end_first_record() {
			if (Rank == 0) {
				write_headers();
			}
			firstpointwritten = true;
			return true;
		}

		void write_headers() {
			FilePathParts fpp(datafilename());
			sort_field_atts();
			if (SaveDFNHeader) {
				std::string aseggdffile = fpp.directory + fpp.stem + ".dfn";
				write_aseggdf_header(aseggdffile);
			}

			if (SaveCSVHeader) {
				std::string csvfile = fpp.directory + fpp.stem + ".csv";
				write_csv_header(csvfile);
			}

			if (SaveHDRHeader) {
				std::string hdrfile = fpp.directory + fpp.stem + ".hdr";
				write_simple_header(hdrfile);
			}

			if (SaveI3Header) {
				std::string i3file = fpp.directory + fpp.stem + ".i3";
				write_i3_header(i3file);
			}
		};


		bool write(const int& val, const spcOutputField& of, const int& pointindex) {
			return write_scalar(val, of, pointindex);
		}

		bool write(const size_t& val, const spcOutputField& of, const int& pointindex) {
			return write_scalar(val, of, pointindex);
		}

		bool write(const float& val, const spcOutputField& of, const int& pointindex) {
			return write_scalar(val, of, pointindex);
		}

		bool write(const double& val, const spcOutputField& of, const int& pointindex) {
			return write_scalar(val, of, pointindex);
		}

		bool write(const char& val, const spcOutputField& of, const int& pointindex) {
			return write_scalar(val, of, pointindex);
		}


		bool write(const std::vector<int>& vals, const spcOutputField& of, const int& pointindex) {
			return write_vector(vals, of, pointindex);
		}

		bool write(const std::vector<size_t>& vals, const spcOutputField& of, const int& pointindex) {
			return write_vector(vals, of, pointindex);
		}

		bool write(const std::vector<float>& vals, const spcOutputField& of, const int& pointindex) {
			return write_vector(vals, of, pointindex);
		}

		bool write(const std::vector<double>& vals, const spcOutputField& of, const int& pointindex) {
			return write_vector(vals, of, pointindex);
		}
		bool write(const std::vector<char>& vals, const spcOutputField& of, const int& pointindex) {
			return write_vector(vals, of, pointindex);
		}

		template <typename T>
		bool write_scalar(const T& val, const spcOutputField& f, const int& pointindex) {
			set_formatflags(f);
			const size_t& w = f->asciicolumn.width;
			const size_t& d = f->asciicolumn.decimals;
			buffer << std::setw(w) << std::setprecision(d) << val;
			return true;
		}

		template <typename T>
		bool write_vector(const std::vector<T>& vals, const spcOutputField& f, const int& pointindex) {
			set_formatflags(f);
			const size_t& w = f->asciicolumn.width;
			const size_t& d = f->asciicolumn.decimals;
			for (size_t i = 0; i < vals.size(); i++) {
				buffer << std::setw(w) << std::setprecision(d) << vals[i];
			}
			return true;
		}

		bool writevrnt(const int& pointindex, const cVrnt& vrnt, const cAsciiColumnField& c) {
			spcOutputField sp = getfield(c.name);
			if (!sp) {
				sp = addfield(c);
			}

			auto WriteVisitor = [&](auto& vals) {
				write(vals, sp, pointindex);
				};

			std::visit(WriteVisitor, vrnt);

			return true;
		}

		void write_aseggdf_header(const std::string pathname) {
			std::ofstream ofs(pathname);
			ofs << "DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A76" << std::endl;
			for (const auto& f : flist) {
				std::string s = f->asciicolumn.aseggdf_header_record();
				ofs << s;
			}
			ofs << "DEFN " << flist.size() + 1 << " ST=RECD,RT=; END DEFN" << std::endl;
		};

		void write_simple_header(const std::string pathname) {
			std::ofstream ofs(pathname);
			for (const auto& f : flist) {
				std::string s = f->asciicolumn.simple_header_record();
				ofs << s;
			}
		};

		void write_i3_header(const std::string pathname) {
			std::ofstream ofs(pathname);
			ofs << "[IMPORT ARCHIVE]" << std::endl;
			ofs << "FILEHEADER\t0" << std::endl;
			ofs << "RECORDFORM\tFIXED" << std::endl;
			for (const auto& f : flist) {
				std::string s = f->asciicolumn.i3_header_record();
				ofs << s;
			}
		};

		void write_csv_header(const std::string pathname) {

			cKeyVecCiStr v = collect_all_att_names();
			const std::vector<std::string> porder = preferred_sort_order();
			v = v.preferred_sort(porder);

			std::ofstream ofs(pathname);
			ofs << "Name,Bands,Format";
			for (size_t j = 0; j < v.size(); j++) {
				ofs << "," << v[j].first;
			}
			ofs << std::endl;

			for (const auto& xf : flist) {
				const cAsciiColumnField& c = xf->asciicolumn;
				ofs << c.name << ",";
				ofs << c.nbands << ",";
				ofs << c.fmtstr_single();
				for (size_t j = 0; j < v.size(); j++) {
					std::string& key = v[j].first;
					std::string val = c.get_att(key);
					ofs << ",";
					if (val.size() > 0) ofs << val;
				}
				ofs << std::endl;
			}
		};

		cKeyVecCiStr collect_all_att_names() {
			cKeyVecCiStr v;
			for (const auto& f : flist) {
				for (const auto& [key, value] : f->atts) {
					v.add(key, std::string());
				}
			}
			return v;
		}
	};

#if defined HAVE_NETCDF
	class cNetCDFOutputManager : public cOutputManager {

	private:
		GeophysicsNetCDF::GFile NC;

	public:

		cNetCDFOutputManager(const cBlock& b, const int& size, const int& rank) {
			cOutputManager::initialise(b, size, rank);
			initialise(b);
		}

		~cNetCDFOutputManager() {};

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
#ifdef ENABLE_MPI
			rank = cMpiEnv::world_rank();
#endif
			if (rank == 0) {
				glog.logmsg(0, "Creating Output NetCDF DataFile %s\n", DataFileName.c_str());
				GFile inncfile(srcfile, NcFile::FileMode::read);
				GFile outncfile(datafilename(), NcFile::FileMode::replace);
				outncfile.subsample(inncfile, subsample, include_varnames, exclude_varnames);
			}

#ifdef ENABLE_MPI
			cMpiEnv::world_barrier();
#endif

			glog.logmsg(0, "Opening Output NetCDF DataFile %s\n", DataFileName.c_str());
			NC.open(datafilename(), NcFile::FileMode::write);
			return true;
		};

		spcOutputField addfield(const cOutputField& f_uninitialised) {
			spcOutputField sp = getfield(f_uninitialised.name);
			if (sp) {
				if (firstpointwritten == false) {//Already exists and on first record so must be duplicate
					glog.errormsg(_SRC_, "Conflict in output field names. Field %s has already been added to the output file\n", f_uninitialised.name.c_str());
				}
			}
			else {
				cOutputField of = f_uninitialised;
				addvar(of);
				sp = std::make_shared<cOutputField>(of);
				flist.push_back(sp);
				return sp;
			}
			return sp;
		};


		spcOutputField addfield(
			const std::string& _name,//name
			const std::string& _description,//description
			const std::string& _units,//units	
			const size_t& _bands,//number of bands
			const BinaryStorageType& _ncstoragetype,//binary storage tyrpe
			const std::string& _ncdimname,//dimension names		
			const char& _fmtchar,//ascii form I, F, E
			const size_t& _width,//ascii width
			const size_t& _decimals//ascii number of decimals places			
		) {
			spcOutputField of = getfield(_name);
			if (!of) {
				cOutputField f(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _decimals);
				addvar(f);
				of = std::make_shared<cOutputField>(f);
				flist.push_back(of);
				return of;
			}
			return of;
		};

		bool writefield(
			const std::string& _name,//name
			const std::string& _description,//description
			const std::string& _units,//units	
			const size_t& _bands,//number of bands
			const BinaryStorageType& _ncstoragetype,//binary storage tyrpe
			const std::string& _ncdimname,//dimension names		
			const char& _fmtchar,//ascii form I, F, E
			const size_t& _width,//ascii width
			const size_t& _decimals,//ascii number of decimals places			
			const int& pointindex,
			const size_t& vals
		) {
			spcOutputField of = getfield(_name);
			if (!of) {
				of = addfield(_name, _description, _units, _bands, _ncstoragetype, _ncdimname, _fmtchar, _width, _decimals);
			}
			write(of, pointindex, vals);
			return true;
		}

		void set_manager(cOutputManager* m) {

		};

		void begin_point_output() {};

		void end_point_output() {};

		bool addvar(cOutputField& of) {
			if (NC.hasVar(of.name) == false) {
				NcDim dim;
				if (of.ncdimname.size()) {
					dim = NC.addDim(of.ncdimname, of.nbands);
				}

				of.var = std::make_shared<GVar>(NC.addgetSampleVar(of.name, of.nctype(), dim));

				for (const auto& [key, value] : of.atts) {
					of.var->add_attribute(key, value);
				}
			}
			else {
				if (NC.isLineVar(NC.getVar(of.name))) {
					GVar v = NC.getLineVar(of.name);
					of.var = std::make_shared<GVar>(v);
				}
				else if (NC.isSampleVar(NC.getVar(of.name))) {
					GVar v = NC.getSampleVar(of.name);
					of.var = std::make_shared<GVar>(v);
				}

			}
			return true;
		}

		template <typename T>
		bool write(spcOutputField of, const int& pointindex, const T& val) {
			return of->var->putRecord(pointindex, val);
		}

		virtual bool write(const int& val, const spcOutputField& of, const int& pointindex) {
			return of->var->putRecord(pointindex, val);
		}
		virtual bool write(const size_t& val, const spcOutputField& of, const int& pointindex) {
			return of->var->putRecord(pointindex, val);
		}
		virtual bool write(const float& val, const spcOutputField& of, const int& pointindex) {
			return of->var->putRecord(pointindex, val);
		}
		virtual bool write(const double& val, const spcOutputField& of, const int& pointindex) {
			return of->var->putRecord(pointindex, val);
		}

		virtual bool write(const std::vector<int>& vals, const spcOutputField& of, const int& pointindex) {
			return of->var->putRecord(pointindex, vals);
		}
		virtual bool write(const std::vector<size_t>& vals, const spcOutputField& of, const int& pointindex) {
			return of->var->putRecord(pointindex, vals);
		}
		virtual bool write(const std::vector<float>& vals, const spcOutputField& of, const int& pointindex) {
			return of->var->putRecord(pointindex, vals);
		}
		virtual bool write(const std::vector<double>& vals, const spcOutputField& of, const int& pointindex) {
			return of->var->putRecord(pointindex, vals);
		}
	};
#endif

};

