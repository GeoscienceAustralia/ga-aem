/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <iostream>
#include <cstdio>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>

#include "gaaem_version.h"
#include "file_utils.h"
#include "tdemsystem.h"
#include "vector_utils.h"
#include "galeisbstdem.h"
#include "stacktrace.h"
#include "logger.h"
#include "streamredirecter.h"

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

#if defined _MPI_ENABLED
	#include "mpi_wrapper.h"
#endif

#if defined _OPENMP	
	//This thread lock must be set when fftw is being initialised
	omp_lock_t fftw_thread_lock;
#endif

void finalise() {
	#if defined _MPI_ENABLED				
		MPI_Finalize();
	#endif		
}

int finaliseandexit(){		
	finalise();	
	return EXIT_FAILURE;
}

int test() {
	try {
		std::string datfile = "test.dat";		
		std::string dfnfile = "test.dfn";
		//std::string datfile = "D:\\projects\\2021_Mundi_MinExCRC_GSNSW\\final\\line_data_bfield\\MinEx_CRC_Mundi_AEM_BField_sub10.dat";
		//std::string dfnfile = "D:\\projects\\2021_Mundi_MinExCRC_GSNSW\\final\\line_data_bfield\\MinEx_CRC_Mundi_AEM_BField.dfn";
		cAsciiColumnFile F;
		F.openfile(datfile);
		F.read_dfn(dfnfile);
		int linefieldindex = F.fieldindexbyname("Line");
		std::cout << "Line field index " << linefieldindex << std::endl;

		cStopWatch w;
		size_t fs = F.file_size();
		w.reportnow();
		std::cout << "File size " << fs << std::endl;

		w.reset();
		size_t rl = F.get_record_length();
		w.reportnow();
		std::cout << "Record length " << rl << std::endl;

		w.reset();
		size_t nr = F.nrecords();
		w.reportnow();
		std::cout << "Number of records " << nr << std::endl;

		w.reset();
		std::string rec = F.get_record(0);
		std::cout << rec << std::endl << std::endl;
		w.reportnow();

		w.reset();
		rec = F.get_record(nr - 1);
		std::cout << rec << std::endl << std::endl;
		w.reportnow();

		std::vector<unsigned int> listart;
		std::vector<unsigned int> licount;
		std::vector<unsigned int> lnum;

		w.reset();
		size_t nrec = F.scan_for_line_index(linefieldindex, listart, licount, lnum);
		w.reportnow();

		w.reset();
		std::vector<bool> gbf = F.scan_for_groupby_fields(licount);
		w.reportnow();

		std::cout << "Number of records " << nrec << std::endl;
		std::cout << "Number of flight lines " << licount.size() << std::endl;

	}
	catch (const std::runtime_error& e) {
		std::cerr << e.what();
	}
	return 0;
}


int main(int argc, char** argv) {	
	
	//return test();
	
	std::string commandline = commandlinestring(argc, argv);
	
	int mpisize = 1;
	int mpirank = 0;
	bool usingopenmp = false;
	int openmpsize = 1;	
	std::string controlfile;
	std::string mpipname = "No MPI - Standalone";

	#if defined _MPI_ENABLED		
		MPI_Init(&argc, &argv);
		mpirank = cMpiEnv::world_rank();
		mpisize = cMpiEnv::world_size();
		mpipname = cMpiEnv::processor_name();		
	#endif

	std::string wlogpath = "warning.log";
	if(mpirank==0) deletefile(wlogpath);
	#if defined _MPI_ENABLED
		cMpiEnv::world_barrier();
	#endif

	std::ofstream log(wlogpath, std::ios::app);
	cStreamRedirecter cerrredirect(log, std::cerr);
	if (mpirank == 0) std::cerr << "Warning log opening " << timestamp() << std::endl;
	
	if (argc < 2) {
		glog.logmsg(0, "Usage: %s control_file_name [number_of_openmp_threads]\n", argv[0]);
		glog.logmsg(0, "       Not enough command line arguments\n");
		return finaliseandexit();
	}
	else if (argc > 3) {
		glog.logmsg(0, "Usage: %s control_file_name [number_of_openmp_threads]\n", argv[0]);
		glog.logmsg(0, "       Too many command line arguments\n");
		return finaliseandexit();
	}
	else if (argc == 2) {
		controlfile = argv[1];
		usingopenmp = false;
	}
	else if (argc == 3 && mpisize > 1) {		
		glog.logmsg(0, "**Error: You may not use OpenMP with MPI\n");
		glog.logmsg(0, "**       Do not use [number_of_openmp_threads] when launched with mpiexec or mpirun\n");		
		return finaliseandexit();
	}
	else if (argc == 3) {			
		usingopenmp = true;
		openmpsize = atoi(argv[2]);
		#if defined _OPENMP			
			glog.set_num_omp_threads(openmpsize);
			int openmpmaxthreads = omp_get_max_threads();
			if (openmpsize > openmpmaxthreads) {
				std::string msg = strprint("**Warning: The number of requested threads (%d) is more than the processors available (%d)\n", openmpsize, openmpmaxthreads);
				std::cerr << msg << std::endl;
				glog.logmsg(0,msg);
			}
			else if(openmpsize < 1) {
				glog.logmsg(0, "**Error: %d is a silly number of threads\n", openmpsize);
				return finaliseandexit();
			}
		#else if 
		    glog.logmsg(0,"Usage: %s control_file_name [number_of_openmp_threads]\n", argv[0]);
			glog.logmsg(0,"       **Error: This executable has not been compiled with OpenMP enabbled\n");
			glog.logmsg(0,"       **Compile with OpenMP or do not specify [number_of_openmp_threads]\n");
			return finaliseandexit();
		#endif		
	}

	controlfile = std::string(argv[1]);	
	if (usingopenmp) {	
		#if defined _OPENMP			
			omp_init_lock(&fftw_thread_lock);		
			#pragma omp parallel num_threads(openmpsize)
			{			
				int openmprank = omp_get_thread_num();			
				cSBSInverter I(controlfile, openmpsize, openmprank, usingopenmp, commandline);
			}
			std::cerr << "Warning log closing " << timestamp() << std::endl;
		#endif
	}
	else {								
		cSBSInverter I(controlfile, mpisize, mpirank, usingopenmp, commandline);		
		#if defined _MPI_ENABLED
			cMpiEnv::world_barrier();
		#endif
		if (mpirank == 0) std::cerr << "Warning log closing " << timestamp() << std::endl;
	}
	
	finalise();
	return EXIT_SUCCESS;	
}