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
#include "cinverter.h"
#include "csbsinverter.h"
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

int main(int argc, char** argv) {	
		
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
		#elif 
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
				//cSBSInverter I(controlfile, openmpsize, openmprank, usingopenmp, commandline);				
			}
			std::cerr << "Warning log closing " << timestamp() << std::endl;
		#endif
	}
	else {								
		//cSBSInverter I(controlfile, mpisize, mpirank, usingopenmp, commandline);				
		std::unique_ptr<cInverter> I = std::make_unique<cSBSInverter>(controlfile, mpisize, mpirank, usingopenmp, commandline);
		//std::unique_ptr<cSBSInverterX> I = std::make_unique<cSBSInverterX>(controlfile, mpisize, mpirank, usingopenmp, commandline);
		#if defined _MPI_ENABLED
			cMpiEnv::world_barrier();
		#endif
		if (mpirank == 0) std::cerr << "Warning log closing " << timestamp() << std::endl;
	}
	
	finalise();
	return EXIT_SUCCESS;	
}
