/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cassert>
#include <iostream>
#include <string>
#include <memory>
#include <filesystem>
#include <complex>

#include "string_print.hpp"
#include "logger.hpp"
#include "file_utils.hpp"
#include "general_utils.hpp"
#include "vector_utils.hpp"
#include "streamredirecter.hpp"
#include "gaaem_version.hpp"
#include "aem_coredefs.hpp"
#include "sbsinverter.hpp"
#include "spectralaemsystem.hpp"
#include "tdemsystem.hpp"
#include "inverter.hpp"

class cLogger glog; //The global instance of the log file manager

#ifdef ENABLE_MPI
	#include "mpi_wrapper.hpp"
#endif

#ifdef _OPENMP
	#include <omp.h>
#endif


using namespace AEM;
using namespace AEM::INVERTER::SBSINVERTER;

void finalise() {
#ifdef ENABLE_MPI
	glog.logmsg(0, "Finalizing MPI\n");
	cMpiEnv::stop();
#endif
};

int finaliseandexit() {
	finalise();
	return EXIT_FAILURE;
};

fs::path get_warning_log_path() {
	std::string s = "warning.log";
	int k = 1;
	do {
		if (fs::exists(s)) {
			std::error_code ec;
			bool status = std::filesystem::remove(s, ec);
			if (status == false) {
				// Path must be being used by another process so increment the suffix
				s = strprint("warning_%d.log", k);
				k++;
			}
			else return fs::path(s);
		}
		else return fs::path(s);
	} while(true);
};

int main(int argc, char** argv) {

	std::string commandline = commandlinestring(argc, argv);
	int size = 1;
	int rank = 0;
	bool usingopenmp = false;
	//int openmpsize = 1;
	fs::path controlfile;
	std::string mpipname = "No MPI - Standalone";

	#ifdef ENABLE_MPI
		cMpiEnv::start(argc, argv);
		rank = cMpiEnv::world_rank();
		size = cMpiEnv::world_size();
		mpipname = cMpiEnv::processor_name();
		//glog.logmsg(0, "%s\n", commandline.c_str());
		//glog.logmsg(0, "%s\n", versionstring(GAAEM_VERSION, __TIME__, __DATE__).c_str());
		//glog.logmsg(0, "MPI Started Processes=%d\tRank=%d\tProcessor name = %s\n", mpisize, mpirank, mpipname.c_str());
	#endif

	fs::path wlogpath;
	if (rank == 0) {
		wlogpath = get_warning_log_path();
	};

	#ifdef ENABLE_MPI
		cMpiEnv::world_barrier();
	#endif

	std::ofstream log(wlogpath, std::ios_base::app);
	cStreamRedirecter cerrredirect(log, std::cerr);
	if (rank == 0) std::cerr << "Warning log opening " << timestamp() << std::endl;

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
		controlfile = fs::path(argv[1]);
		usingopenmp = false;
	}
	else if (argc == 3 && size > 1) {
		glog.logmsg(0, "**Error: You may not use OpenMP with MPI\n");
		glog.logmsg(0, "**       Do not use [number_of_openmp_threads] when launched with mpiexec or mpirun\n");
		return finaliseandexit();
	}
	else if (argc == 3) {
		#if defined _OPENMP
			usingopenmp = true;
			controlfile = fs::path(argv[1]);
			size = atoi(argv[2]);
			glog.set_num_omp_threads(size);
			int openmpmaxthreads = omp_get_max_threads();
			if (size > openmpmaxthreads) {
				std::string msg = strprint("**Warning: The number of requested threads (%d) is more than the processors available (%d).\n", size, openmpmaxthreads);
				std::cerr << msg << std::endl;
				glog.logmsg(0, msg);
			}
			else if (size < 1) {
				glog.logmsg(0, "%d is a silly number of threads.\n", size);
				return finaliseandexit();
			}
		#elif 
			glog.logmsg(0, "Usage: %s control_file_name [number_of_openmp_threads]\n", argv[0]);
			glog.logmsg(0, "       **Error: This executable has not been compiled with OpenMP enabbled\n");
			glog.logmsg(0, "       **Compile with OpenMP or do not specify [number_of_openmp_threads]\n");
			return finaliseandexit();
		#endif
	}

	AEM::SystemType systype = AEM::aem_system_type(SBSINVERTER::get_stmpath(controlfile));

	if (usingopenmp) {
		#if defined _OPENMP
			if (systype == AEM::SystemType::SpectralTimeDomain) {
				#pragma omp parallel num_threads(size)
				{
					rank = omp_get_thread_num();
					SBSInverter<SpectralAEMSystem, std::complex<double>> I(controlfile, size, rank, usingopenmp, commandline);
				}
			}
			else{
				//For OpenMP the FFTW planning cannot be done in parallel
				omp_lock_t fftw_thread_lock;
				omp_init_lock(&fftw_thread_lock);
				#pragma omp parallel num_threads(size)
				{
					int openmprank = omp_get_thread_num();
					SBSInverter<TDEmSystem, double> I(controlfile, size, rank, usingopenmp, commandline, &fftw_thread_lock);
				}
			}
			std::cerr << "Warning log closing " << timestamp() << std::endl;
		#endif
	}
	else {
		if (systype == AEM::SystemType::SpectralTimeDomain) {
			SBSInverter<SpectralAEMSystem, std::complex<double>> I(controlfile, size, rank, usingopenmp, commandline);
		}
		else {
			SBSInverter<TDEmSystem, double> I(controlfile, size, rank, usingopenmp, commandline);
		}
		#ifdef ENABLE_MPI
			cMpiEnv::world_barrier();
		#endif
		if (rank == 0) std::cerr << "Warning log closing " << timestamp() << std::endl;
	}
	finalise();
	return EXIT_SUCCESS;
}

