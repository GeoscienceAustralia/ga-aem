/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include "tdemsystem.h"

#include <cstdio>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>

#include "file_utils.h"
#include "file_formats.h"
#include "lem.h"
#include "tdemsystem.h"
#include "matrix_ops.h"
#include "vector_utils.h"
#include "galeisbstdem.h"
#include "stacktrace.h"
#include "large_loop.h"
#include "logger.h"

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

#if defined _MPI_ENABLED
	#include "mpi_wrapper.h"
#endif

#if defined _OPENMP	
	//This thread lock must be set when fftw is being initialised
	omp_lock_t fftw_thread_lock;
#endif

int main(int argc, char** argv)
{
	int exitstatus;
	int mpisize = 1;
	int mpirank = 0;
	std::string mpipname = "No MPI - Standalone";

#if defined _MPI_ENABLED			
	MPI_Init(&argc, &argv);
	mpirank = cMpiEnv::world_rank();
	mpisize = cMpiEnv::world_size();
	mpipname = cMpiEnv::processor_name();
	if (mpirank == 0)printf("MPI Started Processes=%d\tRank=%d\tProcessor name = %s\n", mpisize, mpirank, mpipname.c_str());
#endif
	if (mpirank == 0) {
		printf("%s\n", commandlinestring(argc, argv).c_str());
		printf("%s\n", versionstring(GAAEM_VERSION, __TIME__, __DATE__).c_str());
	}

	if (argc < 2) {
		printf("Usage: %s control_file_name [number_of_openmp_threads]\n", argv[0]);
		exitstatus = EXIT_FAILURE;
	}
	else if (argc == 2) {
		bool usingopenmp = false;
		std::string controlfile(argv[1]);
		cSBSInverter I(controlfile, mpisize, mpirank, usingopenmp);		
		exitstatus = EXIT_SUCCESS;
	}
	else if (argc == 3) {
		int openmpsize = atoi(argv[2]);
		if (mpisize > 1) {
			if (mpirank == 0) {
				printf("**Error: You may not use OpenMP with MPI\n");
				printf("**       Do not use [number_of_openmp_threads] when launched with mpiexec or mpirun\n");
			}
			exitstatus = EXIT_FAILURE;
		}
		else if (openmpsize < 1) {
			printf("**Error: %d is a silly number of threads\n", openmpsize);
			exitstatus = EXIT_FAILURE;
		}
		else {
#if defined _OPENMP				
			bool usingopenmp = true;
			std::string controlfile = std::string(argv[1]);

			int openmpmaxthreads = omp_get_max_threads();
			if (openmpsize > openmpmaxthreads) {
				printf("**Warning: The number of requested threads (%d) is more than the processors available (%d)\n", openmpsize, openmpmaxthreads);
			}

			printf("OpenMP threading Processes=%d\n", openmpsize);
			omp_init_lock(&fftw_thread_lock);
			glog.set_num_omp_threads(openmpsize);

#pragma omp parallel num_threads(openmpsize)
			{
				int openmprank = omp_get_thread_num();
				cSBSInverter I(controlfile, openmpsize, openmprank, usingopenmp);				
			}
			exitstatus = EXIT_SUCCESS;
#else
			printf("**Error: This executable has not been compiled with OpenMP enabbled\n");
			printf("**       Compile with OpenMP or do not specify [number_of_openmp_threads]\n");
			exitstatus = EXIT_FAILURE;
#endif
		}
	}
	else {
		printf("Usage: %s control_file_name [number_of_openmp_threads]\n", argv[0]);
		printf("Too many command line arguments\n");
		exitstatus = EXIT_FAILURE;
	}

#if defined _MPI_ENABLED
	if (mpirank == 0)printf("Finalizing MPI\n");
	MPI_Finalize();
#endif	
	return exitstatus;
}