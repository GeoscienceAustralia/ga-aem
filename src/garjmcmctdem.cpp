/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <vector>
#include <exception>
using namespace std;

#include "gaaem_version.h"
#include "general_utils.h"
#include "file_utils.h"
#include "rjmcmc1dtdeminverter.h"
#include "logger.h"
#include "stacktrace.h"

#if defined _MPI_ENABLED
	#include "mpi_wrapper.h"
#endif

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

int main(int argc, char* argv[])
{
	int exitstatus;
	int mpisize = 1;
	int mpirank = 0;
	std::string mpipname = "Standalone";
	#if defined _MPI_ENABLED
		MPI_Init(&argc, &argv);
		mpirank = cMpiEnv::world_rank();
		mpisize = cMpiEnv::world_size();
		mpipname = cMpiEnv::processor_name();
		if (mpirank == 0)printf("MPI Started Processes=%d\tRank=%d\tProcessor name = %s\n", mpisize, mpirank, mpipname.c_str());
	#endif

	if (mpirank == 0){
		printf("%s\n", commandlinestring(argc, argv).c_str());
		printf("%s\n", versionstring(GAAEM_VERSION, __TIME__, __DATE__).c_str());
	}

	if(argc!=2){
		printf("Usage: %s control_file_name\n",argv[0]);
		exitstatus = EXIT_FAILURE;
	}
	else{
		glog.logmsg(0, "Executing %s %s\n", argv[0], argv[1]);
		glog.logmsg(0, "Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
		glog.logmsg(0, "Working directory %s\n", getcurrentdirectory().c_str());
		glog.logmsg(0, "This process starting at %s\n", timestamp().c_str());

		std::string executable = string(argv[0]);
		std::string controlfile = string(argv[1]);

		glog.logmsg(0, "Reading control file %s\n", controlfile.c_str());
		try {
			rjmcmc1dTDEmInverter I(executable, controlfile, (size_t)mpisize, (size_t)mpirank);

			glog.logmsg(0,"Starting Inversion\n");
			while (I.readnextrecord_thisprocess()){
				I.parsecurrentrecord();
				I.sample();
				double stime = I.samplingtime;
				double norm_mfit = I.LowestMisfit.get_chi2() / double(I.ndata);
				glog.logmsg("Rec %6lu\t %3lu\t %5lu\t %10lf nmfit=%.1lf stime=%.3lfs\n", I.CurrentRecord, I.flightnumber, I.linenumber, I.fidnumber, norm_mfit, stime);
			}
			glog.logmsg("This process finishing at %s\n", timestamp().c_str());
			exitstatus = EXIT_SUCCESS;
		}
		catch (std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Unable to complete inversion. Exiting..." << std::endl;
			exitstatus = EXIT_FAILURE;
		}
	}

	#if defined _MPI_ENABLED
		MPI_Finalize();
	#endif

	return exitstatus;
}
