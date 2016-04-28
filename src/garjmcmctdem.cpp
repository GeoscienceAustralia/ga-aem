/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <vector>
using namespace std;

#include "general_utils.h"
#include "file_utils.h"
#include "random.h"
#include "tdemsystem.h"
#include "rjmcmc1d.h"
#include "rjmcmc1dtdeminverter.h"

#if defined _MPI_ENABLED
	#include "mpi.h"
#endif


int main(int argc, char* argv[])
{
	int exitstatus;
	int mpisize = 1;
	int mpirank = 0;
	std::string mpipname = "Standalone";
	#if defined _MPI_ENABLED			
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
		mpipname = mpi_processername();
		if (mpirank == 0)printf("MPI Started Processes=%d\tRank=%d\tProcessor name = %s\n", mpisize, mpirank, mpipname.c_str());
	#endif

	if (mpirank == 0){
		printf("%s\n", commandlinestring(argc, argv).c_str());
		printf("%s\n", versionstring(VERSION, __TIME__, __DATE__).c_str());
	}

	if(argc!=2){						
		printf("Usage: %s control_file_name\n",argv[0]);
		exitstatus = EXIT_FAILURE;
	}	
	else{
		rootmessage("Executing %s %s\n", argv[0], argv[1]);
		rootmessage("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		rootmessage("Working directory %s\n", getcurrentdirectory().c_str());
		rootmessage("This process starting at %s\n", timestamp().c_str());

		std::string executable = string(argv[0]);
		std::string controlfile = string(argv[1]);

		rootmessage("Reading control file %s\n", controlfile.c_str());
		rjmcmc1dTDEmInverter I(executable, controlfile, (size_t)mpisize, (size_t)mpirank);

		rootmessage(I.fp_log, "Starting Inversion\n");
		while (I.readnextrecord_thisprocess()){
			I.parsecurrentrecord();
			I.sample();
			double stime = I.samplingtime;
			double norm_mfit = I.mLowestMisfit.misfit() / double(I.ndata);
			message(I.fp_log, "Rec %6lu\t %3lu\t %5lu\t %10lf nmfit=%.1lf stime=%.3lfs\n", I.CurrentRecord, I.flightnumber, I.linenumber, I.fidnumber, norm_mfit, stime);
		}
		message(I.fp_log, "This process finishing at %s\n", timestamp().c_str());
		exitstatus = EXIT_SUCCESS;
	}

	#if defined _MPI_ENABLED
		MPI_Finalize();	
	#endif	

	return exitstatus;
}



