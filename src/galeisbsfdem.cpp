/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cmath>


#if defined _MPI_ENABLED
#include <mpi.h>
#endif
#include <cstring>

#include "gaaem_version.h"
#include "general_utils.h"
#include "general_types.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "vector_utils.h"
#include "fdemsystem.h"
#include "galeisbsfdem.h"
#include "logger.h"

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace


int process(int argc, char** argv, size_t Size, size_t Rank)
{
	if (argc >= 2){
		glog.logmsg(0,"Executing %s %s\n", argv[0], argv[1]);
		glog.logmsg(0,"Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
		glog.logmsg(0,"Working directory %s\n", getcurrentdirectory().c_str());
	}
	else{
		glog.errormsg(_SRC_,"Not enough input arguments\n");
	}

	std::string controlfile(argv[1]);
	
	FDEmSampleInverter I(controlfile,Size,Rank);	
	size_t record = 0;
	while (I.readnextrecord()){		
		record++;
		if ((record - 1) % Size != Rank)continue;	
		I.parserecord();
		I.invert();
	}
	I.finish();
	return 0;
};

int main(int argc, char** argv) {
	std::string stmfile = "C:/Users/rossc/code/repos/ga-aem-eigen/examples/resolve/stmfiles/Resolve-2006.stm";	
	
	std::vector<double> c = { 0.1, 1.0, 0.05 };
	std::vector<double> t = { 10, 5 };
	cEarth1D e(c,t);
	cFDEmSystem F;
	F.readsystemdescriptorfile(stmfile);	
	
	cFDEmGeometry g(30,0,10,0);	
	F.setgeometry(g);
	F.setearth(e);	
	F.setupcomputations();
	cvector fm = F.ppms();

	//cvector dba = F.dppms(cLayeredEarthModeller::CalculationType::DB, 0);
	//double delta = g.Height * 0.001;
	//g.Height -= delta / 2.0;

	int dlayer = 2;
	cvector dba = F.dppms(cLayeredEarthModeller::CalculationType::DC, dlayer);
	double delta = e.conductivity[dlayer] * 0.001;
	e.conductivity[dlayer] -= delta / 2.0;

	/*int dlayer = 2;
	cvector dba = F.dppms(cLayeredEarthModeller::CalculationType::DT, dlayer);
	double delta = e.thickness[dlayer] * 0.001;
	e.thickness[dlayer] -= delta / 2.0;*/
	
	F.setgeometry(g);
	F.setearth(e);
	F.setupcomputations();
	fm = F.ppms();

	//g.Height += delta;
	e.conductivity[dlayer] += delta;
	//e.thickness[dlayer] += delta;
	F.setgeometry(g);
	F.setearth(e);
	F.setupcomputations();
	cvector fm1 = F.ppms();
	
	for (auto i = 0; i < fm.size(); i++) {
		std::cout << fm[i] << " " << fm1[i] << std::endl;
	}		

	cvector dbn = (fm1 - fm) / delta;
	for (auto i = 0; i < fm.size(); i++) {
		std::cout << dbn[i] << " " << dba[i] << std::endl;
	}
}

int main1(int argc, char** argv)
{
#ifdef _MPI_ENABLED
	int Size, Rank;
	if (argc<2){
		printf("Executing %s\n", argv[0]);
		printf("Usage: %s control_file_name\n", argv[0]);
		printf("Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
		return 0;
	}
	glog.logmsg("Starting MPI\n");
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	int len;
	char processor_name[MPI_MAX_PROCESSOR_NAME + 1];
	MPI_Get_processor_name(processor_name, &len);
	processor_name[len] = '\0';
	glog.logmsg("MPI Started Processes=%d\tRank=%d\tProcessor name = %s\n", Size, Rank, processor_name);
	process(argc, argv, (size_t)Size, (size_t)Rank);
	MPI_Finalize();
#elif defined _OPENMP
	int Size;
	if (argc<2){
		printf("Executing %s\n", argv[0]);
		printf("Usage: %s control_file_name number_of_threads\n", argv[0]);
		printf("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		return 0;
	}
	if (argc == 2){
		Size = 1;
	}
	else{
		Size = atoi(argv[2]);
		int maxthreads = omp_get_max_threads();
		if (Size<1){
			printf("Error: %d is a stupid number of threads\n", Size);
			exit(1);
		}
		if (Size >= maxthreads){
			printf("Warning: The number of requested threads (%d) is more than the processors available (%d)\n", Size, maxthreads);
		}
	}

	omp_init_lock(&fftw_thread_lock);
#pragma omp parallel num_threads(Size)
	{
		int Rank = omp_get_thread_num();
		glog.logmsg("OpenMP threading Processes=%d\tRank=%d\n", Size, Rank);
		process(argc, argv, Size, Rank);
	}

#else
	int Size = 1;
	int Rank = 0;
	if (argc<2){
		printf("Executing %s\n", argv[0]);
		printf("Usage: %s control_file_name number_of_threads\n", argv[0]);
		printf("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		return 0;
	}
	glog.logmsg("Stand-alone Processes=%d\tRank=%d\tProcessor name = %s\n", Size, Rank);
	process(argc, argv, Size, Rank);
#endif	
	return 0;
}


