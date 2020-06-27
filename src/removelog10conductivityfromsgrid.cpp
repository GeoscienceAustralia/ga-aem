/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>

#include "gaaem_version.h"
#include "general_utils.h"
#include "file_utils.h"

class cLogger glog; //The global instance of the log file manager
class cStackTrace gtrace; //The global instance of the stacktrace

int main(int argc, char** argv)
{	
	if (argc == 3){
		glog.logmsg("Executing %s %s\n", argv[0], argv[1]);
		glog.logmsg("Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
		glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
	}
	else{
		glog.logmsg("Executing %s\n", argv[0]);
		glog.logmsg("Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
		glog.logmsg("Working directory %s\n", getcurrentdirectory().c_str());
		glog.logmsg("Error: Weong number of input arguments\n");
		glog.logmsg("Usage: %s input_sgrid_dir output_sgrid_dir\n",argv[0]);		
		return 0;
	}

	double t1 = gettime();
	std::string indir  = argv[1];
	std::string outdir = argv[2];
	std::vector<std::string> flist = getfilelist(indir, "sg");
	for (size_t i = 0; i < flist.size(); i++) {
		std::string isg = flist[i];
		std::string osg = outdir + pathseparatorstring() + extractfilename(isg);
		std::string isd = isg + ".data";
		std::string osd = outdir + pathseparatorstring() + extractfilename(isd);		
		std::printf("Converting %zu of %zu: %s\n",i+1,flist.size(),extractfilename(osg).c_str());
		std::string s;
		bool status;

		//Header file
		FILE* isgfp = fileopen(isg, "r");
		FILE* osgfp = fileopen(osg, "w");
		while ((status = filegetline(isgfp, s))) {			
			if(strcmp(s.c_str(),"*painted*variable:Log10Conductivity") == 0){
				std::fprintf(osgfp,"*painted*variable:Conductivity\n");
			}
			else if(strcmp(s.c_str(),"PROPERTY 2 Log10Conductivity") == 0) continue;
			else if(strcmp(s.c_str(),"PROP_UNIT 2 log10S/m") == 0) continue;
			else if(strcmp(s.c_str(),"PROP_NO_DATA_VALUE 2 -999") == 0) continue;						
			else {
				std::fprintf(osgfp, "%s\n", s.c_str());
			}
		}
		fclose(isgfp);
		fclose(osgfp);

		//Data file
		FILE* isdfp = fileopen(isd, "r");
		FILE* osdfp = fileopen(osd, "w");				
		status = filegetline(isdfp, s);
		std::fprintf(osdfp, "%s\n",s.c_str());
		
		status = filegetline(isdfp, s);
		auto t = tokenize(s);		
		if (t.size() != 9) {
			std::printf("Error unknown format: number of items in line 2 of %s is not 9",isd.c_str());
			exit(1);
		}
		if (strcmp(t[5].c_str(),"Log10Conductivity")!=0) {
			std::printf("Error unknown format: item 6 in line 2 of %s must be Log10Conductivity", isd.c_str());
			exit(1);
		}
		std::fprintf(osdfp, "*   X   Y   Z  Conductivity I   J   K\n");

		status = filegetline(isdfp, s);
		std::fprintf(osdfp, "%s\n", s.c_str());

		while ((status = filegetline(isdfp, s))) {
			t = tokenize(s);
			std::fprintf(osdfp, "%s %s %s %s %s %s %s\n",t[0].c_str(), t[1].c_str(), t[2].c_str(), t[3].c_str(), t[5].c_str(), t[6].c_str(), t[7].c_str());				
		}
		fclose(isdfp);
		fclose(osdfp);
	}
	
	double t2 = gettime();
	std::printf("Done ... Elapsed time = %.2lf seconds\n", t2 - t1);
	//prompttocontinue();
	return 0;
}
