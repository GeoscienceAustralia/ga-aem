/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _tdemsysteminfo_H
#define _tdemsysteminfo_H

#include "tdemcomponentinfo.h"
#include "inputmanager.h"

class cTDEmSystemInfo {

public:
	const size_t XCOMP = 0;
	const size_t YCOMP = 1;
	const size_t ZCOMP = 2;
	const size_t XZAMP = 3;

	cTDEmSystem T;
	std::string SystemFile;
	size_t nwindows = 0;
	size_t ncomps = 0;
	size_t nchans = 0;
	cTDEmComponentInfo CompInfo[3];
	std::vector<cTDEmData> predicted;
	std::string units;
	
	bool invertXPlusZ = false;
	bool invertPrimaryPlusSecondary = false;
	bool reconstructPrimary = false;
	
	void initialise(const cBlock& b, const size_t nsoundings) {
		std::string dummy;
		if (b.getvalue("InvertTotalField", dummy)) {
			glog.errormsg(_SRC_, "InvertTotalField is no longer an option, use InvertPrimaryPlusSecondary instead\n");
		};

		std::string stmfile = b.getstringvalue("SystemFile");
		fixseparator(stmfile);

		glog.logmsg(0, "Reading system file %s\n", stmfile.c_str());
		T.readsystemdescriptorfile(stmfile);
		glog.log("==============System file %s\n", stmfile.c_str());
		glog.log(T.STM.get_as_string());
		glog.log("==========================================================================\n");
		nwindows = T.NumberOfWindows;

		invertXPlusZ = b.getboolvalue("InvertXPlusZ");
		invertPrimaryPlusSecondary = b.getboolvalue("InvertPrimaryPlusSecondary");
		reconstructPrimary = b.getboolvalue("ReconstructPrimaryFieldFromInputGeometry");

		CompInfo[XCOMP].initialise(b.findblock("XComponent"), "X", nwindows, nsoundings);
		CompInfo[YCOMP].initialise(b.findblock("YComponent"), "Y", nwindows, nsoundings);
		CompInfo[ZCOMP].initialise(b.findblock("ZComponent"), "Z", nwindows, nsoundings);

		ncomps = 0;
		if (CompInfo[XCOMP].Use) ncomps++;
		if (CompInfo[YCOMP].Use) ncomps++;
		if (CompInfo[ZCOMP].Use) ncomps++;

		if (invertXPlusZ) {
			CompInfo[XCOMP].Use = true;
			CompInfo[ZCOMP].Use = true;
		}
		nchans = nwindows * ncomps;
	}

	void set_units(cInputManager* IM) {
		for (size_t ci = 0; ci < 3; ci++) {
			cFieldDefinition& fd = CompInfo[ci].fdS;
			if (fd.varname.size() > 0) {
				cAsciiColumnField c;
				IM->get_acsiicolumnfield(fd, c);
				std::string u = c.get_att("units");
				if (units.size() > 0){
					if (ciequal(u,units) == false) {
						std::ostringstream msg;
						msg << "Error: units must be the same on all EM system/components. " << units << "does niot match " << u << "." << std::endl;
						glog.errormsg(_SRC_,msg.str().c_str());
					}
				}
				else units = u;
			}
		}
	}
	
};

#endif
