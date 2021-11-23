/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _tdemcomponentinfo_H
#define _tdemcomponentinfo_H


#include <vector>

#include "inputmanager.h"
#include "fielddefinition.h"

class cTDEmComponentInfo {

private:
	size_t nWindows;
	size_t nSoundings;

public:

	class SoundingData {
	public:
		double  P = 0.0;
		std::vector<double> S;
		std::vector<double> E;
	};
	std::vector<SoundingData> data;

	std::string Name;
	bool Use = false;
	cFieldDefinition fdP;
	cFieldDefinition fdS;
	cFieldDefinition fdE;
	bool EstimateNoiseFromModel = false;
	std::vector<double> mn;
	std::vector<double> an;
	
	cTDEmComponentInfo() {};
	void initialise(const cBlock& b, const std::string& name, const size_t& nwindows, const size_t& nsoundings) {
		Name = name;
		if (b.Entries.size() == 0) {
			Use = false;
			return;
		}
		Use = b.getboolvalue("Use");
		if (Use == false)return;

		EstimateNoiseFromModel = b.getboolvalue("EstimateNoiseFromModel");

		if (EstimateNoiseFromModel) {
			mn = b.getdoublevector("MultiplicativeNoise");
			an = b.getdoublevector("AdditiveNoise");
			if (an.size() == 1) {
				an = std::vector<double>(nwindows, an[0]);
			}
			else if (an.size() != nwindows) {
				glog.errormsg(_SRC_, "Must have exactly 1 or nwindows AdditiveNoise values\n");
			};

			if (mn.size() == 1) {
				mn = std::vector<double>(nwindows, mn[0]);
			}
			if (mn.size() != nwindows) {
				glog.errormsg(_SRC_, "Must have exactly 1 or nwindows MultiplicativeNoise values\n");
			}
		}

		fdP.initialise(b, "Primary");
		fdS.initialise(b, "Secondary");
		fdE.initialise(b, "Noise");

		nSoundings = nsoundings;
		nWindows = nwindows;
		data.resize(nSoundings);
		for (size_t si = 0; si < nSoundings; si++) {
			data[si].S.resize(nWindows);
			data[si].E.resize(nWindows);
		}
	}

	const size_t& nw() const
	{
		return nWindows;
	}

	void readdata(const std::unique_ptr<cInputManager>& IM, const size_t& soundingindex)
	{
		const size_t& si = soundingindex;
		if (Use == false) return;
		SoundingData& d = data[si];
		IM->read(fdP, d.P);
		IM->read(fdS, d.S, nw());
		if (EstimateNoiseFromModel) {
			for (size_t wi = 0; wi < nWindows; wi++) {
				const double v = 0.01 * mn[wi] * d.S[wi];
				d.E[wi] = std::hypot(an[wi], v);
			}
		}
		else {
			IM->read(fdE, d.E, nWindows);
		}

	}
};

#endif
