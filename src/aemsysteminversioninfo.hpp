/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include "inputmanager.hpp"
#include "fielddefinition.hpp"
#include "aemsystem.hpp"

namespace AEM {
	class cScaleFactorsStruct {

	public:
		double input;
		double ref;
		double std;
		double min;
		double max;
		double invmodel;
	};

	class AEMComponentInversionInfo {

	private:
		size_t nWindows = 0;
		size_t nSoundings = 0;

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

		cInvertibleFieldDefinition fdSF;
		cScaleFactorsStruct SF;

		AEMComponentInversionInfo() {};
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

			cBlock sfb = b.findblock("ScaleFactor");
			if (sfb.empty() == false) {
				fdSF.initialise(b, "ScaleFactor");
			}

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

	class AEMSystemInversionInfo {

	public:
		std::unique_ptr<AEMSystem> System;
		std::string SystemFile;
		size_t nwindows = 0;
		size_t ncomps = 0;
		size_t nchans = 0;
		AEMComponentInversionInfo CompInfo[NCOMP];
		std::vector<TDEmResponse> predicted;
		std::string units;

		bool InvertXZAmplitude = false;
		bool InvertPrimaryPlusSecondary = false;
		bool ReconstructPrimary = false;

		AEMSystemInversionInfo(cBlock& b, const size_t nsoundings) :
			System(std::make_unique<TDEmSystem>(b.getstringvalue("SystemFile")))
		{
			initialise(b, nsoundings);
		};

		void initialise(const cBlock& b, const size_t nsoundings) {
			std::string stmfile = b.getstringvalue("SystemFile");
			glog.log_to_file(strprint("==============System file %s\n", stmfile.c_str()));
			glog.log_to_file(System->system_descriptor_block().get_as_string());
			glog.log_to_file("==========================================================================\n");
			nwindows = System->nWindows();

			std::string dummy;
			if (b.getvalue("InvertTotalField", dummy)) {
				glog.errormsg(_SRC_, "InvertTotalField is no longer an option, use InvertPrimaryPlusSecondary instead\n");
			};

			//bool status = false;
			if (b.getvalue("InvertXPlusZ", InvertXZAmplitude)) {
				glog.warningmsg("'InvertXPlusZ' is deprecated, please use 'InvertXZAmplitude' instead.");
			}
			else (b.getvalue("InvertXZAmplitude", InvertXZAmplitude));

			InvertPrimaryPlusSecondary = b.getboolvalue("InvertPrimaryPlusSecondary");
			ReconstructPrimary = b.getboolvalue("ReconstructPrimaryFieldFromInputGeometry");

			CompInfo[XCOMP].initialise(b.findblock("XComponent"), "X", nwindows, nsoundings);
			CompInfo[YCOMP].initialise(b.findblock("YComponent"), "Y", nwindows, nsoundings);
			CompInfo[ZCOMP].initialise(b.findblock("ZComponent"), "Z", nwindows, nsoundings);

			ncomps = 0;
			if (CompInfo[XCOMP].Use) ncomps++;
			if (CompInfo[YCOMP].Use) ncomps++;
			if (CompInfo[ZCOMP].Use) ncomps++;

			if (InvertXZAmplitude) {
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
					if (units.size() > 0) {
						if (ciequal(u, units) == false) {
							std::ostringstream msg;
							msg << "Error: units must be the same on all EM system/components. " << units << "does niot match " << u << "." << std::endl;
							glog.errormsg(_SRC_, msg.str().c_str());
						}
					}
					else units = u;
				}
			}
		}

	};
};
