/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include "vector_utils.hpp"
#include "inputmanager.hpp"
#include "fielddefinition.hpp"
#include "aemsystem.hpp"
#include "tdemsystem.hpp"
#include "spectralaemsystem.hpp"

using namespace IOManager;
using namespace AEM::INVERTER;

namespace AEM {


	class cInvertibleFieldDefinition {

	private:
		int  poffset = -1;//offset into paramter index 

	public:

		bool solve = false;//should we solve for this parameter or not
		cFieldDefinition input;// input values
		cFieldDefinition refval;// reference model
		cFieldDefinition refvalstd;// std deviation uncertainty
		cFieldDefinition minval;// minimum bound
		cFieldDefinition maxval;// maximum bound

		cInvertibleFieldDefinition() {};

		cInvertibleFieldDefinition(const cBlock& parent, const std::string& key) {
			initialise(parent, key);
		};

		void set_poffset(const int& _poffset) {
			poffset = _poffset;
		};

		int get_poffset() const {
			if (solve) return poffset;
			return -1;
		};

		bool initialise(const cBlock& parent, const std::string& key) {
			std::string id = parent.findkey(key);
			if (id.compare(undefinedvalue<std::string>()) != 0) {
				return initialise_from_entry(parent, key);
			}
			else {
				cBlock b = parent.findblock(key);
				if (b.empty() == true) {
					//std::string msg = strprint("Could not find control file block: %s.", key.c_str());
					//glog.warningmsg(msg);
					return false;
				}

				if (initialise_from_block(b) == false) {
					std::string msg = strprint("Could not parse control file block: %s.", key.c_str());
					glog.errormsg(_SRC_, msg);
				}
				return true;
			}
		};

		bool initialise_from_entry(const cBlock& b, const std::string& key) {
			poffset = -1;
			input.initialise(b, key);
			return true;
		};

		bool initialise_from_block(const cBlock& b) {
			b.get("solve", solve, false);
			input.initialise(b, "input");
			refval.initialise(b, "ref");
			refvalstd.initialise(b, "std");
			minval.initialise(b, "min");
			maxval.initialise(b, "max");
			return true;
		};

		bool bound() const {

			if (solve && minval.isinitialised() && maxval.isinitialised()) {
				return true;
			}
			else return false;
		}
	};	
	using IFDMap = std::map<std::string, cInvertibleFieldDefinition, caseinsensetiveless<std::string>>;

	class GeometryStore {

	public:
		TDEmGeometry input;
		TDEmGeometry refval;
		TDEmGeometry refvalstd;
		TDEmGeometry minval;
		TDEmGeometry maxval;
		TDEmGeometry pfr;// Primary field reconstruction 
		TDEmGeometry invmodel;// Inversion model
	};

	class EarthStore {

	public:
		Earth1D refval;
		Earth1D refvalstd;
		Earth1D minval;
		Earth1D maxval;
		Earth1D invmodel;

		void sanity_check() const {

			size_t nc = refval.conductivity.size();
			size_t nt = refval.thickness.size();

			std::ostringstream oss;
			if (nc != nt + 1) {
				oss << "The conductivity and/or thickness do not have the correct number of layers\n";
			}

			if (refval.conductivity.size() > 0) {
				if (min(refval.conductivity) <= 0) oss << "The conductivity ref is <= 0 in at least one layer\n";
			}

			if (refvalstd.conductivity.size() > 0) {
				if (min(refvalstd.conductivity) <= 0) oss << "The conductivity std is <= 0\n";
			}

			if (minval.conductivity.size() > 0) {
				if (minval.conductivity.size() != nc) oss << "The conductivity min does not have the correct number of layers\n";
				if (maxval.conductivity.size() != nc) oss << "The conductivity max does not have the correct number of layers\n";
				if (min(minval.conductivity) <= 0) oss << "The conductivity min is <= 0 in at least one layer in at least one layer\n";
				if (min(maxval.conductivity) <= 0) oss << "The conductivity max is <= 0 in at least one layer in at least one layer\n";
				if (min(maxval.conductivity - minval.conductivity) <= 0) oss << "The conductivity max <= min in at least one layer\n";
				if (min(refval.conductivity - minval.conductivity) <= 0) oss << "The conductivity ref <= min in at least one layer\n";
				if (min(maxval.conductivity - refval.conductivity) <= 0) oss << "The conductivity ref >= max in at least one layer\n";
			};

			if (refval.thickness.size() > 0) {
				if (min(refval.thickness) <= 0) oss << "The thickness ref is <= 0 in at least one layer\n";
			};

			if (refvalstd.thickness.size() > 0) {
				if (min(refvalstd.thickness) <= 0) oss << "The thickness std is <= 0 in at least one layer\n";
			};

			if (minval.thickness.size() > 0) {
				if (minval.thickness.size() != nt) oss << "The thickness min does not have the correct number of layers\n";
				if (maxval.thickness.size() != nt) oss << "The thickness max does not have the correct number of layers\n";
				if (min(minval.thickness) <= 0) oss << "The thickness min is <= 0 in at least one layer\n";
				if (min(maxval.thickness) <= 0) oss << "The thickness max is <= 0 in at least one layer\n";
				if (min(maxval.thickness - minval.thickness) <= 0) oss << "The thickness max <= min in at least one layer\n";
				if (min(refval.thickness - minval.thickness) <= 0) oss << "The thickness ref <= min in at least one layer\n";
				if (min(maxval.thickness - refval.thickness) <= 0) oss << "The thickness ref >= max in at least one layer\n";
			};

			if (oss.str().size() > 0) {
				glog.errormsg(_SRC_, oss.str());
			};
		}

	};

	class ScaleFactorsStore {

	public:
		double input;
		double refval;
		double refvalstd;
		double minval;
		double maxval;
		double invmodel;
	};

	//RT is Response type double or std::complex<double>
	template <typename RT>
	class SoundingData {
	public:
		std::vector<RT> T;//Total primary + secondary
		std::vector<RT> P;//Primary
		std::vector<RT> S;//Secondary
		std::vector<RT> E;//Noise std estimate
	};

	//RT is Response type double or std::complex<double>
	template <typename RT>
	class AEMComponentInversionInfo {

	private:
		size_t nWindows = 0;
		size_t nSoundings = 0;

	public:
		std::vector<SoundingData<RT>> data;
		std::string Name;
		bool Use = false;
		FDMap fdMap;
		IFDMap ifdMap;
		ScaleFactorsStore sfStore;

		bool EstimateNoiseFromModel = false;
		std::vector<RT> mn;
		std::vector<RT> an;

		AEMComponentInversionInfo() {};

		bool getvector_ri(const cBlock& b, const std::string& key, std::vector<double>& v) {
			bool status = b.getvalue(key, v);
			return status;
		}

		cInvertibleFieldDefinition& get_ifd(const std::string& key) {
			return ifdMap.at(key);
		};

		const cInvertibleFieldDefinition& get_ifd(const std::string& key) const {
			return ifdMap.at(key);
		};

		bool getvector_ri(const cBlock& b, const std::string& key, std::vector<cdouble>& v) {
			std::vector<double> r;
			std::vector<double> i;
			bool status1 = b.getvalue(key+"Real", r);
			bool status2 = b.getvalue(key+"Imag", i);
			if (status1 && status2) {
				complex_merge(r, i, v);
				return true;
			}
			return false;
		}

		template <typename T>
		void add_fielddefinition(const cBlock& b, const std::string& key) {};

		template <> 
		void add_fielddefinition<double>(const cBlock& b, const std::string& key) {
			fdMap[key] = cFieldDefinition(b, key);
		};

		template <>
		void add_fielddefinition<cdouble>(const cBlock& b, const std::string& key) {
			std::string rkey = key + "Real";
			std::string ikey = key + "Imag";
			fdMap[rkey] = cFieldDefinition(b, rkey);
			fdMap[ikey] = cFieldDefinition(b, ikey);
		};

		void add_invertiblefielddefinition(const cBlock& b, const std::string& key) {
			ifdMap[key] = cInvertibleFieldDefinition(b, key);
		};

		void add_fielddefinitions(const cBlock& b) {
			add_fielddefinition<RT>(b, "Primary");
			add_fielddefinition<RT>(b, "Secondary");
			add_fielddefinition<RT>(b, "Noise");
			add_fielddefinition<RT>(b, "Total");
			add_fielddefinition<double>(b, "GA");
			add_fielddefinition<double>(b, "GGA");
			add_invertiblefielddefinition(b, "ScaleFactor");
		};

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
				bool status1 = getvector_ri(b, "MultiplicativeNoise", mn);
				bool status2 = getvector_ri(b, "AdditiveNoise", an);
				if (an.size() == 1) {
					an = std::vector<RT>(nwindows, an[0]);
				}
				else if (an.size() != nwindows) {
					glog.errormsg(_SRC_, "Must have exactly 1 or nwindows AdditiveNoise values\n");
				};

				if (mn.size() == 1) {
					mn = std::vector<RT>(nwindows, mn[0]);
				}
				if (mn.size() != nwindows) {
					glog.errormsg(_SRC_, "Must have exactly 1 or nwindows MultiplicativeNoise values\n");
				}
			}

			add_fielddefinitions(b);

			nSoundings = nsoundings;
			nWindows = nwindows;
			data.resize(nSoundings);
			for (size_t si = 0; si < nSoundings; si++) {
				data[si].S.resize(nWindows);
				data[si].E.resize(nWindows);
			}
		}

		const size_t& nw() const {
			return nWindows;
		}

		void readdata(const std::unique_ptr<cInputManager>& IM, const size_t& soundingindex) {
			readdata_impl<RT>(IM, soundingindex);
		};

		template <typename RT>
		void readdata_impl(const std::unique_ptr<cInputManager>& IM, const size_t& soundingindex) {};

		template <> 
		void readdata_impl<double>(const std::unique_ptr<cInputManager>& IM, const size_t& soundingindex) {
			const size_t& si = soundingindex;
			if (Use == false) return;
			SoundingData<RT>& d = data[si];

			auto& fdP = fdMap["Primary"];
			auto& fdS = fdMap["Secondary"];
			IM->read(fdP, d.P, nw());
			IM->read(fdS, d.S, nw());
			if (EstimateNoiseFromModel) {
				for (size_t wi = 0; wi < nWindows; wi++) {
					const RT v = 0.01 * AEM::ewise_mul(mn[wi],d.S[wi]);
					d.E[wi] = AEM::hypot(an[wi], v);
				}
			}
			else {
				const auto& fdE = fdMap["Noise"];
				IM->read(fdE, d.E, nWindows);
			}
		}

		template <>
		void readdata_impl<cdouble>(const std::unique_ptr<cInputManager>& IM, const size_t& soundingindex) {
			const size_t& si = soundingindex;
			if (Use == false) return;
			SoundingData<RT>& d = data[si];
			const cFieldDefinition fdTr = fdMap["TotalReal"];
			const cFieldDefinition fdTi = fdMap["TotalImag"];
			const cFieldDefinition fdPr = fdMap["PrimaryReal"];
			const cFieldDefinition fdPi = fdMap["PrimaryImag"];
			const cFieldDefinition fdSr = fdMap["SecondaryReal"];
			const cFieldDefinition fdSi = fdMap["SecondaryImag"];
			const cFieldDefinition fdNr = fdMap["NoiseReal"];
			const cFieldDefinition fdNi = fdMap["NoiseImag"];

			IM->read(fdTr, fdTi, d.T, nw());
			IM->read(fdSr, fdSi, d.S, nw());
			IM->read(fdPr, fdPi, d.P, nw());
			IM->read(fdNr, fdNi, d.E, nw());

			if(EstimateNoiseFromModel){
				if (fdTr.isinitialised()) {
					for (size_t wi = 0; wi < nWindows; wi++) {
						const RT v = 0.01 * AEM::ewise_mul(mn[wi], d.T[wi]);
						d.E[wi] = AEM::hypot(an[wi], v);
					}
				}
				else {
					for (size_t wi = 0; wi < nWindows; wi++) {
						const RT v = 0.01 * AEM::ewise_mul(mn[wi], d.S[wi]);
						d.E[wi] = AEM::hypot(an[wi], v);
					}
				}
			}
		}

		//int scalefactor_poffset() const {
		//	const cInvertibleFieldDefinition& ifd = ifdMap.at("ScaleFactor");
		//	return ifd.get_poffset();
		//};
	};

	//RT is Response type double or std::complex<double>
	template <typename AEMSystemClass, typename RT>
	class AEMSystemInversionInfo {

	public:
		
		std::unique_ptr<AEMSystem<RT>> System;
		
		size_t nwindows = 0;
		size_t ncomps = 0;
		size_t nchans = 0;
		AEMComponentInversionInfo<RT> CompInfo[NCOMP];		 
		std::vector<TDEmResponse<RT>> predicted;
		std::string units;

		bool InvertXZAmplitude  = false;
		bool InvertTotalField   = false;
		bool ReconstructPrimary = false;

		AEMSystemInversionInfo(cBlock& b, const size_t nsoundings){
			fs::path stmfile = b.getstringvalue("SystemFile");
			System = AEMSystemClass::unique_ptr(stmfile);
			initialise(b, nsoundings);
		};

		void initialise(const cBlock& b, const size_t nsoundings) {
			std::string stmfile = b.getstringvalue("SystemFile");
			glog.log_to_file(strprint("==============System file %s\n", stmfile.c_str()));
			glog.log_to_file(System->system_descriptor_block().get_as_string());
			glog.log_to_file("==========================================================================\n");
			nwindows = System->nWindows();

			if (b.getvalue("InvertPrimaryPlusSecondary", InvertTotalField)) {
				glog.warningmsg("'InvertPrimaryPlusSecondary' is deprecated, please use 'InvertTotalField' instead\n");
			}
			else (b.getvalue("InvertTotalField", InvertTotalField));

			ReconstructPrimary = false;
			if (InvertTotalField) {
				ReconstructPrimary = b.getboolvalue("ReconstructPrimaryFieldFromInputGeometry");
			}

			if (b.getvalue("InvertXPlusZ", InvertXZAmplitude)) {
				glog.warningmsg("'InvertXPlusZ' is deprecated, please use 'InvertXZAmplitude' instead.");
			}
			else (b.getvalue("InvertXZAmplitude", InvertXZAmplitude));

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
			//nchans = nwindows * ncomps;
		}

		void set_units(cInputManager* IM) {
			for (size_t ci = 0; ci < NCOMP; ci++) {
				cFieldDefinition& fd = CompInfo[ci].fdMap["Secondary"];
				if (fd.get_varname().size() > 0) {
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
