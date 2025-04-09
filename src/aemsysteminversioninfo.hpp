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

	class Parameter {

		private:
			bool _initialised_status = false;
			bool _solve = false;
			bool _bound = false;
			FDMap fdmap;
			int  poffset = -1;//offset into paramter index 

		public:
			static constexpr const char INPUT[] = "Input";
			static constexpr const char REF[] = "Ref";
			static constexpr const char STD[] = "Std";
			static constexpr const char MIN[] = "Min";
			static constexpr const char MAX[] = "Max";

			Parameter() {};

			Parameter(const cBlock& parentblock, const std::string& key) {
				_initialised_status = initialise(parentblock, key);
			};

			bool initialised() const {
				return _initialised_status;
			};

			bool solve() const {
				return _solve;
			};

			bool bound() const {
				return _bound;
			};

			void set_poffset(const int& _poffset) {
				if (_initialised_status && _solve) {
					poffset = _poffset;
				}
				else glog.errormsg(_SRC_, "Cannot set poffset for a fixed parameter.\n");
			};

			int get_poffset() const {
				if(_initialised_status && _solve) return poffset;
				return -1;
			};

			const cFieldDefinition& get_fd(const std::string & key) const {
				if (has(key)) return fdmap.at(key);
				else return cFieldDefinition();
			};

		private:
			bool initialise(const cBlock& parentblock, const std::string& key) {
				std::string id = parentblock.findkey(key);
				if (id.compare(undefinedvalue<std::string>()) != 0) {
					return initialise_from_entry(parentblock, key);
				}
				else {
					cBlock b = parentblock.findblock(key);
					if (b.empty() == true) {
						return false;
					}

					if (initialise_from_block(b) == false) {
						std::string msg = strprint("Could not parse control file block: %s.", key.c_str());
						glog.errormsg(_SRC_, msg);
						return false;
					}
					return true;
				}
			};

			bool initialise_from_entry(const cBlock& b, const std::string& key) {
				poffset = -1;
				_solve = false;
				_bound = false;
				cFieldDefinition fd(b, key);
				fdmap[INPUT] = fd;
				return true;
			};

			bool has(const std::string& key) const {
				if (fdmap.find(key) != fdmap.end())return true;
				return false;
			};

			bool initialise_from_block(const cBlock& b) {
				_solve = false;
				_bound = false;
				poffset = -1;

				b.get("solve", _solve, false);
				add_fielddefinition(b, INPUT);
				add_fielddefinition(b, REF);

				if (has(REF) == true && has(INPUT) == false) {
					fdmap[INPUT] = fdmap[REF];
				}

				if (has(INPUT) == true && has(REF) == false) {
					fdmap[REF] = fdmap[INPUT];
				}

				if (_solve) {
					add_fielddefinition(b, STD);
					add_fielddefinition(b, MIN);
					add_fielddefinition(b, MAX);
					
					if (has(REF) == false) {
						glog.errormsg(_SRC_, "Parameter %s with 'solve = yes' must have a 'Ref' or 'Imput'.\n", b.Name.c_str());
					}

					if (has(STD) == false) {
						glog.errormsg(_SRC_, "Parameter %s with 'solve = yes' must have an 'Std'.\n", b.Name.c_str());
					}

					if (has(MIN) && has(MAX)) {
						_bound = true;
					}
				}

				if (has(INPUT) == false) {
					glog.errormsg(_SRC_, "Parameter %s must have an 'Input' or 'Ref'.\n", b.Name.c_str());
				}

				
				
				return true;
			};

			bool add_fielddefinition(const cBlock& b, const std::string& key) {
				cFieldDefinition fd(b, key);
				if (fd.isinitialised()) {
					fdmap[key] = fd;
					return true;
				}
				return false;
			}


	};

	/*
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
	*/

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
				if (min(refvalstd.thickness) < 0) oss << "The thickness std is < 0 in at least one layer\n";
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

	class GGAOffsetStore {

	public:
		std::vector<double> input;
		std::vector<double> refval;
		std::vector<double> refvalstd;
		std::vector<double> minval;
		std::vector<double> maxval;
		std::vector<double> invmodel;
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
		double GA;//High altitude coupling ga
		double GGA;//Coupling Ratio g/ga
	};

	//RT is Response type double or std::complex<double>
	template <typename RT>
	class ComponentInversionInfo {

	public:
		std::vector<SoundingData<RT>> data;
		std::string Name;
		bool Use = false;
		FDMap fdMap;
		GGAOffsetStore ggaoffsetStore;
		ScaleFactorsStore sfStore;
		Parameter ggaoffset;
		Parameter scalingfactor;

		bool EstimateNoiseFromModel = false;
		std::vector<RT> multiplicative_noise;
		std::vector<RT> additive_noise;

		ComponentInversionInfo() {
			_nElements = value_size<RT>();
		};


		const std::string& name() const {
			return Name;
		}

		std::string longname() const {
			return Name + "-component";
		}

		const size_t& nSoundings() const { return _nSoundings; };
		const size_t& nWindows() const { return _nWindows; };
		const size_t& nElements() const { return _nElements; };
		const size_t& nChannels() const { return _nWindows * _nElements; };

		double get_ga(const size_t& si) const {
			return data[si].GA;
		};

		double get_gga(const size_t& si) const {
			return data[si].GGA;
		};

		bool getvector_ri(const cBlock& b, const std::string& key, std::vector<double>& v) {
			bool status = b.getvalue(key, v);
			return status;
		}

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
		};

		void add_parameters(const cBlock& b) {
			scalingfactor = Parameter(b, "ScaleFactor");
			ggaoffset = Parameter(b, "GGAOffset");
		}

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
				bool status1 = getvector_ri(b, "MultiplicativeNoise", multiplicative_noise);
				bool status2 = getvector_ri(b, "AdditiveNoise", additive_noise);
				if (additive_noise.size() == 1) {
					additive_noise = std::vector<RT>(nwindows, additive_noise[0]);
				}
				else if (additive_noise.size() != nwindows) {
					glog.errormsg(_SRC_, "Must have exactly 1 or nwindows AdditiveNoise values\n");
				};

				if (multiplicative_noise.size() == 1) {
					multiplicative_noise = std::vector<RT>(nwindows, multiplicative_noise[0]);
				}
				if (multiplicative_noise.size() != nwindows) {
					glog.errormsg(_SRC_, "Must have exactly 1 or nwindows MultiplicativeNoise values\n");
				}
			}

			add_fielddefinitions(b);
			add_parameters(b);

			_nSoundings = nsoundings;
			_nWindows = nwindows;
			_nElements = value_size<RT>();
			_nChannels = _nWindows * _nElements;

			data.resize(nSoundings());
			for (size_t si = 0; si < nSoundings(); si++) {
				data[si].S.resize(nWindows());
				data[si].E.resize(nWindows());
			}

			if (ggaoffset.initialised()) {
				ggaoffsetStore.input.resize(nSoundings());
				if (ggaoffset.solve()) {
					ggaoffsetStore.refval.resize(nSoundings());
					ggaoffsetStore.refvalstd.resize(nSoundings());
					ggaoffsetStore.minval.resize(nSoundings());
					ggaoffsetStore.maxval.resize(nSoundings());
				}
			}
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

			const cFieldDefinition& fdP = fdMap["Primary"];
			const cFieldDefinition& fdS = fdMap["Secondary"];
			IM->read(fdP, d.P, nWindows());
			IM->read(fdS, d.S, nWindows());
			if (EstimateNoiseFromModel) {
				for (size_t wi = 0; wi < nWindows(); wi++) {
					const RT v = 0.01 * AEM::ewise_mul(multiplicative_noise[wi], d.S[wi]);
					d.E[wi] = AEM::hypot(additive_noise[wi], v);
				}
			}
			else {
				const auto& fdE = fdMap["Noise"];
				IM->read(fdE, d.E, nWindows());
			}
		};

		template <>
		void readdata_impl<cdouble>(const std::unique_ptr<cInputManager>& IM, const size_t& soundingindex) {
			
			if (Use == false) return;

			const size_t& si = soundingindex;
			SoundingData<RT>& d = data[si];

			const cFieldDefinition& fdTr = fdMap.at("TotalReal");
			const cFieldDefinition& fdTi = fdMap.at("TotalImag");
			const cFieldDefinition& fdPr = fdMap.at("PrimaryReal");
			const cFieldDefinition& fdPi = fdMap.at("PrimaryImag");
			const cFieldDefinition& fdSr = fdMap.at("SecondaryReal");
			const cFieldDefinition& fdSi = fdMap.at("SecondaryImag");
			const cFieldDefinition& fdNr = fdMap.at("NoiseReal");
			const cFieldDefinition& fdNi = fdMap.at("NoiseImag");

			IM->read(fdTr, fdTi, d.T, nWindows());
			IM->read(fdSr, fdSi, d.S, nWindows());
			IM->read(fdPr, fdPi, d.P, nWindows());
			IM->read(fdNr, fdNi, d.E, nWindows());

			const cFieldDefinition& fdGA = fdMap.at("GA");
			IM->read(fdGA, d.GA, 1);
			const cFieldDefinition& fdGGA = fdMap.at("GGA");
			IM->read(fdGGA, d.GGA, 1);

		};

		void read_ggaoffset_data(const std::unique_ptr<cInputManager>& IM, const size_t& soundingindex) {
			if (Use == false) return;
			const Parameter& p = ggaoffset;
			if (p.initialised()) {
				GGAOffsetStore& s = ggaoffsetStore;
				IM->read(p.get_fd(Parameter::INPUT), s.input[soundingindex], 1);
				if (p.solve()) {
					IM->read(p.get_fd(Parameter::REF), s.refval[soundingindex], 1);
					IM->read(p.get_fd(Parameter::STD), s.refvalstd[soundingindex], 1);
					if (p.bound()) {
						IM->read(p.get_fd(Parameter::MIN), s.minval[soundingindex], 1);
						IM->read(p.get_fd(Parameter::MAX), s.maxval[soundingindex], 1);
					}
				}
			}
		};

		void estimate_noise_from_model(const size_t& soundingindex) {
			if (Use == false) return;
			if (EstimateNoiseFromModel == false) return;
			SoundingData<RT>& d = data[soundingindex];
			bool invertpsi = true;//Todo fix this bookmark
			if (invertpsi) {
				double gga = get_gga(soundingindex);
				for (size_t wi = 0; wi < nWindows(); wi++) {
					RT val = d.T[wi];
					val -= gga;
					const RT mn = 0.01 * AEM::ewise_mul(multiplicative_noise[wi], val);
					d.E[wi] = AEM::hypot(additive_noise[wi], mn);
				}
			}
			else {
				for (size_t wi = 0; wi < nWindows(); wi++) {
					const RT mn = 0.01 * AEM::ewise_mul(multiplicative_noise[wi], d.S[wi]);
					d.E[wi] = AEM::hypot(additive_noise[wi], mn);
				}
			}
		};

		private:
			size_t _nWindows = 0;
			size_t _nElements = 0;
			size_t _nChannels = 0;
			size_t _nSoundings = 0;
	};

	//RT is Response type double or std::complex<double>
	template <typename AEMSystemClass, typename RT>
	class SystemInversionInfo {

	public:
		
		std::unique_ptr<AEMSystem<RT>> System;
		ComponentInversionInfo<RT> CI[NCOMP];
		using CompInfo = ComponentInversionInfo<RT>;
		std::vector<TDEmResponse<RT>> predicted;
		std::string Units;

		bool InvertXZAmplitude  = false;
		bool InvertTotalField   = false;
		bool InvertPSI = false;
		bool ReconstructPrimary = false;

		SystemInversionInfo(cBlock& b, const size_t nsoundings){
			fs::path stmfile = b.getstringvalue("SystemFile");
			System = AEMSystemClass::unique_ptr(stmfile);
			initialise(b, nsoundings);
		};

		void initialise(const cBlock& b, const size_t nsoundings) {
			std::string stmfile = b.getstringvalue("SystemFile");
			glog.log_to_file(strprint("==============System file %s\n", stmfile.c_str()));
			glog.log_to_file(System->system_descriptor_block().get_as_string());
			glog.log_to_file("==========================================================================\n");
			_nWindows = System->nWindows();

			if (b.getvalue("InvertPrimaryPlusSecondary", InvertTotalField)) {
				glog.warningmsg("'InvertPrimaryPlusSecondary' is deprecated, please use 'InvertTotalField' instead\n");
			}
			else (b.getvalue("InvertTotalField", InvertTotalField));


			b.get("InvertPSI", InvertPSI, false);

			ReconstructPrimary = false;
			if (InvertTotalField) {
				ReconstructPrimary = b.getboolvalue("ReconstructPrimaryFieldFromInputGeometry");
			}

			if (b.getvalue("InvertXPlusZ", InvertXZAmplitude)) {
				glog.warningmsg("'InvertXPlusZ' is deprecated, please use 'InvertXZAmplitude' instead.");
			}
			else (b.getvalue("InvertXZAmplitude", InvertXZAmplitude));

			CI[XCOMP].initialise(b.findblock("XComponent"), "X", _nWindows, nsoundings);
			CI[YCOMP].initialise(b.findblock("YComponent"), "Y", _nWindows, nsoundings);
			CI[ZCOMP].initialise(b.findblock("ZComponent"), "Z", _nWindows, nsoundings);

			_nActiveComponents = 0;
			if (CI[XCOMP].Use) _nActiveComponents++;
			if (CI[YCOMP].Use) _nActiveComponents++;
			if (CI[ZCOMP].Use) _nActiveComponents++;

			if (InvertXZAmplitude) {
				CI[XCOMP].Use = true;
				CI[ZCOMP].Use = true;
			}
			size_t vsize = value_size<RT>();
			_nChannels = _nWindows * _nActiveComponents * vsize;
		}

		void set_units(cInputManager* IM) {
			for (size_t ci = 0; ci < NCOMP; ci++) {
				cFieldDefinition& fd = CI[ci].fdMap["Secondary"];
				if (fd.get_varname().size() > 0) {
					cAsciiColumnField c;
					IM->get_acsiicolumnfield(fd, c);
					std::string u = c.get_att("units");
					if (Units.size() > 0) {
						if (ciequal(u, Units) == false) {
							std::ostringstream msg;
							msg << "Error: units must be the same on all EM system/components. " << Units << "does niot match " << u << "." << std::endl;
							glog.errormsg(_SRC_, msg.str().c_str());
						}
					}
					else Units = u;
				}
			}
		}

		AEMSystem<RT>& sys() {
			return *System;
		};

		const AEMSystem<RT>& sys() const {
			return *System;
		};

		const size_t& nWindows() const { return _nWindows; };
		const size_t& nActiveComponents() const { return _nActiveComponents; };
		const size_t& nChannels() const { return _nChannels; };

		Vec3d get_scaled_ga(const size_t& si) const {
			Vec3d ga(1.0, 1.0, 1.0);
			for (size_t ci = 0; ci < NCOMP; ci++) {
				const CompInfo& C = CI[ci];
				if (C.Use) ga[ci] = CI[ci].get_ga(si);
			}
			return (1e-7 * 1e15) * ga;
		};

		Vec3d get_gga(const size_t si) {
			Vec3d gga(1.0, 1.0, 1.0);
			for (size_t ci = 0; ci < NCOMP; ci++) {
				CompInfo& C = CI[ci];
				if (C.Use) {
					gga[ci] = CI[ci].get_gga(si);
				}
			}
			return gga;
		};

		std::vector<RT> get_predicted1(const size_t& si, const size_t& ci) const {
			if (InvertPSI) {
				std::vector<RT> v = predicted[si].total(ci);
				CompInfo& C = CI[ci];
				C.get_ga(si)
			}
			else if (InvertTotalField) {
				return predicted[si].total(ci);
			}
			else return predicted[si].secondary(ci);
		};

		std::vector<RT> get_predicted_xzamp1(const size_t& si) const {
			if(InvertTotalField){
				TDEmVectorResponse<RT> T = predicted[si].totalfield();
				return T.xzamp().storage();
			}
			if (InvertPSI) {
				//const std::vector<RT>& p = predicted[si].primary(ci);
				//const std::vector<RT>& s = predicted[si].secondary(ci);
				//const std::vector<RT> t = tx + sx;
			}
		};

		//std::vector<RT> get_xzamp(const size_t si) {
		//	//if(InvertPSI)
		//	if (InvertXZAmplitude) {
		//		const std::vector<RT> px = predicted[si].primary(ci);
		//		const std::vector<RT> sx = predicted[si].secondary(ci);
		//		const std::vector<RT> tx = px + sx;
		//	}
		//};

	private:
		size_t _nWindows = 0;
		size_t _nActiveComponents = 0;
		size_t _nChannels = 0;
	};
};
