/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <vector>
#include <cassert>

constexpr auto DN_LAYER = "layer";
const std::string DN_NONE;
const std::string UNITLESS;

namespace AEM {
	enum class IPType { NONE, COLECOLE, PELTON };
	
	class Earth1D {

	private:

		IPType ip_type = AEM::IPType::NONE;

	public:

		std::vector<double> thickness;
		std::vector<double> conductivity;
		std::vector<double> chargeability;
		std::vector<double> timeconstant;
		std::vector<double> frequencydependence;

		Earth1D() {};

		Earth1D(const size_t nlayers) {
			thickness.resize(nlayers - 1);
			conductivity.resize(nlayers);
			chargeability.resize(nlayers);
			timeconstant.resize(nlayers);
			frequencydependence.resize(nlayers);
		}

		Earth1D(const std::vector<double>& _conductivity, const std::vector<double>& _thickness) {
			conductivity = _conductivity;
			thickness = _thickness;
		}

		Earth1D(
			const std::vector<double>& _conductivity,
			const std::vector<double>& _thickness,
			const std::vector<double>& _chargeability,
			const std::vector<double>& _timeconstant,
			const std::vector<double>& _frequencydependence)
		{
			conductivity = _conductivity;
			thickness = _thickness;
			chargeability = _chargeability;
			timeconstant = _timeconstant;
			frequencydependence = _frequencydependence;
		}

		Earth1D(const size_t nlayers, const double* _conductivity, const double* _thickness) {
			conductivity = std::vector<double>(_conductivity, _conductivity + nlayers);
			thickness = std::vector<double>(_thickness, _thickness + nlayers - 1);
		}

		Earth1D(const size_t nlayers, const double* _conductivity, const double* _thickness, const double* _chargeability, const double* _timeconstant, const double* _frequencydependence) {
			conductivity = std::vector<double>(_conductivity, _conductivity + nlayers);
			thickness = std::vector<double>(_thickness, _thickness + nlayers - 1);
			chargeability = std::vector<double>(_chargeability, _chargeability + nlayers);
			timeconstant = std::vector<double>(_timeconstant, _timeconstant + nlayers);
			frequencydependence = std::vector<double>(_frequencydependence, _frequencydependence + nlayers);
		}

		const IPType& get_iptype() { return ip_type; }
		void set_iptype(const IPType& _iptype) {
			ip_type = _iptype;
		};

		std::vector<double> layer_top_depth() const {
			const size_t n = thickness.size();
			std::vector<double> dtop(n + 1);
			dtop[0] = 0.0;
			for (size_t i = 1; i <= n; i++) {
				dtop[i] = dtop[i - 1] + thickness[i - 1];
			}
			return dtop;
		}

		std::vector<double> layer_bottom_depth() const {
			const size_t n = thickness.size();
			std::vector<double> dbot(n + 1);
			dbot[0] = thickness[0];
			for (size_t i = 1; i < n; i++) {
				dbot[i] = dbot[i - 1] + thickness[i];
			}
			dbot[n] = dbot[n - 1] + thickness[n - 1];
			return dbot;
		}

		size_t nlayers() const {
			return conductivity.size();
		}

		friend std::ostream& operator<<(std::ostream& os, const Earth1D& e) {
			for (size_t i = 0; i < e.nlayers(); i++) {
				if (i < (e.nlayers() - 1)) os << ixd(4) << i << fxd(10, 6) << e.conductivity[i] << fxd(8, 2) << e.thickness[i] << std::endl;
				else                    os << ixd(4) << i << fxd(10, 6) << e.conductivity[i] << "     inf" << std::endl;
			}
			return os;
		}

		void print() const {
			std::cout << this;
		}

		void write(const fs::path& filepath) const {
			std::ofstream ofs = ofstream_ex(filepath);
			ofs << this;
		}

		std::vector<double> dummy_thickness() const {
			return dummy_thickness(thickness);
		}

		//A dummy thickness vector to have a finite last layer thickness
		static std::vector<double> dummy_thickness(const std::vector<double>& t) {
			const size_t nl = 1 + t.size();
			std::vector<double> tout = t;
			if (nl == 1) {
				tout.push_back(1.0);
			}
			else if (nl == 2) {
				tout.push_back(t[0]);
			}
			else {
				//eg t30 / t29 = t29 / t28;
				tout.push_back(t[nl - 2] * t[nl - 2] / t[nl - 3]);
			}
			return tout;
		}

		double mean_weighted_conductivity() const {
			const size_t nl = nlayers();
			if (nl == 1) return conductivity[0];
			//Returns the thickness weighted mean (in linear space) conductivity (but calculated in linear space)
			double sumc = 0.0;
			double sumt = 0.0;
			for (size_t i = 0; i < nl; i++) {
				if (i < (nl - 1)) {
					sumc += conductivity[i] * thickness[i];
					sumt += thickness[i];
				}
				else {
					sumc += conductivity[i] * sumt; //Make basement layer as thick as sum of all overlying
					sumt += sumt;
				}
			}
			return sumc / sumt;
		}

		double mean_weighted_conductivity_log10_calculation() const {
			const size_t nl = nlayers();
			if (nl == 1) return conductivity[0];

			// Returns the thickness weighted mean (in linear space) conductivity (but calculated in 10g10 space)
			double sumc = 0.0;
			double sumt = 0.0;
			for (size_t i = 0; i < nl; i++) {
				if (i < (nl - 1)) {
					sumc += log10(conductivity[i]) * thickness[i];
					sumt += thickness[i];
				}
				else {
					sumc += log10(conductivity[i]) * sumt; // Make basement layer as thick as sum of all overlying
					sumt += sumt;
				}
			}
			return pow(10.0, sumc / sumt);
		}
	};
};
