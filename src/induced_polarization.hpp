/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include "aem_types.hpp"

#include <complex>
#include <string_utils.hpp>
#include <sstream>
#include <string_print.hpp>

namespace AEM {

	enum class IPType { NONE, COLECOLE, PELTON };
	
	IPType iptype_from_string(const char* iptype_string)
	{
		if (ciequal(iptype_string, "NONE")) return IPType::NONE;
		else if (ciequal(iptype_string, "COLECOLE")) return IPType::COLECOLE;
		else if (ciequal(iptype_string, "PELTON")) return IPType::PELTON;
		else {
			std::string msg = strprint("Invalid iptype_string: %s\n%s", iptype_string, _SRC_CSTR_);
			throw std::runtime_error(msg);
		}
	};

	inline cdouble ip_colecole_conductivity(const double& conductivity, const double& chargeability, const double& timeconstant, const double& frequencydependence, const double& omega)
	{
		//c = c0 - c0*(N / (1 + (1 - N)*(j*omega*T) ^ K));
		if (chargeability == 0.0) {
			return cdouble(conductivity, 0.0);
		}
		else {
			cdouble c = conductivity - conductivity * (chargeability / (1.0 + (1.0 - chargeability) * (std::pow(cdouble(0.0, omega * timeconstant), frequencydependence))));
			return  c;
		}
	};

	inline cdouble ip_pelton_conductivity(const double& conductivity, const double& chargeability, const double& timeconstant, const double& frequencydependence, const double& omega)
	{
		//p = p0[1 - m*(1 - (1 - 1/(1 + (j*omega*T) ^ K));
		if (chargeability == 0.0) {
			return cdouble(conductivity, 0.0);
		}
		else {
			double  rho0 = 1.0 / conductivity;
			cdouble rho = rho0 * (1.0 - chargeability * (1.0 - (1.0 / (1.0 + std::pow(cdouble(0.0, omega * timeconstant), frequencydependence)))));
			return  1.0 / rho;
		}
	};

	inline cdouble ip_complex_conductivity(const IPType& iptype, const double& conductivity, const double& chargeability, const double& timeconstant, const double& frequencydependence, const double& omega)
	{
		cdouble complex_conductivity;
		if (iptype == IPType::NONE) complex_conductivity = conductivity;
		else if (iptype == IPType::COLECOLE) complex_conductivity = ip_colecole_conductivity(conductivity, chargeability, timeconstant, frequencydependence, omega);
		else complex_conductivity = ip_pelton_conductivity(conductivity, chargeability, timeconstant, frequencydependence, omega);
		return complex_conductivity;
	};

};
