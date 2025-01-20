/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once
#include "blocklanguage.hpp"
#include "lem.hpp"
#include "layeredearthmodeller.hpp"
#include "aem_coredefs.hpp"
#include "tdemresponse.hpp"
#include "tdemgeometry.hpp"

namespace AEM {
	template <typename ResponseType>
	class AEMSystem {

	protected:
		std::string SystemName;
		std::string SystemType;
		cBlock STM;
		LEModeller LEM;

	public:
		//using RT = ResponseType;

		AEMSystem() {};

		AEMSystem(const fs::path& descriptorpath) {
			read_system_descriptor_file(descriptorpath);
		};

		const cBlock& system_descriptor_block() const { return STM; };
		LEModeller& lem() { return LEM; };

		virtual void read_system_descriptor_file(const fs::path& systemdescriptorfile) = 0;
		virtual const size_t& nWindows() const = 0;

		//virtual const TDEmVectorResponse& forward_model_primary_field(const TDEmGeometry& G) = 0;
		//virtual const TDEmResponse& forward_model(const Earth1D& E, const TDEmGeometry& G) = 0;
		//virtual TDEmVectorResponse derivative(const CalculationType& calc, const TDEmGeometry& G, const TDEmVectorResponse& forward_model) = 0;
		//virtual const TDEmResponse& derivative(const CalculationType& calc, const TDEmGeometry& G, const TDEmResponse& forward_model) = 0;
		//virtual const TDEmResponse& derivative(const CalculationType& calc) = 0;

		virtual const TDEmVectorResponse<ResponseType>& forward_model_primary_field(const TDEmGeometry& G) = 0;
		virtual const TDEmResponse<ResponseType>& forward_model(const Earth1D& E, const TDEmGeometry& G) = 0;
		virtual TDEmVectorResponse<ResponseType> derivative(const CalculationType& calc, const TDEmGeometry& G, const TDEmVectorResponse<ResponseType>& forward_model) = 0;
		virtual const TDEmResponse<ResponseType>& derivative(const CalculationType& calc, const TDEmGeometry& G, const TDEmResponse<ResponseType>& forward_model) = 0;
		virtual const TDEmResponse<ResponseType>& derivative(const CalculationType& calc) = 0;
	};
};

