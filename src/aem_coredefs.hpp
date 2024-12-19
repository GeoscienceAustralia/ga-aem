/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include "lem.hpp"
#include "layeredearthmodeller.hpp"

namespace AEM {
	using namespace LEM2;
	using CalculationType = CT::CalculationType;
	using CMode = CT::CalculationType::Mode;

	inline constexpr size_t XCOMP = 0;
	inline constexpr size_t YCOMP = 1;
	inline constexpr size_t ZCOMP = 2;
	inline constexpr size_t NCOMP = 3;
};
