/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include "earth1d.hpp"
#include "tdemsystem.hpp"
#include "tdemresponse.hpp"
#include "tdemgeometry.hpp"
#include "calculation_type.hpp"
#include "logger.hpp"

//GATDAEM1D_API_EXPORTS is defined because here we are implmenting/building the library
//	must go before #include "gatdaem1d.h"
#define GATDAEM1D_API_EXPORTS

#include "gatdaem1d.h"
#include "gaaem_version.hpp"

#include <cstring>
#include <cstdio>
#include <sstream>
#include <cstdint>
#include <string>
#include <new>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <cassert>

//using namespace AEM;
class cLogger glog; //The global instance of the log file manager

namespace AEM {
	namespace GATDAEM1D {
		//API system handle
		struct gatdaem1d_system_handle_tag {
			uint64_t magic = 0;
			AEM::TDEmSystem* sys;
		};
		static constexpr uint64_t GATDAEM1D_MAGIC = 0x4741544441454D31ULL; // "GATDAEM1D" ish

		//API Error handling
		static thread_local char GATDAEM1D_LAST_ERROR_STRING[1024] = { 0 };

		const char* get_status_string(const gatdaem1d_status status)
		{
			switch (status) {
			case GATDAEM1D_STATUS_OK:
				return "GATDAEM1D: OK";
			case GATDAEM1D_STATUS_ERRBADHANDLE:
				return "GATDAEM1D: Bad handle";
			case GATDAEM1D_STATUS_ERRBADMAGIC:
				return "GATDAEM1D: Bad magic number";
			case GATDAEM1D_STATUS_ERRNULLSYS:
				return "GATDAEM1D: Null system pointer";
			case GATDAEM1D_STATUS_ERRCPPEXCEPTION:
				return "GATDAEM1D: C++ exception";
			default:
				return "GATDAEM1D: Unknown error code";
			};
		};

		static void set_last_error_string(const char* msg) {
			if (!msg) msg = "GATDAEM1D: No error";
			std::snprintf(GATDAEM1D_LAST_ERROR_STRING, sizeof(GATDAEM1D_LAST_ERROR_STRING), "%s", msg);
		};

		// API helpers - not exported to public C interface
		static inline gatdaem1d_system_handle_tag* handle_tag_ptr_from_handle(gatdaem1d_system_handle handle)
		{
			return reinterpret_cast<gatdaem1d_system_handle_tag*>(static_cast<uintptr_t>(handle));
		};

		static inline gatdaem1d_system_handle handle_from_handle_tag_ptr(gatdaem1d_system_handle_tag* ptr_tag)
		{
			return static_cast<gatdaem1d_system_handle>(reinterpret_cast<uintptr_t>(ptr_tag));
		};

		static inline void validate_handle(const gatdaem1d_system_handle_tag* ptr_tag)
		{
			if (!ptr_tag) throw std::runtime_error("Invalid gatdaem1d handle (null)");
			if (ptr_tag->magic != GATDAEM1D_MAGIC) throw std::runtime_error("Invalid gatdaem1d handle (bad magic)");
			if (!ptr_tag->sys) throw std::runtime_error("Invalid gatdaem1d handle (null sys)");
		};

		static TDEmSystem& tdemsystem_reference(const gatdaem1d_system_handle handle) {
			gatdaem1d_system_handle_tag* ptr_tag = handle_tag_ptr_from_handle(handle);
			validate_handle(ptr_tag);
			TDEmSystem& T = *(ptr_tag->sys);
			return T;
		};

		static bool valid_iptype(const gatdaem1d_iptype iptype) {
			bool status = (iptype == GATDAEM1D_IPTYPE_NONE || iptype == GATDAEM1D_IPTYPE_COLECOLE || iptype == GATDAEM1D_IPTYPE_PELTON);
			return status;
		};

		static void check_response_struct(const TDEmSystem& T, const gatdaem1d_response* sR) {
			if (!sR->PX) throw std::runtime_error("Argument struct ->PX is null");
			if (!sR->PY) throw std::runtime_error("Argument struct ->PY is null");
			if (!sR->PZ) throw std::runtime_error("Argument struct ->PZ is null");
			if (!sR->SX) throw std::runtime_error("Argument struct ->SX is null");
			if (!sR->SY) throw std::runtime_error("Argument struct ->SY is null");
			if (!sR->SZ) throw std::runtime_error("Argument struct ->SZ is null");
			if (sR->nwindows != T.nWindows()) throw std::runtime_error("Mismatch in sR->nwindows and T.nWindows()");
		};

		static TDEmGeometry convert_geometry_struct(const gatdaem1d_system_geometry* sG) {
			TDEmGeometry G(sG->tx_height, sG->tx_roll, sG->tx_pitch, sG->tx_yaw, sG->txrx_dx, sG->txrx_dy, sG->txrx_dz, sG->rx_roll, sG->rx_pitch, sG->rx_yaw);
			return G;
		};

		static Earth1D convert_earth_struct(const gatdaem1d_earth* pE) {
			gatdaem1d_iptype iptype = pE->iptype;
			bool validipmodel = pE->iptype == GATDAEM1D_IPTYPE_NONE || pE->iptype == GATDAEM1D_IPTYPE_COLECOLE || pE->iptype == GATDAEM1D_IPTYPE_PELTON;
			if (valid_iptype(pE->iptype) == false) {
				throw std::runtime_error("Invalid gatdaem1d_iptype");
			}
			if (!pE->conductivity) throw std::runtime_error("Argument pE->conductivity is null");
			if (!pE->thickness) throw std::runtime_error("Argument pE->thickness is null");

			Earth1D E;
			if (pE->iptype == GATDAEM1D_IPTYPE_NONE) {
				E = Earth1D(pE->nlayers, pE->conductivity, pE->thickness);
				E.set_iptype(AEM::IPType::NONE);
			}
			else {
				if (!pE->chargeability) throw std::runtime_error("Argument pE->chargeability is null but should be non-null for pE->iptype is GATDAEM1D_IPTYPE_COLECOLE or GATDAEM1D_IPTYPE_PELTON");
				if (!pE->timeconstant) throw std::runtime_error("Argument pE->timeconstant is null but should be non-null for pE->iptype is GATDAEM1D_IPTYPE_COLECOLE or GATDAEM1D_IPTYPE_PELTON");
				if (!pE->frequencydependence) throw std::runtime_error("Argument pE->frequencydependence is null but should be non-null for pE->iptype is GATDAEM1D_IPTYPE_COLECOLE or GATDAEM1D_IPTYPE_PELTON");
				E = Earth1D(pE->nlayers, pE->conductivity, pE->thickness, pE->chargeability, pE->timeconstant, pE->frequencydependence);
				if (pE->iptype == GATDAEM1D_IPTYPE_COLECOLE) E.set_iptype(AEM::IPType::COLECOLE);
				if (pE->iptype == GATDAEM1D_IPTYPE_PELTON) E.set_iptype(AEM::IPType::PELTON);
			}
			return E;
		};

		static void copy_bytes(const AEM::TDEmResponse<double>& R, double* PX, double* PY, double* PZ, double* SX, double* SY, double* SZ)
		{
			const size_t nbytes = R.nWindows() * sizeof(double);
			std::memcpy(PX, R.primary(XCOMP).data(), nbytes);
			std::memcpy(PY, R.primary(YCOMP).data(), nbytes);
			std::memcpy(PZ, R.primary(ZCOMP).data(), nbytes);
			std::memcpy(SX, R.secondary(XCOMP).data(), nbytes);
			std::memcpy(SY, R.secondary(YCOMP).data(), nbytes);
			std::memcpy(SZ, R.secondary(ZCOMP).data(), nbytes);
		};

		static void copy_response(const AEM::TDEmResponse<double>& R, gatdaem1d_response* pR)
		{
			const size_t nbytes = R.nWindows() * sizeof(double);
			std::memcpy(pR->PX, R.primary(XCOMP).data(), nbytes);
			std::memcpy(pR->PY, R.primary(YCOMP).data(), nbytes);
			std::memcpy(pR->PZ, R.primary(ZCOMP).data(), nbytes);
			std::memcpy(pR->SX, R.secondary(XCOMP).data(), nbytes);
			std::memcpy(pR->SY, R.secondary(YCOMP).data(), nbytes);
			std::memcpy(pR->SZ, R.secondary(ZCOMP).data(), nbytes);
		};

		//Copy bytes to the end of the array and return the pointer advanced to the end of the written data
		static double* copy_bytes_array_end_advance_pointer(const AEM::TDEmResponse<double>& R, double* p)
		{
			assert(p != nullptr);
			const size_t nw = R.nWindows();
			assert(R.primary(XCOMP).size() == nw);
			assert(R.secondary(XCOMP).size() == nw);

			const size_t nbytes = nw * sizeof(double);
			std::memcpy(p, R.primary(XCOMP).data(), nbytes); p += nw;
			std::memcpy(p, R.secondary(XCOMP).data(), nbytes); p += nw;
			std::memcpy(p, R.primary(YCOMP).data(), nbytes); p += nw;
			std::memcpy(p, R.secondary(YCOMP).data(), nbytes); p += nw;
			std::memcpy(p, R.primary(ZCOMP).data(), nbytes); p += nw;
			std::memcpy(p, R.secondary(ZCOMP).data(), nbytes); p += nw;
			return p;
		};

		// API entry-point lambda - C++ code that may throw exceptions go inside the lambda body so that exceptions do not get thrown across ABI boundary
		template<class F>
		gatdaem1d_status api_entry_point(F&& f)
		{
			try {
				f();
				return GATDAEM1D_STATUS_OK;
			}
			catch (const std::exception& e) {
				set_last_error_string(e.what());
				return GATDAEM1D_STATUS_ERRCPPEXCEPTION;
			}
			catch (...) {
				set_last_error_string("GATDAEM1D Unknown C++ exception");
				return GATDAEM1D_STATUS_ERRCPPEXCEPTION;
			}
		};
	};
};

using namespace AEM;
using namespace AEM::GATDAEM1D;

//Exported functions start here

// Returns the ABI version the DLL was built with.
GATDAEM1D_API int32_t GATDAEM1D_CALL gatdaem1d_get_abi_version() 
{
	return GATDAEM1D_ABI_VERSION;
};

GATDAEM1D_API const char* GATDAEM1D_CALL gatdaem1d_get_last_error_string(const gatdaem1d_status error_code) {
	std::ostringstream oss;
	if (error_code != GATDAEM1D_STATUS_OK) {
		oss << get_status_string(error_code) << std::endl;
		if (error_code == GATDAEM1D_STATUS_ERRCPPEXCEPTION) {
			oss << GATDAEM1D_LAST_ERROR_STRING << std::endl;
		}
	}
	set_last_error_string(oss.str().c_str());
	return GATDAEM1D_LAST_ERROR_STRING;
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL createhandle(const char* systemfile, gatdaem1d_system_handle* ptr_handle)
{
	return api_entry_point([&]() -> gatdaem1d_status {
		if (!systemfile)  throw std::runtime_error("argument systemfile is null");
		if (!ptr_handle)  throw std::runtime_error("argument out_handle is null");

		// make idempotent for callers that reuse the same variable
		*ptr_handle = 0;

		gatdaem1d_system_handle_tag* ptr_tag = new gatdaem1d_system_handle_tag{};
		try {
			ptr_tag->sys = new TDEmSystem(systemfile);
			ptr_tag->magic = GATDAEM1D_MAGIC;
		}
		catch (...) {
			// prevent leaks if TDEmSystem constructor throws
			delete ptr_tag->sys;
			ptr_tag->sys = nullptr;
			ptr_tag->magic = 0;
			delete ptr_tag;
			throw;
		}

		*ptr_handle = handle_from_handle_tag_ptr(ptr_tag);
		return GATDAEM1D_STATUS::GATDAEM1D_STATUS_OK;
		});
}

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL deletehandle(gatdaem1d_system_handle* ptr_handle)
{
	return api_entry_point([&]() -> gatdaem1d_status {
		if (!ptr_handle) throw std::runtime_error("argument inout_handle is null");

		// allow double-delete and deleting a default-initialized handle
		if (*ptr_handle == 0) return GATDAEM1D_STATUS_OK;

		gatdaem1d_system_handle_tag* h = handle_tag_ptr_from_handle(*ptr_handle);

		// validate before touching internals
		if (!h) throw std::runtime_error("Invalid gatdaem1d handle (null)");
		if (h->magic != GATDAEM1D_MAGIC) throw std::runtime_error("Invalid gatdaem1d handle (bad magic)");

		delete h->sys;
		h->sys = nullptr;
		h->magic = 0;
		delete h;

		*ptr_handle = 0;
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL nsamplesperwaveform(const gatdaem1d_system_handle handle, int32_t* nsamples)
{
	return api_entry_point([&]()->gatdaem1d_status {
		if (!nsamples) throw std::runtime_error("Argument turns is null");
		TDEmSystem& T = tdemsystem_reference(handle);
		*nsamples = static_cast<int>(T.waveform().NumSamples);
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL waveform(const gatdaem1d_system_handle handle, double* time, double* currentwaveform, double* voltagewaveform)
{
	return api_entry_point([&]()->gatdaem1d_status {
		if (!time) throw std::runtime_error("Argument time is null");
		if (!currentwaveform) throw std::runtime_error("Argument currentwaveform is null");
		if (!voltagewaveform) throw std::runtime_error("Argument voltagewaveform is null");
		TDEmSystem& T = tdemsystem_reference(handle);
		for (size_t i = 0; i < T.waveform().NumSamples; i++) {
			time[i] = T.waveform().Time[i];
			currentwaveform[i] = 0.0;
			voltagewaveform[i] = 0.0;
			if (T.waveform().Type == Waveform::Type::TX) {
				currentwaveform[i] = T.waveform().TD_Waveform[i];
			}
			if (T.waveform().Type == Waveform::Type::RX) {
				voltagewaveform[i] = T.waveform().TD_Waveform[i];
			}
		}
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL nwindows(const gatdaem1d_system_handle handle, int32_t* windows)
{
	return api_entry_point([&]()->gatdaem1d_status {
		TDEmSystem& T = tdemsystem_reference(handle);
		if (!windows) throw std::runtime_error("argument windows is null");
		*windows = static_cast<int>(T.nWindows());
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL nturns(const gatdaem1d_system_handle handle, int32_t* turns) {
	return api_entry_point([&]()->gatdaem1d_status {
		if (!turns) throw std::runtime_error("argument turns is null");
		TDEmSystem& T = tdemsystem_reference(handle);
		*turns = static_cast<int>(T.transmitter().nTurns);
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL peakcurrent(const gatdaem1d_system_handle handle, double* peak_current)
{
	return api_entry_point([&]()->gatdaem1d_status {
		if (!peak_current) throw std::runtime_error("argument turns is null");
		TDEmSystem& T = tdemsystem_reference(handle);
		*peak_current = static_cast<double>(T.transmitter().PeakCurrent);
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL looparea(const gatdaem1d_system_handle handle, double* loop_area)
{
	return api_entry_point([&]()->gatdaem1d_status {
		TDEmSystem& T = tdemsystem_reference(handle);
		if (!loop_area) throw std::runtime_error("Argument loop_area is null");
		*loop_area = static_cast<double>(T.transmitter().LoopArea);
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL basefrequency(const gatdaem1d_system_handle handle, double* base_frequency)
{
	return api_entry_point([&]()->gatdaem1d_status {
		TDEmSystem& T = tdemsystem_reference(handle);
		if (!base_frequency) throw std::runtime_error("Argument base_frequency is null");
		*base_frequency = static_cast<double>(T.waveform().BaseFrequency);
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL windowtimes(const gatdaem1d_system_handle handle, double* wt_low, double* wt_high)
{
	return api_entry_point([&]()->gatdaem1d_status {
		TDEmSystem& T = tdemsystem_reference(handle);
		if (!wt_low) throw std::runtime_error("Argument wt_low is null");
		if (!wt_high) throw std::runtime_error("Argument wt_high is null");
		for (size_t i = 0; i < T.nWindows(); i++) {
			wt_low[i] = T.window(i).Low;
			wt_high[i] = T.window(i).High;
		}
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL forwardmodel(const gatdaem1d_system_handle handle, const gatdaem1d_system_geometry* pG, const gatdaem1d_earth* pE, gatdaem1d_response* pR)
{
	return api_entry_point([&]()->gatdaem1d_status {
		TDEmSystem& T = tdemsystem_reference(handle);
		if (!pE) throw std::runtime_error("Argument pE is null");
		if (!pG) throw std::runtime_error("Argument pG is null");
		if (!pR) throw std::runtime_error("Argument pR is null");
		check_response_struct(T, pR);
		TDEmGeometry G = convert_geometry_struct(pG);
		Earth1D E = convert_earth_struct(pE);
		TDEmResponse R = T.forward_model(E, G);
		copy_response(R, pR);
		return GATDAEM1D_STATUS_OK;
	});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL derivative(const gatdaem1d_system_handle handle, const gatdaem1d_system_geometry* pG, const gatdaem1d_earth* pE, const char* dtype, int32_t dlayer, gatdaem1d_response* pR) {
	return api_entry_point([&]()->gatdaem1d_status {
		TDEmSystem& T = tdemsystem_reference(handle);
		if (!pE) throw std::runtime_error("Argument pE is null");
		if (!pG) throw std::runtime_error("Argument pG is null");
		if (!pR) throw std::runtime_error("Argument pR is null");
		if (!dtype) throw std::runtime_error("Argument dtype is null");
		check_response_struct(T, pR);
		TDEmGeometry G = convert_geometry_struct(pG);

		CMode cmode = CalculationType::get_mode_from_string(dtype);
		if (cmode == CMode::NONE) {
			std::string s = CalculationType::possible_values_message();
			std::ostringstream oss;
			oss << "Bad derivative type: (" << dtype << ")" << std::endl;
			oss << s << std::endl;
			throw std::runtime_error(oss.str());
			return GATDAEM1D_STATUS_ERRCPPEXCEPTION;
		}
		else if (cmode == CMode::DC) {
			if (dlayer < 0) throw std::runtime_error("Argument dtype='DC' and dlayer<0 [ensure 0 <= dlayer <= nLayers-1]");
			if (dlayer >= T.lem().nLayers()) throw std::runtime_error("Argument dtype='DC' and dlayer>=nLayers [ensure 0 <= dlayer <= nLayers-1]");
		}
		else if (cmode == CMode::DT) {
			if (dlayer < 0) throw std::runtime_error("Argument dtype='DT' and dlayer<0 [ensure 0 <= dlayer <= nLayers-2]");
			if (dlayer >= T.lem().nLayers() - 1) throw std::runtime_error("Argument dtype='DT' and dlayer>=nLayers-1 [ensure 0 <= dlayer <= nLayers-2]");
		}

		CalculationType ctype(cmode, dlayer);
		TDEmResponse<double> R = T.derivative(G, ctype);
		copy_response(R, pR);
		return GATDAEM1D_STATUS_OK;
		});
};

//Forward model and log-base-e layer conductivity derivatives for all layers. The primary field derivatives will be zero of course.
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL fm_dlogc(const gatdaem1d_system_handle handle, const gatdaem1d_system_geometry* pG, const gatdaem1d_earth* pE, double* outpur_buffer)
{
	return api_entry_point([&]()->gatdaem1d_status {
		if (!pE) throw std::runtime_error("Argument pE is null");
		if (!pG) throw std::runtime_error("Argument pG is null");
		if (!outpur_buffer) throw std::runtime_error("Argument outpur_buffer is null");

		TDEmSystem& T = tdemsystem_reference(handle);
		TDEmGeometry G = convert_geometry_struct(pG);
		Earth1D E = convert_earth_struct(pE);

		const size_t nw = T.nWindows();
		const size_t nl = pE->nlayers;

		//Innitial p to start of the output buffer
		double* p = outpur_buffer;

		//Run forward model and copy to end of output buffer
		TDEmResponse<double> R = T.forward_model(E, G);
		p = copy_bytes_array_end_advance_pointer(R, p);

		//Run layer derivatives and copy to end of output buffer
		for (size_t k = 0; k < nl; k++) {
			CalculationType ctype(CMode::DC, k);
			R = T.derivative(G, ctype);
			const double& c = pE->conductivity[k];
			R.P *= c; // Must scale by conductivity to get derivatibe w.r.t log-base-e(conductivity)
			R.S *= c;
			p = copy_bytes_array_end_advance_pointer(R, p);
		}
		return GATDAEM1D_STATUS_OK;
		});
};

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL iptype_from_string(const char* iptype_string, gatdaem1d_iptype* iptype) {
	return api_entry_point([&]()->gatdaem1d_status {
		if (!iptype_string) throw std::runtime_error("Argument iptype_string is null");
		if (!iptype) throw std::runtime_error("Argument iptype is null");
		AEM::IPType type = AEM::iptype_from_string(iptype_string);
		if (type == IPType::NONE) *iptype = GATDAEM1D_IPTYPE_NONE;
		else if (type == IPType::COLECOLE) *iptype = GATDAEM1D_IPTYPE_COLECOLE;
		else if (type == IPType::PELTON) *iptype = GATDAEM1D_IPTYPE_PELTON;
		else throw std::runtime_error("Argument iptype_string is not a valid IPType string");
		return GATDAEM1D_STATUS_OK;
	});
};



