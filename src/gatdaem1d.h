#pragma once

//====  gatdaem1d.h - ABI macros ========================= */

//Note: GATDAEM1D_API_EXPORTS must ONLY be defined when building the shared library and NOT when consuming/using the already built shared library 
//      GATDAEM1D_STATIC_LINKAGE should ONLY be defined when linking to the static library

#include <stdint.h>

//Begin MACROS
#if defined(__cplusplus)
	#define GATDAEM1D_EXTERN_C extern "C"
#else
	#define GATDAEM1D_EXTERN_C extern
#endif

/* Sanity: only one of these should be defined */
#if defined(GATDAEM1D_API_EXPORTS) && defined(GATDAEM1D_STATIC_LINKAGE)
	#error "Define only one of GATDAEM1D_API_EXPORTS or GATDAEM1D_STATIC_LINKAGE"
#endif

/* Platform */
#if defined(_WIN32) || defined(__CYGWIN__)
	#define GATDAEM1D_PLATFORM_WINDOWS 1
#else
	#define GATDAEM1D_PLATFORM_WINDOWS 0
#endif

/* Export/import attributes */
#if defined(GATDAEM1D_STATIC_LINKAGE)
	/* Static library: no import/export attributes needed */
	#define GATDAEM1D_EXPORT
#else
	/* Shared library */
	#if GATDAEM1D_PLATFORM_WINDOWS
		#if defined(GATDAEM1D_API_EXPORTS)
			#define GATDAEM1D_EXPORT __declspec(dllexport)
		#else
			#define GATDAEM1D_EXPORT __declspec(dllimport)
		#endif
	#else
		/* Linux/macOS: export marked symbols (useful with -fvisibility=hidden) */
		#if defined(__GNUC__) && __GNUC__ >= 4
			#define GATDAEM1D_EXPORT __attribute__((visibility("default")))
		#else
			#define GATDAEM1D_EXPORT
		#endif
	#endif
#endif

// Public API linkage/export macro
#define GATDAEM1D_API GATDAEM1D_EXTERN_C GATDAEM1D_EXPORT

/* Calling convention */
/* Calling convention (relevant mainly on Windows; Pascal users care).
   Default: cdecl. To force stdcall define:
		#define GATDAEM1D_USE_STDCALL 1
		before including this header.
*/
#if GATDAEM1D_PLATFORM_WINDOWS
	#if defined(GATDAEM1D_USE_STDCALL) && (GATDAEM1D_USE_STDCALL)
		#define GATDAEM1D_CALL __stdcall
	#else
		#define GATDAEM1D_CALL __cdecl
	#endif
#else
	#define GATDAEM1D_CALL
#endif

//End MACROS

 // Opaque system handle value for FFIs
typedef uint64_t gatdaem1d_system_handle;

// ABI version: change this after changing the ABI (struct layouts, calling conventions, etc.). 
//    This is not the same as the version of ga_aem
#define GATDAEM1D_ABI_VERSION 100
// Expected ABI version: for runtime checking DLL matches the external wrapper classes (e.g. gatdaem1d_system_handle.m)
//    Users must not just alter this to match the ABI version of whatever DLL they have.
typedef enum GATDAEM1D_EXPECTED_ABI {
	GATDAEM1D_EXPECTED_ABI_VERSION = GATDAEM1D_ABI_VERSION
} GATDAEM1D_EXPECTED_ABI;

// Status codes
typedef int32_t gatdaem1d_status;
typedef enum GATDAEM1D_STATUS {
	GATDAEM1D_STATUS_OK = 0,
	GATDAEM1D_STATUS_ERRBADHANDLE = 1,
	GATDAEM1D_STATUS_ERRBADMAGIC = 2,
	GATDAEM1D_STATUS_ERRNULLSYS = 3,
	GATDAEM1D_STATUS_ERRCPPEXCEPTION = 4
} GATDAEM1D_STATUS;

typedef int32_t gatdaem1d_iptype;
//These must match the order/value in induced_polarization.hpp {enum class IPType { NONE, COLECOLE, PELTON };}
typedef enum GATDAEM1D_IPTYPE {
	GATDAEM1D_IPTYPE_NONE = 0,
	GATDAEM1D_IPTYPE_COLECOLE = 1,
	GATDAEM1D_IPTYPE_PELTON = 2
} GATDAEM1D_IPTYPE;

typedef struct gatdaem1d_system_geometry {
	double tx_height;
	double tx_roll;
	double tx_pitch;
	double tx_yaw;
	double txrx_dx;
	double txrx_dy;
	double txrx_dz;
	double rx_roll;
	double rx_pitch;
	double rx_yaw;
} gatdaem1d_system_geometry;

typedef struct gatdaem1d_earth {
	gatdaem1d_iptype iptype;
	int32_t nlayers;
	double* thickness;
	double* conductivity;
	double* chargeability;
	double* timeconstant;
	double* frequencydependence;
} gatdaem1d_earth;

typedef struct gatdaem1d_response {
	int32_t nwindows;
	double* PX;
	double* PY;
	double* PZ;
	double* SX;
	double* SY;
	double* SZ;
} gatdaem1d_response;

// Function protypes defined as follows
// <linkage/export-type> <return-type> <calling-convention> function(...)
// GATDAEM1D_API         <return-type> GATDAEM1D_CALL       function(...)

GATDAEM1D_API int32_t GATDAEM1D_CALL gatdaem1d_get_abi_version();
GATDAEM1D_API const char* GATDAEM1D_CALL gatdaem1d_get_last_error_string(const gatdaem1d_status s);

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL createhandle(const char* systemfile, gatdaem1d_system_handle* ptr_handle);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL deletehandle(gatdaem1d_system_handle* ptr_handle);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL nsamplesperwaveform(const gatdaem1d_system_handle handle, int32_t* n_samples);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL waveform(const gatdaem1d_system_handle handle, double* time, double* currentwaveform, double* voltagewaveform);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL nwindows(const gatdaem1d_system_handle handle, int32_t* n_windows);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL nturns(const gatdaem1d_system_handle handle, int32_t* n_turns);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL looparea(const gatdaem1d_system_handle handle, double* loop_area);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL basefrequency(const gatdaem1d_system_handle handle, double* base_frequency);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL peakcurrent(const gatdaem1d_system_handle handle, double* peak_current);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL windowtimes(const gatdaem1d_system_handle handle, double* low, double* high);

GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL forwardmodel(const gatdaem1d_system_handle handle, const gatdaem1d_system_geometry* pG, const gatdaem1d_earth* pE, gatdaem1d_response* pR);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL derivative(const gatdaem1d_system_handle handle, const gatdaem1d_system_geometry* pG, const gatdaem1d_earth* pE, const char* dtype, int32_t dlayer, gatdaem1d_response* pR);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL fm_dlogc(const gatdaem1d_system_handle handle, const gatdaem1d_system_geometry* pG, const gatdaem1d_earth* pE, double* output_buffer);
GATDAEM1D_API gatdaem1d_status GATDAEM1D_CALL iptype_from_string(const char* iptype_string, gatdaem1d_iptype* iptype);


