/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _tdem1Dmodel_H
#define _tdem1Dmodel_H

#include "tdemsystem.h"
////////////////////////////////////////////////////////////////////////////

extern void* tdem_create(const char* systemdescriptorfile);
extern void  tdem_delete(void* S);
extern int   tdem_nwindows(void* S);
extern void  tdem_components(void* S, bool* xyz);
extern void  tdem_1Dmodel(void* S, sTDEmGeometry& G, cEarth1D& E, sTDEmResponse& R);
extern void  tdem_1Dderivative(void* S, eCalculationType ctype, size_t derivate_layerindex, sTDEmResponse& R);

#endif
