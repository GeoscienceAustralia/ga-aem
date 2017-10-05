/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <windows.h>
#include <gdiplus.h>
#include <GdiPlusImageCodec.h>
using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")

#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cstring>
#include "general_types.h"
#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "geometry3d.h"

#define VERSION "1.0"

int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes

	ImageCodecInfo* pImageCodecInfo = NULL;

	GetImageEncodersSize(&num, &size);
	if (size == 0)
		return -1;  // Failure

	pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
	if (pImageCodecInfo == NULL)
		return -1;  // Failure

	GetImageEncoders(num, size, pImageCodecInfo);

	for (UINT j = 0; j < num; ++j)
	{
		if (wcscmp(pImageCodecInfo[j].MimeType, format) == 0)
		{
			*pClsid = pImageCodecInfo[j].Clsid;
			free(pImageCodecInfo);
			return j;  // Success
		}
	}

	free(pImageCodecInfo);
	return -1;  // Failure
}

class cMap {

public:
	std::vector<unsigned char> r;
	std::vector<unsigned char> g;
	std::vector<unsigned char> b;

	cMap(){
		setjet();
	}

	cMap(const std::string name){		
		set(name);		
	}

	void set(const std::string name){

		if(strcasecmp(name, "jet")==0){
			setjet();
		}
		else if(strcasecmp(name, "pseudocolor")==0){
			setpseudocolor();
		}
		else if(strcasecmp(name, "pseudocolour")==0){
			setpseudocolor();
		}
		else if(strcasecmp(name, "geosoft")==0){
			setgeosoft();
		}
		else if(strcasecmp(name, "rainbow1")==0){
			setrainbow1();
		}
		else if(strcasecmp(name, "rainbow2")==0){
			setrainbow2();
		}
		else if(strcasecmp(name, "comtal")==0){
			setcomtal();
		}
		else{
			message("Unknown ColourMap %s\n",name.c_str());			
		}		
	}

	void set(const unsigned char* ar, const unsigned char* ag, const unsigned char* ab, size_t n){
		r = std::vector<unsigned char>(ar,ar+256);
		g = std::vector<unsigned char>(ag,ag+256);
		b = std::vector<unsigned char>(ab,ab+256);
	}

	void setjet(){
		static const unsigned char ar[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 131, 135, 139, 143, 147, 151, 155, 159, 163, 167, 171, 175, 179, 183, 187, 191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235, 239, 243, 247, 251, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 251, 247, 243, 239, 235, 231, 227, 223, 219, 215, 211, 207, 203, 199, 195, 191, 187, 183, 179, 175, 171, 167, 163, 159, 155, 151, 147, 143, 139, 135, 131, 128 };
		static const unsigned char ag[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 131, 135, 139, 143, 147, 151, 155, 159, 163, 167, 171, 175, 179, 183, 187, 191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235, 239, 243, 247, 251, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 251, 247, 243, 239, 235, 231, 227, 223, 219, 215, 211, 207, 203, 199, 195, 191, 187, 183, 179, 175, 171, 167, 163, 159, 155, 151, 147, 143, 139, 135, 131, 128, 124, 120, 116, 112, 108, 104, 100, 96, 92, 88, 84, 80, 76, 72, 68, 64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		static const unsigned char ab[] = { 131, 135, 139, 143, 147, 151, 155, 159, 163, 167, 171, 175, 179, 183, 187, 191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235, 239, 243, 247, 251, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 251, 247, 243, 239, 235, 231, 227, 223, 219, 215, 211, 207, 203, 199, 195, 191, 187, 183, 179, 175, 171, 167, 163, 159, 155, 151, 147, 143, 139, 135, 131, 128, 124, 120, 116, 112, 108, 104, 100, 96, 92, 88, 84, 80, 76, 72, 68, 64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };				
		set(ar,ag,ab,256);	
	}

	void setpseudocolor(){
		static const unsigned char ar[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 128, 131, 134, 137, 140, 143, 146, 149, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191, 194, 197, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233, 236, 239, 242, 245, 248, 251, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };
		static const unsigned char ag[] = { 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 128, 131, 134, 137, 140, 143, 146, 149, 152, 155, 158, 161, 164, 167, 170, 173, 176, 179, 182, 185, 188, 191, 194, 197, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233, 236, 239, 242, 245, 248, 251, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 252, 249, 246, 243, 240, 237, 234, 231, 228, 225, 222, 219, 216, 213, 210, 207, 204, 201, 198, 195, 192, 189, 186, 183, 180, 177, 174, 171, 168, 165, 162, 159, 156, 153, 150, 147, 144, 141, 138, 135, 132, 129, 127, 124, 121, 118, 115, 112, 109, 106, 103, 100, 97, 94, 91, 88, 85, 82, 79, 76, 73, 70, 67, 64, 61, 58, 55, 52, 49, 46, 43, 40, 37, 34, 31, 28, 25, 22, 19, 16, 13, 10, 7, 4, 1 };
		static const unsigned char ab[] = { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 252, 249, 246, 243, 240, 237, 234, 231, 228, 225, 222, 219, 216, 213, 210, 207, 204, 201, 198, 195, 192, 189, 186, 183, 180, 177, 174, 171, 168, 165, 162, 159, 156, 153, 150, 147, 144, 141, 138, 135, 132, 129, 127, 124, 121, 118, 115, 112, 109, 106, 103, 100, 97, 94, 91, 88, 85, 82, 79, 76, 73, 70, 67, 64, 61, 58, 55, 52, 49, 46, 43, 40, 37, 34, 31, 28, 25, 22, 19, 16, 13, 10, 7, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		set(ar,ag,ab,256);	
	}

	void setgeosoft(){
		static const unsigned char ar[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 16, 27, 37, 48, 59, 69, 75, 79, 83, 87, 91, 95, 99, 101, 102, 104, 105, 107, 108, 111, 116, 122, 127, 132, 138, 143, 148, 154, 160, 165, 171, 176, 182, 187, 192, 198, 203, 208, 214, 219, 225, 230, 236, 241, 247, 252, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };
		static const unsigned char ag[] = { 0, 13, 25, 38, 51, 63, 76, 87, 93, 99, 106, 112, 118, 124, 131, 137, 144, 150, 156, 163, 169, 175, 182, 188, 194, 200, 207, 212, 216, 219, 222, 225, 228, 231, 234, 238, 241, 244, 248, 251, 254, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 253, 250, 247, 244, 240, 237, 234, 231, 227, 224, 221, 218, 215, 212, 209, 205, 202, 199, 196, 193, 190, 189, 187, 185, 184, 182, 181, 179, 178, 176, 175, 173, 172, 170, 167, 164, 160, 157, 154, 151, 148, 146, 144, 143, 141, 139, 138, 136, 135, 133, 132, 130, 129, 127, 124, 121, 118, 115, 112, 109, 106, 103, 99, 96, 93, 90, 87, 83, 78, 74, 69, 64, 59, 55, 50, 45, 40, 35, 31, 26, 21, 18, 15, 12, 9, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 5, 7, 8, 10, 18, 34, 51, 67, 83, 100, 116, 125, 131, 136, 142, 148, 153, 159 };
		static const unsigned char ab[] = { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 249, 241, 233, 224, 216, 208, 200, 192, 183, 175, 167, 159, 151, 141, 129, 117, 104, 92, 80, 68, 62, 60, 58, 55, 53, 51, 49, 47, 45, 43, 41, 40, 38, 35, 30, 24, 19, 14, 8, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 10, 19, 27, 35, 43, 51, 59, 67, 76, 84, 92, 100, 108, 118, 129, 140, 151, 162, 173, 183, 188, 193, 199, 204, 210, 215, 220, 226, 231, 237, 242, 248, 253, 255, 255, 255, 255, 255, 255, 255 };
		set(ar,ag,ab,256);	
	}

	void setrainbow1(){
		static const unsigned char ar[] = { 254, 254, 248, 242, 235, 229, 223, 217, 210, 204, 198, 192, 186, 179, 173, 167, 161, 154, 148, 142, 136, 130, 124, 118, 112, 106, 100, 93, 87, 81, 75, 68, 62, 56, 50, 44, 37, 31, 25, 19, 12, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 12, 18, 24, 30, 36, 42, 49, 55, 61, 67, 73, 79, 85, 91, 97, 103, 109, 115, 121, 127, 133, 139, 145, 151, 157, 163, 169, 175, 181, 187, 193, 199, 205, 212, 218, 224, 230, 236, 242, 248, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254 };
		static const unsigned char ag[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 12, 18, 24, 30, 36, 42, 49, 55, 61, 67, 73, 79, 85, 91, 97, 103, 109, 115, 121, 127, 133, 139, 145, 151, 157, 163, 169, 175, 181, 187, 193, 199, 205, 212, 218, 224, 230, 236, 242, 248, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 248, 242, 237, 231, 225, 219, 213, 208, 202, 196, 190, 184, 179, 173, 167, 161, 155, 150, 144, 138, 132, 128, 122, 116, 110, 104, 99, 93, 87, 81, 75, 70, 64, 58, 52, 46, 41, 35, 29, 23, 17, 12, 6, 0, 6, 12, 18, 24, 30, 36, 42, 49, 55, 61, 67, 73, 79, 85, 91, 97, 103, 109, 115, 121, 127, 133, 139, 145, 151, 157, 163, 169, 175, 181, 187, 193, 199, 205, 212, 218, 224, 230, 236, 242, 248, 254 };
		static const unsigned char ab[] = { 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 248, 242, 236, 230, 224, 218, 212, 207, 201, 195, 189, 183, 177, 171, 165, 159, 153, 147, 141, 135, 129, 125, 119, 113, 107, 101, 95, 89, 83, 77, 71, 65, 59, 53, 47, 42, 36, 30, 24, 18, 12, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 12, 18, 24, 30, 36, 42, 49, 55, 61, 67, 73, 79, 85, 91, 97, 103, 109, 115, 121, 127, 133, 139, 145, 151, 157, 163, 169, 175, 181, 187, 193, 199, 205, 212, 218, 224, 230, 236, 242, 248, 254 };
		set(ar,ag,ab,256);	
	}

	void setrainbow2(){
		static const unsigned char ar[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 8, 12, 16, 20, 24, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 78, 82, 86, 90, 94, 98, 102, 106, 110, 114, 118, 122, 126, 128, 132, 136, 140, 144, 148, 152, 156, 160, 164, 168, 172, 176, 179, 183, 187, 191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 230, 234, 238, 242, 246, 250, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254 };
		static const unsigned char ag[] = { 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 130, 134, 138, 142, 146, 150, 154, 158, 162, 166, 170, 174, 178, 182, 186, 190, 194, 198, 202, 206, 210, 214, 218, 222, 226, 230, 234, 238, 242, 246, 250, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 250, 246, 242, 238, 234, 230, 226, 222, 218, 214, 209, 205, 201, 197, 193, 189, 185, 181, 177, 173, 169, 165, 161, 157, 153, 149, 145, 141, 137, 133, 129, 125, 121, 117, 113, 109, 105, 101, 97, 93, 89, 85, 81, 77, 73, 69, 65, 61, 57, 53, 49, 45, 40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0 };
		static const unsigned char ab[] = { 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 250, 246, 242, 238, 234, 230, 226, 222, 218, 214, 210, 206, 202, 198, 194, 190, 186, 182, 178, 174, 170, 166, 162, 158, 154, 150, 146, 142, 138, 134, 130, 128, 124, 120, 116, 112, 108, 104, 100, 96, 92, 88, 84, 80, 76, 72, 68, 64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		set(ar,ag,ab,256);	
	}

	void setcomtal(){
		static const unsigned char ar[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 12, 18, 24, 30, 36, 43, 49, 56, 62, 68, 75, 81, 88, 94, 107, 120, 132, 145, 158, 164, 171, 177, 184, 190, 196, 203, 209, 216, 222, 229, 235, 242, 248, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 249, 242, 236, 229, 223, 217, 210, 204, 197, 191, 185, 178, 172, 165, 159, 153, 146, 140, 133, 127, 123, 120, 116, 113, 109, 105, 102, 98, 95, 91, 94, 97, 101, 104, 107, 111, 115, 119, 123, 127, 131, 135, 140, 144, 148, 153, 159, 164, 170, 175, 180, 185, 189, 194, 199, 205, 210, 216, 221, 227, 233, 238, 244, 249, 255, 239, 223, 207, 191, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 202, 228, 255 };
		static const unsigned char ag[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 36, 55, 73, 91, 98, 105, 113, 120, 127, 133, 140, 146, 153, 159, 165, 172, 178, 185, 191, 198, 205, 212, 219, 226, 232, 238, 243, 249, 255, 255, 255, 255, 255, 255, 251, 247, 244, 240, 236, 229, 222, 215, 208, 201, 191, 182, 172, 163, 153, 148, 143, 139, 134, 129, 125, 121, 117, 113, 109, 113, 116, 120, 123, 127, 130, 133, 137, 140, 143, 146, 149, 152, 155, 158, 164, 171, 177, 184, 190, 196, 203, 209, 216, 222, 229, 235, 242, 248, 255, 249, 242, 236, 229, 223, 217, 210, 204, 197, 191, 185, 178, 172, 165, 159, 159, 159, 159, 159, 159, 165, 172, 178, 185, 191, 197, 204, 210, 217, 223, 229, 236, 242, 249, 255, 248, 242, 235, 229, 222, 216, 209, 203, 196, 190, 184, 178, 172, 166, 160, 153, 146, 140, 133, 126, 119, 112, 105, 98, 91, 73, 55, 36, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 51, 76, 102, 127, 130, 133, 137, 140, 143, 149, 156, 162, 169, 175, 202, 228, 255 };
		static const unsigned char ab[] = { 0, 8, 17, 25, 50, 76, 101, 126, 139, 152, 164, 177, 190, 203, 216, 229, 242, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 242, 229, 217, 204, 191, 186, 181, 175, 170, 165, 159, 152, 146, 139, 133, 126, 119, 111, 104, 97, 92, 87, 83, 78, 73, 65, 56, 48, 39, 31, 28, 25, 21, 18, 15, 12, 9, 6, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13, 25, 38, 50, 63, 63, 63, 63, 63, 63, 72, 81, 89, 98, 107, 111, 115, 119, 123, 127, 131, 136, 140, 145, 149, 154, 159, 165, 170, 175, 180, 185, 189, 194, 199, 205, 210, 216, 221, 227, 233, 238, 244, 249, 255, 239, 223, 207, 191, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 202, 228, 255 };
		set(ar,ag,ab,256);	
	}

};

class cGeorefSection{

private:

	Bitmap* pBitmap;
	cMap m_cmap;
	bool m_drawcbar;
	std::vector<double> m_cbarticks;
	int nhpixels;
	int nvpixels;
	double hlength;
	double vlength;
	double h0,h1,h2,dh;
	double v0,v1,v2,dv;
	double x0,x1,x2,dx;
	double y0,y1,y2,dy;
	double gratdiv;

	std::string outdir;
	std::string prefix;
	std::string suffix;

	int nlayers;
	int nsamples;
	int linenumber;


	double LowClip;
	double HighClip;
	bool   Log10Stretch;

	bool SaveJPG;
	bool SavePNG;
	bool SaveEMF;

	double elevation_median;
	bool autozsectiontop;
	bool autozsectionbot;
	double zsectiontop;
	double zsectionbot;

public:
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> e;
	std::vector<std::vector<double>> c;
	std::vector<std::vector<double>> z;
	std::vector<std::vector<double>> cp10;
	std::vector<std::vector<double>> cp90;
	bool spreadfade;	
	double spreadfadelowclip;
	double spreadfadehighclip;
	double lowspreadfade;
	double highspreadfade;
	
	void getoptions(cBlock b){
		
		spreadfade = b.getboolvalue("SpreadFade");
		if(spreadfade){
			spreadfadelowclip   = b.getdoublevalue("Log10SpreadLowClip");
			spreadfadehighclip  = b.getdoublevalue("Log10SpreadHighClip");
			lowspreadfade       = b.getdoublevalue("LowSpreadFade");
			highspreadfade      = b.getdoublevalue("HighSpreadFade");	
		}


		double v = b.getdoublevalue("ElevationTop");
		if(isundefined(v)){
			autozsectiontop = true;
			zsectiontop = v;
		}
		else{
			autozsectiontop = false;
			zsectiontop = v;
		}

		v = b.getdoublevalue("ElevationBottom");
		if (isundefined(v)){
			autozsectionbot = true;
			zsectionbot = v;
		}
		else{
			autozsectionbot = false;
			zsectionbot = v;
		}

		outdir = b.getstringvalue("OutDir");
		prefix = b.getstringvalue("Prefix");
		suffix = b.getstringvalue("Suffix");

		dh    = b.getdoublevalue("HorizontalResolution");
		double vxg   = b.getdoublevalue("VerticalExaggeration");
		dv = dh/vxg;

		std::string cmapname = b.getstringvalue("ColourMap");
		m_cmap.set(cmapname);

		m_drawcbar  = b.getboolvalue("ColourBar");
		if(m_drawcbar){
			m_cbarticks = b.getdoublevector("ColourBarTicks");
		}

		LowClip  = b.getdoublevalue("LowClip");
		HighClip = b.getdoublevalue("HighClip");
		Log10Stretch  = b.getboolvalue("Log10Stretch");

		gratdiv = b.getdoublevalue("ElevationGridDivision");

		SaveJPG = b.getboolvalue("SaveJPG");
		SavePNG = b.getboolvalue("SavePNG");
		SaveEMF = b.getboolvalue("SaveEMF");
	}

	void readdatafile(const cBlock& input, const std::string filename)
	{
		int subsample = input.getintvalue("Subsample");
		if (isundefined(subsample))subsample = 1;

		std::string lstr = input.getstringvalue("Line");
		std::string xstr = input.getstringvalue("Easting");
		std::string ystr = input.getstringvalue("Northing");
		std::string estr = input.getstringvalue("Elevation");
		
		bool isresistivity = false;
		std::string crstr = input.getstringvalue("Conductivity");		
		if (isundefined(crstr)){
			crstr = input.getstringvalue("Resistivity");
			isresistivity = true;
		}

		std::string cp10str = input.getstringvalue("Conductivity_p10");
		std::string cp90str = input.getstringvalue("Conductivity_p90");								

		double cscale=1.0;
		std::string cunits = input.getstringvalue("InputConductivityUnits");
		if(isundefined(cunits)){
			cscale = 1.0;
		}
		else if(strcasecmp(cunits,"S/m") == 0){
			cscale = 1.0;
		}
		else if(strcasecmp(cunits,"mS/m") == 0){
			cscale = 0.001;
		}		
		else{
			message("Unknown InputConductivityUnits %s\n",cunits.c_str());			
		}

		int lcol, xcol, ycol, ecol;
		int crcol1, crcol2, tcol1, tcol2;
		int cp10col1, cp10col2;
		int cp90col1, cp90col2;
		sscanf(lstr.c_str(), "Column %d", &lcol); lcol--;
		sscanf(xstr.c_str(), "Column %d", &xcol); xcol--;
		sscanf(ystr.c_str(), "Column %d", &ycol); ycol--;
		sscanf(estr.c_str(), "Column %d", &ecol); ecol--;

		sscanf(crstr.c_str(), "Column %d-%d", &crcol1, &crcol2); crcol1--; crcol2--;
		nlayers = crcol2 - crcol1 + 1;

		if(spreadfade){
			sscanf(cp10str.c_str(), "Column %d-%d", &cp10col1, &cp10col2); cp10col1--; cp10col2--;
			sscanf(cp90str.c_str(), "Column %d-%d", &cp90col1, &cp90col2); cp90col1--; cp90col2--;
		}

		bool isconstantthickness = false;
		std::vector<double> constantthickness;
		std::string tstr = input.getstringvalue("Thickness");
		if(sscanf(tstr.c_str(), "Column %d-%d", &tcol1, &tcol2) == 2){
			tcol1--; tcol2--;	
			isconstantthickness = false;
		}
		else{			
			constantthickness = input.getdoublevector("Thickness");
			tcol1 = 0; tcol2 = 0;
			isconstantthickness = true;
			if(constantthickness.size() == 0){
				errormessage("Thickness not set\n");
			}
			else if(constantthickness.size() > 1 && constantthickness.size() < nlayers - 1){
				errormessage("Thickness not set correctly\n");
			}
			else if(constantthickness.size() == 1){
				constantthickness = std::vector<double>(nlayers-1, constantthickness[0]);
			}
			else{
				//all good
			}

		}

		FILE* fp = fileopen(filename,"r");
		std::string str;
		std::vector<std::vector<double>> M;

		int k=0;
		while (filegetline(fp, str)){			
			if(k%subsample == 0){
				M.push_back(getdoublevector(str.c_str(), " "));
			}	
			k++;
		}
		fclose(fp);


		nsamples   = (int)M.size();
		linenumber = (int)M[0][lcol];		

		x.resize(nsamples);
		y.resize(nsamples);
		e.resize(nsamples);
		z.resize(nsamples);
		c.resize(nsamples);
		if(spreadfade){
			cp10.resize(nsamples);
			cp90.resize(nsamples);		
		}
		for (int si = 0; si < nsamples; si++){
			x[si] = M[si][xcol];
			y[si] = M[si][ycol];
			e[si] = M[si][ecol];

			c[si].resize(nlayers);
			for (int li = 0; li < nlayers; li++){
				c[si][li] = M[si][crcol1 + li];				
				if(c[si][li] > 0.0){
					if(isresistivity)c[si][li] = 1.0 / c[si][li];
					c[si][li] *= cscale;
				}
			}

			if(spreadfade){
				cp10[si].resize(nlayers);
				for (int li = 0; li < nlayers; li++){
					cp10[si][li] = M[si][cp10col1 + li];
					if(cp10[si][li] > 0.0){
						cp10[si][li] *= cscale;
					}
				}

				cp90[si].resize(nlayers);
				for (int li = 0; li < nlayers; li++){
					cp90[si][li] = M[si][cp90col1 + li];
					if(cp90[si][li] > 0.0){
						cp90[si][li] *= cscale;
					}
				}
			}

			z[si].resize(nlayers+1);			
			z[si][0] = e[si];
			for (int li = 0; li < nlayers; li++){				

				double t;
				if (li < nlayers - 1){
					if (isconstantthickness == true){
						t = constantthickness[li];
					}
					else{
						t = M[si][tcol1 + li];
					}
				}		
				else{
					if (isconstantthickness == true){
						t = constantthickness[li-1];
					}
					else{
						t = M[si][tcol1 + li - 1];
					}					
				}	

				z[si][li + 1] = z[si][li] - t;

			}
		}		
	}

	void process(){
		calculateextents();
		createbitmap();
		generatesectiondata();
		drawgraticule();
		drawcolorbar();		
		fixtransparenttext();
		save();
		deletebitmap();
	}

	void calculateextents()
	{		
		bestfitlineendpoints(x,y,x1,y1,x2,y2);
		double angle = atan2(y2-y1,x2-x1);
		if(angle > PIONTWO || angle <= -PIONTWO){
			std::swap(x1,x2);
			std::swap(y1,y2);			
		}

		double rawlength = distance(x1,y1,x2,y2);
		double hlength = roundupnearest(rawlength,dh);

		x2 = x1 + (x2-x1)*hlength/rawlength;
		y2 = y1 + (y2-y1)*hlength/rawlength;

		int hmarginpixels = getmarginpixels();
		h0 = 0.0;
		h1 = hmarginpixels * dh;
		h2 = h1+hlength;

		nhpixels  = 1 + (int)((h2-h0)/dh);
		dx = dh*(x2-x1)/(h2-h1);
		dy = dh*(y2-y1)/(h2-h1);		
		x0 = x1 - (x2-x1)*h1/hlength;
		y0 = y1 - (y2-y1)*h1/hlength;

		double zmin =  DBL_MAX;
		double zmax = -DBL_MAX;
		for(int si=0; si<nsamples; si++){
			if(z[si][nlayers] < zmin)zmin = z[si][nlayers];
			if(z[si][0]       > zmax)zmax = z[si][0];
		}

		if(autozsectionbot==false){
			zmin = zsectionbot;
		}

		if(autozsectiontop==false){
			zmax = zsectiontop;
		}

		elevation_median = median(&e[0],e.size());

		double rawheight = zmax-zmin;
		vlength = roundupnearest(rawheight,dv);		
		v0 = zmin;
		v1 = zmin;
		v2 = v1 + vlength;
		nvpixels = 1 + (int)(vlength/dv);		
	}

	int wh2lx(double wh){
		int lx = (int)roundnearest((wh-h0)/dh,1.0);
		return lx;
	}

	int wv2ly(double wv){		
		int ly = (int)roundnearest((v2-wv)/dv,1.0);
		return ly;
	}

	int getmarginpixels(){
		int graticulemargin = 65;
		int colourbarmargin = 95;

		int margin = graticulemargin;
		if(m_drawcbar) margin += colourbarmargin;
		return margin;
	}

	void createbitmap(){
		pBitmap = new Bitmap(nhpixels, nvpixels, PixelFormat32bppARGB);		
		for(int i=0; i<nhpixels; i++){
			for(int j=0; j<nvpixels; j++){				
				pBitmap->SetPixel(i,j,Color(255,255,255,255));				
			}
		}		
	}

	void deletebitmap(){
		delete pBitmap;
	}

	double fadevalue(const double& spread){
		if(spread <= spreadfadelowclip)return  lowspreadfade;
		if(spread >= spreadfadehighclip)return highspreadfade;
		return lowspreadfade + (spread - spreadfadelowclip)/(spreadfadehighclip - spreadfadelowclip) * (highspreadfade - lowspreadfade);		
	}

	void fadecolor(Color& clr, const double& fade){
		//fade=0 gives original color
		//fade=1 gives white
		double alpha = 1.0-fade;
		unsigned char w = (unsigned char)(255.0*fade);		
		unsigned char r = w + (unsigned char)(clr.GetRed()   * alpha);
		unsigned char g = w + (unsigned char)(clr.GetGreen() * alpha);
		unsigned char b = w + (unsigned char)(clr.GetBlue()  * alpha);		
		clr = Color(255,r,g,b);				
	}

	void generatesectiondata(){

		Color BkgColor(255,128,128,128);
		Color AirColor(0,255,255,255);
		Color NullsColor(255,128,128,128);

		int hp1 = wh2lx(h1);
		int hp2 = wh2lx(h2);
		int vp1 = wv2ly(v1);
		int vp2 = wv2ly(v2);
		for(int i=hp1; i<=hp2; i++){
			double xp = x0 + (double)i*dx;
			double yp = y0 + (double)i*dy;

			double mind=DBL_MAX;
			int mini;
			for(int si=0; si<nsamples; si++){
				double d = distance(0.0,0.0,xp-x[si],yp-y[si]);
				if(d<mind){
					mini = si;
					mind = d;
				}
			}


			for(int j=vp2; j<=vp1; j++){				
				pBitmap->SetPixel(i,j,BkgColor);

				double zp = v2 - (double)j*dv;
				if(zp>e[mini]){
					pBitmap->SetPixel(i,j,AirColor);					
				}
				else{					
					for(int li=0; li<nlayers; li++){
						if(zp < z[mini][li] && zp >= z[mini][li+1]){
							double conductivity = c[mini][li];							
							if(conductivity < 0.0){
								pBitmap->SetPixel(i,j,NullsColor);
							}
							else{							
								int ind;
								if(Log10Stretch) ind = log10stretch(conductivity,LowClip,HighClip);							
								else ind = linearstretch(conductivity,LowClip,HighClip);

								Color clr(255,m_cmap.r[ind],m_cmap.g[ind],m_cmap.b[ind]);

								if(spreadfade){
									double conductivity_p10 = cp10[mini][li];
									double conductivity_p90 = cp90[mini][li];
									double spread = log10(conductivity_p90) - log10(conductivity_p10);
									double fade   = fadevalue(spread);
									fadecolor(clr,fade);
								}								
								pBitmap->SetPixel(i,j,clr);
							}
							break;
						}						
					}										
				}//layer loop
			}//v pixel loop
		}//h pixel loop		
	}

	REAL getfontsize(){
		//REAL fontsize = (REAL)abs((wv2ly(0.33*gratdiv) - wv2ly(0)));
		REAL fontsize = 8.5;
		return fontsize;
	}

	void drawgraticule(){

		Graphics gr(pBitmap);
		gr.SetCompositingMode(CompositingModeSourceOver);						
		gr.SetTextRenderingHint(TextRenderingHintAntiAlias);		

		int zg1=(int)roundnearest(v1,(int)gratdiv);
		int zg2=(int)roundnearest(v2,(int)gratdiv);

		Pen blackpen(Color::Black,0);		
		SolidBrush blackbrush(Color::Black);
		Font font(L"Arial", getfontsize(), FontStyleRegular, UnitPoint);				

		StringFormat textformat;
		textformat.SetAlignment(StringAlignmentFar);
		textformat.SetLineAlignment(StringAlignmentCenter);
		SizeF layoutsize(32767,32767);
		SizeF textsize;				
		gr.MeasureString(L"-0000 m",-1,&font,layoutsize,&textformat,&textsize);

		for(int zg=zg1; zg<=zg2; zg+=(int)gratdiv){
			gr.DrawLine(&blackpen,wh2lx(h1),wv2ly(zg),wh2lx(h2),wv2ly(zg));

			wchar_t s[20];				
			swprintf(s,20,L"%5d m",zg);						
			gr.MeasureString(s,-1,&font,layoutsize,&textformat,&textsize);

			int tx = wh2lx(h1);
			int ty = wv2ly(zg);
			if(ty > nvpixels)continue;
			if(ty-textsize.Height < 0)continue;	
			PointF txpos((REAL)tx,(REAL)ty);
			gr.DrawString(s, -1, &font, txpos, &textformat, &blackbrush);
		}		    				
	}

	void drawcolorbar(){

		if(m_drawcbar==false)return;

		Pen blackpen(Color::Black, 0);
		SolidBrush blackbrush(Color::Black);

		Graphics gr(pBitmap);

		TextRenderingHint hint = gr.GetTextRenderingHint();		
		gr.SetTextRenderingHint(TextRenderingHintAntiAlias);

		int ph1 = wv2ly(v1);
		int ph2 = wv2ly(v2);
		ph1 = 25; ph2 = 40;

		int pvtop = wv2ly(v1+(v2-v1)*0.95);
		int pvbot = wv2ly(v1+(v2-v1)*0.05);
		for (int j = pvtop; j <= pvbot; j++){			
			int ind = 255*(j - pvbot) / (pvtop - pvbot);
			for (int i = ph1; i <= ph2; i++){
				pBitmap->SetPixel(i, j, Color(255, m_cmap.r[ind], m_cmap.g[ind], m_cmap.b[ind]));
			}
		}
		gr.DrawRectangle(&blackpen, ph1, pvtop, ph2-ph1, pvbot-pvtop);

		Font font(L"Arial", getfontsize(), FontStyleRegular, UnitPoint);
		StringFormat textformat;		
		textformat.SetAlignment(StringAlignmentNear);
		textformat.SetLineAlignment(StringAlignmentCenter);
		SizeF layoutsize(32767, 32767);
		SizeF textsize;
		gr.MeasureString(L"0.0001 m", -1, &font, layoutsize, &textformat, &textsize);


		for (int i = 0; i<m_cbarticks.size(); i++){			
			if(m_cbarticks[i] < LowClip)continue;
			if(m_cbarticks[i] > HighClip)continue;

			int tickv;
			if (Log10Stretch){
				int ind = log10stretch(m_cbarticks[i], LowClip, HighClip);
				tickv = pvbot + ind*(pvtop - pvbot) / 255;
			}
			else{
				int ind = linearstretch(m_cbarticks[i], LowClip, HighClip);
				tickv = pvbot + ind*(pvtop - pvbot) / 255;
			}

			wchar_t s[20];
			swprintf(s, 20,L"%5.3lf", m_cbarticks[i]);

			PointF p;
			p.X = (REAL)ph2+3;
			p.Y = (REAL)tickv;
			gr.DrawString(s, -1, &font, p, &textformat, &blackbrush);
			gr.DrawLine(&blackpen,ph2-2,tickv,ph2+2,tickv);
		}		

		PointF p;
		p.X = (REAL)ph1-2;
		p.Y = (REAL)(pvtop + pvbot) / 2;

		textformat.SetAlignment(StringAlignmentCenter);
		textformat.SetLineAlignment(StringAlignmentFar);
		gr.TranslateTransform(p.X,p.Y);
		gr.RotateTransform(-90);		
		gr.DrawString(L"Conductivity (S/m)", -1, &font, PointF(0,0), &textformat, &blackbrush);

	}

	void fixtransparenttext(){	
		int hp0 = wh2lx(h0);
		int hp1 = wh2lx(h1);
		int vp1 = wv2ly(v1);
		int vp2 = wv2ly(v2);
		for(int i=hp0; i<=hp1; i++){		
			for(int j=vp2; j<=vp1; j++){
				Color c;
				pBitmap->GetPixel(i,j,&c);				
				pBitmap->SetPixel(i, j, Color(255, c.GetR(), c.GetG(), c.GetB()));
			}
		}		
	}

	void save(){
		std::string basename = prefix + strprint("%d", linenumber) + suffix;
		makedirectorydeep(outdir);

		if(SavePNG){		
			std::string imagepath = outdir + basename + ".png";						
			wchar_t wcimagepath[200];
			size_t len;			
			mbstowcs_s(&len, wcimagepath, 200, imagepath.c_str(), _TRUNCATE);

			CLSID pngClsid;
			int ret = GetEncoderClsid(L"image/png", &pngClsid);			
			Status result = pBitmap->Save(wcimagepath, &pngClsid, NULL);

			std::string worldfilepath = outdir + basename + ".pngw";
			saveworldfile(worldfilepath);
		}

		if(SaveJPG){			
			std::string imagepath = outdir + basename + ".jpg";
			wchar_t wcimagepath[500];
			size_t len;
			mbstowcs_s(&len, wcimagepath, 500, imagepath.c_str(), _TRUNCATE);

			CLSID jpgClsid;
			GetEncoderClsid(L"image/jpeg", &jpgClsid);
			Status result = pBitmap->Save(wcimagepath, &jpgClsid, NULL);

			std::string worldfilepath = outdir + basename + ".jgw";
			saveworldfile(worldfilepath);
		}

	}

	void saveworldfile(std::string& worldfilepath){

		//Line 1: A: x component of the pixel width (x-scale)
		//Line 2: D: y component of the pixel width (y-skew)
		//Line 3: B: x component of the pixel height (x-skew)
		//Line 4: E: y component of the pixel height (y-scale), almost always negative
		//Line 5: C: x-coordinate of the center of the upper left pixel
		//Line 6: F: y-coordinate of the center of the upper left pixel

		double angle  = atan2(y2-y1,x2-x1);
		double vshift = (v2 - elevation_median);
		double ix0    = x0 - (dh/dv)*vshift*sin(angle);
		double iy0    = y0 + (dh/dv)*vshift*cos(angle);

		FILE* fp=fileopen(worldfilepath,"w");
		fprintf(fp,"%lf\n",dh*cos(angle));
		fprintf(fp,"%lf\n",dh*sin(angle));
		fprintf(fp,"%lf\n",dh*sin(angle));
		fprintf(fp,"%lf\n",-dh*cos(angle));
		fprintf(fp,"%lf\n",ix0);
		fprintf(fp,"%lf\n",iy0);	
		fclose(fp);
	}

};


int main(int argc, char** argv)
{	
	GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);

	if (argc >= 2){
		message("Executing %s %s\n", argv[0], argv[1]);
		message("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		message("Working directory %s\n", getcurrentdirectory().c_str());
	}
	else{
		message("Executing %s\n", argv[0]);
		message("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		message("Working directory %s\n", getcurrentdirectory().c_str());
		message("Error: Not enough input arguments\n");
		message("Usage: %s controlfilename\n",argv[0]);		
		return 0;
	}

	cBlock b(argv[1]);
	cBlock inputblock   = b.findblock("Input");	
	cBlock sectionblock = b.findblock("Section");	
	std::string infiles = inputblock.getstringvalue("DataFiles");

	std::vector<std::string> filelist =  cDirectoryAccess::getfilelist(infiles);
	double t1 = gettime();	
	for (int i = 0; i < filelist.size(); i++){	
		printf("Processing file %s %3d of %3d\n", filelist[i].c_str(),i+1,filelist.size());
		cGeorefSection S;		
		S.getoptions(sectionblock);		
		S.readdatafile(inputblock,filelist[i].c_str());		
		S.process();		
	}
	double t2 = gettime();
	printf("Done ... Elapsed time = %.2lf seconds\n", t2 - t1);
	
	GdiplusShutdown(gdiplusToken);
	return 0;
}



