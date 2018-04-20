/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _large_loop_H
#define _large_loop_H


#include <vector>
#include "general_constants.h"
#include "geometry3d.h"

class cLargeLoop {

	private:				
		std::vector<cPnt> V; //loop vertices polygon coords of loop vertices - must be a closed loop

	public:

		cLargeLoop(const std::vector<double>& loopvertices_xyz){			
			size_t nv = loopvertices_xyz.size() / 3;
			V.resize(nv);
			for (size_t i = 0; i < nv; i++){
				V[i].x = loopvertices_xyz[i * 3];
				V[i].y = loopvertices_xyz[i * 3 + 1];
				V[i].z = loopvertices_xyz[i * 3 + 2];
			}			
		};

		cLargeLoop(const std::vector<cPnt>& loopvertices){
			V = loopvertices;			
		};

		size_t nvertex() const {return V.size();}

		double area(){ 
			// Assumes the z cords are in the same z=constant plane
			double a = 0;
			for (size_t i = 0; i < V.size()-1; i++){
				a += V[i].x*V[i + 1].y - V[i + 1].x*V[i].y;
			}
			a = abs(a/2.0);
			return a;
		}


		cVec HField(const cVec& r) const 
		{
			cVec H(0.0,0.0,0.0);
			const size_t nv = nvertex();
			for (size_t i = 0; i < nv-1; i++){
				const cVec p1 = V[i];
				const cVec p2 = V[i+1];
				const cVec l  = p2 - p1;
				const double len = l.length();
				
				double bot = l.dot(l);
				double t   = l.dot(r - p1)/bot;
				cVec   p0  = p1 + t*l;
				double h   = (p0 - r).length();
				
				cVec v1 = r - p1;
				cVec v2 = r - p2;
				cVec un = cross(v1, v2).unit();				
				double cosa1 = dot(l, v1) / len / v1.length();
				double cosa2 = dot(l, v2) / len / v2.length();
				cVec dH = un*((cosa1 - cosa2) / h);
				H += dH;
			}
			H /= FOURPI;
			return H;
		}

		cVec BField(const cPnt& r) const 
		{
			return MUZERO*HField(r);
		}
		
};

//void test_large_loop() {
//	std::vector<double> txloop = { 6.622, 0.624, 0.000, 5.114, 4.164, 0.000, 4.259, 7.918, 0.000, 2.821, 8.232, 0.000, 0.954, 7.893, 0.000, -3.924, 4.379, 0.000, -7.082, 1.622, 0.000, -8.675, -0.166, 0.000, -7.101, -1.845, 0.000, -3.908, -4.572, 0.000, 1.015, -7.949, 0.000, 2.828, -8.248, 0.000, 4.294, -7.896, 0.000, 5.108, -4.244, 0.000, 6.625, -0.716, 0.000, 6.622, 0.624, 0.000 };
//	cLargeLoop L(txloop);
//	double area = L.area();
//	cPnt r(-120, -10, -40);
//	cVec B = L.BField(r);
//	printf("%lg %lg %lg\n", B.x, B.y, B.z);
//};


#endif
