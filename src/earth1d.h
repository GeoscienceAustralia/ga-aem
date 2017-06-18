/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _earth1d_H
#define _earth1d_H

#include <vector>
#include <cassert>

class cEarth1D{

public:
	std::vector<double> thickness;
	std::vector<double> conductivity;
	std::vector<double> chargeability;
	std::vector<double> timeconstant;
	std::vector<double> frequencydependence;
	
	cEarth1D(){};

	cEarth1D(const size_t nlayers){
		conductivity.resize(nlayers);
		chargeability.resize(nlayers);
		timeconstant.resize(nlayers);
		frequencydependence.resize(nlayers);
		thickness.resize(nlayers - 1);
	}

	cEarth1D(const std::vector<double>& _conductivity, const std::vector<double>& _thickness){
		conductivity = _conductivity;
		thickness    = _thickness;
	}

	cEarth1D(
		const std::vector<double>& _conductivity,
		const std::vector<double>& _thickness,
		const std::vector<double>& _chargeability,
		const std::vector<double>& _timeconstant,
		const std::vector<double>& _frequencydependence)
	{
		conductivity        = _conductivity;
		thickness           = _thickness;
		chargeability       = _chargeability;
		timeconstant        = _timeconstant;
		frequencydependence = _frequencydependence;		
	}

	cEarth1D(const size_t nlayers, const double* _conductivity, const double* _thickness)
	{
		conductivity = std::vector<double>(_conductivity, _conductivity + nlayers);
		thickness = std::vector<double>(_thickness, _thickness + nlayers - 1);
	}

	cEarth1D(const size_t nlayers, const double* _conductivity, const double* _thickness, const double* _chargeability, const double* _timeconstant, const double* _frequencydependence)
	{
		conductivity  = std::vector<double>(_conductivity, _conductivity+nlayers);
		thickness     = std::vector<double>(_thickness, _thickness + nlayers-1);
		chargeability = std::vector<double>(_chargeability, _chargeability + nlayers);
		timeconstant  = std::vector<double>(_timeconstant, _timeconstant + nlayers);
		frequencydependence = std::vector<double>(_frequencydependence, _frequencydependence + nlayers);		
	}
	
	size_t nlayers(){ return conductivity.size(); }

	void print(){
		for (size_t i = 0; i<nlayers() - 1; i++){
			printf("%d\t%8.6lf\t%6.2lf\n", (int)i, conductivity[i], thickness[i]);
		}
		printf("%d\t%8.6lf\n\n", (int)nlayers(), conductivity[nlayers() - 1]);
	}
};

#endif
