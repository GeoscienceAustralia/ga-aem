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

constexpr auto DN_LAYER = "layer";
const std::string DN_NONE;
const std::string UNITLESS;


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

	std::vector<double> layer_top_depth() const
	{
		const size_t n = thickness.size();
		std::vector<double> dtop(n+1);
		dtop[0] = 0.0;
		for (size_t i = 1; i <= n; i++) {
			dtop[i] = dtop[i-1] + thickness[i-1];
		}		
		return dtop;
	}

	std::vector<double> layer_bottom_depth() const 
	{
		const size_t n = thickness.size();
		std::vector<double> dbot(n+1);		
		dbot[0] = thickness[0];
		for (size_t i = 1; i < n; i++) {
			dbot[i] = dbot[i-1] + thickness[i];
		}
		dbot[n] = dbot[n-1] + thickness[n-1];
		return dbot;
	}
	
	size_t nlayers() const
	{ 
		return conductivity.size();
	}

	void print() const	
	{
		for (size_t i = 0; i<nlayers() - 1; i++){
			printf("%d\t%8.6lf\t%6.2lf\n", (int)i, conductivity[i], thickness[i]);
		}
		printf("%d\t%8.6lf\n\n", (int)nlayers(), conductivity[nlayers() - 1]);
	}

	void write(const std::string& filepath) const 
	{
		FILE* fp = fileopen(filepath, "w");
		size_t nl = conductivity.size();
		for (size_t i = 0; i < nl; i++) {
			if (i < thickness.size()) {
				fprintf(fp, "%e\t%e\n", conductivity[i], thickness[i]);
			}
			else {
				fprintf(fp, "%e\tInf\n", conductivity[i]);
			}
		}
		fclose(fp);
	}
};

#endif
