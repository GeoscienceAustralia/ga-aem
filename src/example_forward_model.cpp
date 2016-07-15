/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

/* Example driver program for simple forward model*/

#include <cstring>

#include "general_utils.h"
#include "file_utils.h"
#include "vector_utils.h"
#include "blocklanguage.h"
#include "le.h"
#include "tdemsystem.h"


int main(int argc, char* argv[])
{	
	//Load the AEM system specification files for the Skytem moments
	//only do this once
	//LM = Low moment pulse and 
	cTDEmSystem LM("..\\..\\examples\\bhmar-skytem\\stmfiles\\Skytem-LM.stm");
	//HM = high moment pulse	
	cTDEmSystem HM("..\\..\\examples\\bhmar-skytem\\stmfiles\\Skytem-HM.stm");

	//Load the system geometry (same for both moments)
	//This changes every fiducial/station
	cTDEmGeometry G;	
	G.tx_height = 30;
	G.tx_roll = 0;       G.tx_pitch = 0; G.tx_yaw = 0;
	G.txrx_dx = -12.62;  G.txrx_dy  = 0;  G.txrx_dz = +2.16;
	G.rx_roll = 0;       G.rx_pitch = 0; G.rx_yaw = 0;

	//Create the earth structure
	//This changes every fiducial/station
	cEarth1D E(3);	
	E.conductivity[0] = 0.010;
	E.conductivity[1] = 0.100;
	E.conductivity[2] = 0.001;
	E.thickness[0] = 20;
	E.thickness[1] = 40;
	//bottom layer is infinite thickness and not set
	
	//Create a response object for each moment (they have different numbers of windwos)
	cTDEmResponse LMR;
	cTDEmResponse HMR;

	//Run the forward model for each moment
	LM.forwardmodel(G, E, LMR);		
	HM.forwardmodel(G, E, HMR);

	//Merge the secondary field vertical (Z) components data into one data vector
	std::vector<double> data = concaternate(LMR.SZ, HMR.SZ);
	
	for (size_t i = 0; i < data.size(); i++){
		printf("%d %g\n", i, data[i]);
	}
}

