/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <math.h>
#include <string.h>

#include "general_utils.h"
#include "blocklanguage.h"
#include "rollpitchyaw.h"
#include "fdemsystem.h"

////////////////////////////////////////////////////////////////////////////////
cFDEmSystem::cFDEmSystem()
{
	
};
////////////////////////////////////////////////////////////////////////////////
cFDEmSystem::~cFDEmSystem()
{  
	
};
///////////////////////////////////////////////////////////////////////////
void cFDEmSystem::readsystemdescriptorfile(std::string systemdescriptorfile)
{	
	STM.loadfromfile(systemdescriptorfile);			
	SystemName = STM.getstringvalue("System.Name");
	std::vector<double> frequencies = STM.getdoublevector("System.Frequencies");
	std::vector<double> separations = STM.getdoublevector("System.Separations");
	std::vector<std::string> orientations = STM.getstringvector("System.Orientations");			

	//NominalHeight = control.getdoublevalue("System.NominalHeight");
	std::vector<COILSETTYPE> orientation;

	//Check system details	
	size_t nc = frequencies.size();
	if(separations.size() != nc || orientations.size() != nc){
		glog.logmsg("Error in system descriptor: Number of frequencies does not match number of seperations and/or orientations\n");
		exit(1);        
	}

	for(size_t i=0; i<nc; i++){						
		if(strcasecmp("HCP",orientations[i]) == 0)orientation.push_back(eHCP);		
		else if(strcasecmp("VCX",orientations[i]) == 0)orientation.push_back(eVCX);		
		else if(strcasecmp("VCP",orientations[i]) == 0)orientation.push_back(eVCP);
		else if(strcasecmp("PER",orientations[i]) == 0)orientation.push_back(ePER);
		else{		   
			glog.logmsg("Error in system descriptor: %s is an unknown orientation\n",orientations[i].c_str());		   
			exit(1);           
		}		
	}	
	
	//Add the coil sets		
	CoilSets.resize(nc);
	for(size_t i=0; i<nc; i++){		
		initialisecoilset(i,frequencies[i],orientation[i],separations[i]);		
	}	
};
////////////////////////////////////////////////////////////////////////////////
void cFDEmSystem::initialisecoilset(size_t cindex, double frequency, COILSETTYPE orientation, double separation)
{       
	sFDEmCoilSet& C = CoilSets[cindex];

	C.Frequency  = frequency;  
	C.L.setfrequency(frequency);

	C.Separation = separation;
	C.tx.npos = Vector3<double>(C.Separation/2.0,0,0);
	C.rx.npos = Vector3<double>(-C.Separation/2.0,0,0);

	C.Orientation = orientation;

	if(orientation==eHCP){
		C.tx.naxis = Vector3<double>(0,0,1);
		C.rx.naxis = Vector3<double>(0,0,1);      
	}
	else if(orientation==eVCX){
		C.tx.naxis = Vector3<double>(1,0,0);
		C.rx.naxis = Vector3<double>(1,0,0);
	}
	else if(orientation==ePER){
		C.tx.naxis = Vector3<double>(1,0,0);
		C.rx.naxis = Vector3<double>(0,1,0);
	}
	else if(orientation==eVCP){
		C.tx.naxis = Vector3<double>(0,1,0);
		C.rx.naxis = Vector3<double>(0,1,0);
	}
	else{
		printf("cFDEmSystem::addcoilset() unknown orientation\n");
		exit(1);
	}	
}
////////////////////////////////////////////////////////////////////////////
void cFDEmSystem::setrollpitchyaw(const double roll, const double pitch, const double yaw)
{
	Matrix33<double> IR  = inverse_rollpitchyaw_matrix(roll,pitch,yaw);
	
	for(size_t i=0;i<NumberOfCoilSets();i++){
		sFDEmCoilSet& C = CoilSets[i];	

		C.tx.rpy   = Vector3<double>(roll,pitch,yaw);		    	
		C.tx.IR    = IR;		
		C.tx.axis  = IR*C.tx.naxis;

		C.rx.rpy   = Vector3<double>(roll,pitch,yaw);		    	
		C.rx.IR    = IR;		
		C.rx.axis  = IR*C.rx.naxis;

		C.tx.pos = C.tx.IR*C.tx.npos;
		C.rx.pos = C.rx.IR*C.rx.npos;

		C.vtxrx  = C.tx.pos - C.rx.pos;

		C.L.setx(-C.vtxrx.e1);
		C.L.sety(-C.vtxrx.e2);		
	}	    
}
////////////////////////////////////////////////////////////////
void cFDEmSystem::setheight(const double height)
{		
	for(size_t i=0;i<NumberOfCoilSets();i++){
		sFDEmCoilSet& C = CoilSets[i];							
		C.L.setz(height + C.rx.pos.e3); 				
		C.L.seth(height + C.tx.pos.e3);
	}	    
}
////////////////////////////////////////////////////////////////
void cFDEmSystem::setgeometry(const sFDEmGeometry& g)
{
	setrollpitchyaw(g.birdroll,g.birdpitch,g.birdyaw);
	setheight(g.birdheight);
}
////////////////////////////////////////////////////////////////////////////////
/*
void cFDEmSystem::setearth(const size_t nlayers, const double* conductivity, const double* thickness)
{	
	for(size_t i=0; i<NumberOfCoilSets(); i++){   
		sFDEmCoilSet& C = CoilSets[i];
		C.L.setconductivitythickness(nlayers,conductivity,thickness);
	}
}
*/
////////////////////////////////////////////////////////////////////////////
void cFDEmSystem::setearth(const std::vector<double>& conductivity, const std::vector<double>& thickness)
{	
	for(size_t i=0; i<NumberOfCoilSets(); i++){   
		sFDEmCoilSet& C = CoilSets[i];
		C.L.setconductivitythickness(conductivity,thickness);
	} 
}
////////////////////////////////////////////////////////////////
void cFDEmSystem::setearth(const cEarth1D& e)
{	
	for(size_t i=0; i<NumberOfCoilSets(); i++){   
		sFDEmCoilSet& C = CoilSets[i];
		C.L.setconductivitythickness(e.conductivity,e.thickness);
	} 
}
////////////////////////////////////////////////////////////////////////////////
void cFDEmSystem::setupcomputations()
{		
	for(size_t i=0; i<NumberOfCoilSets(); i++){   
		sFDEmCoilSet& C = CoilSets[i];	 				 
		C.L.setupcomputations();
	}
}
////////////////////////////////////////////////////////////////
std::vector<double> cFDEmSystem::p()
{		
	std::vector<double> v(NumberOfCoilSets());
	for(size_t i=0; i<NumberOfCoilSets(); i++){
		sFDEmCoilSet& C = CoilSets[i];		 	 		
		v[i] = C.L.p(C.tx.axis,C.rx.axis);
		if(C.Orientation == eVCX) v[i] *= -1.0;
	}
	return v;
}
////////////////////////////////////////////////////////////////
cvector cFDEmSystem::s()
{		
	cvector v(NumberOfCoilSets());		
	for(size_t i=0; i<NumberOfCoilSets(); i++){
		sFDEmCoilSet& C = CoilSets[i];		 	 		
		v[i] = C.L.s(C.tx.axis,C.rx.axis);
		if(C.Orientation == eVCX) v[i] *= -1.0;
	}
	return v;
}
////////////////////////////////////////////////////////////////
cvector cFDEmSystem::ppms()
{		
	cvector v(NumberOfCoilSets());		
	for(size_t i=0; i<NumberOfCoilSets(); i++){
		sFDEmCoilSet& C = CoilSets[i];		 	 		
		v[i] = C.L.ppm(C.tx.axis,C.rx.axis);
		if(C.Orientation == eVCX) v[i] *= -1.0;
	}
	return v;
}
////////////////////////////////////////////////////////////////
cvector cFDEmSystem::dppms(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer)
{			
	cvector v(NumberOfCoilSets());		
	for(size_t i=0; i<NumberOfCoilSets(); i++){
		sFDEmCoilSet& C = CoilSets[i];	
		if(calculationtype == eDB){
			v[i]  = 2.0 * C.L.dppm(eDH,derivativelayer,C.tx.axis,C.rx.axis);		
		}
		else v[i] = C.L.dppm(calculationtype,derivativelayer,C.tx.axis,C.rx.axis);

		if(C.Orientation == eVCX) v[i] *= -1.0;
	}
	return v;
}
////////////////////////////////////////////////////////////////////////////////
cvector cFDEmSystem::noiseestimates(const cvector& response, const sFDEmNoiseModel& noisemodel)
{	
	double anr,mnr,ani,mni;
	cvector v(NumberOfCoilSets()); 	
	for(size_t i=0; i<NumberOfCoilSets(); i++){		
		anr = noisemodel.floor[i].real();
		mnr = 0.01 * response[i].real() * noisemodel.percentage[i].real();			
		ani = noisemodel.floor[i].imag();
		mni = 0.01 * response[i].imag() * noisemodel.percentage[i].imag();
		if(noisemodel.type==MAXPERCENTFLOOR){
			v[i] = cdouble(std::max(anr,mnr),std::max(ani,mni));			
		}
		else if(noisemodel.type==COMBINEPERCENTFLOOR){
			v[i] = cdouble(sqrt(anr*anr+mnr*mnr),sqrt(ani*ani+mni*mni));
		}
		else{	
			glog.logmsg("cFDEmSystem::noiseestimates unknown noise type\n");
			exit(1);        
		}
	}  
	return v;	
}
///////////////////////////////////////////////////////////////////////////
std::vector<double> cFDEmSystem::cv2dv(const cvector& cv)
{
	std::vector<double> dv(NumberOfCoilSets()*2);
	for(size_t i=0; i<NumberOfCoilSets(); i++){
		dv[i*2]   = cv[i].real();
		dv[i*2+1] = cv[i].imag();
	}
	return dv;
}
cvector cFDEmSystem::dv2cv(std::vector<double>& dv)
{
	cvector cv(NumberOfCoilSets());	
	for(size_t i=0; i<NumberOfCoilSets(); i++){	
		cv[i] = cdouble(dv[i*2] , dv[i*2+1]);
	}
	return cv;
}





