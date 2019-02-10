/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include "general_utils.h"
#include "layeredearthmodeller.h"


//Formulation mainly from the book 
//Geo-Electromagnetism, Wait James, R. Academic Press 1982.
 
////////////////////////////////////////////////////////////////////////////////
LayeredEarthModeller::LayeredEarthModeller(){
	initialise();  
};
////////////////////////////////////////////////////////////////////////////////
LayeredEarthModeller::~LayeredEarthModeller()
{	

};
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::initialise()
{      
  NumLayers     = 0;  
  NumAbscissa   = 17;

  LowerFractionalWidth  = 4.44;
  UpperFractionalWidth  = 1.84;


  //Ground EM 31 and 38
  //NumLayers     = 0;  
  //NumAbscissa   = 1581;
  //LowerFractionalWidth  = 4.44;
  //UpperFractionalWidth  = 10.84;

  setnumabscissa(NumAbscissa);
};
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setupcomputations()
{      
	if(earthchanged  || geometrychanged){
		setfrequencyabscissalayers();
		earthchanged=false;
		geometrychanged=false;
	}
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setnumabscissa(size_t nabscissa)
{
	NumAbscissa = nabscissa;	
	Abscissa.resize(NumAbscissa);	
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setnumlayers(size_t nlayers)
{   			
	//must set number of abscissa first
	NumLayers = nlayers;
	Conductivity.resize(NumLayers);
	Thickness.resize(NumLayers-1);
	for(size_t ai=0; ai<NumAbscissa; ai++){
		Abscissa[ai].Layer.resize(NumLayers);
	}
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setconductivitythickness(const size_t nlayers, const double* conductivity, const double* thickness)
{   			
	setnumlayers(nlayers);    
	for(size_t i=0; i<NumLayers; i++)Conductivity[i] = conductivity[i];  	
	for(size_t i=0; i<NumLayers-1; i++)Thickness[i] = thickness[i];	
	setmeanconductivity();		
	setmeanlog10conductivity();	
	earthchanged = true;
};
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness)
{   
	size_t nlayers = conductivity.size();
	setnumlayers(nlayers);	       
	Conductivity = conductivity;  
	Thickness    = thickness;
	setmeanconductivity();	
	setmeanlog10conductivity();	
	earthchanged = true;
};
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::printearth()
{	
	for(size_t i=0; i<NumLayers-1; i++){
       printf("Layer %02zu:\t%10lf mS/m\t%10lf m\n",i+1,1000.0*Conductivity[i],Thickness[i]);	   	   
	}
	printf("Layer %02zu:\t%10lf mS/m\n",NumLayers,1000.0*Conductivity[NumLayers-1]);	   
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setmeanlog10conductivity()
{ 
	//Returns the thickness weighted mean (in linear space) conductivity (but calculated in 10g10 space)
	double sumc = 0.0;
	double sumt = 0.0;
    for(size_t i=0; i<NumLayers; i++){
        if(i<(NumLayers-1)){			
			sumc += log10(Conductivity[i])*Thickness[i];
			sumt += Thickness[i];
		}
		else{
            sumc += log10(Conductivity[i])*sumt; //make basement layer as thick as sum of all overlying
			sumt += sumt;
		}
	}	
	meanlog10conductivity = pow(10.0,sumc/sumt);
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setmeanconductivity()
{ 
	//Returns the thickness weighted mean (in linear space) conductivity (but calculated in linear space)
	double sumc = 0.0;
	double sumt = 0.0;
    for(size_t i=0; i<NumLayers; i++){
        if(i<(NumLayers-1)){			
			sumc += Conductivity[i]*Thickness[i];
			sumt += Thickness[i];
		}
		else{
            sumc += Conductivity[i]*sumt; //make basement layer as thick as sum of all overlying
			sumt += sumt;
		}
	}	
	meanconductivity = sumc/sumt;
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setx(double x)
{
	if(x==X)return;
	X=x; X2=X*X; X4=X2*X2; setr();	
	geometrychanged = true;  
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::sety(double y)
{
	if(y==Y)return;
	Y=y; Y2=Y*Y; Y4=Y2*Y2; setr();
	geometrychanged = true;  
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setz(double z)
{
	if(z==Z)return;
	Z=z; ZH=Z-H; ZH2=ZH*ZH; setR();
	geometrychanged = true;  
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::seth(double h)
{
	if(h==H)return;
	H=h; ZH=Z-H; ZH2=ZH*ZH; setR();
	geometrychanged = true;  
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setr()
{	
	r2=X2+Y2; r=sqrt(r2); r3=r2*r; r4=r2*r2; r5=r4*r; setR();
	geometrychanged = true;  
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setR()
{		
	R2=(X2+Y2+ZH2); R=sqrt(R2); R5=R2*R2*R; R7=R5*R2;
	geometrychanged = true;  
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setgeometry(double x, double y, double z, double h)
{  
	setx(x); sety(y); setz(z); seth(h);  
	geometrychanged = true;  
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setfrequency(const double frequency)
{        	
	Frequency = frequency;
	Omega = TWOPI*frequency;
	MuZeroOmega = MUZERO*Omega;		
	iMuZeroOmega = cdouble(0.0,MuZeroOmega);			
};
/////////////////////////////////////////////////////////
void LayeredEarthModeller::setfrequencyabscissalayers()
{  		
	setintegrationnodes();
	for(size_t ai=0; ai<NumAbscissa; ai++){							
		for(size_t li=0; li<NumLayers; li++){
			double gamma2 = Conductivity[li]*MuZeroOmega;
			cdouble u = sqrt(cdouble(Abscissa[ai].Lambda2,gamma2));
			Abscissa[ai].Layer[li].U = u;
			if(li<NumLayers-1)Abscissa[ai].Layer[li].Exp2UT = exp(-2.0*u*Thickness[li]);
		}
		setlayermatrices(ai); 
		setpmatrix(ai);
	}
};
////////////////////////////////////////////////////////////////////////////////
double LayeredEarthModeller::approximatehalfspace()
{ 	
	if(NumLayers==1)return Conductivity[0];
	
	//peak lambda
	double   peaklambda = sqrt(MuZeroOmega * meanlog10conductivity / 4.0);			
    cdouble  rz = rzero_recursive(peaklambda);	
	cdouble  v  = (1.0+rz);
		
	v = iMuZeroOmega*v*v;	
	cdouble a = (-4.0*peaklambda*peaklambda*rz/v);
	return a.real();
}
////////////////////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::rzero_recursive(double lambda)
{
  //Wait's recursive formulation

  double lambda2 = lambda*lambda;
  double muzeroomega = MuZeroOmega;
  cdouble imuzeroomega(0.0,muzeroomega);

  double gamma2 = muzeroomega*Conductivity[NumLayers-1];
  cdouble u = sqrt(cdouble(lambda2,gamma2));  
  cdouble y = u/imuzeroomega;
  cdouble Nn,tanhv,v1,v2,v4,top,bot;
    
  for(size_t i=NumLayers-2; i>=0; i--){
	gamma2 = muzeroomega*Conductivity[i];
	u = sqrt(cdouble(lambda2,gamma2));
    Nn = u/imuzeroomega;

    v1 = u*Thickness[i];		
	
	//Expand - unstable    	
	//tanh(v) = (1.0 - v4)/(1.0 + v4 + 2.0*v2);
	
	v2 = exp(-2.0*v1);	
	v4 = exp(-4.0*v1);    	
	top = 1.0-v4;
	bot = 1.0+v4+2.0*v2;
    tanhv = top/bot;
						
	y = Nn*(y+Nn*tanhv)/(Nn+y*tanhv);

	if(i==0)break;//can't decrement size_t variable below 0
  }

  cdouble N0(0.0,-lambda/(muzeroomega));  // minus because dividing by i.muzero.omega
  return (N0-y)/(N0+y);

};
////////////////////////////////////////////////////////////////////////////////
inline cdouble LayeredEarthModeller::rzero_propogationmatrix(size_t ai)
{
  //Oldenberg's propogation matrix formulation
  setlayermatrices(ai);
  setpmatrix(ai);
  return Abscissa[ai].P21onP11;
};

////////////////////////////////////////////////////////////////////////////////
inline void LayeredEarthModeller::multiplymatrices(const sPropogationMatrix& m1, const sPropogationMatrix& m2, sPropogationMatrix& m)
{  
	 m.e11 = m1.e11 * m2.e11 + m1.e12 * m2.e21;
	 m.e12 = m1.e11 * m2.e12 + m1.e12 * m2.e22;
	 m.e21 = m1.e21 * m2.e11 + m1.e22 * m2.e21;
	 m.e22 = m1.e21 * m2.e12 + m1.e22 * m2.e22;	 	 
}
////////////////////////////////////////////////////////////////////////////////
inline void LayeredEarthModeller::addmatrices(const sPropogationMatrix& m1, const sPropogationMatrix& m2, sPropogationMatrix& m)
{  
	 m.e11 = m1.e11 + m2.e11;
	 m.e12 = m1.e12 + m2.e12;
	 m.e21 = m1.e21 + m2.e21;
	 m.e22 = m1.e22 + m2.e22;
}
////////////////////////////////////////////////////////////////////////////////
inline void LayeredEarthModeller::setlayermatrices(size_t ai)
{  
  cdouble e,eh,e1,e2;

  sAbscissaNode& A=Abscissa[ai];

  //M1  
  e  = A.Layer[0].U/A.Lambda;
  eh = e/2.0;	 
  e1 = 0.5 + eh;
  e2 = 0.5 - eh;
  
  sPropogationMatrix M;
  M.e11 = e1;
  M.e12 = e2;
  M.e21 = e2;
  M.e22 = e1;
  A.Layer[0].LayerMatrix = M;
    
  for(size_t li=1; li<NumLayers; li++){		 
	 //assumes all pearmabilities are muzero
	 e  = A.Layer[li].U / A.Layer[li-1].U; 	 
	 eh = e/2.0;	 
	 e1 = 0.5 + eh;
	 e2 = 0.5 - eh;
	 M.e11 = e1;
	 M.e12 = e2;
	 M.e21 = e2 * A.Layer[li-1].Exp2UT;
	 M.e22 = e1 * A.Layer[li-1].Exp2UT;
		 	 	 
	 A.Layer[li].LayerMatrix = M;	 
  }
  
  //Set Prematrices - prematrix for layer 0 does not apply
  if(NumLayers>1)A.Layer[1].LayerPreMatrix = A.Layer[0].LayerMatrix;	
  for(size_t li=2; li<NumLayers; li++){			 		  	
	multiplymatrices(A.Layer[li-1].LayerPreMatrix , A.Layer[li-1].LayerMatrix,A.Layer[li].LayerPreMatrix);	
  }

  //Set Postmatrices - postmatrix for layer N-1 does not apply  
  if(NumLayers>1){
	  A.Layer[NumLayers-2].LayerPostMatrix = A.Layer[NumLayers-1].LayerMatrix;	
  }
  for(size_t li=NumLayers-3; li>=0; li--){			 		 
 	multiplymatrices(A.Layer[li+1].LayerMatrix , A.Layer[li+1].LayerPostMatrix,A.Layer[li].LayerPostMatrix);
	if(li==0)break;//cannot decrement size_t variable below 0
  }
  
};
////////////////////////////////////////////////////////////////////////////////
inline void LayeredEarthModeller::setpmatrix(size_t ai)
{  	 	
	sAbscissaNode& A=Abscissa[ai];
	//Set Full matrix
	if(NumLayers==1){
		A.P_Full = A.Layer[0].LayerMatrix;
	}
	else multiplymatrices(Abscissa[ai].Layer[NumLayers-1].LayerPreMatrix , Abscissa[ai].Layer[NumLayers-1].LayerMatrix, A.P_Full);
	A.P21onP11 = A.P_Full.e21/A.P_Full.e11;				
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline sPropogationMatrix LayeredEarthModeller::dMjdCj(size_t ai, size_t li)
{  	 	
	sPropogationMatrix m;
		
    if(li==0){				
		cdouble a=iMuZeroOmega/(4.0*Abscissa[ai].Lambda*Abscissa[ai].Layer[li].U);
		m.e11 = a;
		m.e12 = -a;
		m.e21 = -a;
		m.e22 = a;
		return m;
    }
    else{
		cdouble a=iMuZeroOmega/(4.0*Abscissa[ai].Layer[li-1].U * Abscissa[ai].Layer[li].U);
		cdouble ae=a*Abscissa[ai].Layer[li-1].Exp2UT;        		
		m.e11 = a;
		m.e12 = -a;
		m.e21 = -ae;
		m.e22 = ae;
		return m;
    }   
}
////////////////////////////////////////////////////////////////////////////////
inline sPropogationMatrix LayeredEarthModeller::dMjplus1dCj(size_t ai, size_t li)
{  	 		
	cdouble duds= iMuZeroOmega/(2.0*Abscissa[ai].Layer[li].U);
	cdouble y   = Abscissa[ai].Layer[li+1].U/Abscissa[ai].Layer[li].U;
	
	cdouble dydu = -y/Abscissa[ai].Layer[li].U;
	cdouble dyds = dydu*duds;

	cdouble v = Abscissa[ai].Layer[li].Exp2UT;
	cdouble dvdu=-2.0*Thickness[li] * Abscissa[ai].Layer[li].Exp2UT;

	cdouble dvds=dvdu*duds;
	
	cdouble ydvds=y*dvds;
	cdouble vdyds=v*dyds;


    sPropogationMatrix M;
	M.e11 = 0.5*dyds;
	M.e12 = -M.e11;
	M.e21 = 0.5*(dvds - (ydvds + vdyds));
	M.e22 = 0.5*(dvds + (ydvds + vdyds));
	
	return M;		
	
}
////////////////////////////////////////////////////////////////////////////////
inline sPropogationMatrix LayeredEarthModeller::dPdCj(size_t ai, size_t li)
{  	
	sAbscissaNode& A=Abscissa[ai];

	sPropogationMatrix tmp;
	//One layer case
	if(NumLayers==1)return tmp = dMjdCj(ai,li);

	//Not last layer
	if(li<NumLayers-1){	
		
		tmp = dMjdCj(ai,li);
		sPropogationMatrix a;
		multiplymatrices(tmp , A.Layer[li+1].LayerMatrix, a);
		
        tmp = dMjplus1dCj(ai,li);
		sPropogationMatrix b;
		multiplymatrices(A.Layer[li].LayerMatrix , tmp, b);

		sPropogationMatrix c;
		addmatrices(a,b,c);
		//First layer case
		if(li==0){												
			if(NumLayers==2) return c;			
			sPropogationMatrix result;
			multiplymatrices(c,A.Layer[li+1].LayerPostMatrix, result);
			return result;			
		}				
		else if(li==NumLayers-2){															
			sPropogationMatrix result;
			multiplymatrices(A.Layer[li].LayerPreMatrix , c, result);
			return result;		
		}		
		else{				
			sPropogationMatrix d;
			multiplymatrices(A.Layer[li].LayerPreMatrix , c, d);						
			sPropogationMatrix result;
			multiplymatrices(d,A.Layer[li+1].LayerPostMatrix, result);
			return result;
		}		
	}
	//Last layer
	else{
		tmp = dMjdCj(ai,li);
		sPropogationMatrix result;
        multiplymatrices(A.Layer[li].LayerPreMatrix , tmp, result);		
		return result;
	}
	    		
}
////////////////////////////////////////////////////////////////////////////////
inline cdouble LayeredEarthModeller::dP21onP11dCj(size_t ai, size_t li)
{  	 
	sPropogationMatrix m = dPdCj(ai,li);
	return m.e21/Abscissa[ai].P_Full.e11 - m.e11*Abscissa[ai].P21onP11/Abscissa[ai].P_Full.e11;	    		
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline sPropogationMatrix LayeredEarthModeller::dMjplus1dTj(size_t ai, size_t li)
{  	 			
	cdouble y = Abscissa[ai].Layer[li+1].U / Abscissa[ai].Layer[li].U;			
	cdouble dvdt=-2.0*Abscissa[ai].Layer[li].U*Abscissa[ai].Layer[li].Exp2UT;
	
    sPropogationMatrix m;
	m.e11 = 0.0;
	m.e12 = 0.0;
	m.e21 = 0.5*(1.0-y)*dvdt;
	m.e22 = 0.5*(1.0+y)*dvdt;		
	return m;		
	
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline sPropogationMatrix LayeredEarthModeller::dPdTj(size_t ai, size_t li)
{  	    
	sAbscissaNode& A=Abscissa[ai];

    sPropogationMatrix tmp;
	tmp.e11 = tmp.e12 = tmp.e21 = tmp.e22 =0;

	//One layer case
	if(NumLayers==1) return tmp;

	tmp = dMjplus1dTj(ai,li);
	sPropogationMatrix b;
	multiplymatrices(A.Layer[li].LayerMatrix , tmp, b);

	//First layer case
	if(li==0){
		if(NumLayers==2)return b;
		sPropogationMatrix result;
		multiplymatrices(b,A.Layer[li+1].LayerPostMatrix, result);
		return result;
	}
	else if(li>0 && li<NumLayers-1){				
		sPropogationMatrix d;
		multiplymatrices(A.Layer[li].LayerPreMatrix,b, d);		
		if(li==NumLayers-2){												
			return d;				
		}
		sPropogationMatrix result;
		multiplymatrices(d,A.Layer[li+1].LayerPostMatrix, result);				
		return result;
	}
	//Last layer case
	else{
		message("Warning: LayeredEarthModeller::dPdtj Zero thickness derivative for halfspace layer\n");
		tmp.e11 = tmp.e12 = tmp.e21 = tmp.e22 =0;
		return tmp;
	}
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline cdouble LayeredEarthModeller::dP21onP11dTj(size_t ai, size_t li)
{  	 
	sPropogationMatrix m = dPdTj(ai,li);
	return m.e21/Abscissa[ai].P_Full.e11 - m.e11*Abscissa[ai].P21onP11/Abscissa[ai].P_Full.e11;	    		
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::setintegrationnodes()
{			
	double peak_exp2=2.0/(Z+H);	
	double peak_exp3=3.0/(Z+H);		

	ApproximateHalfspace = approximatehalfspace();
	PeakLambda = sqrt(MuZeroOmega * ApproximateHalfspace/4.0);

	double lp = log(std::min(PeakLambda,peak_exp2));
	double up = log(std::max(PeakLambda,peak_exp3));	
						
	LowerBound = lp - LowerFractionalWidth;
	UpperBound = up + UpperFractionalWidth;		
		
	AbscissaSpacing = (UpperBound - LowerBound)/(double)(NumAbscissa-1);
		
	double lambda;
	double loglambda = LowerBound;
	for(size_t ai=0; ai<NumAbscissa; ai++){
		sAbscissaNode& A=Abscissa[ai];
		lambda = exp(loglambda);
		A.Lambda = lambda;		
		A.Lambda2 = A.Lambda * lambda;
		A.Lambda3 = A.Lambda2 * lambda;
		A.Lambda4 = A.Lambda3 * lambda;
		A.Lambda_r = A.Lambda * r;						
				
		#if defined _WIN32
		A.j0Lambda_r = _j0(A.Lambda_r);		
		A.j1Lambda_r = _j1(A.Lambda_r);				
		#else
		A.j0Lambda_r = j0(A.Lambda_r);		
		A.j1Lambda_r = j1(A.Lambda_r);				
		#endif
		
		loglambda += AbscissaSpacing;
	}	
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dointegrals()
{		
	trapezoid();//the results go into the variable trapezoid_result			
	
	switch(CalculationType)
	{
	case eFM:
		I0.FM=trapezoid_result[0];
		I1.FM=trapezoid_result[1];
		I2.FM=trapezoid_result[2];
		break;
	case eDC:	
		I0.dC=trapezoid_result[0];
		I1.dC=trapezoid_result[1];
		I2.dC=trapezoid_result[2];
		break;
	case eDT:		
		I0.dT=trapezoid_result[0];
		I1.dT=trapezoid_result[1];
		I2.dT=trapezoid_result[2];
		break;
	case eDZ:	
		I0.dZ=trapezoid_result[0];
		I1.dZ=trapezoid_result[1];
		I2.dZ=trapezoid_result[2];
		break;
	case eDH:
		I0.dH=trapezoid_result[0];
		I1.dH=trapezoid_result[1];
		I2.dH=trapezoid_result[2];	
		break;
	case eDX:
		I0.dX=trapezoid_result[0];
		I1.dX=trapezoid_result[1];
		I2.dX=trapezoid_result[2];
		break;
	case eDY:
		I0.dY=trapezoid_result[0];
		I1.dY=trapezoid_result[1];
		I2.dY=trapezoid_result[2];
		break;	
	default:			
		message("Error: LayeredEarthModeller::dointegrals() unknown calculation type %c\n",CalculationType);		
		exit(1);
	}
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void LayeredEarthModeller::trapezoid()
{			
	trapezoid_result[0]=cdouble(0.0,0.0);
	trapezoid_result[1]=cdouble(0.0,0.0);
	trapezoid_result[2]=cdouble(0.0,0.0);

	//First and abscissa
	compute_integrand(0);		
	trapezoid_result[0] += integrand_result[0];
	trapezoid_result[1] += integrand_result[1];
	trapezoid_result[2] += integrand_result[2];
	
	compute_integrand(NumAbscissa-1);		
	trapezoid_result[0] += integrand_result[0];
	trapezoid_result[1] += integrand_result[1];
	trapezoid_result[2] += integrand_result[2];
		
	trapezoid_result[0] *= 0.5;
	trapezoid_result[1] *= 0.5;
	trapezoid_result[2] *= 0.5;
		

	//Cenral Abscissas
	for(size_t ai=1; ai<NumAbscissa-1; ai++){
		compute_integrand(ai);
		trapezoid_result[0] += integrand_result[0];
		trapezoid_result[1] += integrand_result[1];
		trapezoid_result[2] += integrand_result[2];		
	}
		
	trapezoid_result[0] *= AbscissaSpacing;
	trapezoid_result[1] *= AbscissaSpacing;
	trapezoid_result[2] *= AbscissaSpacing;	
}
////////////////////////////////////////////////////////////////////////////////
inline void LayeredEarthModeller::compute_integrand(size_t ai)
{		
	sAbscissaNode& A = Abscissa[ai];	
		
	double lambdar = A.Lambda_r;
	double j0=A.j0Lambda_r;
	double j1=A.j1Lambda_r;
	double e=exp(-(Z+H)*A.Lambda);
	double l2e=A.Lambda2*e;
	double l3e=A.Lambda3*e;
	double l4e=A.Lambda4*e;	
	

	double k0,k1,k2;
	
	cdouble earthkernel;

	switch(CalculationType)
	{
	case eFM:
		earthkernel = -A.P21onP11;		
		k0 = l3e*j0;
	    k1 = l3e*j1;
	    k2 = l2e*j1;
		break;
	case eDX:
		earthkernel = -A.P21onP11;
		k0 = -l4e*j1*X/r;
	    k1 = l4e*(j0-j1/lambdar)*X/r;
	    k2 = l3e*(j0-j1/lambdar)*X/r;
		break;
	case eDY:
		earthkernel = -A.P21onP11;
		k0 = -l4e*j1*Y/r;
	    k1 = l4e*(j0-j1/lambdar)*Y/r;
	    k2 = l3e*(j0-j1/lambdar)*Y/r;
		break;
	case eDZ:
		earthkernel = -A.P21onP11;
		k0 = -l4e*j0;
	    k1 = -l4e*j1;
	    k2 = -l3e*j1;
		break;
	case eDH:
		earthkernel = -A.P21onP11;;
		k0 = -l4e*j0;
	    k1 = -l4e*j1;
	    k2 = -l3e*j1;
		break;
	case eDC:
		earthkernel = -dP21onP11dCj(ai,DerivativeLayer);		
		k0 = l3e*j0;
	    k1 = l3e*j1;
	    k2 = l2e*j1;
		break;
	case eDT:
		earthkernel = -dP21onP11dTj(ai,DerivativeLayer);
		k0 = l3e*j0;
	    k1 = l3e*j1;
	    k2 = l2e*j1;
		break;
	default:
		message("Error: LayeredEarthModeller::compute_integrand() unknown calculation type %c\n",CalculationType);		
		exit(1);			
	}	
					
	integrand_result[0]=earthkernel*k0;	
	integrand_result[1]=earthkernel*k1;	
	integrand_result[2]=earthkernel*k2;	
}
////////////////////////////////////////////////////////////////////////////////
void LayeredEarthModeller::PTFM(Matrix33<double>& T)
{	
	T.e11 = (3.0*X2 - R2)/R5; 
	T.e12 = 3.0*X*Y/R5;
	T.e13 = 3.0*X*ZH/R5;

	//Note error in Fitterman and Yin paper should no be minus sign at element 2,1
	T.e21 = T.e12;
	T.e22 = (3.0*Y2 - R2)/R5;
	T.e23 = 3.0*Y*ZH/R5;

	T.e31 = T.e13;
	T.e32 = T.e23;
	T.e33 = (3.0*ZH2 - R2)/R5;	
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dPTdX(Matrix33<double>& T)
{
	T.e11 = -3.0*X*(2.0*X2-3.0*Y2-3.0*ZH2)/R7;
	T.e12 = -3.0*Y*(4.0*X2-Y2-ZH2)/R7;
	T.e13 = -3.0*ZH*(4.0*X2-Y2-ZH2)/R7;
	T.e21 = T.e12;
	T.e22 = 3.0*X*(X2-4.0*Y2+ZH2)/R7;
	T.e23 = -15.0*Y*ZH/R7*X;
	T.e31 = T.e13;
	T.e32 = T.e23;
	T.e33 = 3.0*X*(X2+Y2-4.0*ZH2)/R7;
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dPTdY(Matrix33<double>& T)
{
	T.e11 = -3.0*Y*(4.0*X2-Y2-ZH2)/R7;
	T.e12 = 3.0*X*(X2-4.0*Y2+ZH2)/R7;
	T.e13 = -15.0*X*ZH/R7*Y;
	T.e21 = T.e12;
	T.e22 = 3.0*Y*(3.0*X2-2.0*Y2+3.0*ZH2)/R7;
	T.e23 = 3.0*ZH*(X2-4.0*Y2+ZH2)/R7;
	T.e31 = T.e13;
	T.e32 = T.e23;
	T.e33 = 3.0*Y*(X2+Y2-4.0*ZH2)/R7;
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dPTdZ(Matrix33<double>& T)
{	
	T.e11 = -3.0*ZH*(4.0*X2-Y2-ZH2)/R7;
	T.e12 = -15.0*X*Y/R7*ZH;
	T.e13 = 3.0*X*(X2+Y2-4.0*ZH2)/R7;
	T.e21 = T.e12;
	T.e22 = 3.0*ZH*(X2-4.0*Y2+ZH2)/R7;
	T.e23 = 3.0*Y*(X2+Y2-4.0*ZH2)/R7;
	T.e31 = T.e13;
	T.e32 = T.e23;
	T.e33 = 3.0*ZH*(3.0*X2+3.0*Y2-2.0*ZH2)/R7;	
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dPTdH(Matrix33<double>& T)
{
	//of course this is just minus d/dZ		
	T.e11 = 3.0*ZH*(4.0*X2-Y2-ZH2)/R7;
	T.e12 = 15.0*X*Y/R7*ZH;
	T.e13 = -3.0*X*(X2+Y2-4.0*ZH2)/R7;
	T.e21 = T.e12;
	T.e22 = -3.0*ZH*(X2-4.0*Y2+ZH2)/R7;
	T.e23 = -3.0*Y*(X2+Y2-4.0*ZH2)/R7;
	T.e31 = T.e13;
	T.e32 = T.e23;
	T.e33 = -3.0*ZH*(3.0*X2+3.0*Y2-2.0*ZH2)/R7;
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::STFM(Matrix33<cdouble>& T)
{				
	T.e11 = ((X2/r2 - Y2/r2)*I2.FM/r - I0.FM*X2/r2);
	T.e12 = (X*Y/r2)*(2.0*I2.FM/r - I0.FM);
	T.e13 = (-X/r)*I1.FM;	

	//Note error in Fitterman and Yin paper should not be minus sign at element 2,1
	T.e21 = T.e12;	
	T.e22 = ((Y2/r2 - X2/r2)*I2.FM/r - I0.FM*Y2/r2);
	T.e23 = (-Y/r)*I1.FM;	

	T.e31 = -T.e13;
	T.e32 = -T.e23;
	T.e33 = -I0.FM;	
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dSTdC(Matrix33<cdouble>& T)
{		
	T.e11 = ((X2/r2 - Y2/r2)*I2.dC/r - I0.dC*X2/r2);
	T.e12 = (X*Y/r2)*(2.0*I2.dC/r - I0.dC);
	T.e13 = (-X/r)*I1.dC;

	T.e21 = T.e12;	
	T.e22 = ((Y2/r2 - X2/r2)*I2.dC/r - I0.dC*Y2/r2);
	T.e23 = (-Y/r)*I1.dC;

	T.e31 = -T.e13;
	T.e32 = -T.e23;
	T.e33 = -I0.dC;	
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dSTdT(Matrix33<cdouble>& T)
{		
	T.e11 = ((X2/r2 - Y2/r2)*I2.dT/r - I0.dT*X2/r2);
	T.e12 = (X*Y/r2)*(2.0*I2.dT/r - I0.dT);
	T.e13 = (-X/r)*I1.dT;

	T.e21 = T.e12;	
	T.e22 = ((Y2/r2 - X2/r2)*I2.dT/r - I0.dT*Y2/r2);
	T.e23 = (-Y/r)*I1.dT;

	T.e31 = -T.e13;
	T.e32 = -T.e23;
	T.e33 = -I0.dT;
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dSTdX(Matrix33<cdouble>& T)
{	
	T.e11 = (X4*I2.dX-I0.dX*X4*r-I2.FM*X2*X-I0.dX*X2*r*Y2-2.0*I0.FM*X*r*Y2+5.0*X*Y2*I2.FM-Y4*I2.dX)/r5;
	T.e12 = 2.0*Y/r3*I2.FM-Y/r2*I0.FM-6.0*X2*Y/r5*I2.FM+2.0*X2*Y/r4*I0.FM+2.0*X*Y/r3*I2.dX-X*Y/r2*I0.dX;
	T.e13 = -(I1.FM*Y2+X2*X*I1.dX+X*I1.dX*Y2)/r3;

	T.e21 = T.e12;
	T.e22 = -(I2.dX*X4-I2.FM*X2*X+I0.dX*Y2*r*X2-2.0*I0.FM*Y2*X*r+5.0*X*Y2*I2.FM+I0.dX*Y4*r-Y4*I2.dX)/r5;
	T.e23 = Y/r3*I1.FM*X-Y/r*I1.dX;

	T.e31 = -T.e13;
	T.e32 = -T.e23;
	T.e33 = -I0.dX;	
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dSTdY(Matrix33<cdouble>& T)
{	
	T.e11 = (X4*I2.dY-I0.dY*X4*r+2.0*I0.FM*X2*Y*r-I0.dY*X2*r*Y2-5.0*X2*Y*I2.FM+Y2*Y*I2.FM-Y4*I2.dY)/r5;
	T.e12 = 2.0*X/r3*I2.FM-X/r2*I0.FM-6.0*X*Y2/r5*I2.FM+2.0*X*Y2/r4*I0.FM+2.0*X*Y/r3*I2.dY-X*Y/r2*I0.dY;
	T.e13 = X/r3*I1.FM*Y-X/r*I1.dY;

	T.e21 = T.e12;
	T.e22 = -(I2.dY*X4+I0.dY*Y2*r*X2+2.0*I0.FM*Y*r*X2-5.0*X2*Y*I2.FM+I0.dY*Y4*r+Y2*Y*I2.FM-Y4*I2.dY)/r5;
	T.e23 = -(I1.FM*X2+Y*I1.dY*X2+Y2*Y*I1.dY)/r3;

	T.e31 = -T.e13;
	T.e32 = -T.e23;
	T.e33 = -I0.dY;	
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dSTdZ(Matrix33<cdouble>& T)
{
	T.e11 = ((X2/r2-Y2/r2)*I2.dZ/r-I0.dZ*X2/r2);
	T.e12 = X*Y/r2*(2.0*I2.dZ/r-I0.dZ);
	T.e13 = -X/r*I1.dZ;

	T.e21 = T.e12;
	T.e22 = ((Y2/r2-X2/r2)*I2.dZ/r-I0.dZ*Y2/r2);
	T.e23 = -Y/r*I1.dZ;

	T.e31 = -T.e13;
	T.e32 = -T.e23;
	T.e33 = -I0.dZ;
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::dSTdH(Matrix33<cdouble>& T)
{
	T.e11 = ((X2/r2-Y2/r2)*I2.dH/r-I0.dH*X2/r2);
	T.e12 = X*Y/r2*(2.0*I2.dH/r-I0.dH);
	T.e13 = -X/r*I1.dH;

	T.e21 = T.e12;
	T.e22 = ((Y2/r2-X2/r2)*I2.dH/r-I0.dH*Y2/r2);
	T.e23 = -Y/r*I1.dH;

	T.e31 = -T.e13;
	T.e32 = -T.e23;
	T.e33 = -I0.dH;	
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::PT(const CALCULATIONTYPE& calculationtype, Matrix33<double>& T)
{		
	CalculationType = calculationtype;
	switch(calculationtype)
	{
	case eFM: PTFM(T);break;
	case eDX: dPTdX(T);break;
	case eDY: dPTdY(T);break;
	case eDZ: dPTdZ(T);break;
	case eDH: dPTdH(T);break;
	case eDC: break;
	case eDT: break;
	default:
		message("LayeredEarthModeller::pt() unknown calculation type %c\n",calculationtype);		
		exit(1);			
	}	
}
////////////////////////////////////////////////////////////////
void LayeredEarthModeller::ST(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer, Matrix33<cdouble>& T)
{	
	CalculationType = calculationtype;
	DerivativeLayer = derivativelayer;
	dointegrals();

	switch(calculationtype)
	{
	case eFM: STFM(T);  break;
	case eDX: dSTdX(T); break;
	case eDY: dSTdY(T); break;
	case eDZ: dSTdZ(T); break;
	case eDH: dSTdH(T); break;
	case eDC: dSTdC(T); break;
	case eDT: dSTdT(T); break;
	default:
		message("LayeredEarthModeller::sfield_if unknown calculation type %c\n",calculationtype);		
		exit(1);			
	}		
}
////////////////////////////////////////////////////////////////
Vector3<double> LayeredEarthModeller::pfield_if(const CALCULATIONTYPE& calculationtype, const Vector3<double>& txdir)
{	
	Matrix33<double> T;
	PT(calculationtype,T);	
	return T*txdir;
}
////////////////////////////////////////////////////////////////
Vector3<cdouble> LayeredEarthModeller::sfield_if(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer, const Vector3<double>& txdir)
{	
	Matrix33<cdouble> T;
	ST(calculationtype,derivativelayer,T);	
	return T*txdir;	
}
////////////////////////////////////////////////////////////////
double LayeredEarthModeller::p(const Vector3<double>& txdir, const Vector3<double>& rxdir)
{	
	Vector3<double> v = pfield_if(eFM,txdir);	
	return v.dot(rxdir);
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::s(const Vector3<double>& txdir, const Vector3<double>& rxdir)
{	
	Vector3<cdouble> v = sfield_if(eFM,0,txdir);	
	return v.dot(rxdir);			
}
////////////////////////////////////////////////////////////////
double  LayeredEarthModeller::dp(const CALCULATIONTYPE& calculationtype, const Vector3<double>& txdir,  const Vector3<double>& rxdir)
{
	Vector3<double> v = pfield_if(calculationtype,txdir);	
	return v.dot(rxdir);			
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::ds(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer, const Vector3<double>& txdir,  const Vector3<double>& rxdir)
{	
	Vector3<cdouble> v = sfield_if(calculationtype,derivativelayer,txdir);	
	return v.dot(rxdir);			
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::dsdx(const Vector3<double>& txdir, const Vector3<double>& rxdir)
{			 		
	Vector3<cdouble> sf = sfield_if(eDX,0,txdir);
	return sf.dot(rxdir);			
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::dsdy(const Vector3<double>& txdir, const Vector3<double>& rxdir)
{	
	Vector3<cdouble> sf = sfield_if(eDY,0,txdir);	
	return sf.dot(rxdir);		
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::dsdz(const Vector3<double>& txdir, const Vector3<double>& rxdir)
{	
	Vector3<cdouble> sf = sfield_if(eDZ,0,txdir);	
	return sf.dot(rxdir);		
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::dsdh(const Vector3<double>& txdir, const Vector3<double>& rxdir)
{	
	Vector3<cdouble> sf = sfield_if(eDH,0,txdir);	
	return sf.dot(rxdir);		
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::dsdc(const size_t dlayer,const Vector3<double>& txdir, const Vector3<double>& rxdir)
{	
	Vector3<cdouble> sf = sfield_if(eDC,dlayer,txdir);	
	return sf.dot(rxdir);		
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::dsdt(const size_t dlayer,const Vector3<double>& txdir, const Vector3<double>& rxdir)
{	
	Vector3<cdouble> sf = sfield_if(eDT,dlayer,txdir);	
	return sf.dot(rxdir);		
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::ppm(const Vector3<double>& txdir, const Vector3<double>& rxdir)
{
	double  pf = p(txdir,rxdir);
	cdouble sf = s(txdir,rxdir);	
	return 1.0e6*(sf/pf);	
}
////////////////////////////////////////////////////////////////
cdouble LayeredEarthModeller::dppm(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer, const Vector3<double>& txdir, const Vector3<double>& rxdir)
{
	double  pf = p(txdir,rxdir);
	cdouble dsf = ds(calculationtype,derivativelayer,txdir,rxdir);	
	return 1.0e6*(dsf/pf);	
}
////////////////////////////////////////////////////////////////
