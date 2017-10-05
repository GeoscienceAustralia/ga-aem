/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cmath>

#define VERSION "1.0"

#if defined _MPI_ENABLED
#include <mpi.h>
#endif
#include <cstring>
#include "general_utils.h"
#include "general_types.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "vector_utils.h"
#include "matrix_ops.h"
#include "fdemsystem.h"
#include "galeisbsfdem.h"



FDEmSampleInverter::FDEmSampleInverter(std::string controlfilename, size_t size, size_t rank)
{	
	mSize = size;
	mRank = rank;
	mControlFilename = controlfilename;
	message("Loading control file %s\n", controlfilename.c_str());
	mControl.loadfromfile(mControlFilename);
	initialise();
}

void FDEmSampleInverter::initialise()
{ 		
	Logfile = mControl.getstringvalue("LogFile");
	if(mSize != 1){
		std::string tmp=stringvalue(mRank,".%04d");
		Logfile = insert_after_filename(Logfile,tmp);		
	}	

	message("Opening log file %s\n",Logfile.c_str());
	fp_log = fileopen(Logfile,"w");
	message(fp_log,"Logfile opened on %s\n",timestamp().c_str());
	message(fp_log,"Control file %s\n",mControlFilename.c_str());
	message(fp_log,"Version %s Compiled at %s on %s\n",VERSION,__TIME__,__DATE__);
	message(fp_log,"Working directory %s\n",getcurrentdirectory().c_str());
	message(fp_log,"Processes=%d\tRank=%d\n",mSize,mRank);
	mControl.write(fp_log);
	fflush(fp_log);

	//Load control file			
	Dump = mControl.getboolvalue("Dump");
	if(Dump){
		DumpPath = mControl.getstringvalue("DumpPath");
		fixseparator(DumpPath);
		if(DumpPath[DumpPath.length()-1] != pathseparator()){
			DumpPath.append(pathseparatorstring());
		}
		makedirectory(DumpPath.c_str());
	}

	initialise_emsystems();
	getoptions();
	getinputs();
	parsecolumns();
	getoutputs();

	setupdata();
	setupparameters();
	resizematrices();		
}
void FDEmSampleInverter::initialise_emsystems()
{
	cBlock b   = mControl.findblock("EMSystem");	
	std::string SystemFile = b.getstringvalue("SystemFile");
	EMSystem.readsystemdescriptorfile(SystemFile);	
	nCoilsets = EMSystem.NumberOfCoilSets();
	nChannels = nCoilsets * 2;	
	EMSystem.STM.write(fp_log);
}
void FDEmSampleInverter::getoptions()
{ 	
	cBlock b = mControl.findblock("Earth");
	nlayers  = (size_t)b.getintvalue("NumberOfLayers");

	//Options
	b = mControl.findblock("Options");
	SolveConductivity = b.getboolvalue("SolveConductivity");
	SolveThickness    = b.getboolvalue("SolveThickness");
	SolveBirdHeight   = b.getboolvalue("SolveBirdHeight");	

	AlphaC = b.getdoublevalue("AlphaConductivity");
	AlphaT = b.getdoublevalue("AlphaThickness");
	AlphaG = b.getdoublevalue("AlphaGeometry");
	AlphaS = b.getdoublevalue("AlphaSmoothness");

	thickness_weight_referenceconductivity = true;
	thickness_weight_smoothness = true;

	MaxIterations      = (size_t)b.getintvalue("MaximumIterations");
	MinimumPhiD        = b.getdoublevalue("MinimumPhiD");
	MinimumImprovement = b.getdoublevalue("MinimumPercentageImprovement");
}
void FDEmSampleInverter::getinputs()
{	
	//Input
	cBlock b = mControl.findblock("Input");
	
	DataFileName = b.getstringvalue("DataFile");
	
	DataFileHeaderLines = b.getsizetvalue("Headerlines");
	if (isundefined(DataFileHeaderLines)){
		DataFileHeaderLines = 0;
	}

	DataFileSubsample = b.getsizetvalue("Subsample");
	if (isundefined(DataFileSubsample)){
		DataFileSubsample = 1;
	}

	rootmessage(fp_log, "Opening Input DataFile %s\n", DataFileName.c_str());
	DataFilePointer = fileopen(DataFileName, "r");
	DataFileRecord = 0;

	rootmessage(fp_log, "Opening Output DataFile %s\n", OutputDataFile.c_str());
	ofp = fileopen(OutputDataFile, "w");
	Outputrecord = 1;
		
}
void FDEmSampleInverter::getoutputs()	
{
	//Output
	cBlock b = mControl.findblock("Output");
	OutputDataFile   = b.getstringvalue("DataFile");
	OutputHeaderFile = b.getstringvalue("HeaderFile");	
	OutputConductivity = b.getboolvalue("Conductivity");
	OutputThickness = b.getboolvalue("Thickness");
	OutputPositiveLayerTopDepths = b.getboolvalue("PositiveLayerTopDepths");
	OutputNegativeLayerTopDepths = b.getboolvalue("NegativeLayerTopDepths");
	OutputLayerTopElevations = b.getboolvalue("LayerTopElevations");
	OutputObservedData = b.getboolvalue("ObservedData");
	OutputPredictedData = b.getboolvalue("PredictedData");
	OutputPredictedBirdHeight = b.getboolvalue("PredictedBirdHeight");

	if (mSize>1){
		std::string tmp = stringvalue(mRank, ".%04d");
		OutputDataFile   = insert_after_filename(OutputDataFile, tmp);
		OutputHeaderFile = insert_after_filename(OutputHeaderFile, tmp);
	}
	message("Opening OutputDataFile file %s\n", OutputDataFile.c_str());
	ofp = fileopen(OutputDataFile, "w");

	message("Opening OutputHeaderFile file %s\n", OutputHeaderFile.c_str());
	hfp = fileopen(OutputHeaderFile, "w");
	Outputrecord = 1;

	

}
void FDEmSampleInverter::parsecolumns()
{	
	cBlock b = mControl.findblock("Input.Columns");
	sn.set(b,"SurveyNumber");
	dn.set(b,"DateNumber");
	fn.set(b,"FlightNumber");
	ln.set(b,"LineNumber");
	fidn.set(b,"FidNumber");
	xord.set(b,"Easting");
	yord.set(b,"Northing");
	elevation.set(b,"GroundElevation");
	altimeter.set(b,"Altimeter");

	birdheight.set(b,"BirdHeight");
	birdroll.set(b,"BirdRoll");
	birdpitch.set(b,"BirdPitch");
	birdyaw.set(b,"BirdYaw");
	
	emchannels.set(b,"EMChannels");
	std_emchannels.set(b,"StdDevEMChannels");

	cBlock rm = b.findblock("ReferenceModel");	
	ref_conductivity.set(rm,"Conductivity");
	ref_thickness.set(rm,"Thickness");
	ref_birdheight.set(rm,"BirdHeight");

	cBlock sd = b.findblock("StdDevReferenceModel");
	std_conductivity.set(sd,"Conductivity");
	std_thickness.set(sd,"Thickness");
	std_birdheight.set(sd,"BirdHeight");
}

void FDEmSampleInverter::setupdata()
{  	
	size_t dindex=0;
	ndata=0;	
	
	ndata += nChannels;  
	emIndex = dindex;
	dindex += nChannels;
					
	vObs  = std::vector<double>(ndata);
	vErr  = std::vector<double>(ndata);	
	vPred = std::vector<double>(ndata);	
}
void FDEmSampleInverter::setupparameters()
{
	nparam = 0;
	ngeomparam = 0;

	size_t pindex = 0;	
	if(SolveConductivity){	
		cIndex=pindex;
		nparam += nlayers;
		pindex += nlayers;
	}
	
	if(SolveThickness){
		tIndex=pindex;
		nparam += nlayers-1;
		pindex += nlayers-1;
	}
	
	if(SolveBirdHeight){		
		gIndex=pindex;
		birdheightIndex=pindex;
		pindex++;
		nparam++;
		ngeomparam++;
	}
			
	vParam       = std::vector<double>(nparam);	
	vRefParam    = std::vector<double>(nparam);	
	vRefParamStd = std::vector<double>(nparam);			
}
void FDEmSampleInverter::resizematrices()
{			
	Wc.newsize(nparam,nparam);
	Wt.newsize(nparam,nparam);
	Wg.newsize(nparam,nparam);
	Wm.newsize(nparam,nparam);

	Wd.newsize(ndata,ndata);
	J.newsize(ndata,nparam);

	JtWd.newsize(nparam,ndata);
	JtWdJ.newsize(nparam,nparam);
	A.newsize(nparam,nparam);
	
	x.resize(nparam);
	b.resize(nparam);
}

bool FDEmSampleInverter::readnextrecord(){
	if (DataFileRecord == 0){
		//Skip header lines
		for (size_t i = 0; i < DataFileHeaderLines; i++){
			bool status = filegetline(DataFilePointer, DataFileRecordString);
			if (status == false)return status;
			DataFileRecord++;
		}
	}
	else{
		//Skip lines for subsampling
		for (size_t i = 0; i < DataFileSubsample - 1; i++){
			bool status = filegetline(DataFilePointer, DataFileRecordString);
			if (status == false)return status;
			DataFileRecord++;
		}
	}

	bool status = filegetline(DataFilePointer, DataFileRecordString);
	if (status == false)return status;
	DataFileRecord++;
	return true;
}

bool FDEmSampleInverter::parserecord()
{ 			
	DataFileFieldStrings = fieldparsestring(DataFileRecordString.c_str(), " ,\t\r\n");
	if (DataFileFieldStrings.size() <= 1)return false;
		
	Id.uniqueid     = DataFileRecord;
	Id.surveynumber = (size_t)intvalue(sn);		
	Id.daynumber    = (size_t)intvalue(dn);
	Id.flightnumber = (size_t)intvalue(fn);
	Id.linenumber   = (size_t)intvalue(ln);
	Id.fidnumber    = doublevalue(fidn);
	Location.x      = doublevalue(xord); 	
	Location.y      = doublevalue(yord); 	
	Location.groundelevation = doublevalue(elevation); 
	Location.z      = doublevalue(altimeter); 
	
	GI.birdheight = doublevalue(birdheight); 
	GI.birdroll   = doublevalue(birdroll); 
	GI.birdpitch  = doublevalue(birdpitch); 
	GI.birdyaw    = doublevalue(birdyaw); 
			
	std::vector<double> chan    = doublevector(emchannels,nChannels);
	std::vector<double> stdchan = doublevector(std_emchannels,nChannels);
	
	DI.inphase.resize(nCoilsets);
	DI.quadrature.resize(nCoilsets);
	DE.inphase.resize(nCoilsets);
	DE.quadrature.resize(nCoilsets);
	for(size_t i=0; i<nCoilsets; i++){
		DI.inphase[i]    = chan[i*2];
		DI.quadrature[i] = chan[i*2+1];
		DE.inphase[i]    = stdchan[i*2];
		DE.quadrature[i] = stdchan[i*2+1];
	}
			
	ER.conductivity = doublevector(ref_conductivity,nlayers);	
	if(SolveConductivity){
		ES.conductivity = doublevector(std_conductivity,nlayers);
	}

	ER.thickness = doublevector(ref_thickness,nlayers-1);	
	if(SolveThickness){
		ES.thickness    = doublevector(std_thickness,nlayers-1);	
	}
	
	GR = GI;
	if(SolveBirdHeight){
		GR.birdheight = doublevalue(ref_birdheight);	
		GS.birdheight = doublevalue(std_birdheight);	
	}

	return true;
}

void FDEmSampleInverter::writeresult()
{ 	
	Outputcolumn=1;		
	std::string buf;
	std::string hdr;

	//Id	
	addhdrstring(hdr,"uniqueid");
	buf += strprint("%7d ",Id.uniqueid);

	addhdrstring(hdr,"survey");
	buf += strprint("%4d ",Id.surveynumber);

	addhdrstring(hdr,"date");
	buf += strprint("%8d ",Id.daynumber);

	addhdrstring(hdr,"flight");
	buf += strprint("%4d ",Id.flightnumber);

	addhdrstring(hdr,"line");
	buf += strprint("%10d ",Id.linenumber);

	addhdrstring(hdr,"fid");
	buf += strprint("%10.2lf ",Id.fidnumber);

	//Location
	addhdrstring(hdr,"easting");
	buf += strprint("%8.1lf ",Location.x);

	addhdrstring(hdr,"northing");
	buf += strprint("%9.1lf ",Location.y);

	addhdrstring(hdr,"elevation");
	buf += strprint("%7.2lf ",Location.groundelevation);

	addhdrstring(hdr,"altimeter");
	buf += strprint("%7.2lf ",Location.z);
	
	//Observed Geometry
	addhdrstring(hdr,"birdheight");	
	buf += strprint("%7.2lf ",GI.birdheight);
	addhdrstring(hdr,"birdroll");	
	buf += strprint("%7.2lf ",GI.birdroll);
	addhdrstring(hdr,"birdpitch");	
	buf += strprint("%7.2lf ",GI.birdpitch);
	addhdrstring(hdr,"birdyaw");	
	buf += strprint("%7.2lf ",GI.birdyaw);
	
	//Predicted TX height
	if(OutputPredictedBirdHeight){
		addhdrstring(hdr,"predicted_birdheight");	
		buf += strprint("%7.2lf ",GM.birdheight);
	}

	//Earth
	addhdrstring(hdr,"nlayers");
	buf += strprint("%3d ",nlayers);

	addhdrstring(hdr,"conductivity",nlayers);
	for(size_t i=0;i<nlayers;i++){		
		buf += strprint("%14.6le ",EM.conductivity[i]);
	}

	addhdrstring(hdr,"thickness",nlayers-1);
	for(size_t i=0;i<nlayers-1;i++){		
		buf += strprint("%8.2lf ",EM.thickness[i]);
	}

	if(OutputPositiveLayerTopDepths){
		addhdrstring(hdr,"depth_top",nlayers-1);
		double tsum=0.0;
		for(size_t i=0;i<nlayers-1;i++){		
			buf += strprint("%8.2lf ",tsum);
			tsum += EM.thickness[i];			
		}
	}

	if(OutputNegativeLayerTopDepths){
		addhdrstring(hdr,"depth_top_negative",nlayers-1);
		double tsum=0.0;
		for(size_t i=0;i<nlayers-1;i++){		
			buf += strprint("%8.2lf ",-tsum);
			tsum += EM.thickness[i];			
		}
	}

	if(OutputLayerTopElevations){
		addhdrstring(hdr,"elevation_top",nlayers-1);
		double etop = Location.groundelevation;
		for(size_t i=0;i<nlayers-1;i++){		
			buf += strprint("%8.2lf ",etop);
			etop -= EM.thickness[i];			
		}
	}
	
	//Observed
	if(OutputObservedData){
		addhdrstring(hdr,"observed_data",nChannels);
		for(size_t i=0;i<nCoilsets;i++){			
			buf += strprint("%8.2lf ",DI.inphase[i]);
			buf += strprint("%8.2lf ",DI.quadrature[i]);
		}		
	}

	if(OutputObservedData){
		addhdrstring(hdr,"observed_noise",nChannels);
		for(size_t i=0;i<nCoilsets;i++){			
			buf += strprint("%8.2lf ",DE.inphase[i]);
			buf += strprint("%8.2lf ",DE.quadrature[i]);
		}		
	}

	//Predicted
	if(OutputPredictedData){
		addhdrstring(hdr,"predicted_data",nChannels);
		for(size_t i=0;i<nCoilsets;i++){			
			buf += strprint("%8.2lf ",DM.inphase[i]);
			buf += strprint("%8.2lf ",DM.quadrature[i]);
		}		
	}

	addhdrstring(hdr,"AlphaC");
	buf += strprint("%14.6le ",AlphaC);

	addhdrstring(hdr,"AlphaT");
	buf += strprint("%14.6le ",AlphaT);

	addhdrstring(hdr,"AlphaG");
	buf += strprint("%14.6le ",AlphaG);

	addhdrstring(hdr,"AlphaS");
	buf += strprint("%14.6le ",AlphaS);

	addhdrstring(hdr,"PhiD");
	buf += strprint("%14.6le ",LastPhiD);

	addhdrstring(hdr,"PhiM");
	buf += strprint("%14.6le ",LastPhiM);

	addhdrstring(hdr,"PhiC");
	buf += strprint("%14.6le ",LastPhiC);

	addhdrstring(hdr,"PhiT");
	buf += strprint("%14.6le ",LastPhiT);

	addhdrstring(hdr,"PhiG");
	buf += strprint("%14.6le ",LastPhiG);

	addhdrstring(hdr,"PhiS");
	buf += strprint("%14.6le ",LastPhiS);

	addhdrstring(hdr,"Lambda");
	buf += strprint("%14.6le ",LastLambda);

	addhdrstring(hdr,"Iterations");
	buf += strprint("%3d",LastIteration);

	//Carriage return		
	buf += strprint("\n");
	fprintf(ofp,buf.c_str());
	fflush(ofp);

	if(Outputrecord==1){
		fprintf(hfp,hdr.c_str());
		fclose(hfp);
	}
	Outputrecord++;

}
void FDEmSampleInverter::addhdrstring(std::string& hstr, const std::string s, size_t nband)
{
	addhdrstring(hstr,s.c_str(),nband);
}
void FDEmSampleInverter::addhdrstring(std::string& hstr, const char* s, size_t nband)
{			
	if(nband==1){
		hstr += strprint("%d\t%s\n",Outputcolumn,s);
		Outputcolumn++;
	}
	else{
		hstr += strprint("%d-%d\t%s\n",Outputcolumn,Outputcolumn+nband-1,s);
		Outputcolumn += nband;
	}			
}

void FDEmSampleInverter::fillsample()
{	
	filldata();
	fillparameters();
	fillWd();	
	fillWc();	
	fillWt();	 
	fillWg();	
	fillLtL();

	if (Dump){
		dumptofile(GI, "geometry_in.dat");
		dumptofile(ER, "earth_in.dat");
	}

	if(Dump){	
		/*
		dumptofile(GR,"geometry_start.dat");		
		dumptofile(ER,"earth_start.dat");

		dumptofile(GR,"geometry_ref.dat");			
		dumptofile(ER,"earth_ref.dat");

		dumptofile(ES,"earth_std.dat");

		FILE* fp = fileopen(DumpPath+"Id.dat","w");
		fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf",Id.uniqueid,Id.surveynumber,Id.daynumber,Id.flightnumber,Id.linenumber,Id.fidnumber,Location.x,Location.y,Location.groundelevation,Location.z);
		fclose(fp);
		*/
	}	
}  
void FDEmSampleInverter::filldata()
{		
	//Data  	
	for(size_t i=0; i<nCoilsets; i++){	  
		vObs[emIndex+i*2]=DI.inphase[i];
		vObs[emIndex+i*2+1]=DI.quadrature[i];
		vErr[emIndex+i*2]=DE.inphase[i];
		vErr[emIndex+i*2+1]=DE.quadrature[i];
	}
	
	if(Dump){		
		dumptofile(vObs,"observed.dat");		
		dumptofile(vErr,"observed_std.dat");	
	}

	if(Dump){		
		dumptofile(vObs,"observed.dat");		
		dumptofile(vErr,"observed_std.dat");
	}
}
void FDEmSampleInverter::fillparameters()
{
	if(SolveConductivity){
		for(size_t i=0; i<nlayers; i++){
			vRefParam[i+cIndex]=log10(ER.conductivity[i]);
			vRefParamStd[i+cIndex]=ES.conductivity[i];
		}
	}
	
	if(SolveThickness){
		for(size_t i=0; i<nlayers-1; i++){
			vRefParam[i+tIndex]=log10(ER.thickness[i]);
			vRefParamStd[i+tIndex]=ES.thickness[i];
		}
	}
	
	if(SolveBirdHeight){
		vRefParam[birdheightIndex]=GI.birdheight;
		vRefParamStd[birdheightIndex]=GS.birdheight;
	}	
}
void FDEmSampleInverter::fillWd()
{
	Wd = 0.0;	
	for(size_t i=0; i<ndata; i++){			
		Wd[i][i] = 1.0/(ndata*vErr[i]*vErr[i]);
	}

	//if(Dump) writetofile(Wd,DumpPath+"Wd.dat");
}
void FDEmSampleInverter::fillWc()
{
	Wc=0;	
	if(SolveConductivity==false)return;

	std::vector<double> t(nlayers);
	if(nlayers==1){
		t[0]=1;
	}
	else if(nlayers==2){
		t[0]=ER.thickness[0];
		t[1]=ER.thickness[0];
	}
	else{
		for(size_t i=0;i<(nlayers-1);i++){	
			t[i]=ER.thickness[i];
		}
		t[nlayers-1]=(t[nlayers-2]/t[nlayers-3])*t[nlayers-2];
	}


	double tsum=0;
	for(size_t i=0; i<nlayers; i++)tsum += t[i];

	for(size_t i=0; i<nlayers; i++){						
		size_t p = i + cIndex;		
		if(nlayers>1 && thickness_weight_referenceconductivity){
			Wc[p][p]  = (double)nlayers * (t[i]/tsum) / (vRefParamStd[p]*vRefParamStd[p]);			
		}		
		else{
			Wc[p][p]  = 1.0/((double)nlayers*vRefParamStd[p]*vRefParamStd[p]);
		}
	}	
	//if(Dump)writetofile(Wc,DumpPath+"Wc.dat");		
}
void FDEmSampleInverter::fillWt()
{
	Wt=0;	
	if(SolveThickness==false)return;

	for(size_t i=0; i<nlayers-1; i++){						
		size_t p = i + tIndex;
		Wt[p][p] = 1.0/((double)(nlayers-1)*vRefParamStd[p]*vRefParamStd[p]);
	}
	//if(Dump)writetofile(Wt,DumpPath+"Wt.dat");
}
void FDEmSampleInverter::fillWg()
{
	Wg=0;
	if(ngeomparam<=0)return;

	for(size_t i=0; i<ngeomparam; i++){						
		size_t p = i + gIndex;		
		Wg[p][p] = 1.0/((double)ngeomparam*vRefParamStd[p]*vRefParamStd[p]);
	}	
	//if(Dump)writetofile(Wg,DumpPath+"Wg.dat");	
}
void FDEmSampleInverter::fillLtL()
{ 	
	if(AlphaS == 0 || nlayers<3) return;
	L = Matrix<double>(nlayers - 2, nparam);
	LtL = Matrix<double>(nparam, nparam);

	L=0.0;	

	std::vector<double> t(nlayers);
	if(nlayers==1){
		t[0]=1.0;
	}
	else{
		for(size_t i=0;i<(nlayers-1);i++){			
			t[i]=ER.thickness[i];
		}
		t[nlayers-1]=(t[nlayers-2]/t[nlayers-3])*t[nlayers-2];
	}	

	if(nlayers>2 && SolveConductivity){
		size_t neqn=0;
		if(thickness_weight_smoothness){									
			for(size_t li=1; li<nlayers-1; li++){		
				double t1,t2,t3;
				size_t pindex = cIndex + li;
				t1 = t[li-1];
				t2 = t[li];
				t3 = t[li+1];												
				L[neqn][pindex-1] =  1.0/(0.5*(t1+t2));
				L[neqn][pindex  ] = -1.0/(0.5*(t1+t2)) - 1.0/(0.5*(t2+t3));
				L[neqn][pindex+1] =  1.0/(0.5*(t2+t3));				
				neqn++;
			}							
		}
		else{
			for(size_t li=1; li<nlayers-1; li++){			
				size_t pindex = cIndex + li;
				L[neqn][pindex-1]=  1.0;
				L[neqn][pindex  ]= -2.0;
				L[neqn][pindex+1]=  1.0;
				neqn++;
			}							
		}
	}	
	const double s = (1.0/(double)(nlayers-2));
	LtL = s*transpose_mult(L,L);

	//if(Dump){
	//writetofile(L,DumpPath+"L.dat");		
	//writetofile(LtL,DumpPath+"LtL.dat");
	//}
}

void FDEmSampleInverter::dumptofile(const std::vector<double>& v, std::string path) 
{		
	FILE* fp = fileopen(DumpPath+path,"w");	
	for(size_t c=0; c<nCoilsets; c++){	  
		fprintf(fp,"%lf\t",v[emIndex+c*2]);		
		fprintf(fp,"%lf\t",v[emIndex+c*2+1]);
		fprintf(fp,"\n");
	}	
	fclose(fp);	
}

void FDEmSampleInverter::dumptofile(const cEarth1D& e, std::string path) 
{	
	FILE* fp = fileopen(DumpPath+path,"w");
	size_t nlayers = e.conductivity.size();
	for(size_t i=0; i<nlayers; i++){		
		if(i<nlayers-1)fprintf(fp,"%lf\t%lf\n",e.conductivity[i],e.thickness[i]);	
		else fprintf(fp,"%lf\t200.0\n",e.conductivity[i]);
	}
	fclose(fp);
}

void FDEmSampleInverter::dumptofile(const sFDEmGeometry& g, std::string path) 
{		
	FILE* fp = fileopen(DumpPath+path,"w");	
	fprintf(fp,"height\t%lf\n",g.birdheight);				
	fclose(fp);	
}

double FDEmSampleInverter::phiData(const std::vector<double>& g) 
{		
	double val;
	double sum=0.0;
	for(size_t i=0; i<ndata; i++){
		val =  (vObs[i]-g[i])/vErr[i];
		sum += (val*val);
	}
	return (sum/ndata);		
}

double FDEmSampleInverter::phiModel(const std::vector<double>& p) 
{		
	double val;
	double sum=0.0;
	for(size_t i=0; i<nparam; i++){
		val = (vRefParam[i]-p[i])/vRefParamStd[i];
		sum += (val*val);
	}
	return (sum/nparam);		
}

double FDEmSampleInverter::phiModel(const std::vector<double>& p, double& phic, double& phit, double& phig, double& phis) 
{		
	phic = phiC(p);
	phit = phiT(p);
	phig = phiG(p);
	phis = phiS(p);

	double v=0;
	v += AlphaC * phic;
	v += AlphaT * phit;
	v += AlphaG * phig;
	v += AlphaS * phis;
	return v;	
}
double FDEmSampleInverter::phiC(const std::vector<double>& p) 
{	
	if(SolveConductivity==false)return 0.0;
	std::vector<double> v = p - vRefParam;	
	return mtDm(v,Wc);	
}
double FDEmSampleInverter::phiT(const std::vector<double>& p) 
{	
	if(SolveThickness==false)return 0.0;
	std::vector<double> v= p - vRefParam;
	return mtDm(v,Wt);	
}
double FDEmSampleInverter::phiG(const std::vector<double>& p) 
{	
	if(SolveBirdHeight==false)return 0.0;
	std::vector<double> v = p - vRefParam;
	return mtDm(v,Wg);
}
double FDEmSampleInverter::phiS(const std::vector<double>& p) 
{				
	if(AlphaS == 0)return 0.0;	
	else return mtAm(p,LtL);			
}

cEarth1D FDEmSampleInverter::get_earth(const std::vector<double>& parameters)
{ 		
	cEarth1D e = ER;
	if(SolveConductivity){
		for(size_t li=0; li<nlayers; li++){
			e.conductivity[li] = pow10(parameters[li+cIndex]);
		}
	}
	
	if(SolveThickness){
		for(size_t li=0; li<nlayers-1; li++){ 	
			e.thickness[li] = pow10(parameters[li+tIndex]);
		}  
	}
	return e;			
}
sFDEmGeometry FDEmSampleInverter::get_geometry(const std::vector<double>& parameters)
{ 	
	sFDEmGeometry g = GR;		
	if(SolveBirdHeight){
		g.birdheight = parameters[birdheightIndex];
	}
	return g;
}
sFDEmData FDEmSampleInverter::get_data(const std::vector<double>& data)
{ 		
	sFDEmData d;	
	d.inphase.resize(nCoilsets);
	d.quadrature.resize(nCoilsets);	
	for(size_t i=0; i<nCoilsets; i++){
		d.inphase[i]    = data[emIndex+i*2];
		d.quadrature[i] = data[emIndex+i*2+1];
	}	
	return d;
}

void FDEmSampleInverter::forwardmodel(const std::vector<double>& parameters, std::vector<double>& predicted, bool computederivatives)
{ 	
	cEarth1D e = get_earth(parameters);
	sFDEmGeometry g = get_geometry(parameters);
			
	EMSystem.setearth(e);
	EMSystem.setrollpitchyaw(g.birdroll,g.birdpitch,g.birdyaw);	
	EMSystem.setheight(g.birdheight);
	EMSystem.setupcomputations();

	//Forwardmodel
	cvector fm;
	fm = EMSystem.ppms();		

	for(size_t i=0; i<nCoilsets; i++){
		predicted[emIndex+i*2] = fm[i].real();
		predicted[emIndex+i*2+1] = fm[i].imag();
	}	
	
	if(computederivatives){
		//Derivatives		
		cvector deriv;

		J=0.0;//initialise to zero
		
		if(SolveConductivity){
			for(size_t li=0; li<nlayers; li++){			
				deriv = EMSystem.dppms(eDC,li);			
				size_t pindex = li+cIndex;
				for(size_t i=0; i<nCoilsets; i++){
					//multiply by natural log(10) as parameters are in logbase10 units
					J[emIndex+i*2][pindex] = log(10.0)*e.conductivity[li]*deriv[i].real();
					J[emIndex+i*2+1][pindex] = log(10.0)*e.conductivity[li]*deriv[i].imag();					
				}
			}
		}

		if(SolveThickness){
			for(size_t li=0; li<nlayers-1; li++){			
				deriv = EMSystem.dppms(eDT,li);	
				size_t pindex = li+tIndex;
				for(size_t i=0; i<nCoilsets; i++){
					//multiply by natural log(10) as parameters are in logbase10 units
					J[emIndex+i*2][pindex] = log(10.0)*e.thickness[li]*deriv[i].real();
					J[emIndex+i*2+1][pindex] = log(10.0)*e.thickness[li]*deriv[i].imag();
				}			
			}
		}

		if(SolveBirdHeight){			
			deriv = EMSystem.dppms(eDB,0);		
			size_t pindex = birdheightIndex;			
			for(size_t i=0; i<nCoilsets; i++){
				J[emIndex+i*2][pindex]   = deriv[i].real();
				J[emIndex+i*2+1][pindex] = deriv[i].imag();
			}
		}		
		
		JtWd  = transpose_mult(J,Wd);
		JtWdJ = mult(JtWd,J);
	
		if(Dump){		
			writetofile(J,DumpPath + "J.dat");	
			writetofile(JtWd,DumpPath + "JtWd.dat");			
			writetofile(JtWdJ,DumpPath + "JtWdJ.dat");
		}
	}	
}

void FDEmSampleInverter::fillAb(const double lambda)
{ 	
	//Ax = b
	//A = [J'WdJ + lambda*Wm + lanbda*LtL]
	//x = m(n+1)
	//b = J'Wd(d - g(m) + Jm) + labda*Wm*m0
	//dm = m(n+1) - m = x - m

	const std::vector<double>& m  = vParam;
	const std::vector<double>& d  = vObs;
	const std::vector<double>& g  = vPred;
	const std::vector<double>& m0 = vRefParam;	
	
	double zc=lambda*AlphaC;
	double zt=lambda*AlphaT;
	double zg=lambda*AlphaG;
	double zs=lambda*AlphaS;

	Wm = 0.0;
	if(zc>0) Wm += zc*Wc;
	if(zt>0) Wm += zt*Wt;
	if(zg>0) Wm += zg*Wg;

	b = JtWd*(d - g + J*m) + Wm*m0;			
	if(zs>0){
		A = JtWdJ + Wm + zs*LtL;
	}
	else{
		A = JtWdJ + Wm;	
	}

	if(Dump){		
		writetofile(d,DumpPath+"d.dat");		
		writetofile(g,DumpPath+"g.dat");		
		writetofile(m,DumpPath+"m.dat");		
		writetofile(m0,DumpPath+"m0.dat");		
		writetofile(J,DumpPath+"J.dat");		
		writetofile(Wc,DumpPath+"Wc.dat");		
		writetofile(Wt,DumpPath+"Wt.dat");		
		writetofile(A,DumpPath+"A.dat");		
		writetofile(b,DumpPath+"b.dat");		
		writetofile(A,DumpPath+"A.dat");			
	}
}

void FDEmSampleInverter::solveAxb()
{ 					
	Matrix<double> pinvA = pseudoinverse(A);
	x = pinvA*b;

	//writetofile(A,DumpPath+"A.dat");	
	//writetofile(pinvA,DumpPath+"pinvA.dat");	
	//prompttocontinue();
}

void FDEmSampleInverter::invert()
{
	if (Dump){
		FILE* fp = fileopen(DumpPath + "record.dat", "w");
		fprintf(fp, "Record\t%lu", DataFileRecord);
		fclose(fp);
	}

	double t1 = gettime();
	fillsample();
	iterate();
	writeresult();
	double t2 = gettime();

	double etime = t2 - t1;
	message(fp_log,"Rec %6lu\t %3lu\t %5lu\t %10lf ...",DataFileRecord,Id.flightnumber,Id.linenumber,Id.fidnumber);
	message(fp_log,"Its=%3lu\tPhiD=%6.2lf\t%s time=%.1lfs\n",LastIteration,LastPhiD,TerminationReason.c_str(),etime);
}
void FDEmSampleInverter::iterate()
{				
	std::vector<double> dm(nparam);		
	std::vector<double> gtemp(ndata);
	std::vector<double> mtemp(nparam);	
	double phidtemp,phimtemp,phictemp,phittemp,phigtemp,phistemp;

	size_t iteration = 0;
	vParam = vRefParam;
	forwardmodel(vParam,vPred,false);	
	LastPhiD   = phiData(vPred);
	LastPhiM   = phiModel(vParam,phictemp,phittemp,phigtemp,phistemp);	
	LastPhiC   = phictemp;
	LastPhiT   = phittemp;
	LastPhiG   = phigtemp;
	LastPhiS   = phistemp;		
	LastLambda = 1e8;	
	TargetPhiD = LastPhiD * 0.7;
	LastIteration = 0;


	TerminationReason = "Has not terminated";	
	cEarth1D earth = get_earth(vParam);
	sFDEmGeometry geometry = get_geometry(vParam);	

	bool keepiterating=true;
	if(LastPhiD < MinimumPhiD){
		keepiterating=false;		
		TerminationReason = "Reached minimum";	
	}

	while(keepiterating==true){		
		iteration++;
		TargetPhiD = LastPhiD * 0.7;

		mtemp  = vParam;
		gtemp  = vPred;	

		if(TargetPhiD < MinimumPhiD)TargetPhiD=MinimumPhiD;
		forwardmodel(mtemp,gtemp,true);

		cTrial t = targetsearch(LastLambda,TargetPhiD);

		dm = parameterchange(t.lambda);		
		mtemp = vParam + (t.stepfactor*dm);
		forwardmodel(mtemp,gtemp,false);
		phidtemp = phiData(gtemp);
		phimtemp = phiModel(mtemp,phictemp,phittemp,phigtemp,phistemp);				
		double percentchange = 100.0*(LastPhiD-phidtemp)/(LastPhiD);					

		if(phidtemp<LastPhiD){						
			vParam = mtemp;
			vPred  = gtemp;	
			 			
			DM = get_data(vPred);
			EM = get_earth(vParam);
			GM = get_geometry(vParam);	
					
			LastPhiD   = phidtemp;
			LastPhiM   = phimtemp;
			LastPhiC   = phictemp;
			LastPhiT   = phittemp;
			LastPhiG   = phigtemp;
			LastPhiS   = phistemp;
			LastLambda = t.lambda;
			LastIteration = iteration;

			if(LastPhiD <= MinimumPhiD){
				keepiterating = false;
				TerminationReason = "Reached minimum";					
			}
			else if(percentchange < MinimumImprovement){
				keepiterating = false;
				TerminationReason = "Small % improvement";					
			}					
		}		
		else{
			keepiterating = false;
			TerminationReason = "No improvement";
		}	

		if(iteration >= MaxIterations){
			keepiterating=false;
			TerminationReason = "Too many iterations";
		}		
	}
	
	//forwardmodel(vParam,gtemp,true);	
	//LayerSensitivity1 = compute_doi_sensitivity1();	
	//LayerSensitivity2 = compute_doi_sensitivity2();	

	if(Dump){ 		
		dumptofile(earth,"earth_inv.dat");		
		dumptofile(geometry,"geometry_inv.dat");		
		dumptofile(vPred,"predicted.dat");					
		FILE* fp=fileopen(DumpPath + "iteration.dat","w");
		fprintf(fp,"Iteration\t%d\n",(int)LastIteration);
		fprintf(fp,"TargetPhiD\t%lf\n",TargetPhiD); 	
		fprintf(fp,"PhiD\t%lf\n",LastPhiD); 	
		fprintf(fp,"Lambda\t%lf\n",LastLambda); 	
		fprintf(fp,"\n");		
		fclose(fp); 					
	}
}
void FDEmSampleInverter::finish()
{
	fclose(ofp);

	message(fp_log,"Logfile closing at %s\n", timestamp().c_str());
	fclose(fp_log);
}

std::vector<double> FDEmSampleInverter::parameterchange(const double lambda)
{ 	
	fillAb(lambda);
	solveAxb();			
	std::vector<double> dm = x - vParam;	

	if(SolveConductivity){
		for(size_t li=0; li<nlayers; li++){			
			size_t pindex = li+cIndex;
			const double maxcond=10;
			const double mincond=1e-6;
			if(vParam[pindex] + dm[pindex] > log10(maxcond)){
				//printf("upper limit li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] =  log10(maxcond) - vParam[pindex];
			}
			else if(vParam[pindex] + dm[pindex] < log10(mincond)){
				//printf("lower limit li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] =  log10(mincond) - vParam[pindex];
			}
		}			
	}

	if(SolveThickness){
		for(size_t li=0; li<nlayers-1; li++){			
			size_t pindex = li+tIndex;
			if(dm[pindex] > 0.5){
				//printf("li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] =  0.5;				
			}
			else if(dm[pindex] < -0.5){				
				//printf("li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = -0.5;
			}
		}			
	}

	if(SolveBirdHeight){
		size_t pindex = birdheightIndex;
		if(dm[pindex] > 0.5){			
			//printf("li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
			dm[pindex] =  0.5;				
		}
		else if(dm[pindex] < -0.5){				
			//printf("li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
			dm[pindex] = -0.5;
		}

		if(vParam[pindex] + dm[pindex] > 1000){				
			dm[pindex] =  1000 - vParam[pindex];				
		}
		else if(vParam[pindex] + dm[pindex] < 10){				
			dm[pindex] =  10 - vParam[pindex];				
		}
	}

	//writetofile(x,DumpPath+"x.dat");		
	//writetofile(dm,DumpPath+"dm.dat");		
	return dm;		
}
cTrial FDEmSampleInverter::targetsearch(const double currentlambda, const double target)
{				
	cTrialCache T;
	T.target = target;		
	bracket_t b = brackettarget(T,target,currentlambda);
	//printtrials(T);
	//prompttocontinue();
	
	if(b==BRACKETED){
		//bracketed target - find with Brents Method
		double newphid=DBL_MIN;
		double lambda = brentsmethod(T,target,newphid);	
		//printtrials(T);
		//prompttocontinue();
		return T.findlambda(lambda);		
	}	
	else if(b==MINBRACKETED){
		//bracketed minimum but above target - take smallest phid
		return T.minphidtrial();		
	} 
	else if(b==ALLBELOW){
		//all below target	- take the largest lambda			
		return T.maxlambdatrial();		
	}
	else if(b==ALLABOVE){
		//all above target - take smallest phid
		return T.minphidtrial();		
	}
	else{
		errormessage("Error unknown value %d returned from target brackettarget()\n",b);
	}	
	return T.minphidtrial();
}
bool FDEmSampleInverter::istargetbraketed(cTrialCache& T)
{				
	double target = T.target;
	T.sort_lambda();		
	for (size_t i = 0; i<T.trial.size()-1; i++){
		if(T.trial[i].phid >= target && T.trial[i+1].phid<=target){			
			return true;
		}
		if(T.trial[i].phid <=target && T.trial[i+1].phid>=target){			
			return true;
		}
	}		
	return false;	
}
bool FDEmSampleInverter::isminbraketed(cTrialCache& T)
{		
	size_t index = T.minphidindex();	
	if(index==0 || index==T.trial.size()-1){		
		return false;
	}

	double fa=T.trial[index-1].phid;
	double fb=T.trial[index].phid;
	double fc=T.trial[index+1].phid;	
	if ((fb<fa) && (fb<fc)){					
		return true;
	}		
	return false;	
}
bracket_t FDEmSampleInverter::brackettarget(cTrialCache& T, const double target, const double currentlambda)
{		
	double startx = log10(currentlambda);	
	if(LastIteration == 0){		
		std::vector<double> x;		
		x.push_back(8);x.push_back(6);
		x.push_back(4);x.push_back(2);
		x.push_back(1);x.push_back(0);		
		for(size_t k=0;k<x.size();k++){
			trialfunction(T,pow10(x[k]));
			bool tarbrak = istargetbraketed(T);			
			if(tarbrak){											
				return BRACKETED;//target bracketed		
			}					
		}

		double minv=DBL_MAX;
		for(size_t k=0;k<T.trial.size();k++){
			if(fabs(T.trial[k].phid-target) < minv){
				minv = fabs(T.trial[k].phid-target);
				startx = log10(T.trial[k].lambda);
			}
		}		
	}
	else{		
		trialfunction(T,pow10(startx));
	}

	std::vector<double> x;		
	x.push_back(+1); x.push_back(-1);
	x.push_back(+2); x.push_back(-2);
	x.push_back(+3); x.push_back(-3);
	for(size_t k=0; k<x.size(); k++){		
		trialfunction(T,pow10(startx + x[k]));				
		bool tarbrak = istargetbraketed(T);
		if(tarbrak){				
			return BRACKETED;//target bracketed		
		}
	}
	//printtrials(T);
	//prompttocontinue();

	if(T.maxphid() < target){
		return ALLBELOW;//all below target	
	}

	bool minbrak = isminbraketed(T);		
	if(minbrak)return MINBRACKETED;//min bracketed											

	return ALLABOVE;//all above target
}
double FDEmSampleInverter::brentsmethod(cTrialCache& T, const double target, double& newphid)
{
	T.sort_lambda();
	size_t index=T.trial.size()-1;
	for(size_t i=T.trial.size()-1;i>=1;i--){
		double f1=T.trial[i].phid   - target;
		double f2=T.trial[i-1].phid - target;
		if(f1 * f2 <= 0.0){
			index=i;
			break;
		}			
	}

	//Adapted from http://en.wikipedia.org/wiki/Brent's_method
	double xerrorTol=0.01;
	double yerrorTol=target*0.1;//10% accuracy is good enough

	double a = log10(T.trial[index-1].lambda);
	double b = log10(T.trial[index].lambda);
	double fa = T.trial[index-1].phid - target;
	double fb = T.trial[index].phid   - target;	
	if(fa * fb >= 0.0){
		warningmessage("Target must be bracketed for cSBSInverter::brentsmethod()\n");
	}	

	double c = 0;
	double d = DBL_MAX;
	double fc = 0;
	double s = 0;
	double fs = 0;

	// if f(a) f(b) >= 0 then error-exit
	if (fa * fb >= 0)
	{
		if (fa < fb){
			newphid=fa+target;
			return pow10(a);
		}
		else{
			newphid=fb+target;
			return pow10(b);
		}
	}

	// if |f(a)| < |f(b)| then swap (a,b) end if
	if (fabs(fa) < fabs(fb)){
		double tmp;
		tmp =  a;   a = b;  b = tmp;
		tmp = fa; fa = fb; fb = tmp;
	}

	c = a;
	fc = fa;
	bool mflag = true;
	int i = 0;

	while (!(fb==0) && (fabs(a-b) > xerrorTol))
	{
		if ((fa != fc) && (fb != fc))
			// Inverse quadratic interpolation
			s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
		else
			// Secant Rule
			s = b - fb * (b - a) / (fb - fa);

		double tmp2 = (3 * a + b) / 4;
		if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) || (mflag && (fabs(s - b) >= (fabs(b - c) / 2))) || (!mflag && (fabs(s - b) >= (fabs(c - d) / 2))))
		{
			s = (a + b) / 2;
			mflag = true;
		}
		else
		{
			if ((mflag && (fabs(b - c) < xerrorTol)) || (!mflag && (fabs(c - d) < xerrorTol)))
			{
				s = (a + b) / 2;
				mflag = true;
			}
			else
				mflag = false;
		}
		fs = trialfunction(T,pow10(s)) - target;		
		if(fabs(fs)<yerrorTol){
			newphid =fs+target;
			return pow10(s);
		}

		d = c;
		c = b;
		fc = fb;
		if (fa * fs < 0) { b = s; fb = fs; }
		else { a = s; fa = fs; }

		// if |f(a)| < |f(b)| then swap (a,b) end if
		if (fabs(fa) < fabs(fb))
		{ double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp; }
		i++;
		if (i > 20){
			printtrials(T);
			warningmessage("Too many bisections in cSBSInverter::brentsmethod()\n");	
			newphid=fb+target;
			return pow10(b);
		}
	}
	if(fb<fa){
		newphid=fb+target;
		return pow10(b);
	}
	else{
		newphid=fa+target;
		return pow10(a);
	}
}
void FDEmSampleInverter::printtrials(cTrialCache T)
{	
	T.sort_lambda();
	printf("\n");
	printf("CurrentLambda = %lf CurrentPhid = %lf    Target = %lf\n",LastLambda,LastPhiD,T.target);
	printf("N    Stepfactor       Lambda          Phid\n");
	for(size_t i=0;i<T.trial.size(); i++){
		printf("%2d %12g %12g %12g\n",(int)T.trial[i].order,T.trial[i].stepfactor,T.trial[i].lambda,T.trial[i].phid);
	}
	printf("\n");	
}
double FDEmSampleInverter::goldensearch(double a, double b, double c, double xtol, const double lambda, const std::vector<double>& m, const std::vector<double>& dm, std::vector<double>& g, cTrialCache& cache)
{
	//adapted from http://en.wikipedia.org/wiki/Golden_section_search	
	double resphi = 2 - ((1 + sqrt(5.0))/2.0);
	double x;
	if(c - b > b - a){
		x = b + resphi * (c - b);
	}
	else{
		x = b - resphi * (b - a);
	}

	//if(fabs(c - a) < tau * (fabs(b) + fabs(x))){
	//  return (c + a) / 2; 
	//}
	if(fabs(c-a)<xtol){
		return (c + a) / 2; 
	}

	double fx = cache.sfsearch(x);	
	if(fx<0){			
		cTrial t;
		std::vector<double> p = m + x*dm;
		forwardmodel(p, g, false);
		fx = phiData(g);				
		t.stepfactor=x;
		t.phid=fx;
		t.phim=phiModel(p);
		t.order=cache.trial.size();
		t.lambda=lambda;
		cache.trial.push_back(t);
		//printf("%lf %lf\n",x,fx);
	}

	double fb = cache.sfsearch(b);	
	if(fb<0){		
		cTrial t;
		std::vector<double> p = m + b*dm;
		forwardmodel(p, g, false);
		fb = phiData(g);		
		t.stepfactor=b;
		t.phid=fb;
		t.phim=phiModel(p);
		t.order=cache.trial.size();
		t.lambda=lambda;
		cache.trial.push_back(t);
		//printf("%lf %lf\n",b,fb);
	}	

	assert(fx != fb);
	if(fx < fb){
		if(c - b > b - a){
			return goldensearch(b, x, c, xtol, lambda, m, dm, g, cache);
		}
		else{
			return goldensearch(a, x, b, xtol, lambda, m, dm, g, cache);
		}
	}
	else{
		if(c - b > b - a){
			return goldensearch(a, b, x, xtol, lambda, m, dm, g, cache);
		}
		else{
			return goldensearch(x, b, c, xtol, lambda, m, dm, g, cache);
		}
	}
}
double FDEmSampleInverter::trialfunction(cTrialCache& T, const double triallambda)
{			
	std::vector<double> dm(nparam);	
	std::vector<double> p(nparam);
	std::vector<double> g(ndata);	
	dm = parameterchange(triallambda);
	cTrialCache cache;
	cTrial t0;		
	t0.phid  = LastPhiD;
	t0.phim  = LastPhiM;
	t0.stepfactor = 0.0;	
	t0.lambda = triallambda;	
	t0.order  = cache.trial.size();
	cache.trial.push_back(t0);

	cTrial t1;		
	p  = vParam + dm;
	forwardmodel(p,g,false);
	t1.phid = phiData(g);	
	t1.phim = phiModel(p);
	t1.stepfactor = 1.0;	
	t1.lambda = triallambda;	
	t1.order  = cache.trial.size();
	cache.trial.push_back(t1);

	double pcdiff = 100 * (t1.phid - t0.phid)/t0.phid;
	if(pcdiff > 0.0 || pcdiff<-1.0){
		//ie dont do not do golden search
		//if only tiny improvement				
		double xtol = 0.1; 
		double gsf = goldensearch(0.0,0.38196601125010510,1.0,xtol,triallambda,vParam,dm,g,cache);
		cTrial t3;		
		p  = vParam + gsf*dm;
		forwardmodel(p,g,false);
		t3.phid = phiData(g);	
		t3.phim = phiModel(p);
		t3.stepfactor = gsf;	
		t3.lambda = triallambda;	
		t3.order  = cache.trial.size();
		cache.trial.push_back(t3);
		//double gsf = goldensearch(0.0,0.5,1.0,tau,vParam,dm,g,cache);
		//printf("gsf=%lf\n",gsf);				
	}	
	//printtrials(cache);
	//prompttocontinue();
	size_t minindex = cache.minphidindex();	

	cTrial t = cache.trial[minindex];	
	t.order  = T.trial.size();
	T.trial.push_back(t);				
	return t.phid;
}

int FDEmSampleInverter::intvalue(const FieldDefinition& cd)
{ 	
	return cd.intvalue(DataFileFieldStrings);
}
double FDEmSampleInverter::doublevalue(const FieldDefinition& cd)
{ 	
	return cd.doublevalue(DataFileFieldStrings);
}
std::vector<int> FDEmSampleInverter::intvector(const FieldDefinition& cd, const size_t& n)
{
	return cd.intvector(DataFileFieldStrings, n);
}
std::vector<double> FDEmSampleInverter::doublevector(const FieldDefinition& cd, const size_t& n)
{
	return cd.doublevector(DataFileFieldStrings, n);
}

int process(int argc, char** argv, size_t Size, size_t Rank)
{
	if (argc >= 2){
		message("Executing %s %s\n", argv[0], argv[1]);
		message("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		message("Working directory %s\n", getcurrentdirectory().c_str());
	}
	else{
		errormessage("Not enough input arguments\n");
	}

	std::string controlfile(argv[1]);
	
	FDEmSampleInverter I(controlfile,Size,Rank);	
	size_t record = 0;
	while (I.readnextrecord()){		
		record++;
		if ((record - 1) % Size != Rank)continue;	
		I.parserecord();
		I.invert();
	}
	I.finish();
	return 0;
};

int main(int argc, char** argv)
{
#ifdef _MPI_ENABLED
	int Size, Rank;
	if (argc<2){
		printf("Executing %s\n", argv[0]);
		printf("Usage: %s control_file_name\n", argv[0]);
		printf("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		return 0;
	}
	message("Starting MPI\n");
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	int len;
	char processor_name[MPI_MAX_PROCESSOR_NAME + 1];
	MPI_Get_processor_name(processor_name, &len);
	processor_name[len] = '\0';
	message("MPI Started Processes=%d\tRank=%d\tProcessor name = %s\n", Size, Rank, processor_name);
	process(argc, argv, (size_t)Size, (size_t)Rank);
	MPI_Finalize();
#elif defined _OPENMP
	int Size;
	if (argc<2){
		printf("Executing %s\n", argv[0]);
		printf("Usage: %s control_file_name number_of_threads\n", argv[0]);
		printf("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		return 0;
	}
	if (argc == 2){
		Size = 1;
	}
	else{
		Size = atoi(argv[2]);
		int maxthreads = omp_get_max_threads();
		if (Size<1){
			printf("Error: %d is a stupid number of threads\n", Size);
			exit(1);
		}
		if (Size >= maxthreads){
			printf("Warning: The number of requested threads (%d) is more than the processors available (%d)\n", Size, maxthreads);
		}
	}

	omp_init_lock(&fftw_thread_lock);
#pragma omp parallel num_threads(Size)
	{
		int Rank = omp_get_thread_num();
		message("OpenMP threading Processes=%d\tRank=%d\n", Size, Rank);
		process(argc, argv, Size, Rank);
	}

#else
	int Size = 1;
	int Rank = 0;
	if (argc<2){
		printf("Executing %s\n", argv[0]);
		printf("Usage: %s control_file_name number_of_threads\n", argv[0]);
		printf("Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		return 0;
	}
	message("Stand-alone Processes=%d\tRank=%d\tProcessor name = %s\n", Size, Rank);
	process(argc, argv, Size, Rank);
#endif	
	return 0;
}


