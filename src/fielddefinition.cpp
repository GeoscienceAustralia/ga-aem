/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include "float.h"
#include "stdlib.h"
#include "general_utils.h"
#include "fielddefinition.h"

FieldDefinition::FieldDefinition()
{
	firstcolumn=1;
}
void FieldDefinition::set(const cBlock& b, const std::string& fieldname)
{ 	
	int col;	
	std::string rhs = b.getstringvalue(fieldname);
	if(strncasecmp(rhs,"Column",6)==0){
		sscanf(&(rhs.c_str()[6]),"%d",&col);
		coldeftype = COLDEF_COLUMNVALUE;
		column     = (size_t)col;		
	}
	else if(strncasecmp(rhs,"-Column",7)==0){
		sscanf(&(rhs.c_str()[7]),"%d",&col);
		coldeftype = COLDEF_NEGATIVECOLUMNVALUE;
		column = (size_t)col;		
	}
	else if(strncasecmp(rhs,"Unavailable",11)==0){
		coldeftype = COLDEF_UNAVAILABLE;
		column     = 0;		
	}
	else{		
		//sscanf(rhs.c_str(),"%lf",&val);
		coldeftype = COLDEF_DEFAULTVALUE;
		column = (size_t)0;		
		defaultvector = b.getdoublevector(fieldname);
	}	
}
int FieldDefinition::intvalue(const std::vector<std::string>& fields) const
{ 	
	int v;
	if(coldeftype == COLDEF_COLUMNVALUE){	
		sscanf(fields[column-firstcolumn].c_str(),"%d",&v);
		return v;	    
	}
	else if(coldeftype == COLDEF_NEGATIVECOLUMNVALUE){		
		sscanf(fields[column-firstcolumn].c_str(),"%d",&v);
		return -v;
	}
	else if(coldeftype == COLDEF_DEFAULTVALUE){				
		if(defaultvector.size()==0){
			return cBlock::ud_int();
		}
		return (int)defaultvector[0];
	}
	else if(coldeftype == COLDEF_UNAVAILABLE){
		return cBlock::ud_int();
	}
	else{
		warningmessage("FieldDefinition::intvalue() unknown column definition\n");
	}
	return cBlock::ud_int();
	
}
double FieldDefinition::doublevalue(const std::vector<std::string>& fields) const 
{ 	
	double v;
	if(coldeftype == COLDEF_COLUMNVALUE){	
		sscanf(fields[column-firstcolumn].c_str(),"%lf",&v);
		return v;	    
	}
	else if(coldeftype == COLDEF_NEGATIVECOLUMNVALUE){		
		sscanf(fields[column-firstcolumn].c_str(),"%lf",&v);
		return -v;
	}
	else if(coldeftype == COLDEF_DEFAULTVALUE){				
		if(defaultvector.size()==0){
			return cBlock::ud_double();
		}
		return defaultvector[0];
	}
	else if(coldeftype == COLDEF_UNAVAILABLE){
		return cBlock::ud_double();
	}
	else{
		warningmessage("FieldDefinition::doublevalue() unknown column definition\n");
	}
	return cBlock::ud_double();
}
std::vector<int> FieldDefinition::intvector(const std::vector<std::string>& fields, const size_t& n) const
{ 
	int v;
	if(coldeftype == COLDEF_COLUMNVALUE){		
		std::vector<int> vec(n);
		for(size_t i=0; i<n; i++){		   
		   sscanf(fields[i+column-firstcolumn].c_str(),"%d",&v);
           vec[i] = v;	    
		}
		return vec;
	}	
	else if(coldeftype == COLDEF_NEGATIVECOLUMNVALUE){		
		std::vector<int> vec(n);
		for(size_t i=0; i<n; i++){			
		   sscanf(fields[i+column-firstcolumn].c_str(),"%d",&v);
           vec[i] = -v;	    
		}
		return vec;
	}
	else if(coldeftype == COLDEF_DEFAULTVALUE){				
		size_t deflen = defaultvector.size();				
		std::vector<int> vec(n);
		for(size_t i=0; i<n; i++){
			if(deflen==1) vec[i] = (int)defaultvector[0];
			else vec[i] = (int)defaultvector[i];
		}		
		return vec;
	}
	else if(coldeftype == COLDEF_UNAVAILABLE){		
		std::vector<int> vec(0);
		return vec;
	}
	else{
		warningmessage("FieldDefinition::intvector() unknown column definition type\n");		
	}
	std::vector<int> vec(0);
	return vec;
}
std::vector<double> FieldDefinition::doublevector(const std::vector<std::string>& fields, const size_t& n) const
{ 	
	double v;	
	if(coldeftype == COLDEF_COLUMNVALUE){		
		std::vector<double> vec(n);
		for(size_t i=0; i<n; i++){		   
		   sscanf(fields[i+column-firstcolumn].c_str(),"%lf",&v);
           vec[i] = v;	    
		}
		return vec;
	}	
	else if(coldeftype == COLDEF_NEGATIVECOLUMNVALUE){		
		std::vector<double> vec(n);
		for(size_t i=0; i<n; i++){			
		   sscanf(fields[i+column-firstcolumn].c_str(),"%lf",&v);
           vec[i] = -v;	    
		}
		return vec;
	}
	else if(coldeftype == COLDEF_DEFAULTVALUE){				
		std::vector<double> vec(n);
		int deflen = (int)defaultvector.size();				
		for(size_t i=0; i<n; i++){
			if(deflen==1)vec[i] = defaultvector[0];
			else vec[i] = defaultvector[i];
		}		
		return vec;
	}
	else if(coldeftype == COLDEF_UNAVAILABLE){		
		std::vector<double> vec(0);
		return vec;
	}
	else{
		warningmessage("FieldDefinition::doublevector() unknown column definition type\n");
	}		
	std::vector<double> vec(0);
	return vec;
}



