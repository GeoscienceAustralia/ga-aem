/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

//---------------------------------------------------------------------------
#ifndef _blocklanguage_H
#define _blocklanguage_H

#include <cstring>
#include <vector>
#include <climits>
#include <cfloat>
#include <stdint.h>

////////////////////////////////////////////////////////////////////////////////
class cBlock{

private:
	std::string delimiters = " ,\t";

public:

cBlock();
cBlock(const std::string& filename);

std::string Filename;
std::string Name;
std::vector<std::string> Entries;
std::vector<cBlock> Blocks;

void loadfromfile(const std::string& filemname);
void loadfromfilepointer(FILE* fp, bool rootlevel);

void print(size_t n=0) const;
void write(FILE* fp, size_t n=0) const;
std::string identifier(std::string entry) const;
void printidentifiers() const;
std::string value(std::string entry) const;
void printvalues() const;
std::string getentry(std::string id) const;
cBlock findblock(std::string name) const;
std::vector<cBlock> findblocks(std::string name) const;
std::string findidentifer(std::string id) const;

static std::string ud_string(){return "ENTRY NOT FOUND";}
static short  ud_short(){return SHRT_MIN;}
static size_t ud_size_t(){ return UINT64_MAX;}
static int    ud_int(){return INT_MIN;}
static float  ud_float(){return  -FLT_MAX;}
static double ud_double(){return -DBL_MAX;}



std::string getstringvalue(std::string id) const;
bool getboolvalue(std::string id) const;
short getshortvalue(std::string id) const;
int getintvalue(std::string id) const;
size_t getsizetvalue(std::string id) const;
float getfloatvalue(std::string id) const;
double getdoublevalue(std::string id) const;

std::vector<int> getintvector(std::string id) const;
std::vector<double> getdoublevector(std::string id) const;
std::vector<std::string> getstringvector(std::string id) const;
std::vector< std::vector<double> > getdoublematrix(std::string id) const;
std::vector<std::string> getblockstrings(std::string id) const;

//Get many std::strings from different lines with numeric end
std::vector<std::string> getmultiplestrings(std::string id) const;
std::vector<double> getmultipledoubles(std::string id) const;
std::vector<int> getmultipleints(std::string id) const;

bool getvalue(std::string id, bool& value) const;
bool getvalue(std::string id, short& value) const;
bool getvalue(std::string id, int& value) const;
bool getvalue(std::string id, size_t& value) const;
bool getvalue(std::string id, float& value) const;
bool getvalue(std::string id, double& value) const;
bool getvalue(std::string id, std::string& value) const;

std::vector<std::vector<std::string>> getblockleftright(const std::string id) const;

};
////////////////////////////////////////////////////////////////////////////////
#endif
