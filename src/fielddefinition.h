/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

//---------------------------------------------------------------------------
#ifndef fielddefinitionH
#define fielddefinitionH

#include <cstring>
#include <vector>

#include "blocklanguage.h"

#define COLDEF_COLUMNVALUE  0
#define COLDEF_NEGATIVECOLUMNVALUE  1
#define COLDEF_DEFAULTVALUE 2
#define COLDEF_UNAVAILABLE  3


////////////////////////////////////////////////////////////////////////////////
class FieldDefinition{

public:

	FieldDefinition();
	void set(const cBlock& b, const std::string& fieldname);
	int intvalue(const std::vector<std::string>& fields) const;
	double doublevalue(const std::vector<std::string>& fields) const;	
	std::vector<int> intvector(const std::vector<std::string>& fields, const size_t& n) const;
	std::vector<double> doublevector(const std::vector<std::string>& fields, const size_t& n) const;
	
	int coldeftype;
	size_t column;	
	std::vector<double> defaultvector;
	size_t firstcolumn;

};
////////////////////////////////////////////////////////////////////////////////
#endif
