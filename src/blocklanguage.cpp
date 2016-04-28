/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <iostream>
#include <cstdlib>

#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"

cBlock::cBlock()
{
	
}

cBlock::cBlock(const std::string& filename)
{
	loadfromfile(filename);
}

void cBlock::loadfromfile(const std::string& filename)
{
	Filename = filename;
	FILE* fp = fileopen(filename, "r");
	if (fp == NULL){
		printf("cBlock::loadfromfile - Could not open file: %s\n", filename.c_str());
		throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
	}
	loadfromfilepointer(fp, true);
	fclose(fp);
}

void cBlock::loadfromfilepointer(FILE* fp, bool rootlevel)
{
	std::string linestr;
	while (filegetline(fp, linestr)){
		std::string s = trim(linestr);
		if (s.length() == 0)continue;

		std::vector<std::string> vs = tokenize(s);
		if (vs.size() >= 2){
			if (strcasecmp(vs[1], "End") == 0){
				break;
			}
			else if (strcasecmp(vs[1], "Begin") == 0){
				if (rootlevel == true){
					Name = vs[0];
					rootlevel = false;
				}
				else{
					cBlock newblock;
					newblock.Name = vs[0];
					newblock.loadfromfilepointer(fp, false);
					Blocks.push_back(newblock);
				}
				continue;
			}
			Entries.push_back(s);
		}
		else{
			Entries.push_back(s);
		}
	}
}

void cBlock::print(size_t n) const
{
	for (size_t j = 0; j < n; j++) std::cout << "  ";
	 std::cout << Name.c_str() << " Begin" << std::endl;

	for (size_t i = 0; i < Entries.size(); i++){
		for (size_t j = 0; j < n + 1; j++) std::cout << "  ";
		 std::cout << Entries[i].c_str() << std::endl;
	}

	for (size_t i = 0; i < Blocks.size(); i++){
		Blocks.at(i).print(n + 1);
	}

	for (size_t j = 0; j < n; j++) std::cout << "  ";
	 std::cout << Name.c_str() << " End" << std::endl;
	return;
}

void cBlock::write(FILE* fp, size_t n) const
{
	for (size_t j = 0; j < n; j++)fprintf(fp, "\t");
	fprintf(fp, "%s Begin\n", Name.c_str());

	for (size_t i = 0; i < Entries.size(); i++){
		for (size_t j = 0; j < n + 1; j++)fprintf(fp, "\t");
		fprintf(fp, "%s\n", Entries[i].c_str());
	}

	for (size_t i = 0; i < Blocks.size(); i++){
		Blocks[i].write(fp, n + 1);
	}

	for (size_t j = 0; j < n; j++)fprintf(fp, "\t");
	fprintf(fp, "%s End\n", Name.c_str());
	return;
}
std::string cBlock::identifier(std::string entry) const
{
	size_t index = entry.find("=");
	std::string id = entry.substr(0, index - 1);
	return trim(id);
}
void cBlock::printidentifiers() const
{
	for (size_t i = 0; i < Entries.size(); i++){
		 std::cout << identifier(Entries.at(i)).c_str() << std::endl;
	}

	for (size_t i = 0; i < Blocks.size(); i++){
		Blocks.at(i).printidentifiers();
	}

}
std::string cBlock::value(std::string entry) const
{
	size_t index = entry.find("=");
	size_t len = entry.size() - index - 1;
	if (len == 0)return std::string("");
	std::string s = entry.substr(index + 1, len);
	return trim(s);
}
void cBlock::printvalues()  const
{
	for (size_t i = 0; i < Entries.size(); i++){
		 std::cout << value(Entries.at(i)).c_str() << std::endl;
	}

	for (size_t i = 0; i < Blocks.size(); i++){
		Blocks.at(i).printvalues();
	}
}
std::string cBlock::getentry(std::string id) const
{
	size_t index = id.find(".");
	if (index != std::string::npos){
		std::string start = id.substr(0, index);
		std::string end = id.substr(index + 1, id.size() - index - 1);
		if (strcasecmp(start, Name) == 0){
			return getentry(end);
		}
		else{
			cBlock b = findblock(start);
			return b.getentry(end);
		}
	}
	else{
		return findidentifer(id);
	}
}
cBlock cBlock::findblock(std::string name) const
{
	size_t index = name.find(".");

	if (index != std::string::npos){
		std::string start = name.substr(0, index);
		std::string end = name.substr(index + 1, name.size() - index - 1);
		cBlock b = findblock(start);
		return b.findblock(end);
	}

	for (size_t i = 0; i < Blocks.size(); i++){
		if (strcasecmp(Blocks[i].Name, name) == 0){
			return Blocks[i];
		}
	}
	return cBlock();

}
std::vector<cBlock> cBlock::findblocks(std::string name) const
{
	std::vector<cBlock> v;
	size_t index = name.find(".");
	if (index != std::string::npos){
		std::string start = name.substr(0, index);
		std::string end = name.substr(index + 1, name.size() - index - 1);
		cBlock b = findblock(start);
		return b.findblocks(end);
	}

	for (size_t i = 0; i < Blocks.size(); i++){
		if (strcasecmp(Blocks[i].Name, name) == 0){
			v.push_back(Blocks[i]);
		}
	}
	return v;
};
std::string cBlock::findidentifer(std::string id) const
{
	for (size_t i = 0; i < Entries.size(); i++){
		if (strcasecmp(identifier(Entries[i]), id) == 0){
			return Entries[i];
		}
	}
	return ud_string();
}
std::string cBlock::getstringvalue(std::string id) const
{
	return value(getentry(id));
}

short cBlock::getshortvalue(std::string id) const
{
	short v;
	int status;
	std::string entry = getentry(id);
	if (strcasecmp(entry, ud_string()) == 0){
		return ud_short();
	}
	else{
		status = sscanf(value(entry).c_str(), "%hd", &v);
	}

	if (status == 1) return v;
	else return ud_short();
}
int cBlock::getintvalue(std::string id) const
{
	int v;
	int status;
	std::string entry = getentry(id);
	if (strcasecmp(entry, ud_string()) == 0){
		return ud_int();
	}
	else{
		status = sscanf(value(entry).c_str(), "%d", &v);
	}

	if (status == 1) return v;
	else return ud_int();
}
size_t cBlock::getsizetvalue(std::string id) const
{
	size_t v;
	int status;
	std::string entry = getentry(id);
	if (strcasecmp(entry, ud_string()) == 0){
		return ud_size_t();
	}
	else{		
		long tmp;		
		status = std::sscanf(value(entry).c_str(), "%ld", &tmp);
		v = (size_t)tmp;
	}

	if (status == 1) return v;
	else return ud_size_t();
}
float cBlock::getfloatvalue(std::string id) const
{
	float v;
	int status;
	std::string entry = getentry(id);
	if (strcasecmp(entry, ud_string()) == 0){
		return ud_float();
	}
	else{
		status = sscanf(value(entry).c_str(), "%f", &v);
	}

	if (status == 1) return v;
	else return ud_float();
}
double cBlock::getdoublevalue(std::string id) const
{
	double v;
	int status;
	std::string entry = getentry(id);
	if (strcasecmp(entry, ud_string()) == 0){
		return ud_double();
	}
	else{
		status = sscanf(value(entry).c_str(), "%lf", &v);
	}

	if (status == 1) return v;
	else return ud_double();
}

std::vector<int> cBlock::getintvector(std::string id) const
{
	int v;
	std::vector<int> vec;
	std::string s = getstringvalue(id);

	if (strcasecmp(s, ud_string()) == 0){
		return vec;
	}

	std::vector<std::string> fields = fieldparsestring(s.c_str(), delimiters.c_str());
	for (size_t i = 0; i < fields.size(); i++){
		sscanf(fields[i].c_str(), "%d", &v);
		vec.push_back(v);
	}
	return vec;
}
std::vector<double> cBlock::getdoublevector(std::string id) const
{
	double v;
	std::vector<double> vec;
	std::string s = getstringvalue(id);

	if (strcasecmp(s, ud_string()) == 0){
		return vec;
	}

	std::vector<std::string> fields = fieldparsestring(s.c_str(), delimiters.c_str());
	for (size_t i = 0; i < fields.size(); i++){
		sscanf(fields[i].c_str(), "%lf", &v);
		vec.push_back(v);
	}
	return vec;
}
std::vector<std::string> cBlock::getstringvector(std::string id) const
{
	std::vector<std::string> fields;
	std::string s = getstringvalue(id);

	if (strcasecmp(s, ud_string()) == 0){
		return fields;
	}

	fields = fieldparsestring(s.c_str(), delimiters.c_str());
	for (size_t i = 0; i < fields.size(); i++){
		fields[i] = trim(fields[i]);
	}
	return fields;
}
std::vector<std::vector<double>> cBlock::getdoublematrix(std::string id) const
{
	std::vector<std::vector<double>> matrix;
	cBlock b = findblock(id);
	for (size_t i = 0; i < b.Entries.size(); i++){
		double v;
		std::vector<double> vec;

		std::string s = b.Entries[i];
		std::vector<std::string> fields = fieldparsestring(s.c_str(), delimiters.c_str());
		for (size_t j = 0; j < fields.size(); j++){
			sscanf(fields[j].c_str(), "%lf", &v);
			vec.push_back(v);
		}
		matrix.push_back(vec);
	}
	return matrix;
}
bool cBlock::getboolvalue(std::string id) const
{
	std::string s = getstringvalue(id);
	size_t k = s.find_first_of(" \t\r\n");
	std::string value = s.substr(0, k);
	if (strcasecmp(value, "yes") == 0)return true;
	else if (strcasecmp(value, "no") == 0)return false;
	else if (strcasecmp(value, "true") == 0)return true;
	else if (strcasecmp(value, "false") == 0)return false;
	else if (strcasecmp(value, "1") == 0)return true;
	else if (strcasecmp(value, "0") == 0)return false;
	else if (strcasecmp(value, "off") == 0)return false;
	else if (strcasecmp(value, "on") == 0)return true;
	else return false;
}

bool cBlock::getvalue(std::string id, bool& value) const
{
	if (getentry(id).compare(ud_string()) == 0){
		//value = false;
		return false;
	}
	value = getboolvalue(id);
	return true;
}
bool cBlock::getvalue(std::string id, short& value) const
{
	if (getentry(id).compare(ud_string()) == 0){
		//value = ud_short();
		return false;
	}
	value = getshortvalue(id);
	return true;
}
bool cBlock::getvalue(std::string id, int& value) const
{
	if (getentry(id).compare(ud_string()) == 0){
		//value = ud_int();
		return false;
	}
	value = getintvalue(id);
	return true;
}
bool cBlock::getvalue(std::string id, size_t& value) const
{
	if (getentry(id).compare(ud_string()) == 0){
		//value = ud_size_t();
		return false;
	}
	value = getsizetvalue(id);
	return true;
}
bool cBlock::getvalue(std::string id, float& value) const
{
	if (getentry(id).compare(ud_string()) == 0){
		//value = ud_float();
		return false;
	}
	value = getfloatvalue(id);
	return true;
}
bool cBlock::getvalue(std::string id, double& value) const
{
	if (getentry(id).compare(ud_string()) == 0){
		//value = ud_double();
		return false;
	}
	value = getdoublevalue(id);
	return true;
}
bool cBlock::getvalue(std::string id, std::string& value) const
{
	if (getentry(id).compare(ud_string()) == 0){
		//value = ud_string();
		return false;
	}
	value = getstringvalue(id);
	return true;
}


std::vector<int>    cBlock::getmultipleints(std::string id) const
{
	std::vector<std::string> str = getmultiplestrings(id);
	std::vector<int> result;

	int val;
	for (size_t i = 0; i < str.size(); i++){
		sscanf(str[i].c_str(), "%d", &val);
		result.push_back(val);
	}
	return result;
}
std::vector<double> cBlock::getmultipledoubles(std::string id) const
{
	std::vector<std::string> str = getmultiplestrings(id);
	std::vector<double> result;

	double val;
	for (size_t i = 0; i < str.size(); i++){
		sscanf(str[i].c_str(), "%lf", &val);
		result.push_back(val);
	}
	return result;
}
std::vector<std::string> cBlock::getmultiplestrings(std::string id) const
{
	int i;
	std::vector<std::string> result;

	std::string str = getstringvalue(id);
	if (str.length() > 0 && strcasecmp(str, ud_string()) != 0){
		result.push_back(str);
	}
	for (i = 0; i < 100; i++){
		std::string token = id;
		token += strprint("%d", i);

		str = getstringvalue(token);
		if (strcasecmp(str, ud_string()) != 0){
			result.push_back(str);
		}
	}
	return result;
}
std::vector<std::string> cBlock::getblockstrings(std::string id) const
{
	std::vector<std::string> result;
	cBlock b = findblock(id);
	for (size_t i = 0; i < b.Entries.size(); i++){
		std::string s = b.Entries[i];
		result.push_back(s);
	}
	return result;
}

std::vector<std::vector<std::string>> cBlock::getblockleftright(const std::string id) const
{
	std::vector<std::string> s = getblockstrings(id);
	std::vector<std::vector<std::string>> v(s.size());
	for (size_t i = 0; i < s.size(); i++){
		v[i] = split(s[i], '=');
		if(v[i].size() == 1)v[i].push_back("");		
		v[i][0] = trim(v[i][0]);
		v[i][1] = trim(v[i][1]);
	}
	return v;
}



