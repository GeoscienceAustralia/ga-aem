/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _file_utils_H
#define _file_utils_H

#include <cstdint>
#include <cstring>
#include <vector>


struct sFilePathParts{
	std::string directory;
	std::string prefix;
	std::string extension;	
};

char pathseparator();
std::string pathseparatorstring();
void fixseparator(std::string& path);
void removetrailingseparator(std::string& path);
void addtrailingseparator(std::string& path);
bool exists(std::string path);
bool isdirectory(std::string path);
bool isfile(std::string path);
bool isabsolutepath(std::string path);
FILE* fileopen(const std::string filepath, const std::string mode);
std::string getcurrentdirectory();
bool makedirectory(std::string dirname);
bool makedirectorydeep(std::string dirname);
int copyfile(std::string src, std::string dest);
int deletefile(std::string src);

sFilePathParts getfilepathparts(const std::string& path);
std::string extractfiledirectory(const std::string& pathname);
std::string extractfilename(const std::string& pathname);
std::string extractfilename_noextension(const std::string& pathname);
std::string extractfileextension(const std::string& pathname);

std::string insert_before_filename(std::string input, std::string insertion);
std::string insert_after_filename(std::string input, std::string insertion);
std::string insert_after_extension(std::string input, std::string insertion);

int64_t filesize(const std::string& path);
std::vector<std::string> getfilelist(const std::string& path, const std::string& extension);
std::vector<std::string> sortfilelistbysize(std::vector<std::string>& filelist,int sortupordown);
void recursivefilelist(const std::string& path, const std::string& extension, std::vector<std::string>& list);
void recursivefilelist(const std::string& path, const std::string& extension, FILE* outfile);
std::vector<std::string> directoryheirachy(std::string dirname);

size_t countlines(const std::string filename);
size_t countlines1(const std::string filename);

class cDirectoryAccess  
{
	private:
		static std::vector<std::string> getfilelist_single(const std::string& searchpattern);
public:
    char pathseperator;	
	
	cDirectoryAccess();
	~cDirectoryAccess();

	static std::vector<std::string> getfilelist(const std::string& searchpattern);
	static bool wildcmp(std::string& wildpattern, std::string& stringpattern);

private:
	
};

#endif

