/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _csamplebunch_H
#define _csamplebunch_H

#include <vector>
#include <algorithm>

class cSampleBunch {

//public:
	std::vector<size_t> indices;
	size_t master;//index of the "central" sample (not necessarily middle of array
public:
	cSampleBunch() {};
	cSampleBunch(const std::vector<size_t>& _indices, const size_t& masterindex) {
		indices = _indices;
		for (size_t i = 0; i < indices.size(); i++) {
			if (indices[i] == masterindex) {
				master = i;
				return;
			}
		}
		glog.errormsg("Error in cSampleBunch()");
		return;
	}

	const size_t& master_index() const {
		return master;
	}

	const size_t& master_record() const {
		return indices[master];
	}

	const size_t& record(const size_t& bunchindex) const {
		return indices[bunchindex];
	}

	size_t size() const {
		return indices.size();
	}
};

#endif
