#pragma once

#include <string>
#include <iostream>
#include "logger.hpp"

namespace AEM {
	class CalculationType {

	public: enum class Mode { FM, DC, DT, DH, DR, DX, DY, DZ, DTX_HEIGHT, DTX_ROLL, DTX_PITCH, DTX_YAW, DRX_ROLL, DRX_PITCH, DRX_YAW, NONE };
	private: inline static std::vector<std::string> ModeNames{ "FM", "DC", "DT", "DH", "DR", "DX", "DY", "DZ", "NONE" };
	private: inline static std::vector<std::string> ModeDescriptions{ "Forward Model", "Layer conductivity derivative", "Layer thickness derivative", "Tx Height derivative", "Tx-Rx radial distance derivative", "Tx-Rx Dx horizontal distance derivative", "Tx-Rx Dy horizontal distance derivative", "Tx-Rx Dz horizontal distance derivative", "No calculation" };

	public:
		CalculationType() {};

		CalculationType(Mode _mode, size_t _layer) : mode(_mode), layer(_layer) {
			if (layer > 2000) {//Sanity check on layer number 
				glog.errormsg(_SRC_, "Sorry but %zu is a ridiculous derivative layer number.\n", layer);
			}
			if (mode == Mode::DC || mode == Mode::DT) {
				int dummy = 0;
			}
		};

		CalculationType(Mode _mode) : mode(_mode) {
			if (mode == Mode::DC || mode == Mode::DT) {
				glog.errormsg(_SRC_, "Must set a layer number.\n");
			}
			layer = std::numeric_limits<size_t>::max();
		};

		const Mode& get_mode() const {
			return mode;
		};

		const size_t& get_layer() const {
			return layer;
		};

		std::string mode_name() const {
			size_t index = get_index_from_mode(mode);
			return ModeNames[index];
		}

		std::string string() const {
			std::ostringstream oss;
			oss << mode_name();
			if (mode == Mode::DC || mode == Mode::DT) {
				oss << "(Layer: " << layer << ")";
			}
			return oss.str();
		};

		static Mode lookup_mode(size_t index) {
			if (index > last_index()) {
				glog.errormsg(_SRC_, "Bad index (%zu) for calculation mode lookup.\n%s", index, possible_values_message().c_str());
			}
			return get_mode_from_index(index);
		};

	private:
		Mode mode = Mode::NONE;
		size_t layer = std::numeric_limits<size_t>::max();

		static Mode get_mode_from_index(const int& _index) {
			return static_cast<Mode>(_index);
		};

		static size_t get_index_from_mode(const Mode& _mode) {
			return static_cast<int>(_mode);
		};

		static size_t last_index() {
			return get_index_from_mode(Mode::NONE);
		};

		static std::string possible_values_message() {
			std::ostringstream oss;
			oss << "Possible values for calculation mode index are" << std::endl;
			for (int i = 0; i < last_index(); i++) {
				oss << i << " (" << ModeNames[i] << ") " << ModeDescriptions[i] << std::endl;
			}
			return oss.str();
		}
	};
}
