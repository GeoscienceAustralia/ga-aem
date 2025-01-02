#pragma once
#include "fftw3.h"

namespace AEM {

	class FFTWPlanWrapper {

	private:
		fftw_plan Plan = nullptr;

	public:

		// Default constructor
		FFTWPlanWrapper() {
			Plan = nullptr;
		}

		FFTWPlanWrapper(const fftw_plan plan) {
			setplan(plan);
		}

		// Move constructor
		FFTWPlanWrapper(FFTWPlanWrapper&& other) noexcept
			: Plan(other.Plan)
		{
			other.Plan = nullptr;
		};

		// Copy assignment operator
		FFTWPlanWrapper& operator=(FFTWPlanWrapper& other) noexcept {
			setplan(other.Plan);
			other.Plan = nullptr;
			return *this;
		};

		~FFTWPlanWrapper() {
			destroy();
		};

		void setplan(const fftw_plan plan) {
			Plan = plan;
		}

		void destroy() const {
			if (Plan) {
				fftw_destroy_plan(Plan);
			}
		}

		void execute() const {
			fftw_execute(Plan);
		}

		void print() const {
			fftw_print_plan(Plan);
		}

	};
};

