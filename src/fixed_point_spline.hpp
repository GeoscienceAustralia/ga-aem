#pragma once

// Fixed point spline for repetitive re-evaluation at fixed interpolation points
template<typename T>
class FixedPointSpline {

private:
	size_t nNodes = 0;
	size_t nInterp = 0;
	std::vector<T> xnodes;
	std::vector<T> xinterp;

	std::vector<size_t> klo; //Work vectors for particular set of nodes and interpolation points
	std::vector<size_t> khi;
	std::vector<T> h2; 
	std::vector<T> a;
	std::vector<T> b;
	std::vector<T> a3ma;
	std::vector<T> b3mb;

	std::vector<T> u; // Work vector for a particular set of yvals
	std::vector<T> coefficients; // Coefficients for a particular set of yvals
	std::vector<T> interpolation_results; // Interpolation result for a particular set of yvals

	void resize(const size_t numnodes, const size_t numinterp) {
		nNodes = numnodes;
		nInterp = numinterp;

		h2.resize(numinterp);
		a.resize(numinterp);
		b.resize(numinterp);
		klo.resize(numinterp);
		khi.resize(numinterp);
		a3ma.resize(numinterp);
		b3mb.resize(numinterp);
		u.resize(nNodes-1);
		coefficients.resize(nNodes);
		interpolation_results.resize(numinterp);
	};

	void compute_coefficients(const std::vector<T>& yvals, T yp1 = FixedPointSpline<T>::UNDEFINED, T ypn = FixedPointSpline<T>::UNDEFINED) {
		const size_t& n = nNodes;
		if (yp1 == FixedPointSpline<T>::UNDEFINED)
			coefficients[0] = u[0] = 0.0;
		else {
			coefficients[0] = -0.5;
			u[0] = (3.0 / (xnodes[1] - xnodes[0])) * ((yvals[1] - yvals[0]) / (xnodes[1] - xnodes[0]) - yp1);
		}

		for (size_t i = 1; i < n - 1; i++) {
			const T sig = (xnodes[i] - xnodes[i - 1]) / (xnodes[i + 1] - xnodes[i - 1]);
			const T p = sig * coefficients[i - 1] + 2.0;
			coefficients[i] = (sig - 1.0) / p;
			u[i] = (yvals[i + 1] - yvals[i]) / (xnodes[i + 1] - xnodes[i]) - (yvals[i] - yvals[i - 1]) / (xnodes[i] - xnodes[i - 1]);
			u[i] = (6.0 * u[i] / (xnodes[i + 1] - xnodes[i - 1]) - sig * u[i - 1]) / p;
		}

		T qn;
		T un;
		if (ypn == FixedPointSpline<T>::UNDEFINED)
			qn = un = 0.0;
		else {
			qn = 0.5;
			un = (3.0 / (xnodes[n - 1] - xnodes[n - 2])) * (ypn - (yvals[n - 1] - yvals[n - 2]) / (xnodes[n - 1] - xnodes[n - 2]));
		}
		coefficients[n - 1] = (un - qn * u[n - 2]) / (qn * coefficients[n - 2] + 1.0);

		for (size_t k = n - 1; k-- > 0;) {
			//Note the unusual syntax for decrement of unsigned variable
			coefficients[k] = coefficients[k] * coefficients[k + 1] + u[k];
		}
	};

public:

	inline static const T UNDEFINED = std::numeric_limits<T>::min();

	FixedPointSpline() {}

	void initialise(const std::vector<T>& _xnodes, const std::vector<T>& _xinterp) {
		xnodes = _xnodes;
		xinterp = _xinterp;
		resize(xnodes.size(), xinterp.size());
		for (size_t i = 0; i < nInterp; i++) {
			klo[i] = 0;
			khi[i] = nNodes - 1;
			while ((khi[i] - klo[i]) > 1) {
				size_t k = (klo[i] + khi[i]) >> 1;
				if (xnodes[k] > xinterp[i]) khi[i] = k;
				else klo[i] = k;
			}
			const T h = (xnodes[khi[i]] - xnodes[klo[i]]);
			if (h == 0.0) glog.errormsg(_SRC_, "Bad xi input to routine\n");
			a[i] = (xnodes[khi[i]] - xinterp[i]) / h;
			b[i] = (xinterp[i] - xnodes[klo[i]]) / h;
			h2[i] = h * h;
			a3ma[i] = a[i] * a[i] * a[i] - a[i];
			b3mb[i] = b[i] * b[i] * b[i] - b[i];
		}
	}

	void compute_interpolation(const std::vector<T>& yvals) {
		compute_coefficients(yvals);
		for (size_t k = 0; k < nInterp; k++) {
			const T scale = (h2[k]) / 6.0;
			interpolation_results[k] = (a[k] * yvals[klo[k]] + b[k] * yvals[khi[k]] + scale * (a3ma[k] * coefficients[klo[k]] + b3mb[k] * coefficients[khi[k]]));
		}
	}

	const std::vector<T>& interpolated_values() const { return interpolation_results; }
};
