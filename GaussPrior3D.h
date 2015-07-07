#ifndef _GaussPrior3D_
#define _GaussPrior3D_

#include <vector>
#include "Distributions/Distribution.h"

class GaussPrior3D:public Distribution
{
	private:
		// Data limits
		double t_min, t_max, t_range;

		// The means
		double mean_logA, mean_logDuration, mean_logSkew;

		// The coefficients
		double co_ADuration, co_ASkew, co_durationSkew;

		// The standard deviations
		// (conditional! for the latter two anyway)
		double sig_logA, sig_logDuration, sig_logSkew;

		double perturb_parameters();

	public:
		GaussPrior3D(double t_min, double t_max);

		void fromPrior();

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;

		static const int weight_parameter = 1;
};

#endif

