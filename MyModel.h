#ifndef _MyModel_
#define _MyModel_

#include "Model.h"
#include "Data.h"
#include "RJObject.h"
#include "GaussPrior3D.h"
#include <vector>

class MyModel:public DNest3::Model
{
	private:
		// Reference to the data
		static const Data& data;

		// A flat background level
		double background;

		// The spikes
		RJObject<GaussPrior3D> spikes;

		// Time delay
		double time_delay;

		// Magnification ratio
		double mag_ratio;

		// Extra white noise on the poisson rate
		std::vector<double> noise_normals;
		double noise_sigma, noise_L;

		// Poisson mean
		std::vector<long double> mu;

		// Calculates mu from scratch
		void calculate_mu();

	public:
		MyModel();

		// Generate the point from the prior
		void fromPrior();

		// Metropolis-Hastings proposals
		double perturb();

		// Likelihood function
		double logLikelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
};

#endif

