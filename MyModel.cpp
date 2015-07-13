/*
* Copyright (c) 2009, 2010, 2011, 2012 Brendon J. Brewer.
*
* This file is part of DNest3.
*
* DNest3 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DNest3 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DNest3. If not, see <http://www.gnu.org/licenses/>.
*/

#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace DNest3;

const Data& MyModel::data = Data::get_instance();

MyModel::MyModel()
//:spikes(4, 100, false, ClassicMassInf1D(data.get_t_min(), data.get_t_max(),
//				1E-3*data.get_y_mean(), 1E3*data.get_y_mean()))
:spikes(4, 100, false, GaussPrior3D(data.get_t_min() - 0.1*data.get_t_range(), data.get_t_max()))
,noise_normals(data.get_t().size())
,mu(data.get_t().size())
{

}

void MyModel::calculate_mu()
{
	const vector<double>& t_left = data.get_t_left();
	const vector<double>& t_right = data.get_t_right();

	// Update or from scratch?
	// NEVER UPDATE BECAUSE COMPONENT IS A LOG-AMPLITUDE, NOT AN AMPLITUDE!
	bool update = false;

	// Get the components
	const vector< vector<double> >& components = (update)?(spikes.get_added()):
				(spikes.get_components());

	// Set the background level
	if(!update)
		mu.assign(mu.size(), background);

	double amplitude, duration, skew, tc;
	double rise, fall;

	// First, the non-delayed signal
	for(size_t j=0; j<components.size(); j++)
	{
		tc = components[j][0];
		amplitude = exp(components[j][1]);
		duration = exp(components[j][2]);
		skew = exp(components[j][3]);

		rise = duration/(1. + skew);
		fall = rise*skew;

		for(size_t i=0; i<mu.size(); i++)
		{
			if(tc <= t_left[i])
			{
				// Bin to the right of peak
				mu[i] += amplitude*fall/data.get_dt()*
						(exp((tc - t_left[i])/fall) -
						 exp((tc - t_right[i])/fall));
			}
			else if(tc >= t_right[i])
			{
				// Bin to the left of peak
				mu[i] += -amplitude*rise/data.get_dt()*
						(exp((t_left[i] - tc)/rise) -
						 exp((t_right[i] - tc)/rise));
			}
			else
			{
				// Part to the left
				mu[i] += -amplitude*rise/data.get_dt()*
						(exp((t_left[i] - tc)/rise) -
						 1.);

				// Part to the right
				mu[i] += amplitude*fall/data.get_dt()*
						(1. -
						 exp((tc - t_right[i])/fall));
			}
		}
	}

	// Now, the delayed signal
	for(size_t j=0; j<components.size(); j++)
	{
		tc = components[j][0] + time_delay;
		amplitude = mag_ratio*exp(components[j][1]);
		duration = exp(components[j][2]);
		skew = exp(components[j][3]);

		rise = duration/(1. + skew);
		fall = rise*skew;

		for(size_t i=0; i<mu.size(); i++)
		{
			if(tc <= t_left[i])
			{
				// Bin to the right of peak
				mu[i] += amplitude*fall/data.get_dt()*
						(exp((tc - t_left[i])/fall) -
						 exp((tc - t_right[i])/fall));
			}
			else if(tc >= t_right[i])
			{
				// Bin to the left of peak
				mu[i] += -amplitude*rise/data.get_dt()*
						(exp((t_left[i] - tc)/rise) -
						 exp((t_right[i] - tc)/rise));
			}
			else
			{
				// Part to the left
				mu[i] += -amplitude*rise/data.get_dt()*
						(exp((t_left[i] - tc)/rise) -
						 1.);

				// Part to the right
				mu[i] += amplitude*fall/data.get_dt()*
						(1. -
						 exp((tc - t_right[i])/fall));
			}
		}
	}



//	// Compute the OU process
//	vector<double> y(mu.size());
//	double alpha = exp(-1./noise_L);
//	// y[i+1] = a*y[i] + sigma*n[i]
//	// S^2 = a^2 S^2 + sigma^2
//	// S^2 = sigma^2/(1 - a^2)
//	// S = sigma/sqrt(1 - a^2)

//	for(size_t i=0; i<mu.size(); i++)
//	{
//		if(i==0)
//			y[i] = noise_sigma/sqrt(1. - alpha*alpha)*noise_normals[i];
//		else
//			y[i] = alpha*y[i-1] + noise_sigma*noise_normals[i];
//		mu[i] *= exp(y[i]);
//	}
}

void MyModel::fromPrior()
{
	background = tan(M_PI*(0.97*randomU() - 0.485));
	background = exp(background);

	time_delay = exp(log(1E-3*data.get_t_range()) + log(1E3)*randomU());
	mag_ratio = exp(3.*randn());

	spikes.fromPrior();

//	noise_sigma = exp(log(1E-3) + log(1E3)*randomU());
//	noise_L = exp(log(1E-2*data.get_t_range())
//			+ log(1E3)*randomU());
	calculate_mu();
}

double MyModel::perturb()
{
	double logH = 0.;

	int which;

	// Single parameters
	if(randomU() <= 0.5)
	{
		which = randInt(3);
		if(which == 0)
		{
			background = log(background);
			background = (atan(background)/M_PI + 0.485)/0.97;
			background += pow(10., 1.5 - 6.*randomU())*randn();
			background = mod(background, 1.);
			background = tan(M_PI*(0.97*background - 0.485));
			background = exp(background);
		}
		else if(which == 1)
		{
			time_delay = log(time_delay);
			time_delay += log(1E3)*randh();
			wrap(time_delay,
				log(1E-3*data.get_t_range()),
				log(data.get_t_range()));
			time_delay = exp(time_delay);
		}
		else if(which == 2)
		{
			mag_ratio = log(mag_ratio);
			logH -= -0.5*pow(mag_ratio/3., 2);
			mag_ratio += 3.*randh();
			logH += -0.5*pow(mag_ratio/3., 2);
			mag_ratio = exp(mag_ratio);
		}
//		else if(which == 3)
//		{
//			noise_sigma = log(noise_sigma);
//			noise_sigma += log(1E3)*randh();
//			wrap(noise_sigma, log(1E-3), log(1.));
//			noise_sigma = exp(noise_sigma);
//		}
//		else
//		{
//			noise_L = log(noise_L);
//			noise_L += log(1E3)*randh();
//			wrap(noise_L, log(1E-2*data.get_t_range()), log(10.*data.get_t_range()));
//			noise_L = exp(noise_L);
//		}
	}
	else
	{
		// Compound parameters
		which = 0;//randInt(2);
		if(which == 0)
			logH += spikes.perturb();
//		else
//		{
//			int num = exp(log((double)noise_normals.size())*randomU());
//			for(int i=0; i<num; i++)
//			{
//				int k = randInt(noise_normals.size());
//				noise_normals[k] = randn();
//			}
//		}
	}

	calculate_mu();

	return logH;
}

double MyModel::logLikelihood() const
{
        const vector<double>& t = data.get_t();
	const vector<int>& y = data.get_y();

	double logl = 0.;
	for(size_t i=0; i<t.size(); i++)
		logl += -mu[i] + y[i]*log(mu[i]) - gsl_sf_lngamma(y[i] + 1.);
 
	return logl;
}

void MyModel::print(std::ostream& out) const
{
	out<<background<<' '<<time_delay<<' '<<mag_ratio<<' ';
	spikes.print(out);

	for(size_t i=0; i<mu.size(); i++)
		out<<mu[i]<<' ';
}

string MyModel::description() const
{
	return string("");
}

