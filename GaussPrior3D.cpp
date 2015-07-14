#include "GaussPrior3D.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>
#include <gsl/gsl_cdf.h>

using namespace DNest3;

GaussPrior3D::GaussPrior3D(double t_min, double t_max)
:t_min(t_min)
,t_max(t_max)
,t_range(t_max - t_min)
{

}

void GaussPrior3D::fromPrior()
{
	// The three means
	mean_logA = exp(tan(M_PI*(0.97*randomU() - 0.485)));
	mean_logDuration = -10. +  20.*randomU();
	mean_logSkew = -10. + 20.*randomU();

	co_ADuration = -10. + 20.*randomU();
	co_ASkew = -10. + 20.*randomU();
	co_durationSkew = -10. + 20.*randomU();

	sig_logA = exp(randn());
	sig_logDuration = exp(randn());
	sig_logSkew = exp(randn());
}

double GaussPrior3D::perturb_parameters()
{
	double logH = 0.;

	int which = randInt(9);

	if(which == 0)
	{
		mean_logA = log(mean_logA);
		mean_logA = (atan(mean_logA)/M_PI + 0.485)/0.97;
		mean_logA += randh();
		wrap(mean_logA, 0., 1.);
		mean_logA = tan(M_PI*(0.97*mean_logA - 0.485));
		mean_logA = exp(mean_logA);
	}
	if(which == 1)
	{
		mean_logDuration += 20.*randh();
		wrap(mean_logDuration, -10., 10.);
	}
	if(which == 2)
	{
		mean_logSkew += 20.*randh();
		wrap(mean_logSkew, -10., 10.);
	}
	if(which == 3)
	{
		co_ADuration += 20.*randh();
		wrap(co_ADuration, -10., 10.);
	}
	if(which == 4)
	{
		co_ASkew += 20.*randh();
		wrap(co_ASkew, -10., 10.);
	}
	if(which == 5)
	{
		co_durationSkew += 20.*randh();
		wrap(co_durationSkew, -10., 10.);
	}
	if(which == 6)
	{
		sig_logA = log(sig_logA);
		logH -= -0.5*pow(sig_logA, 2);
		sig_logA += randh();
		logH += -0.5*pow(sig_logA, 2);
		sig_logA = exp(sig_logA);
	}
	if(which == 7)
	{
		sig_logDuration = log(sig_logDuration);
		logH -= -0.5*pow(sig_logDuration, 2);
		sig_logDuration += randh();
		logH += -0.5*pow(sig_logDuration, 2);
		sig_logDuration = exp(sig_logDuration);
	}
	if(which == 8)
	{
		sig_logSkew = log(sig_logSkew);
		logH -= -0.5*pow(sig_logSkew, 2);
		sig_logSkew += randh();
		logH += -0.5*pow(sig_logSkew, 2);
		sig_logSkew = exp(sig_logSkew);
	}

	return logH;
}

/*
* vec[0] = spike time
* vec[1] = log(amplitude)
* vec[2] = log(duration)
* vec[3] = log(skew)
*/

double GaussPrior3D::log_pdf(const std::vector<double>& vec) const
{
	double logP = 0.;

	if(vec[0] < t_min || vec[0] > t_max)
		return -1E250;

	double n1 = (vec[1] - mean_logA)/sig_logA;
	double n2 = (vec[2] - mean_logDuration - co_ADuration*n1)/sig_logDuration;

	logP += -log(sig_logA) - 0.5*pow((vec[1] - mean_logA)/sig_logA, 2);
	logP += -log(sig_logDuration) - 0.5*pow((vec[2] - mean_logDuration - co_ADuration*n1)/sig_logDuration, 2);
	logP += -log(sig_logSkew) - 0.5*pow((vec[3] -  mean_logSkew - co_ASkew*n1 - co_durationSkew*n2)/sig_logSkew, 2);

	return logP;
}

void GaussPrior3D::from_uniform(std::vector<double>& vec) const
{
	double n1 = gsl_cdf_ugaussian_Pinv(vec[1]);
	double n2 = gsl_cdf_ugaussian_Pinv(vec[2]);

	vec[0] = t_min + (t_max - t_min)*vec[0];
	vec[1] = mean_logA + sig_logA*n1;
	vec[2] = mean_logDuration + co_ADuration*n1 + sig_logDuration*n2;
	vec[3] = mean_logSkew + co_ASkew*n1 + co_durationSkew*n2
			+ sig_logSkew*gsl_cdf_ugaussian_Pinv(vec[3]);
}

void GaussPrior3D::to_uniform(std::vector<double>& vec) const
{
	vec[0] = (vec[0] - t_min)/(t_max - t_min);

	double n1 = (vec[1] - mean_logA)/sig_logA;
	double n2 = (vec[2] - mean_logDuration - co_ADuration*n1)/sig_logDuration;
	double n3 = (vec[3] - mean_logSkew - co_ASkew*n1 - co_durationSkew*n2)/sig_logSkew;

	vec[1] = gsl_cdf_ugaussian_P(n1);
	vec[2] = gsl_cdf_ugaussian_P(n2);
	vec[3] = gsl_cdf_ugaussian_P(n3);
}

void GaussPrior3D::print(std::ostream& out) const
{
	out<<mean_logA<<" "<<mean_logDuration<<" "<<mean_logSkew<<" ";
	out<<co_ADuration<<" "<<co_ASkew<<" "<<co_durationSkew<<" ";
	out<<sig_logA<<" "<<sig_logDuration<<" "<<sig_logSkew<<" ";
}

