// Copyright (C) [2020] [Jonathan Schmalfuß]
//
// This file is part of the modern Generic Scientific Library (mGSL).
// mGSL is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// mGSL is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License 
// along with mGSL; if not, see <https://www.gnu.org/licenses>.

#include "mgsl_rstat.h"

#include <cmath>

using namespace mgsl_rstat::rstat;

void workspace::add(const double x)noexcept
{
	/* update mean and variance */
	const double delta = x - mean;
	double d_n = static_cast<double>(++(n));
	const double delta_n = delta / d_n;
	const double delta_nsq = delta_n * delta_n;
	const double term1 = delta * delta_n * (d_n - 1.0);
	mean += delta_n;
	M4 += term1 * delta_nsq * (d_n * d_n - 3.0 * d_n + 3.0) +
	      6.0 * delta_nsq * M2 - 4.0 * delta_n * M3;
	M3 += term1 * delta_n * (d_n - 2.0) - 3.0 * delta_n * M2;
	M2 += term1;

	/* update median and min and max*/
	median_workspace.add(x);
}//void workspace::add(const double x)noexcept

size_t workspace::get_n()const noexcept
{
	return n;
}//size_t get_n()const noexcept

double workspace::get_min()const noexcept
{
	return median_workspace.get_min();
}//double workspace::get_min()const noexcept

double workspace::get_max()const noexcept
{
	return median_workspace.get_max();
}//double workspace::get_max()const noexcept

double workspace::get_mean()const noexcept
{
	return mean;
}//double workspace::get_mean()const noexcept

double workspace::get_variance()const noexcept
{
	if (n > 1)
	{
		const double d_n = static_cast<double>(n);
		return (M2 / (d_n - 1.0));
	}
	else
		return 0.0;
}//double workspace::get_variance()const noexcept

double workspace::get_sd()const noexcept
{
	return std::sqrt(get_variance());
}//double workspace::get_sd()const noexcept

double workspace::get_rms()const noexcept
{
	double rms = 0.0;
	if (n > 0)
	{
		const double sigma = get_sd();
		const double d_n = static_cast<double>(n);
		const double a = sqrt((n - 1.0) / d_n);
		rms = std::hypot(mean, a * sigma);
	}
	return rms;
}//double workspace::get_rms()const noexcept

double workspace::get_sd_mean()const noexcept
{
	if (n > 0)
	{
		return get_sd() / std::sqrt(static_cast<double>(n));
	}
	else
		return 0.0;
}//double workspace::get_sd_mean()const noexcept

double workspace::get_median()const noexcept
{
	return median_workspace.get_p_quantile();
}//double workspace::get_median()const noexcept

double workspace::get_skew()const noexcept
{
	if ( n > 0)
	{
		const double d_n = static_cast<double>(n);
		const double fac = std::pow(d_n - 1.0, 1.5) / d_n;
		return ((fac * M3) / std::pow(M2, 1.5));
	}
	else
		return 0.0;
}//double workspace::get_skew()const noexcept

double workspace::get_kurtosis()const noexcept
{
	if (n > 0)
	{
		const double d_n = static_cast<double>(n);
		const double fac = ((d_n - 1.0) / d_n) * (d_n - 1.0);
		return ((fac * M4) / (M2 * M2) - 3.0);
	}
	else
		return 0.0;
}//double workspace::get_kurtosis()const noexcept
