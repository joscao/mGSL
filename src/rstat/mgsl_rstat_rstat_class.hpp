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

#ifndef MGSL_RSTAT_RSTAT_CLASS_HPP
#define MGSL_RSTAT_RSTAT_CLASS_HPP
#pragma once

#include <cmath>

#include "mgsl_rstat_quantile_class.hpp"

namespace mgsl {
namespace rstat {

class rstat {
public:
	void add(const double) noexcept;

	size_t get_n()const noexcept;
	double get_min()const noexcept;
	double get_max()const noexcept;
	double get_mean()const noexcept;
	double get_variance()const noexcept;
	double get_sd()const noexcept;
	double get_rms()const noexcept;
	double get_sd_mean()const noexcept;
	double get_median()const noexcept;
	double get_skew()const noexcept;
	double get_kurtosis()const noexcept;
private:
	double mean{0.0};   /* current mean */
	double M2{0.0};     /* M_k = sum_{i=1..n} [ x_i - mean_n ]^k */
	double M3{0.0};
	double M4{0.0};
	size_t n{0};        /* number of data points added */
	quantile median_workspace{0.5};
};

//inline member functions


inline void rstat::add(const double x)noexcept
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
}//void rstat::add(const double x)noexcept

inline size_t rstat::get_n()const noexcept
{
	return n;
}//size_t rstat::get_n()const noexcept

inline double rstat::get_min()const noexcept
{
	return median_workspace.get_min();
}//double rstat::get_min()const noexcept

inline double rstat::get_max()const noexcept
{
	return median_workspace.get_max();
}//double rstat::get_max()const noexcept

inline double rstat::get_mean()const noexcept
{
	return mean;
}//double rstat::get_mean()const noexcept

inline double rstat::get_variance()const noexcept
{
	if (n > 1)
	{
		const double d_n = static_cast<double>(n);
		return (M2 / (d_n - 1.0));
	}
	else
		return 0.0;
}//double rstat::get_variance()const noexcept

inline double rstat::get_sd()const noexcept
{
	return std::sqrt(get_variance());
}//double rstat::get_sd()const noexcept

inline double rstat::get_rms()const noexcept
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
}//double rstat::get_rms()const noexcept

inline double rstat::get_sd_mean()const noexcept
{
	if (n > 0)
	{
		return get_sd() / std::sqrt(static_cast<double>(n));
	}
	else
		return 0.0;
}//double rstat::get_sd_mean()const noexcept

inline double rstat::get_median()const noexcept
{
	return median_workspace.get_p_quantile();
}//double rstat::get_median()const noexcept

inline double rstat::get_skew()const noexcept
{
	if ( n > 0)
	{
		const double d_n = static_cast<double>(n);
		const double fac = std::pow(d_n - 1.0, 1.5) / d_n;
		return ((fac * M3) / std::pow(M2, 1.5));
	}
	else
		return 0.0;
}//double rstat::get_skew()const noexcept

inline double rstat::get_kurtosis()const noexcept
{
	if (n > 0)
	{
		const double d_n = static_cast<double>(n);
		const double fac = ((d_n - 1.0) / d_n) * (d_n - 1.0);
		return ((fac * M4) / (M2 * M2) - 3.0);
	}
	else
		return 0.0;
}//double rstat::get_kurtosis()const noexcept



}// END NAMEPSACE rstat
}//END NAMESPACE mgsl

#endif //MGSL_RSTAT_RSTAT_CLASS_HPP