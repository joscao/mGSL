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
#include "mgsl_rstat_simple_rstat_class.hpp"

namespace mgsl {
namespace rstat {

class rstat{
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
	simple_rstat accumulator{};
	quantile median{0.5};
};

//inline member functions

inline void rstat::add(const double x)noexcept
{
	//update simple_rstat accumulator
	accumulator.add(x);

	/* update median and min and max*/
	median.add(x);
}//void rstat::add(const double x)noexcept

inline size_t rstat::get_n()const noexcept
{
	return accumulator.get_n();
}//size_t rstat::get_n()const noexcept

inline double rstat::get_min()const noexcept
{
	return median.get_min();
}//double rstat::get_min()const noexcept

inline double rstat::get_max()const noexcept
{
	return median.get_max();
}//double rstat::get_max()const noexcept

inline double rstat::get_mean()const noexcept
{
	return accumulator.get_mean();
}//double rstat::get_mean()const noexcept

inline double rstat::get_variance()const noexcept
{
	return accumulator.get_variance();
}//double rstat::get_variance()const noexcept

inline double rstat::get_sd()const noexcept
{
	return std::sqrt(get_variance());
}//double rstat::get_sd()const noexcept

inline double rstat::get_rms()const noexcept
{
	return accumulator.get_rms();
}//double rstat::get_rms()const noexcept

inline double rstat::get_sd_mean()const noexcept
{
	return accumulator.get_sd_mean();
}//double rstat::get_sd_mean()const noexcept

inline double rstat::get_median()const noexcept
{
	return median.get_p_quantile();
}//double rstat::get_median()const noexcept

inline double rstat::get_skew()const noexcept
{
	return accumulator.get_skew();
}//double rstat::get_skew()const noexcept

inline double rstat::get_kurtosis()const noexcept
{
	return accumulator.get_kurtosis();
}//double rstat::get_kurtosis()const noexcept



}// END NAMEPSACE rstat
}//END NAMESPACE mgsl

#endif //MGSL_RSTAT_RSTAT_CLASS_HPP