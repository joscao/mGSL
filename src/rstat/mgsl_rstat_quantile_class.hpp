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

#pragma once

#include <array>

namespace mgsl{
namespace rstat {

class quantile {
public:
	quantile(const double) noexcept;

	void add(const double);

	double get_min()const noexcept;
	double get_max()const noexcept;
	double get_p_quantile()const noexcept;
private:
	double p;       			/* p-quantile */
	std::array<double, 5> q;    /* heights q_i */
	std::array<int, 5> npos;    /* positions n_i */
	std::array<double, 5> np;   /* desired positions n_i' */
	std::array<double, 5> dnp;  /* increments dn_i' */
	std::size_t n;        		/* number of data added */
};


}// END NAMEPSACE rstat
}//END NAMESPACE mgsl

