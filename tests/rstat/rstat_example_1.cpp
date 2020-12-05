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

#include <iostream>
#include <array>
#include <algorithm>
#include "../../src/rstat/mgsl_rstat.h"

int main()
{
	std::array<double, 5> data{17.2, 18.1, 16.5, 18.3, 12.6};

	mgsl::rstat::rstat accumulator;

	//add elements to running statistics accumulator: 2 liner
	auto add_elem = [&accumulator](const double & elem) {accumulator.add(elem);};
	std::for_each(data.begin(), data.end(), add_elem);

	//alternative more classical way: 3 liner
	//for (const auto & point:data){
	//	accumulator.add(point);
	//}

	std::cout << "The dataset is [" << data[0] << ", " << data[1] << ", " << data[2] << ", " << data[3] << ", " << data[4] << "]\n";
	//The dataset is [17.2, 18.1, 16.5, 18.3, 12.6]
	std::cout << "The mean is " << accumulator.get_mean() << "\n";
	//The mean is 16.54
	std::cout << "The variance is " << accumulator.get_variance() << "\n";
	//The variance is 5.373
	std::cout << "The largest value is " << accumulator.get_max() << "\n";
	//The largest value is 18.3
	std::cout << "The smallest value is " << accumulator.get_min() << "\n";
	//The smallest value is 12.6
	std::cout << "The median is " << accumulator.get_median() << "\n";
	//The median is 17.2
	std::cout << "The standard deviation is " << accumulator.get_sd() << "\n";
	//The standard deviation is 2.31797
	std::cout << "The standard deviation of the mean is " << accumulator.get_sd_mean() << "\n";
	//The standard deviation of the mean is 1.03663
	std::cout << "The skew is " << accumulator.get_skew() << "\n";
	//The skew is -0.829058
	std::cout << "The root mean square is " << accumulator.get_rms() << "\n";
	//The root mean square is 16.6694
	std::cout << "The kurtosis is " << accumulator.get_kurtosis() << "\n";
	//The kurtosis is -1.2217
	std::cout << "The number of items in the accumulator is " << accumulator.get_n() << "\n";
	//The number of items in the accumulator is 5

	//resetting accumulator through copy assignment constructor
	accumulator = mgsl::rstat::rstat();
	std::cout << "The number of items in the accumulator is " << accumulator.get_n() << "\n";
	//The number of items in the accumulator is 0
}