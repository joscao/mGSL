#pragma once

#include <cstddef>
#include <array>

namespace mgsl {

namespace simple_rstat {
class workspace {
public:
	void add(const double) noexcept;

	size_t get_n()const noexcept;
	double get_mean()const noexcept;
	double get_variance()const noexcept;
	double get_sd()const noexcept;
	double get_rms()const noexcept;
	double get_sd_mean()const noexcept;
	double get_skew()const noexcept;
	double get_kurtosis()const noexcept;
private:
	double mean{0.0};     /* current mean */
	double M2{0.0};       /* M_k = sum_{i=1..n} [ x_i - mean_n ]^k */
	double M3{0.0};
	double M4{0.0};
	size_t n{0};        /* number of data points added */
};
}// END NAMEPSACE simple_rstat
}//END NAMESPACE msl