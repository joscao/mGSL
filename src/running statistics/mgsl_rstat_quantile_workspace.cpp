#include "mgsl_rstat.h"

#include <stdexcept>
#include <algorithm>

using namespace mgsl_rstat::quantile;

double calc_psq
(
    const double qp1, const double q, const double qm1,
    const double d, const double np1, const double n, const double nm1
)
{
	double outer = d / (np1 - nm1);
	double inner_left = (n - nm1 + d) * (qp1 - q) / (np1 - n);
	double inner_right = (np1 - n - d) * (q - qm1) / (n - nm1);

	return q + outer * (inner_left + inner_right);
} // double calc_psq()

double workspace::get_min()const noexcept
{
	if (n > 5)
	{
		return *q.begin();
	}
	else
	{
		/* not yet initialized */
		auto temp_q = q; //to keep function being const
		std::nth_element(temp_q.begin(), temp_q.begin(), temp_q.begin() + n);
		return temp_q[0];
	}
}//double workspace::get_min()const noexcept

double workspace::get_max()const noexcept
{
	if (n > 5)
	{
		return *(q.end() - 1);
	}
	else
	{
		/* not yet initialized */
		auto temp_q = q; //to keep function being const
		std::nth_element(temp_q.begin(), temp_q.begin(), temp_q.begin() + n, std::greater<double>());
		return temp_q[0];
	}
}//double workspace::get_max()const noexcept

double workspace::get_p_quantile()const noexcept
{
	if (n > 5)
	{
		return q[2];
	}
	else
	{
		/* not yet initialized */
		auto temp_q = q; //to keep function being const
		std::sort(temp_q.begin(), temp_q.begin() + n);
		//gsl_stats_quantile_from_sorted_data extracted
		if (n == 0) return 0.0;
		const double index = p * (n - 1);
		const size_t lhs = static_cast<size_t>(index);
		const double delta = index - lhs;
		double result;
		if (lhs == n - 1)
		{
			result = temp_q[lhs * 1] ;
		}
		else
		{
			result = (1 - delta) * temp_q[lhs * 1] + delta * temp_q[(lhs + 1) * 1] ;
		}
		return result ;
	}
}//double workspace::get_p_quantile()const noexcept

workspace::workspace(const double p_)noexcept:
p{p_},
q{},
npos{1, 2, 3, 4, 5},
np{1.0, 1.0 + 2.0 * p_, 1.0 + 4.0 * p, 3.0 + 2.0 * p, 5.0},
dnp{0.0, 0.5 * p, p, 0.5 * (1.0 + p), 1.0},
n{0}
{}//workspace::workspace(const double p_) noexcept

void workspace::add(const double x)
{
	if (n < 5)
	{
		q[n] = x;
	}
	else
	{
		int k = -1;

		if (n == 5)
		{
			/* initialization: sort the first five heights */
			std::sort(q.begin(), q.end());
		}

		/* step B1: find k such that q_k <= x < q_{k+1} */
		if (x < q[0])
		{
			q[0] = x;
			k = 0;
		}
		else if (x >= q[4])
		{
			q[4] = x;
			k = 3;
		}
		else
		{
			for (int i = 0; i <= 3; ++i)
			{
				if (q[i] <= x && x < q[i + 1])
				{
					k = i;
					break;
				}
			}
		}

		if (k < 0)
		{
			/* we could get here if x is nan */
			throw std::invalid_argument("invalid input argument x: void workspace::add(const double x)");
		}

		/* step B2(a): update n_i */
		for (int i = k + 1; i <= 4; ++i)
			++(npos[i]);

		/* step B2(b): update n_i' */
		for (int i = 0; i < 5; ++i)
			np[i] += dnp[i];

		/* step B3: update heights */
		for (int i = 1; i <= 3; ++i)
		{
			double ni = (double) npos[i];
			double d = np[i] - ni;

			if ((d >= 1.0 && (npos[i + 1] - npos[i] > 1)) ||
			        (d <= -1.0 && (npos[i - 1] - npos[i] < -1)))
			{
				int dsign = (d > 0.0) ? 1 : -1;
				double qp1 = q[i + 1];
				double qi = q[i];
				double qm1 = q[i - 1];
				double np1 = (double) npos[i + 1];
				double nm1 = (double) npos[i - 1];
				double qp = calc_psq(qp1, qi, qm1, (double) dsign,
				                     np1, ni, nm1);

				if (qm1 < qp && qp < qp1)
					q[i] = qp;
				else
				{
					/* use linear formula */
					q[i] += dsign * (q[i + dsign] - qi) / ((double) npos[i + dsign] - ni);
				}

				npos[i] += dsign;
			}
		}
	}

	++n;

}//void workspace::add(const double x)

