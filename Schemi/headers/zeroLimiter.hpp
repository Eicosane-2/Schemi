/*
 * zeroLimiter.hpp
 *
 *  Created on: 2020/02/19
 *      Author: Maxim Boldyrev
 *
 *      Zero slope limiter class (piecewise interpolation).
 */

#ifndef ZEROLIMITER_HPP_
#define ZEROLIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class zeroLimiter: public abstractLimiter
{
public:
	virtual vector calculate(const vector&, const vector&) const noexcept
			override;

	virtual tensor calculate(const tensor&, const tensor&) const noexcept
			override;

	virtual tensor3 calculate(const tensor3&, const tensor3&) const noexcept
			override;

	virtual vector calculateNoRSLimit(const vector&,
			const vector&) const noexcept override;

	virtual tensor calculateNoRSLimit(const tensor&,
			const tensor&) const noexcept override;

	virtual tensor3 calculateNoRSLimit(const tensor3&,
			const tensor3&) const noexcept override;
};
}  // namespace schemi

#endif /* ZEROLIMITER_HPP_ */
