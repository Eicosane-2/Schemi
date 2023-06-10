/*
 * minmod2Limiter.hpp
 *
 *  Created on: 2021/11/25
 *      Author: Maxim Boldyrev
 *
 *      Minmod ver. 2 slope limiter class.
 */

#ifndef MINMOD2LIMITER_HPP_
#define MINMOD2LIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class minmod2Limiter: public abstractLimiter
{
	scalar minmod2LimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar minmod2LimiterCalculation(const scalar r) const noexcept;
public:
	vector calculate(const vector & r, const vector & gradientC) const noexcept
			override;

	tensor calculate(const tensor & r, const tensor & gradientC) const noexcept
			override;

	tensor3 calculate(const tensor3 & r,
			const tensor3 & gradientC) const noexcept override;

	vector calculateNoRightLimit(const vector & r,
			const vector & gradientC) const noexcept override;

	tensor calculateNoRightLimit(const tensor & r,
			const tensor & gradientC) const noexcept override;

	tensor3 calculateNoRightLimit(const tensor3 & r,
			const tensor3 & gradientC) const noexcept override;
};
}  // namespace schemi

#endif /* MINMOD2LIMITER_HPP_ */
