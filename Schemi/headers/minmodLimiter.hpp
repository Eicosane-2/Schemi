/*
 * minmodLimiter.hpp
 *
 *  Created on: 2019/12/29
 *      Author: Maxim Boldyrev
 *
 *      Minmod slope limiter class.
 */

#ifndef MINMODLIMITER_HPP_
#define MINMODLIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class minmodLimiter: public abstractLimiter
{
	scalar minmodLimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar minmodLimiterCalculation(const scalar r) const noexcept;
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

#endif /* MINMODLIMITER_HPP_ */
