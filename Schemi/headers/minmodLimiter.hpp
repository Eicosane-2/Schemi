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
	vector calculate(const vector & r, const vector & gradient) const noexcept
			override;

	tensor calculate(const tensor & r, const tensor & gradient) const noexcept
			override;

	tensor3 calculate(const tensor3 & r,
			const tensor3 & gradient) const noexcept override;

	vector calculate3OLimit(const vector & r,
			const vector & gradient) const noexcept override;

	tensor calculate3OLimit(const tensor & r,
			const tensor & gradient) const noexcept override;

	tensor3 calculate3OLimit(const tensor3 & r,
			const tensor3 & gradient) const noexcept override;
};
}  // namespace schemi

#endif /* MINMODLIMITER_HPP_ */
