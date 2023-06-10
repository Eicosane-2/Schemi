/*
 * superbeeLimiter.hpp
 *
 *  Created on: 2021/02/08
 *      Author: Maxim Boldyrev
 *
 *       Superbee slope limiter class.
 */

#ifndef SUPERBEELIMITER_HPP_
#define SUPERBEELIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class superbeeLimiter: public abstractLimiter
{
	scalar superbeeLimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar superbeeLimiterCalculation(const scalar r) const noexcept;
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

#endif /* SUPERBEELIMITER_HPP_ */
