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
	vector calculate(const vector & r, const vector & gradient) const noexcept
			override;

	tensor calculate(const tensor & r, const tensor & gradient) const noexcept
			override;

	tensor3 calculate(const tensor3 & r,
			const tensor3 & gradient) const noexcept override;

	vector calculateNoRSLimit(const vector & r,
			const vector & gradient) const noexcept override;

	tensor calculateNoRSLimit(const tensor & r,
			const tensor & gradient) const noexcept override;

	tensor3 calculateNoRSLimit(const tensor3 & r,
			const tensor3 & gradient) const noexcept override;
};
}  // namespace schemi

#endif /* SUPERBEELIMITER_HPP_ */
