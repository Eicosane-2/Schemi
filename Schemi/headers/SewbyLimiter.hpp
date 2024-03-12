/*
 * SewbyLimiter.hpp
 *
 *  Created on: 2023/11/05
 *      Author: Maxim Boldyrev
 */

#ifndef SEWBYLIMITER_HPP_
#define SEWBYLIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class SewbyLimiter: public abstractLimiter
{
	scalar SewbyLimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar SewbyLimiterCalculation(const scalar r) const noexcept;
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

#endif /* SEWBYLIMITER_HPP_ */
