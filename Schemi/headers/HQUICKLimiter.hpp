/*
 * HQUICKLimiter.hpp
 *
 *  Created on: 2020/11/25
 *      Author: Maxim Boldyrev
 *
 *      HQUICK slope limiter class (not TVD).
 */

#ifndef HQUICKLIMITER_HPP_
#define HQUICKLIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class HQUICKLimiter: public abstractLimiter
{
	scalar HQUICKLimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar HQUICKLimiterCalculation(const scalar r) const noexcept;
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

#endif /* HQUICKLIMITER_HPP_ */
