/*
 * KorenLimiter.hpp
 *
 *  Created on: 2025/05/21
 *      Author: Maxim Boldyrev
 */

#ifndef KORENLIMITER_HPP_
#define KORENLIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class KorenLimiter: public abstractLimiter
{
	scalar KorenLimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar KorenLimiterCalculation(const scalar r) const noexcept;
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

#endif /* KORENLIMITER_HPP_ */
