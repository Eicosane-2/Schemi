/*
 * vanLeerLimiter.hpp
 *
 *  Created on: 2019/12/15
 *      Author: Maxim Boldyrev
 *
 *      van Leer slope limiter class.
 */

#ifndef VANLEERLIMITER_HPP_
#define VANLEERLIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class vanLeerLimiter: public abstractLimiter
{
	scalar vanLeerLimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar vanLeerLimiterCalculation(const scalar r) const noexcept;
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

#endif /* VANLEERLIMITER_HPP_ */
