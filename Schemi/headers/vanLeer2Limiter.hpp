/*
 * vanLeer2.hpp
 *
 *  Created on: 2020/11/29
 *      Author: Maxim Boldyrev
 *
 *      van Leer ver. 2 slope limiter class.
 */

#ifndef VANLEER2LIMITER_HPP_
#define VANLEER2LIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class vanLeer2Limiter: public abstractLimiter
{
	scalar vanLeer2LimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar vanLeer2LimiterCalculation(const scalar r) const noexcept;
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

#endif /* VANLEER2LIMITER_HPP_ */
