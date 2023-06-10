/*
 * vanAlbada2Limiter.hpp
 *
 *  Created on: 2021/11/24
 *      Author: Maxim Boldyrev
 *
 *      van Albada ver. 2 slope limiter class.
 */

#ifndef VANALBADA2LIMITER_HPP_
#define VANALBADA2LIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class vanAlbada2Limiter: public abstractLimiter
{
	scalar vanAlbada2LimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar vanAlbada2LimiterCalculation(const scalar r) const noexcept;
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

#endif /* VANALBADA2LIMITER_HPP_ */
