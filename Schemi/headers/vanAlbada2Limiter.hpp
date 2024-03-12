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

#endif /* VANALBADA2LIMITER_HPP_ */
