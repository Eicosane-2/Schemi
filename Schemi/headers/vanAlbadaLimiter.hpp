/*
 * vanAlbadaLimiter.hpp
 *
 *  Created on: 2020/11/22
 *      Author: Maxim Boldyrev
 *
 *      van Albada slope limiter class.
 */

#ifndef VANALBADALIMITER_HPP_
#define VANALBADALIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class vanAlbadaLimiter: public abstractLimiter
{
	scalar vanAlbadaLimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar vanAlbadaLimiterCalculation(const scalar r) const noexcept;
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

#endif /* VANALBADALIMITER_HPP_ */
