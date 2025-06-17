/*
 * SchmidtmannLimiter.hpp
 *
 *  Created on: 2025/05/24
 *      Author: Maxim Boldyrev
 */

#ifndef SCHMIDTMANNLIMITER_HPP_
#define SCHMIDTMANNLIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class SchmidtmannLimiter: public abstractLimiter
{
	const scalar alpha { 1 }, beta { 2 }, gamma { 1.5 };

	scalar SchmidtmannLimiterCalculation(const scalar r,
			const scalar xiR) const noexcept;

	scalar SchmidtmannLimiterCalculation(const scalar r) const noexcept;
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

#endif /* SCHMIDTMANNLIMITER_HPP_ */
