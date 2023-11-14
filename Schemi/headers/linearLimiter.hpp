/*
 * linearLimiter.hpp
 *
 *  Created on: 2020/02/19
 *      Author: Maxim Boldyrev
 *
 *      Linear slope limiter class (no limitation actually, thus not TVD).
 */

#ifndef LINEARLIMITER_HPP_
#define LINEARLIMITER_HPP_

#include "abstractLimiter.hpp"

namespace schemi
{
class linearLimiter: public abstractLimiter
{
public:
	vector calculate(const vector&, const vector & gradient) const noexcept
			override;

	tensor calculate(const tensor&, const tensor & gradient) const noexcept
			override;

	tensor3 calculate(const tensor3&, const tensor3 & gradient) const noexcept
			override;

	vector calculate3OLimit(const vector&,
			const vector & gradient) const noexcept override;

	tensor calculate3OLimit(const tensor&,
			const tensor & gradient) const noexcept override;

	tensor3 calculate3OLimit(const tensor3&,
			const tensor3 & gradient) const noexcept override;
};
}  // namespace schemi

#endif /* LINEARLIMITER_HPP_ */
