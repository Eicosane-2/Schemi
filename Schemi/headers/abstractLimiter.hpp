/*
 * abstractLimiter.hpp
 *
 *  Created on: 2020/02/18
 *      Author: Maxim Boldyrev
 *
 *      Interface class for TVD limiters.
 */

#ifndef ABSTRACTTVDLIMITER_HPP_
#define ABSTRACTTVDLIMITER_HPP_

#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"

namespace schemi
{
class abstractLimiter
{
public:
	virtual ~abstractLimiter() noexcept =0;

	virtual vector calculate(const vector& /*r*/,
			const vector& /*gradientC*/) const noexcept =0;

	virtual tensor calculate(const tensor& /*r*/,
			const tensor& /*gradientC*/) const noexcept =0;

	virtual tensor3 calculate(const tensor3& /*r*/,
			const tensor3& /*gradientC*/) const noexcept =0;

	virtual vector calculateNoRightLimit(const vector& /*r*/,
			const vector& /*gradientC*/) const noexcept =0;

	virtual tensor calculateNoRightLimit(const tensor& /*r*/,
			const tensor& /*gradientC*/) const noexcept =0;

	virtual tensor3 calculateNoRightLimit(const tensor3& /*r*/,
			const tensor3& /*gradientC*/) const noexcept =0;
};
}  // namespace schemi

#endif /* ABSTRACTTVDLIMITER_HPP_ */
