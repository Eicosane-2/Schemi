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

#include <memory>
#include <string>

#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"

namespace schemi
{
class abstractLimiter
{
public:
	virtual ~abstractLimiter() noexcept =0;

	static std::unique_ptr<abstractLimiter> createLimiter(
			const std::string & name);

	virtual vector calculate(const vector& /*r*/,
			const vector& /*gradient*/) const noexcept =0;

	virtual tensor calculate(const tensor& /*r*/,
			const tensor& /*gradient*/) const noexcept =0;

	virtual tensor3 calculate(const tensor3& /*r*/,
			const tensor3& /*gradient*/) const noexcept =0;

	virtual vector calculateNoRSLimit(const vector& /*r*/,
			const vector& /*gradient*/) const noexcept =0;

	virtual tensor calculateNoRSLimit(const tensor& /*r*/,
			const tensor& /*gradient*/) const noexcept =0;

	virtual tensor3 calculateNoRSLimit(const tensor3& /*r*/,
			const tensor3& /*gradient*/) const noexcept =0;
};
}  // namespace schemi

#endif /* ABSTRACTTVDLIMITER_HPP_ */
