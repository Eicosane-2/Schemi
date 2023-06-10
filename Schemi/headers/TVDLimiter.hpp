/*
 * TVDLimiter.hpp
 *
 *  Created on: 2019/12/15
 *      Author: Maxim Boldyrev
 *
 *      Functions for TVD-limited gradients calculation.
 */

#ifndef TVDLIMITER_HPP_
#define TVDLIMITER_HPP_

#include "abstractLimiter.hpp"
#include "boundaryConditionValue.hpp"
#include "vector.hpp"
#include "volumeField.hpp"

namespace schemi
{
volumeField<vector> TVDLimiter(const volumeField<vector> & gradient,
		const volumeField<scalar> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder);

volumeField<tensor> TVDLimiter(const volumeField<tensor> & gradient,
		const volumeField<vector> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder);

volumeField<tensor3> TVDLimiter(const volumeField<tensor3> & gradient,
		const volumeField<tensor> & value,
		const abstractLimiter & limiterObjectP,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder);
}  // namespace schemi

#endif /* TVDLIMITER_HPP_ */
