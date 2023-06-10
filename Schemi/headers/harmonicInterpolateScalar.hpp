/*
 * harmonicInterpolateScalar.hpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#ifndef HARMONICINTERPOLATESCALAR_HPP_
#define HARMONICINTERPOLATESCALAR_HPP_

#include "boundaryConditionValue.hpp"
#include "globalConstants.hpp"
#include "scalar.hpp"
#include "surfaceField.hpp"
#include "volumeField.hpp"

namespace schemi
{
surfaceField<scalar> harmonicInterpolate(const volumeField<scalar> & inField,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder);
}  // namespace schemi

#endif /* HARMONICINTERPOLATESCALAR_HPP_ */
