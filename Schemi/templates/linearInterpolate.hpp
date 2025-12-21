/*
 * linearInterpolate.hpp
 *
 *  Created on: 2019/11/23
 *      Author: Maxim Boldyrev
 *
 *      Linear interpolation functions.
 */

#ifndef LINEARINTERPOLATE_HPP_
#define LINEARINTERPOLATE_HPP_

#include "returnTypeDivergence.hpp"
#include "returnTypeGradient.hpp"
#include "volumeField.hpp"

namespace schemi
{
template<typename typeOfValue>
volumeField<typeOfValue> linearInterpolate(
		const surfaceField<typeOfValue> & inField) noexcept
{
	auto & mesh_ { inField.meshRef() };

	volumeField<typeOfValue> retVolField { mesh_, typeOfValue { 0 } };

	for (std::size_t i = 0; i < retVolField.size(); ++i)
	{
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfIndex { mesh_.surfacesOfCells()[i][j] };

			retVolField.val()[i] += inField.cval()[surfIndex]
					* (1
							- mesh_.cellSurfaceDistances()[i].second[j]
									/ mesh_.cellSurfaceDistances()[i].first);
		}
	}

	return retVolField;
}

template<typename typeOfValue>
surfaceField<typeOfValue> linearInterpolate(
		const volumeField<typeOfValue> & inField,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder)
{
	auto & mesh_ { inField.meshRef() };

	surfaceField<typeOfValue> retSurfField { mesh_, typeOfValue { 0 },
			inField.boundCond() };

	for (std::size_t i = 0; i < retSurfField.size(); ++i)
	{
		switch (inField.boundCond()[i].first)
		{
		case boundaryConditionType::innerSurface:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };
			const std::size_t neiIndex { mesh_.surfaceNeighbour()[i] };

			retSurfField.val()[i] = inField.cval()[ownIndex]
					* mesh_.surfOwnW()[i]
					+ inField.cval()[neiIndex] * mesh_.surfNeiW()[i];
		}
			break;
		[[unlikely]] default:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };

			retSurfField.val()[i] = bncCalc.boundaryConditionValueSurface(
					inField.cval()[ownIndex], inField.boundCond()[i], ownIndex,
					i, compt);
		}
			break;
		}
	}
	return retSurfField;
}
} // namespace schemi

#endif /* LINEARINTERPOLATE_HPP_ */
