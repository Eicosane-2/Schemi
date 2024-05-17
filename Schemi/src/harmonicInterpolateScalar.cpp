/*
 * harmonicInterpolateScalar.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "harmonicInterpolateScalar.hpp"

schemi::surfaceField<schemi::scalar> schemi::harmonicInterpolate(
		const volumeField<scalar> & inField,
		const boundaryConditionValue & bncCalc, const std::size_t compt)
{
	auto & mesh_ { inField.meshRef() };

	surfaceField<scalar> retSurfField { mesh_, 0., inField.boundCond() };

	for (std::size_t i = 0; i < retSurfField.size(); ++i)
	{
		switch (inField.boundCond()[i].first)
		{
		case boundaryConditionType::innerSurface:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };
			const std::size_t neiIndex { mesh_.surfaceNeighbour()[i] };

			const scalar reversedValue = mesh_.surfOwnW()[i]
					/ inField()[ownIndex]
					+ mesh_.surfNeiW()[i] / inField()[neiIndex];

			retSurfField.r()[i] = 1. / reversedValue;
		}
			break;
		default:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };

			retSurfField.r()[i] = bncCalc.boundaryConditionValueSurface(
					inField()[ownIndex], inField.boundCond()[i], ownIndex, i,
					compt);
		}
			break;
		}
	}
	return retSurfField;
}
