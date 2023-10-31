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
			const scalar surfOwnR { (mesh_.surfaces()[i].rC()
					- mesh_.cells()[ownIndex].rC()).mag() };
			const scalar surfNeiR { (mesh_.surfaces()[i].rC()
					- mesh_.cells()[neiIndex].rC()).mag() };

			const scalar reversedValue = (surfNeiR / inField()[ownIndex]
					+ surfOwnR / inField()[neiIndex]) / (surfOwnR + surfNeiR);

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
