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
	auto & mesh { inField.meshRef() };

	surfaceField<scalar> retSurfField { mesh, 0., inField.boundCond() };

	for (std::size_t i = 0; i < retSurfField.size(); ++i)
	{
		switch (inField.boundCond()[i].first)
		{
		case boundaryConditionType::innerSurface:
		{
			const std::size_t ownIndex { mesh.surfaceOwner()[i] };
			const std::size_t neiIndex { mesh.surfaceNeighbour()[i] };
			const scalar surfOwnR { (mesh.surfaces()[i].rC()
					- mesh.cells()[ownIndex].rC()).mag() };
			const scalar surfNeiR { (mesh.surfaces()[i].rC()
					- mesh.cells()[neiIndex].rC()).mag() };

			const scalar reversedValue = (surfNeiR / inField.ref()[ownIndex]
					+ surfOwnR / inField.ref()[neiIndex])
					/ (surfOwnR + surfNeiR);

			retSurfField.ref_r()[i] = 1. / reversedValue;
		}
			break;
		default:
		{
			const std::size_t ownIndex { mesh.surfaceOwner()[i] };

			retSurfField.ref_r()[i] = bncCalc.boundaryConditionValueSurface(
					inField.ref()[ownIndex], inField.boundCond()[i], ownIndex,
					i, compt);
		}
			break;
		}
	}
	return retSurfField;
}
