/*
 * pressureStarClass.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "pressureStarClass.hpp"

#include "vectorVectorDotProduct.hpp"

schemi::scalar schemi::pressureStarClass::pressureStar(
		const abstractMixtureThermodynamics & mix,
		const std::valarray<scalar> & densityState,
		const vector & momentumState, const scalar totalEnergyState,
		const scalar rhokState) const noexcept
{
	scalar pressureOutput;

	const scalar kinEnState { 0.5 * (momentumState & momentumState)
			/ densityState[0] };

	const scalar internalEnergyState { totalEnergyState - kinEnState - rhokState };

	std::valarray<scalar> concState(scalar(0), densityState.size());

	for (std::size_t k = 1; k < densityState.size(); ++k)
	{
		concState[k] = densityState[k] / mix.Mv()[k - 1];
		concState[0] += concState[k];
	}

	pressureOutput = mix.pFromUv(concState, internalEnergyState);

	return pressureOutput;
}

schemi::pressureStarClass::~pressureStarClass() noexcept
{
}
