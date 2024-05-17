/*
 * harmonicInterpolate.hpp
 *
 *  Created on: 2020/06/17
 *      Author: Maxim Boldyrev
 *
 *      Harmonic interpolation functions.
 */

#ifndef HARMONICINTERPOLATE_HPP_
#define HARMONICINTERPOLATE_HPP_

#include "surfaceField.hpp"
#include "volumeField.hpp"

namespace schemi
{
template<typename typeOfValue>
volumeField<typeOfValue> harmonicInterpolate(
		const surfaceField<typeOfValue> & inField) noexcept
{
	auto & mesh_ { inField.meshP() };

	volumeField<typeOfValue> retVolField { mesh_, typeOfValue { 0 } };

	for (std::size_t i = 0; i < retVolField.size(); ++i)
	{
		const vector & cellR { mesh_.cells()[i].rC() };

		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfIndex { mesh_.surfacesOfCells()[i][j] };

			retVolField.r()[i] += inField()[surfIndex]
					* mesh_.cellSurfaceDistancesR()[i].second[j];
		}

		retVolField.r()[i] /= mesh_.cellSurfaceDistancesR()[i].first;
	}

	return retVolField;
}
}  // namespace schemi

#endif /* HARMONICINTERPOLATE_HPP_ */
