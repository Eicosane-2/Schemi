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
	auto & mesh { inField.meshP() };

	volumeField<typeOfValue> retVolField { mesh, typeOfValue { 0 } };

	for (std::size_t i = 0; i < retVolField.size(); ++i)
	{
		const vector & cellR { mesh.cells()[i].rC() };
		scalar sumOfVecMag { 0 };

		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfIndex { mesh.surfacesOfCells()[i][j] };
			const vector & surfaceR { mesh.surfaces()[surfIndex].rC() };
			const scalar deltaRMag1 { 1. / (surfaceR - cellR).mag() };

			retVolField.ref_r()[i] += inField.ref()[surfIndex] * deltaRMag1;
			sumOfVecMag += deltaRMag1;
		}

		retVolField.ref_r()[i] /= sumOfVecMag;
	}

	return retVolField;
}
}  // namespace schemi

#endif /* HARMONICINTERPOLATE_HPP_ */
