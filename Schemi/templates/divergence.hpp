/*
 * divergence.hpp
 *
 *  Created on: 2019/11/30
 *      Author: Maxim Boldyrev
 *
 *      Functions for divergence calculation.
 */

#ifndef DIVERGENCE_HPP_
#define DIVERGENCE_HPP_

#include "returnTypeDivergence.hpp"

namespace schemi
{
template<typename Type>
volumeField<returnTypeDivergence<Type>> divergence(
		const volumeField<Type> & inField,
		const boundaryConditionValue & bncCalc, const std::size_t compt = 0)
{
	auto & mesh { inField.meshRef() };

	const surfaceField<Type> interpolatedField { linearInterpolate(inField,
			bncCalc, compt) };

	volumeField<returnTypeDivergence<Type>> divergence { mesh,
			returnTypeDivergence<Type> { 0 } };

	for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
	{
		returnTypeDivergence<Type> cellDivValue { 0 };

		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh.surfacesOfCells()[i][j] };
			vector normalVector { 0 };

			if (mesh.surfaceOwner()[surfaceIndex] == i)
				normalVector = mesh.surfaces()[surfaceIndex].N();
			else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				normalVector = mesh.surfaces()[surfaceIndex].N() * (-1);
			else
				throw exception("Couldn't choose normal's orientation.",
						errorsEnum::systemError);

			cellDivValue += (interpolatedField.ref()[surfaceIndex]
					& normalVector) * mesh.surfaces()[surfaceIndex].S();
		}

		divergence.ref_r()[i] = cellDivValue / mesh.cells()[i].V();
	}

	return divergence;
}

template<typename Type>
volumeField<returnTypeDivergence<Type>> divergence(
		const surfaceField<Type> & inField)
{
	auto & mesh { inField.meshRef() };

	volumeField<returnTypeDivergence<Type>> divergence { mesh,
			returnTypeDivergence<Type> { 0 } };

	for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
	{
		returnTypeDivergence<Type> cellDivValue { 0 };

		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh.surfacesOfCells()[i][j] };
			vector normalVector { 0 };

			if (mesh.surfaceOwner()[surfaceIndex] == i)
				normalVector = mesh.surfaces()[surfaceIndex].N();
			else if (mesh.surfaceNeighbour()[surfaceIndex] == i)
				normalVector = mesh.surfaces()[surfaceIndex].N() * (-1);
			else
				throw exception("Couldn't choose normal's orientation.",
						errorsEnum::systemError);

			cellDivValue += (inField.ref()[surfaceIndex] & normalVector)
					* mesh.surfaces()[surfaceIndex].S();
		}

		divergence.ref_r()[i] = cellDivValue / mesh.cells()[i].V();
	}

	return divergence;
}
}  // namespace schemi

#endif /* DIVERGENCE_HPP_ */
