/*
 * gradient.hpp
 *
 *  Created on: 2019/11/28
 *      Author: Maxim Boldyrev
 *
 *      Funcitons for gradient calculation.
 */

#ifndef GRADIENT_HPP_
#define GRADIENT_HPP_

#include "linearInterpolate.hpp"
#include "returnTypeGradient.hpp"

namespace schemi
{
template<typename Type>
volumeField<returnTypeGradient<Type>> grad(const volumeField<Type> & inField,
		const boundaryConditionValue & bncCalc, const std::size_t compt = 0)
{
	auto & mesh { inField.meshRef() };

	const surfaceField<Type> interpolatedField { linearInterpolate(inField,
			bncCalc, compt) };

	volumeField<returnTypeGradient<Type>> gradient { mesh, returnTypeGradient<
			Type> { 0 } };

	for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
	{
		returnTypeGradient<Type> cellGradValue { 0 };

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
						errors::systemError);

			cellGradValue += (interpolatedField.ref()[surfaceIndex]
					* normalVector) * mesh.surfaces()[surfaceIndex].S();
		}

		gradient.ref_r()[i] = cellGradValue / mesh.cells()[i].V();
	}

	return gradient;
}

template<typename Type>
volumeField<returnTypeGradient<Type>> grad(const surfaceField<Type> & inField)
{
	auto & mesh { inField.meshRef() };

	volumeField<returnTypeGradient<Type>> gradient { mesh, returnTypeGradient<
			Type> { 0 } };

	for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
	{
		returnTypeGradient<Type> cellGradValue { 0 };

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
						errors::systemError);

			cellGradValue += (inField.ref()[surfaceIndex] * normalVector)
					* mesh.surfaces()[surfaceIndex].S();
		}

		gradient.ref_r()[i] = cellGradValue / mesh.cells()[i].V();
	}

	return gradient;
}
}  // namespace schemi

#endif /* GRADIENT_HPP_ */
