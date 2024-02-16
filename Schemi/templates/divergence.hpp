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

#include "linearInterpolate.hpp"
#include "returnTypeDivergence.hpp"
#include "intExpPow.hpp"

namespace schemi
{
template<typename Type>
volumeField<returnTypeDivergence<Type>> divergence(
		const volumeField<Type> & inField,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder)
{
	auto & mesh_ { inField.meshRef() };

	const surfaceField<Type> interpolatedField { linearInterpolate(inField,
			bncCalc, compt) };

	volumeField<returnTypeDivergence<Type>> divergence { mesh_,
			returnTypeDivergence<Type> { 0 } };

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		returnTypeDivergence<Type> cellDivValue { 0 };

		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };
			vector normalVector { 0 };

			if (mesh_.surfaceOwner()[surfaceIndex] == i)
				normalVector = mesh_.surfaces()[surfaceIndex].N();
			else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				normalVector = mesh_.surfaces()[surfaceIndex].N() * (-1);
			else
				throw exception("Couldn't choose normal's orientation.",
						errors::systemError);

			cellDivValue += (interpolatedField()[surfaceIndex] & normalVector)
					* mesh_.surfaces()[surfaceIndex].S();
		}

		divergence.r()[i] = cellDivValue / mesh_.cells()[i].V();
	}

	return divergence;
}

template<typename Type>
volumeField<returnTypeDivergence<Type>> divergence(
		const surfaceField<Type> & inField)
{
	auto & mesh_ { inField.meshRef() };

	volumeField<returnTypeDivergence<Type>> divergence { mesh_,
			returnTypeDivergence<Type> { 0 } };

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		returnTypeDivergence<Type> cellDivValue { 0 };

		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfaceIndex { mesh_.surfacesOfCells()[i][j] };
			vector normalVector { 0 };

			if (mesh_.surfaceOwner()[surfaceIndex] == i)
				normalVector = mesh_.surfaces()[surfaceIndex].N();
			else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				normalVector = mesh_.surfaces()[surfaceIndex].N() * (-1);
			else
				throw exception("Couldn't choose normal's orientation.",
						errors::systemError);

			cellDivValue += (inField()[surfaceIndex] & normalVector)
					* mesh_.surfaces()[surfaceIndex].S();
		}

		divergence.r()[i] = cellDivValue / mesh_.cells()[i].V();
	}

	return divergence;
}

template<typename Type>
surfaceField<returnTypeDivergence<Type>> surfDivergence(
		const volumeField<Type> & inField,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder)
{
	auto & mesh_ { inField.meshRef() };

	surfaceField<returnTypeDivergence<Type>> divergence { mesh_,
			returnTypeDivergence<Type> { 0 } };

	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
		switch (inField.boundCond()[i].first)
		{
		case boundaryConditionType::innerSurface:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };
			const std::size_t neiIndex { mesh_.surfaceNeighbour()[i] };

			const vector deltaVec { mesh_.cells()[neiIndex].rC()
					- mesh_.cells()[ownIndex].rC() };

			divergence.r()[i] = (inField()[neiIndex] - inField()[ownIndex])
					/ pow<scalar, 2>(deltaVec.mag()) & deltaVec;
		}
			break;
		case boundaryConditionType::calculatedParallelBoundary:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };

			const Type outerCellValue { bncCalc.boundaryConditionValueCell(
					inField()[ownIndex], inField.boundCond()[i], ownIndex, i,
					compt) };

			const vector deltaR(
					(mesh_.surfaces()[i].rC() - mesh_.cells()[ownIndex].rC())
							- bncCalc.parallelism.cSdR().boundCond()[i].second);

			const scalar deltaRMag { deltaR.mag() };

			const vector deltaRNorm(deltaR / deltaRMag);

			const Type deltaV { outerCellValue - inField()[ownIndex] };

			divergence.r()[i] = (deltaV / deltaRMag) & deltaRNorm;
		}
			break;
		default:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };

			const Type outerCellValue { bncCalc.boundaryConditionValueCell(
					inField()[ownIndex], inField.boundCond()[i], ownIndex, i,
					compt) };

			const vector deltaR(
					(mesh_.surfaces()[i].rC() - mesh_.cells()[ownIndex].rC())
							* 2.);

			const scalar deltaRMag { deltaR.mag() };

			const vector deltaRNorm(deltaR / deltaRMag);

			const Type deltaV { outerCellValue - inField()[ownIndex] };

			divergence.r()[i] = (deltaV / deltaRMag) & deltaRNorm;
		}
			break;
		}

	return divergence;
}
}  // namespace schemi

#endif /* DIVERGENCE_HPP_ */
