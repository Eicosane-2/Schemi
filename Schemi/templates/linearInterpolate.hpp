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

	volumeField<std::pair<scalar, std::vector<scalar>> > cellDistances { mesh_,
			std::pair<scalar, std::vector<scalar>>(0., 0) };

	for (std::size_t i = 0; i < retVolField.size(); ++i)
	{
		const vector & cellR { mesh_.cells()[i].rC() };

		cellDistances.r()[i].second.resize(mesh_.surfacesOfCells()[i].size());

		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfIndex { mesh_.surfacesOfCells()[i][j] };
			const vector & surfaceR { mesh_.surfaces()[surfIndex].rC() };
			const scalar deltaRMag { (surfaceR - cellR).mag() };

			cellDistances.r()[i].second[j] = deltaRMag;
			cellDistances.r()[i].first += deltaRMag;
		}
	}

	for (std::size_t i = 0; i < retVolField.size(); ++i)
	{
		for (std::size_t j = 0; j < mesh_.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfIndex { mesh_.surfacesOfCells()[i][j] };

			retVolField.r()[i] += inField()[surfIndex]
					* (1
							- cellDistances()[i].second[j]
									/ cellDistances()[i].first);
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
			const scalar surfOwnR { (mesh_.surfaces()[i].rC()
					- mesh_.cells()[ownIndex].rC()).mag() };
			const scalar surfNeiR { (mesh_.surfaces()[i].rC()
					- mesh_.cells()[neiIndex].rC()).mag() };
			retSurfField.r()[i] = (inField()[ownIndex] * surfNeiR
					+ inField()[neiIndex] * surfOwnR) / (surfOwnR + surfNeiR);
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

template<typename T>
surfaceField<T> divergenceLinearInterpolate(const volumeField<T> & inField,
		const volumeField<returnTypeGradient<T>> & parentField,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder)
{
	auto & mesh_ { inField.meshRef() };

	surfaceField<T> retSurfField { mesh_, T { 0 }, inField.boundCond() };

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
			retSurfField.r()[i] = (inField()[ownIndex] * surfNeiR
					+ inField()[neiIndex] * surfOwnR) / (surfOwnR + surfNeiR);
		}
			break;
		case boundaryConditionType::calculatedParallelBoundary:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };

			const returnTypeGradient<T> outerCellValue(
					bncCalc.boundaryConditionValueCell(parentField()[ownIndex],
							parentField.boundCond()[i], ownIndex, i, compt));

			const vector deltaR(
					(mesh_.surfaces()[i].rC() - mesh_.cells()[ownIndex].rC())
							- bncCalc.parallelism.cSdR().boundCond()[i].second);

			const scalar deltaRMag { deltaR.mag() };

			const vector deltaRNorm(deltaR / deltaRMag);

			const returnTypeGradient<T> deltaV(
					outerCellValue - parentField()[ownIndex]);

			retSurfField.r()[i] = (deltaV / deltaRMag) & deltaRNorm;
		}
			break;
		default:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };

			const returnTypeGradient<T> outerCellValue(
					bncCalc.boundaryConditionValueCell(parentField()[ownIndex],
							parentField.boundCond()[i], ownIndex, i, compt));

			const vector deltaR(
					(mesh_.surfaces()[i].rC() - mesh_.cells()[ownIndex].rC())
							* 2.);

			const scalar deltaRMag { deltaR.mag() };

			const vector deltaRNorm(deltaR / deltaRMag);

			const returnTypeGradient<T> deltaV(
					outerCellValue - parentField()[ownIndex]);

			retSurfField.r()[i] = (deltaV / deltaRMag) & deltaRNorm;
		}
			break;
		}
	}
	return retSurfField;
}

template<typename T>
surfaceField<T> gradientLinearInterpolate(const volumeField<T> & inField,
		const volumeField<returnTypeDivergence<T>> & parentField,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder)
{
	auto & mesh_ { inField.meshRef() };

	surfaceField<T> retSurfField { mesh_, T { 0 }, inField.boundCond() };

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
			retSurfField.r()[i] = (inField()[ownIndex] * surfNeiR
					+ inField()[neiIndex] * surfOwnR) / (surfOwnR + surfNeiR);
		}
			break;
		case boundaryConditionType::calculatedParallelBoundary:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };

			const returnTypeDivergence<T> outerCellValue {
					bncCalc.boundaryConditionValueCell(parentField()[ownIndex],
							parentField.boundCond()[i], ownIndex, i, compt) };

			const vector deltaR { (mesh_.surfaces()[i].rC()
					- mesh_.cells()[ownIndex].rC())
					- bncCalc.parallelism.cSdR().boundCond()[i].second };

			const scalar deltaRMag { deltaR.mag() };

			const vector deltaRNorm(deltaR / deltaRMag);

			const returnTypeDivergence<T> deltaV { outerCellValue
					- parentField()[ownIndex] };

			retSurfField.r()[i] = (deltaV / deltaRMag) * deltaRNorm;
		}
			break;
		default:
		{
			const std::size_t ownIndex { mesh_.surfaceOwner()[i] };

			const returnTypeDivergence<T> outerCellValue {
					bncCalc.boundaryConditionValueCell(parentField()[ownIndex],
							parentField.boundCond()[i], ownIndex, i, compt) };

			const vector deltaR(
					(mesh_.surfaces()[i].rC() - mesh_.cells()[ownIndex].rC())
							* 2.);

			const scalar deltaRMag { deltaR.mag() };

			const vector deltaRNorm(deltaR / deltaRMag);

			const returnTypeDivergence<T> deltaV { outerCellValue
					- parentField()[ownIndex] };

			retSurfField.r()[i] = (deltaV / deltaRMag) * deltaRNorm;
		}
			break;
		}
	}
	return retSurfField;
}
}  // namespace schemi

#endif /* LINEARINTERPOLATE_HPP_ */
