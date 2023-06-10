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
	auto & mesh { inField.meshRef() };

	volumeField<typeOfValue> retVolField { mesh, typeOfValue { 0 } };

	volumeField<std::pair<scalar, std::vector<scalar>> > cellDistances { mesh,
			std::pair<scalar, std::vector<scalar>>(0., 0) };

	for (std::size_t i = 0; i < retVolField.size(); ++i)
	{
		const vector & cellR { mesh.cells()[i].rC() };

		cellDistances.ref_r()[i].second.resize(
				mesh.surfacesOfCells()[i].size());

		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfIndex { mesh.surfacesOfCells()[i][j] };
			const vector & surfaceR { mesh.surfaces()[surfIndex].rC() };
			const scalar deltaRMag { (surfaceR - cellR).mag() };

			cellDistances.ref_r()[i].second[j] = deltaRMag;
			cellDistances.ref_r()[i].first += deltaRMag;
		}
	}

	for (std::size_t i = 0; i < retVolField.size(); ++i)
	{
		for (std::size_t j = 0; j < mesh.surfacesOfCells()[i].size(); ++j)
		{
			const std::size_t surfIndex { mesh.surfacesOfCells()[i][j] };

			retVolField.ref_r()[i] += inField.ref()[surfIndex]
					* (1
							- cellDistances.ref()[i].second[j]
									/ cellDistances.ref()[i].first);
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
	auto & mesh { inField.meshRef() };

	surfaceField<typeOfValue> retSurfField { mesh, typeOfValue { 0 },
			inField.boundCond() };

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
			retSurfField.ref_r()[i] = (inField.ref()[ownIndex] * surfNeiR
					+ inField.ref()[neiIndex] * surfOwnR)
					/ (surfOwnR + surfNeiR);
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

template<typename T>
surfaceField<T> gradientLinearInterpolate(const volumeField<T> & inField,
		const volumeField<returnTypeDivergence<T>> & parentField,
		const boundaryConditionValue & bncCalc, const std::size_t compt =
				componentPlaceholder)
{
	auto & mesh { inField.meshRef() };

	surfaceField<T> retSurfField { mesh, T { 0 }, inField.boundCond() };

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
			retSurfField.ref_r()[i] = (inField.ref()[ownIndex] * surfNeiR
					+ inField.ref()[neiIndex] * surfOwnR)
					/ (surfOwnR + surfNeiR);
		}
			break;
		default:
		{
			const std::size_t ownIndex { mesh.surfaceOwner()[i] };

			const returnTypeDivergence<T> outerCellValue {
					bncCalc.boundaryConditionValueCell(
							parentField.ref()[ownIndex],
							parentField.boundCond()[i], ownIndex, i, compt) };

			const vector deltaR(
					(mesh.surfaces()[i].rC() - mesh.cells()[ownIndex].rC())
							* 2.);

			const scalar deltaRMag { deltaR.mag() };

			const vector deltaRNorm(deltaR / deltaRMag);

			const returnTypeDivergence<T> deltaV { outerCellValue
					- parentField.ref()[ownIndex] };

			retSurfField.ref_r()[i] = (deltaV / deltaRMag) * deltaRNorm;
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
	auto & mesh { inField.meshRef() };

	surfaceField<T> retSurfField { mesh, T { 0 }, inField.boundCond() };

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
			retSurfField.ref_r()[i] = (inField.ref()[ownIndex] * surfNeiR
					+ inField.ref()[neiIndex] * surfOwnR)
					/ (surfOwnR + surfNeiR);
		}
			break;
		default:
		{
			const std::size_t ownIndex { mesh.surfaceOwner()[i] };

			const returnTypeGradient<T> outerCellValue(
					bncCalc.boundaryConditionValueCell(
							parentField.ref()[ownIndex],
							parentField.boundCond()[i], ownIndex, i, compt));

			const vector deltaR(
					(mesh.surfaces()[i].rC() - mesh.cells()[ownIndex].rC())
							* 2.);

			const scalar deltaRMag { deltaR.mag() };

			const vector deltaRNorm(deltaR / deltaRMag);

			const returnTypeGradient<T> deltaV(
					outerCellValue - parentField.ref()[ownIndex]);

			retSurfField.ref_r()[i] = (deltaV / deltaRMag) & deltaRNorm;
		}
			break;
		}
	}
	return retSurfField;
}
}  // namespace schemi

#endif /* LINEARINTERPOLATE_HPP_ */
