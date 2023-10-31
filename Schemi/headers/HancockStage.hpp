/*
 * HancockStage.hpp
 *
 *  Created on: 2020/02/23
 *      Author: Maxim Boldyrev
 *
 *      Function for semi-timestep time integration in Hancock method for second order in time.
 */

#ifndef HANCOCKSTAGE_HPP_
#define HANCOCKSTAGE_HPP_

#include "boundaryConditionValue.hpp"
#include "MPIHandler.hpp"
#include "quadraticSurface.hpp"
#include "homogeneousPhase.hpp"
#include "surfaceField.hpp"
#include "volumeField.hpp"

namespace schemi
{
template<typename T>
volumeField<T> HancockDivergence(
		const surfaceField<vector> & surfaceOwnerSideVelocity,
		const surfaceField<vector> & surfaceNeighbourSideVelocity,
		const surfaceField<T> & surfaceOwnerSideT,
		const surfaceField<T> & surfaceNeighbourSideT)
{
	auto & mesh_ { surfaceOwnerSideVelocity.meshRef() };

	volumeField<T> divTVHancock { mesh_, T(0) };

	for (std::size_t i = 0; i < divTVHancock.size(); ++i)
	{
		const std::vector<std::size_t> & surfacesOfCell_i {
				mesh_.surfacesOfCells()[i] };

		for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
		{
			const std::size_t surfaceIndex { surfacesOfCell_i[j] };

			if (i == mesh_.surfaceOwner()[surfaceIndex])
				divTVHancock.r()[i] += ((surfaceOwnerSideT()[surfaceIndex]
						* surfaceOwnerSideVelocity()[surfaceIndex])
						& mesh_.surfaces()[surfaceIndex].N())
						* mesh_.surfaces()[surfaceIndex].S();
			else if (i == mesh_.surfaceNeighbour()[surfaceIndex])
				divTVHancock.r()[i] += ((surfaceNeighbourSideT()[surfaceIndex]
						* surfaceNeighbourSideVelocity()[surfaceIndex])
						& (mesh_.surfaces()[surfaceIndex].N() * (-1)))
						* mesh_.surfaces()[surfaceIndex].S();
			else
				throw exception(
						"Cell is neither owner, nor neighbour to surface.",
						errors::systemError);
		}
		divTVHancock.r()[i] /= mesh_.cells()[i].V();
	}

	return divTVHancock;
}

template<typename T>
void HancockTimeIntegration(const volumeField<T> & flowDivergence,
		surfaceField<T> & surfaceOwnerSideT,
		surfaceField<T> & surfaceNeighbourSideT, const scalar halfTimestep)
{
	auto & mesh_ { surfaceOwnerSideT.meshRef() };

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const std::vector<std::size_t> & surfacesOfCell_i {
				mesh_.surfacesOfCells()[i] };

		for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
		{
			const std::size_t surfaceIndex { surfacesOfCell_i[j] };

			if (i == mesh_.surfaceOwner()[surfaceIndex])
				surfaceOwnerSideT.r()[surfaceIndex] -= flowDivergence()[i]
						* halfTimestep;
			else if (i == mesh_.surfaceNeighbour()[surfaceIndex])
				surfaceNeighbourSideT.r()[surfaceIndex] -= flowDivergence()[i]
						* halfTimestep;
			else
				throw exception(
						"Cell is neither owner, nor neighbour to surface.",
						errors::systemError);
		}
	}
}

void HancockStage(homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
		homogeneousPhase<quadraticSurface> & surfaceNeighbourSide,
		const boundaryConditionValue & boundaryConditionValueCalc,
		const MPIHandler & parallelism);
}  // namespace schemi

#endif /* HANCOCKSTAGE_HPP_ */
