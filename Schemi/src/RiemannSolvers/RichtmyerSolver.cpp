/*
 * RichtmyerSolver.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "RichtmyerSolver.hpp"

#include "vector.hpp"

std::tuple<schemi::conservativeFlows, schemi::starFields> schemi::RichtmyerSolver::calculateFlows(
		const homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
		const homogeneousPhase<quadraticSurface> & surfaceNeighbourSide) const
{
	auto & mesh_ { surfaceOwnerSide.pressure.meshRef() };

	conservativeFlows numFluxes(mesh_,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());
	starFields starValues(mesh_,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());

	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
	{
		scalar deltaMag;

		switch (mesh_.bndType()[i])
		{
		case boundaryConditionType::innerSurface:
		{
			const auto ownIndex = mesh_.surfaceOwner()[i];
			const auto neiIndex = mesh_.surfaceNeighbour()[i];

			deltaMag = (mesh_.cells()[ownIndex].rC()
					- mesh_.cells()[neiIndex].rC()).mag();
		}
			break;
		case boundaryConditionType::calculatedParallelBoundary:
		{
			const auto ownIndex = mesh_.surfaceOwner()[i];

			const auto outerDelta =
					parallelism.cSdR().boundCond()[i].second.mag();

			const auto innerDelta = (mesh_.surfaces()[i].rC()
					- mesh_.cells()[ownIndex].rC()).mag();

			deltaMag = outerDelta + innerDelta;
		}
			break;
		default:
		{
			const auto ownIndex = mesh_.surfaceOwner()[i];

			deltaMag =
					(mesh_.surfaces()[i].rC() - mesh_.cells()[ownIndex].rC()).mag()
							* 2;
		}
			break;
		}

		const scalar turbulentPressureOwner { twothirds
				* surfaceOwnerSide.rhokTurb.cval()[i] };
		const scalar turbulentPressureNeighbour { twothirds
				* surfaceNeighbourSide.rhokTurb.cval()[i] };

		std::valarray<scalar> densityState(numFluxes.density.size());
		vector velocityState;
		scalar pressureState;
		vector aState;
		scalar bState;

		for (std::size_t k = 0; k < densityState.size(); ++k)
			densityState[k] =
					(surfaceNeighbourSide.density[k].cval()[i]
							+ surfaceOwnerSide.density[k].cval()[i]) / 2
							+ (0.5 * mesh_.timestep() / deltaMag
									* (surfaceOwnerSide.density[k].cval()[i]
											* surfaceOwnerSide.velocity.cval()[i]
											- surfaceNeighbourSide.density[k].cval()[i]
													* surfaceNeighbourSide.velocity.cval()[i])
									& mesh_.surfaces()[i].N());

		const vector momentumState =
				(surfaceNeighbourSide.momentum.cval()[i]
						+ surfaceOwnerSide.momentum.cval()[i]) / 2
						+ (0.5 * mesh_.timestep() / deltaMag
								* (surfaceOwnerSide.momentum.cval()[i]
										* surfaceOwnerSide.velocity.cval()[i]
										+ tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
												* (turbulentPressureOwner
														+ surfaceOwnerSide.pressure.cval()[i])
										- surfaceNeighbourSide.momentum.cval()[i]
												* surfaceNeighbourSide.velocity.cval()[i]
										- tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
												* (turbulentPressureNeighbour
														+ surfaceNeighbourSide.pressure.cval()[i]))
								& mesh_.surfaces()[i].N());
		const scalar totalEnergyState =
				(surfaceNeighbourSide.totalEnergy.cval()[i]
						+ surfaceOwnerSide.totalEnergy.cval()[i]) / 2
						+ (0.5 * mesh_.timestep() / deltaMag
								* ((surfaceOwnerSide.totalEnergy.cval()[i]
										+ (turbulentPressureOwner
												+ surfaceOwnerSide.pressure.cval()[i]))
										* surfaceOwnerSide.velocity.cval()[i]
										- (surfaceNeighbourSide.totalEnergy.cval()[i]
												+ (turbulentPressureNeighbour
														+ surfaceNeighbourSide.pressure.cval()[i]))
												* surfaceNeighbourSide.velocity.cval()[i])
								& mesh_.surfaces()[i].N());
		const scalar rhokState =
				(surfaceNeighbourSide.rhokTurb.cval()[i]
						+ surfaceOwnerSide.rhokTurb.cval()[i]) / 2
						+ (0.5 * mesh_.timestep() / deltaMag
								* (surfaceOwnerSide.rhokTurb.cval()[i]
										* surfaceOwnerSide.velocity.cval()[i]
										- surfaceNeighbourSide.rhokTurb.cval()[i]
												* surfaceNeighbourSide.velocity.cval()[i])
								& mesh_.surfaces()[i].N());
		const scalar rhoepsState =
				(surfaceNeighbourSide.rhoepsTurb.cval()[i]
						+ surfaceOwnerSide.rhoepsTurb.cval()[i]) / 2
						+ (0.5 * mesh_.timestep() / deltaMag
								* (surfaceOwnerSide.rhoepsTurb.cval()[i]
										* surfaceOwnerSide.velocity.cval()[i]
										- surfaceNeighbourSide.rhoepsTurb.cval()[i]
												* surfaceNeighbourSide.velocity.cval()[i])
								& mesh_.surfaces()[i].N());
		const vector rhoaState =
				(surfaceNeighbourSide.rhoaTurb.cval()[i]
						+ surfaceOwnerSide.rhoaTurb.cval()[i]) / 2
						+ (0.5 * mesh_.timestep() / deltaMag
								* (surfaceOwnerSide.rhoaTurb.cval()[i]
										* surfaceOwnerSide.velocity.cval()[i]
										- surfaceNeighbourSide.rhoaTurb.cval()[i]
												* surfaceNeighbourSide.velocity.cval()[i])
								& mesh_.surfaces()[i].N());
		const scalar rhobState =
				(surfaceNeighbourSide.rhobTurb.cval()[i]
						+ surfaceOwnerSide.rhobTurb.cval()[i]) / 2
						+ (0.5 * mesh_.timestep() / deltaMag
								* (surfaceOwnerSide.rhobTurb.cval()[i]
										* surfaceOwnerSide.velocity.cval()[i]
										- surfaceNeighbourSide.rhobTurb.cval()[i]
												* surfaceNeighbourSide.velocity.cval()[i])
								& mesh_.surfaces()[i].N());

		velocityState = momentumState / densityState[0];
		pressureState = pressureStar(*(surfaceOwnerSide.phaseThermodynamics),
				densityState, momentumState, totalEnergyState, rhokState);
		aState = rhoaState / densityState[0];
		bState = rhobState / densityState[0];

		for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
			numFluxes.density[k].val()[i] = densityState[k] * velocityState;

		numFluxes.momentum.val()[i] = momentumState * velocityState
				+ (rhokState * twothirds + pressureState)
						* tensor(1, 0, 0, 0, 1, 0, 0, 0, 1);

		numFluxes.totalEnergy.val()[i] = (totalEnergyState
				+ (rhokState * twothirds + pressureState)) * velocityState;

		numFluxes.rhokTurb.val()[i] = rhokState * velocityState;
		numFluxes.rhoepsTurb.val()[i] = rhoepsState * velocityState;
		numFluxes.rhoaTurb.val()[i] = rhoaState * velocityState;
		numFluxes.rhobTurb.val()[i] = rhobState * velocityState;

		starValues.c.v[0].val()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c.v[k].val()[i] = densityState[k]
					/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			starValues.c.v[0].val()[i] += starValues.c.v[k].cval()[i];
		}
		starValues.rho.val()[i] = densityState[0];
		starValues.v.val()[i] = velocityState;
		starValues.p.val()[i] = pressureState;
		starValues.a.val()[i] = aState;
		starValues.b.val()[i] = bState;
	}

	return std::make_tuple(numFluxes, starValues);
}

