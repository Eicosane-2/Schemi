/*
 * RichtmyerSolver.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "RichtmyerSolver.hpp"

#include "vectorVectorDotProduct.hpp"

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
				* surfaceOwnerSide.rhokTurb()[i] };
		const scalar turbulentPressureNeighbour { twothirds
				* surfaceNeighbourSide.rhokTurb()[i] };

		std::valarray<scalar> densityState(numFluxes.density.size());
		vector velocityState;
		scalar pressureState;
		vector aState;
		scalar bState;

		for (std::size_t k = 0; k < densityState.size(); ++k)
			densityState[k] =
					(surfaceNeighbourSide.density[k]()[i]
							+ surfaceOwnerSide.density[k]()[i]) / 2
							+ (0.5 * mesh_.timestep() / deltaMag
									* (surfaceOwnerSide.density[k]()[i]
											* surfaceOwnerSide.velocity()[i]
											- surfaceNeighbourSide.density[k]()[i]
													* surfaceNeighbourSide.velocity()[i])
									& mesh_.surfaces()[i].N());

		const vector momentumState =
				(surfaceNeighbourSide.momentum()[i]
						+ surfaceOwnerSide.momentum()[i]) / 2
						+ (0.5 * mesh_.timestep() / deltaMag
								* (surfaceOwnerSide.momentum()[i]
										* surfaceOwnerSide.velocity()[i]
										+ tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
												* (turbulentPressureOwner
														+ surfaceOwnerSide.pressure()[i])
										- surfaceNeighbourSide.momentum()[i]
												* surfaceNeighbourSide.velocity()[i]
										- tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
												* (turbulentPressureNeighbour
														+ surfaceNeighbourSide.pressure()[i]))
								& mesh_.surfaces()[i].N());
		const scalar totalEnergyState =
				(surfaceNeighbourSide.totalEnergy()[i]
						+ surfaceOwnerSide.totalEnergy()[i]) / 2
						+ (0.5 * mesh_.timestep() / deltaMag
								* ((surfaceOwnerSide.totalEnergy()[i]
										+ (turbulentPressureOwner
												+ surfaceOwnerSide.pressure()[i]))
										* surfaceOwnerSide.velocity()[i]
										- (surfaceNeighbourSide.totalEnergy()[i]
												+ (turbulentPressureNeighbour
														+ surfaceNeighbourSide.pressure()[i]))
												* surfaceNeighbourSide.velocity()[i])
								& mesh_.surfaces()[i].N());
		const scalar rhokState = (surfaceNeighbourSide.rhokTurb()[i]
				+ surfaceOwnerSide.rhokTurb()[i]) / 2
				+ (0.5 * mesh_.timestep() / deltaMag
						* (surfaceOwnerSide.rhokTurb()[i]
								* surfaceOwnerSide.velocity()[i]
								- surfaceNeighbourSide.rhokTurb()[i]
										* surfaceNeighbourSide.velocity()[i])
						& mesh_.surfaces()[i].N());
		const scalar rhoepsState = (surfaceNeighbourSide.rhoepsTurb()[i]
				+ surfaceOwnerSide.rhoepsTurb()[i]) / 2
				+ (0.5 * mesh_.timestep() / deltaMag
						* (surfaceOwnerSide.rhoepsTurb()[i]
								* surfaceOwnerSide.velocity()[i]
								- surfaceNeighbourSide.rhoepsTurb()[i]
										* surfaceNeighbourSide.velocity()[i])
						& mesh_.surfaces()[i].N());
		const vector rhoaState = (surfaceNeighbourSide.rhoaTurb()[i]
				+ surfaceOwnerSide.rhoaTurb()[i]) / 2
				+ (0.5 * mesh_.timestep() / deltaMag
						* (surfaceOwnerSide.rhoaTurb()[i]
								* surfaceOwnerSide.velocity()[i]
								- surfaceNeighbourSide.rhoaTurb()[i]
										* surfaceNeighbourSide.velocity()[i])
						& mesh_.surfaces()[i].N());
		const scalar rhobState = (surfaceNeighbourSide.rhobTurb()[i]
				+ surfaceOwnerSide.rhobTurb()[i]) / 2
				+ (0.5 * mesh_.timestep() / deltaMag
						* (surfaceOwnerSide.rhobTurb()[i]
								* surfaceOwnerSide.velocity()[i]
								- surfaceNeighbourSide.rhobTurb()[i]
										* surfaceNeighbourSide.velocity()[i])
						& mesh_.surfaces()[i].N());

		velocityState = momentumState / densityState[0];
		pressureState = pressureStar(*(surfaceOwnerSide.phaseThermodynamics),
				densityState, momentumState, totalEnergyState, rhokState);
		aState = rhoaState / densityState[0];
		bState = rhobState / densityState[0];

		for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
			numFluxes.density[k].r()[i] = densityState[k] * velocityState;

		numFluxes.momentum.r()[i] = momentumState * velocityState
				+ (rhokState * twothirds + pressureState)
						* tensor(1, 0, 0, 0, 1, 0, 0, 0, 1);

		numFluxes.totalEnergy.r()[i] = (totalEnergyState
				+ (rhokState * twothirds + pressureState)) * velocityState;

		numFluxes.rhokTurb.r()[i] = rhokState * velocityState;
		numFluxes.rhoepsTurb.r()[i] = rhoepsState * velocityState;
		numFluxes.rhoaTurb.r()[i] = rhoaState * velocityState;
		numFluxes.rhobTurb.r()[i] = rhobState * velocityState;

		starValues.c.v[0].r()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c.v[k].r()[i] = densityState[k]
					/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			starValues.c.v[0].r()[i] += starValues.c.v[k]()[i];
		}
		starValues.rho.r()[i] = densityState[0];
		starValues.v.r()[i] = velocityState;
		starValues.p.r()[i] = pressureState;
		starValues.a.r()[i] = aState;
		starValues.b.r()[i] = bState;
	}

	return std::make_tuple(numFluxes, starValues);
}

