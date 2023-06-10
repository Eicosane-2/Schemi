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
	auto & mesh { surfaceOwnerSide.pressure.meshRef() };

	conservativeFlows numFluxes(mesh,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());
	starFields starValues(mesh,
			surfaceOwnerSide.phaseThermodynamics->Mv().size());

	for (std::size_t i = 0; i < mesh.surfacesSize(); ++i)
	{
		scalar deltaMag;

		if (mesh.bndType()[i] == boundaryConditionType::innerSurface)
		{
			const auto ownIndex = mesh.surfaceOwner()[i];
			const auto neiIndex = mesh.surfaceNeighbour()[i];

			deltaMag = (mesh.cells()[ownIndex].rC()
					- mesh.cells()[neiIndex].rC()).mag();
		}
		else
		{
			const auto ownIndex = mesh.surfaceOwner()[i];

			deltaMag = mesh.cells()[ownIndex].rC().mag() * 2;
		}

		const scalar turbulentPressureOwner { twothirds
				* surfaceOwnerSide.rhokTurb.ref()[i] };
		const scalar turbulentPressureNeighbour { twothirds
				* surfaceNeighbourSide.rhokTurb.ref()[i] };

		std::valarray<scalar> densityState(numFluxes.density.size());
		vector velocityState;
		scalar pressureState;
		vector aState;
		scalar bState;

		for (std::size_t k = 0; k < densityState.size(); ++k)
			densityState[k] =
					(surfaceNeighbourSide.density[k].ref()[i]
							+ surfaceOwnerSide.density[k].ref()[i]) / 2
							+ (0.5 * mesh.timestep() / deltaMag
									* (surfaceOwnerSide.density[k].ref()[i]
											* surfaceOwnerSide.velocity.ref()[i]
											- surfaceNeighbourSide.density[k].ref()[i]
													* surfaceNeighbourSide.velocity.ref()[i])
									& mesh.surfaces()[i].N());

		const vector momentumState =
				(surfaceNeighbourSide.momentum.ref()[i]
						+ surfaceOwnerSide.momentum.ref()[i]) / 2
						+ (0.5 * mesh.timestep() / deltaMag
								* (surfaceOwnerSide.momentum.ref()[i]
										* surfaceOwnerSide.velocity.ref()[i]
										+ tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
												* (turbulentPressureOwner
														+ surfaceOwnerSide.pressure.ref()[i])
										- surfaceNeighbourSide.momentum.ref()[i]
												* surfaceNeighbourSide.velocity.ref()[i]
										- tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
												* (turbulentPressureNeighbour
														+ surfaceNeighbourSide.pressure.ref()[i]))
								& mesh.surfaces()[i].N());
		const scalar totalEnergyState =
				(surfaceNeighbourSide.totalEnergy.ref()[i]
						+ surfaceOwnerSide.totalEnergy.ref()[i]) / 2
						+ (0.5 * mesh.timestep() / deltaMag
								* ((surfaceOwnerSide.totalEnergy.ref()[i]
										+ (turbulentPressureOwner
												+ surfaceOwnerSide.pressure.ref()[i]))
										* surfaceOwnerSide.velocity.ref()[i]
										- (surfaceNeighbourSide.totalEnergy.ref()[i]
												+ (turbulentPressureNeighbour
														+ surfaceNeighbourSide.pressure.ref()[i]))
												* surfaceNeighbourSide.velocity.ref()[i])
								& mesh.surfaces()[i].N());
		const scalar rhokState =
				(surfaceNeighbourSide.rhokTurb.ref()[i]
						+ surfaceOwnerSide.rhokTurb.ref()[i]) / 2
						+ (0.5 * mesh.timestep() / deltaMag
								* (surfaceOwnerSide.rhokTurb.ref()[i]
										* surfaceOwnerSide.velocity.ref()[i]
										- surfaceNeighbourSide.rhokTurb.ref()[i]
												* surfaceNeighbourSide.velocity.ref()[i])
								& mesh.surfaces()[i].N());
		const scalar rhoepsState =
				(surfaceNeighbourSide.rhoepsTurb.ref()[i]
						+ surfaceOwnerSide.rhoepsTurb.ref()[i]) / 2
						+ (0.5 * mesh.timestep() / deltaMag
								* (surfaceOwnerSide.rhoepsTurb.ref()[i]
										* surfaceOwnerSide.velocity.ref()[i]
										- surfaceNeighbourSide.rhoepsTurb.ref()[i]
												* surfaceNeighbourSide.velocity.ref()[i])
								& mesh.surfaces()[i].N());
		const vector rhoaState =
				(surfaceNeighbourSide.rhoaTurb.ref()[i]
						+ surfaceOwnerSide.rhoaTurb.ref()[i]) / 2
						+ (0.5 * mesh.timestep() / deltaMag
								* (surfaceOwnerSide.rhoaTurb.ref()[i]
										* surfaceOwnerSide.velocity.ref()[i]
										- surfaceNeighbourSide.rhoaTurb.ref()[i]
												* surfaceNeighbourSide.velocity.ref()[i])
								& mesh.surfaces()[i].N());
		const scalar rhobState =
				(surfaceNeighbourSide.rhobTurb.ref()[i]
						+ surfaceOwnerSide.rhobTurb.ref()[i]) / 2
						+ (0.5 * mesh.timestep() / deltaMag
								* (surfaceOwnerSide.rhobTurb.ref()[i]
										* surfaceOwnerSide.velocity.ref()[i]
										- surfaceNeighbourSide.rhobTurb.ref()[i]
												* surfaceNeighbourSide.velocity.ref()[i])
								& mesh.surfaces()[i].N());

		velocityState = momentumState / densityState[0];
		pressureState = pressureStar(*(surfaceOwnerSide.phaseThermodynamics),
				densityState, momentumState, totalEnergyState, rhokState);
		aState = rhoaState / densityState[0];
		bState = rhobState / densityState[0];

		for (std::size_t k = 0; k < numFluxes.density.size(); ++k)
			numFluxes.density[k].ref_r()[i] = densityState[k] * velocityState;

		numFluxes.momentum.ref_r()[i] = momentumState * velocityState
				+ (rhokState * twothirds + pressureState)
						* tensor(1, 0, 0, 0, 1, 0, 0, 0, 1);

		numFluxes.totalEnergy.ref_r()[i] = (totalEnergyState
				+ (rhokState * twothirds + pressureState)) * velocityState;

		numFluxes.rhokTurb.ref_r()[i] = rhokState * velocityState;
		numFluxes.rhoepsTurb.ref_r()[i] = rhoepsState * velocityState;
		numFluxes.rhoaTurb.ref_r()[i] = rhoaState * velocityState;
		numFluxes.rhobTurb.ref_r()[i] = rhobState * velocityState;

		starValues.c[0].ref_r()[i] = 0;
		for (std::size_t k = 1; k < densityState.size(); ++k)
		{
			starValues.c[k].ref_r()[i] = densityState[k]
					/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			starValues.c[0].ref_r()[i] += starValues.c[k].ref()[i];
		}
		starValues.rho.ref_r()[i] = densityState[0];
		starValues.v.ref_r()[i] = velocityState;
		starValues.p.ref_r()[i] = pressureState;
		starValues.a.ref_r()[i] = aState;
		starValues.b.ref_r()[i] = bState;
	}

	return std::make_tuple(numFluxes, starValues);
}

