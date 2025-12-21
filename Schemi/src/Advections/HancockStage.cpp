/*
 * HancockStage.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "HancockStage.hpp"

void schemi::HancockStage(homogeneousPhase<quadraticSurface> & surfaceOwnerSide,
		homogeneousPhase<quadraticSurface> & surfaceNeighbourSide,
		const boundaryConditionValue & boundaryConditionValueCalc,
		const MPIHandler & parallelism)
{
	auto & mesh_ { surfaceOwnerSide.pressure.meshRef() };

	/*Creating Hancock divergences.*/
	std::vector<volumeField<scalar>> divRhoVHancock {
			surfaceOwnerSide.phaseThermodynamics->Mv().size() + 1, volumeField<
					scalar> { mesh_, 0 } };

	volumeField<vector> divRhoVVHancock(
			volumeField<vector> { mesh_, vector(0) });
	volumeField<scalar> divRhoEVHancock { mesh_, 0 };
	volumeField<scalar> divRhokVHancock { mesh_, 0 };
	volumeField<scalar> divRhoEpsVHancock { mesh_, 0 };
	volumeField<vector> divRhoaVHancock { mesh_, vector(0) };
	volumeField<scalar> divRhobVHancock { mesh_, 0 };

	/*Hancock divergence calculation.*/
	/*Rho.*/
	for (std::size_t k = 0; k < divRhoVHancock.size(); ++k)
		divRhoVHancock[k] = HancockDivergence(surfaceOwnerSide.velocity,
				surfaceNeighbourSide.velocity, surfaceOwnerSide.density[k],
				surfaceNeighbourSide.density[k]);
	/*Momentum.*/
	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const std::vector<std::size_t> & surfacesOfCell_i {
				mesh_.surfacesOfCells()[i] };

		for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
		{
			const std::size_t surfaceIndex { surfacesOfCell_i[j] };

			if (i == mesh_.surfaceOwner()[surfaceIndex])
				divRhoVVHancock.val()[i] +=
						(((surfaceOwnerSide.momentum.cval()[surfaceIndex]
								* surfaceOwnerSide.velocity.cval()[surfaceIndex])
								& mesh_.surfaces()[surfaceIndex].N())
								+ (mesh_.surfaces()[surfaceIndex].N()
										* (surfaceOwnerSide.pressure.cval()[surfaceIndex]
												+ twothirds
														* surfaceOwnerSide.rhokTurb.cval()[surfaceIndex])))
								* mesh_.surfaces()[surfaceIndex].S();
			else if (i == mesh_.surfaceNeighbour()[surfaceIndex])
				divRhoVVHancock.val()[i] +=
						(((surfaceNeighbourSide.momentum.cval()[surfaceIndex]
								* surfaceNeighbourSide.velocity.cval()[surfaceIndex])
								& (mesh_.surfaces()[surfaceIndex].N() * -1))
								+ ((mesh_.surfaces()[surfaceIndex].N() * -1)
										* (surfaceNeighbourSide.pressure.cval()[surfaceIndex]
												+ twothirds
														* surfaceNeighbourSide.rhokTurb.cval()[surfaceIndex])))
								* mesh_.surfaces()[surfaceIndex].S();
			else
				[[unlikely]]
				throw exception(
						"Cell is neither owner, nor neighbour to surface.",
						errors::systemError);
		}
		divRhoVVHancock.val()[i] /= mesh_.cells()[i].V();
	}
	/*Total energy.*/
	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const std::vector<std::size_t> & surfacesOfCell_i {
				mesh_.surfacesOfCells()[i] };

		for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
		{
			const std::size_t surfaceIndex { surfacesOfCell_i[j] };

			if (i == mesh_.surfaceOwner()[surfaceIndex])
				divRhoEVHancock.val()[i] +=
						(((surfaceOwnerSide.velocity.cval()[surfaceIndex]
								* surfaceOwnerSide.totalEnergy.cval()[surfaceIndex])
								& mesh_.surfaces()[surfaceIndex].N())
								+ ((surfaceOwnerSide.velocity.cval()[surfaceIndex]
										* (surfaceOwnerSide.pressure.cval()[surfaceIndex]
												+ twothirds
														* surfaceOwnerSide.rhokTurb.cval()[surfaceIndex]))
										& mesh_.surfaces()[surfaceIndex].N()))
								* mesh_.surfaces()[surfaceIndex].S();
			else if (i == mesh_.surfaceNeighbour()[surfaceIndex])
				divRhoEVHancock.val()[i] +=
						(((surfaceNeighbourSide.velocity.cval()[surfaceIndex]
								* surfaceNeighbourSide.totalEnergy.cval()[surfaceIndex])
								& (mesh_.surfaces()[surfaceIndex].N() * -1))
								+ ((surfaceNeighbourSide.velocity.cval()[surfaceIndex]
										* (surfaceNeighbourSide.pressure.cval()[surfaceIndex]
												+ twothirds
														* surfaceNeighbourSide.rhokTurb.cval()[surfaceIndex]))
										& (mesh_.surfaces()[surfaceIndex].N()
												* -1)))
								* mesh_.surfaces()[surfaceIndex].S();
			else
				[[unlikely]]
				throw exception(
						"Cell is neither owner, nor neighbour to surface.",
						errors::systemError);
		}
		divRhoEVHancock.val()[i] /= mesh_.cells()[i].V();
	}
	if (surfaceOwnerSide.turbulence->turbulence())
	{
		/*Rhok.*/
		divRhokVHancock = HancockDivergence(surfaceOwnerSide.velocity,
				surfaceNeighbourSide.velocity, surfaceOwnerSide.rhokTurb,
				surfaceNeighbourSide.rhokTurb);
		/*RhoEpsilon.*/
		divRhoEpsVHancock = HancockDivergence(surfaceOwnerSide.velocity,
				surfaceNeighbourSide.velocity, surfaceOwnerSide.rhoepsTurb,
				surfaceNeighbourSide.rhoepsTurb);
		if (surfaceOwnerSide.turbulence->aField())
		{
			/*Rhoa.*/
			divRhoaVHancock = HancockDivergence(surfaceOwnerSide.velocity,
					surfaceNeighbourSide.velocity, surfaceOwnerSide.rhoaTurb,
					surfaceNeighbourSide.rhoaTurb);
			if (surfaceOwnerSide.turbulence->bField())
				/*Rhob.*/
				divRhobVHancock = HancockDivergence(surfaceOwnerSide.velocity,
						surfaceNeighbourSide.velocity,
						surfaceOwnerSide.rhobTurb,
						surfaceNeighbourSide.rhobTurb);
		}
	}

	/*Hancock time integration.*/
	const scalar halfTimestep { 0.5 * mesh_.timestep() };

	/*Rho.*/
	for (std::size_t k = 0; k < divRhoVHancock.size(); ++k)
		HancockTimeIntegration(divRhoVHancock[k], surfaceOwnerSide.density[k],
				surfaceNeighbourSide.density[k], halfTimestep);
	/*Momentum.*/
	HancockTimeIntegration(divRhoVVHancock, surfaceOwnerSide.momentum,
			surfaceNeighbourSide.momentum, halfTimestep);
	/*Total energy.*/
	HancockTimeIntegration(divRhoEVHancock, surfaceOwnerSide.totalEnergy,
			surfaceNeighbourSide.totalEnergy, halfTimestep);
	if (surfaceOwnerSide.turbulence->turbulence())
	{
		/*Rhok.*/
		HancockTimeIntegration(divRhokVHancock, surfaceOwnerSide.rhokTurb,
				surfaceNeighbourSide.rhokTurb, halfTimestep);
		/*RhoEpsilon.*/
		HancockTimeIntegration(divRhoEpsVHancock, surfaceOwnerSide.rhoepsTurb,
				surfaceNeighbourSide.rhoepsTurb, halfTimestep);
		if (surfaceOwnerSide.turbulence->aField())
		{
			/*Rhoa.*/
			HancockTimeIntegration(divRhoaVHancock, surfaceOwnerSide.rhoaTurb,
					surfaceNeighbourSide.rhoaTurb, halfTimestep);
			if (surfaceOwnerSide.turbulence->bField())
				/*Rhob.*/
				HancockTimeIntegration(divRhobVHancock,
						surfaceOwnerSide.rhobTurb,
						surfaceNeighbourSide.rhobTurb, halfTimestep);
		}
	}

	/*Recalculation of other fields.*/
	{
		/*Owner side.*/
		surfaceOwnerSide.concentration.v[0].val() = 0;
		for (std::size_t k = 1; k < surfaceOwnerSide.concentration.v.size();
				++k)
		{
			surfaceOwnerSide.concentration.v[k].val() =
					(surfaceOwnerSide.density[k]
							/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1]).cval();

			surfaceOwnerSide.concentration.v[0] +=
					surfaceOwnerSide.concentration.v[k];
		}

		if (surfaceOwnerSide.turbulence->turbulence())
		{
			surfaceOwnerSide.kTurb.val() = (surfaceOwnerSide.rhokTurb
					/ surfaceOwnerSide.density[0]).cval();

			surfaceOwnerSide.epsTurb.val() = (surfaceOwnerSide.rhoepsTurb
					/ surfaceOwnerSide.density[0]).cval();

			if (surfaceOwnerSide.turbulence->aField())
			{
				surfaceOwnerSide.aTurb.val() = (surfaceOwnerSide.rhoaTurb
						/ surfaceOwnerSide.density[0]).cval();

				if (surfaceOwnerSide.turbulence->bField())
					surfaceOwnerSide.bTurb.val() = (surfaceOwnerSide.rhobTurb
							/ surfaceOwnerSide.density[0]).cval();
			}
		}

		surfaceOwnerSide.velocity.val() = (surfaceOwnerSide.momentum
				/ surfaceOwnerSide.density[0]).cval();

		{
			const auto v2 = surfaceOwnerSide.velocity
					& surfaceOwnerSide.velocity;

			surfaceOwnerSide.internalEnergy.val() =
					(surfaceOwnerSide.totalEnergy
							- surfaceOwnerSide.density[0] * v2 * 0.5
							- surfaceOwnerSide.rhokTurb).cval();
		}

		surfaceOwnerSide.pressure.val() =
				surfaceOwnerSide.phaseThermodynamics->pFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy.cval());

		surfaceOwnerSide.temperature.val() =
				surfaceOwnerSide.phaseThermodynamics->TFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy.cval());

		surfaceOwnerSide.HelmholtzEnergy.val() =
				surfaceOwnerSide.phaseThermodynamics->Fv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature.cval());

		surfaceOwnerSide.entropy.val() =
				surfaceOwnerSide.phaseThermodynamics->Sv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature.cval());

		/*Update parallel boundary conditions*/
		parallelism.correctBoundaryValues(surfaceOwnerSide);

		/*Neighbour side.*/
		for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
			if (mesh_.surfaceNeighbour()[i] != mesh_.nonexistCell())
			{
				surfaceNeighbourSide.concentration.v[0].val()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[k].val()[i] =
							surfaceNeighbourSide.density[k].cval()[i]
									/ surfaceNeighbourSide.phaseThermodynamics->Mv()[k
											- 1];

					surfaceNeighbourSide.concentration.v[0].val()[i] +=
							surfaceNeighbourSide.concentration.v[k].cval()[i];
				}

				surfaceNeighbourSide.velocity.val()[i] =
						surfaceNeighbourSide.momentum.cval()[i]
								/ surfaceNeighbourSide.density[0].cval()[i];

				if (surfaceNeighbourSide.turbulence->turbulence())
				{
					surfaceNeighbourSide.kTurb.val()[i] =
							surfaceNeighbourSide.rhokTurb.cval()[i]
									/ surfaceNeighbourSide.density[0].cval()[i];

					surfaceNeighbourSide.epsTurb.val()[i] =
							surfaceNeighbourSide.rhoepsTurb.cval()[i]
									/ surfaceNeighbourSide.density[0].cval()[i];

					if (surfaceNeighbourSide.turbulence->aField())
					{
						surfaceNeighbourSide.aTurb.val()[i] =
								surfaceNeighbourSide.rhoaTurb.cval()[i]
										/ surfaceNeighbourSide.density[0].cval()[i];

						if (surfaceNeighbourSide.turbulence->bField())
							surfaceNeighbourSide.bTurb.val()[i] =
									surfaceNeighbourSide.rhobTurb.cval()[i]
											/ surfaceNeighbourSide.density[0].cval()[i];
					}
				}

				const scalar v2 { surfaceNeighbourSide.velocity.cval()[i]
						& surfaceNeighbourSide.velocity.cval()[i] };

				surfaceNeighbourSide.internalEnergy.val()[i] =
						surfaceNeighbourSide.totalEnergy.cval()[i]
								- surfaceNeighbourSide.density[0].cval()[i] * v2
										* 0.5
								- surfaceNeighbourSide.rhokTurb.cval()[i];

				std::valarray<scalar> concentrations_i(
						surfaceNeighbourSide.concentration.v.size());
				for (std::size_t k = 0; k < concentrations_i.size(); ++k)
					concentrations_i[k] =
							surfaceNeighbourSide.concentration.v[k].cval()[i];

				surfaceNeighbourSide.pressure.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->pFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy.cval()[i]);

				surfaceNeighbourSide.temperature.val() =
						surfaceNeighbourSide.phaseThermodynamics->TFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy.cval()[i]);

				surfaceNeighbourSide.HelmholtzEnergy.val() =
						surfaceNeighbourSide.phaseThermodynamics->Fv(
								concentrations_i,
								surfaceNeighbourSide.temperature.cval()[i]);

				surfaceNeighbourSide.entropy.val() =
						surfaceNeighbourSide.phaseThermodynamics->Sv(
								concentrations_i,
								surfaceNeighbourSide.temperature.cval()[i]);
			}
			else
			{
				/*Coping from surface owner side.*/
				surfaceNeighbourSide.concentration.v[0].val()[i] = 0;
				surfaceNeighbourSide.density[0].val()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[k].val()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.concentration.v[k].cval()[i],
									surfaceOwnerSide.concentration.v[k].boundCond()[i],
									i, i);

					surfaceNeighbourSide.density[k].val()[i] =
							surfaceNeighbourSide.concentration.v[k].cval()[i]
									* surfaceNeighbourSide.phaseThermodynamics->Mv()[k
											- 1];

					surfaceNeighbourSide.concentration.v[0].val()[i] +=
							surfaceNeighbourSide.concentration.v[k].cval()[i];

					surfaceNeighbourSide.density[0].val()[i] +=
							surfaceNeighbourSide.density[k].cval()[i];
				}

				surfaceNeighbourSide.velocity.val()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.velocity.cval()[i],
								surfaceOwnerSide.velocity.boundCond()[i], i, i);

				surfaceNeighbourSide.momentum.val()[i] =
						surfaceNeighbourSide.velocity.val()[i]
								* surfaceNeighbourSide.density[0].cval()[i];

				surfaceNeighbourSide.pressure.val()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.pressure.cval()[i],
								surfaceOwnerSide.pressure.boundCond()[i], i, i);

				if (surfaceNeighbourSide.turbulence->turbulence())
				{
					surfaceNeighbourSide.kTurb.val()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.kTurb.cval()[i],
									surfaceOwnerSide.kTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhokTurb.val()[i] =
							surfaceNeighbourSide.density[0].cval()[i]
									* surfaceNeighbourSide.kTurb.cval()[i];

					surfaceNeighbourSide.epsTurb.val()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.epsTurb.cval()[i],
									surfaceOwnerSide.epsTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhoepsTurb.val()[i] =
							surfaceNeighbourSide.density[0].cval()[i]
									* surfaceNeighbourSide.epsTurb.cval()[i];

					if (surfaceNeighbourSide.turbulence->aField())
					{
						surfaceNeighbourSide.aTurb.val()[i] =
								boundaryConditionValueCalc.boundaryConditionValueCell(
										surfaceOwnerSide.aTurb.cval()[i],
										surfaceOwnerSide.aTurb.boundCond()[i],
										i, i);

						surfaceNeighbourSide.rhoaTurb.val()[i] =
								surfaceNeighbourSide.aTurb.cval()[i]
										* surfaceNeighbourSide.density[0].cval()[i];

						if (surfaceNeighbourSide.turbulence->bField())
						{
							surfaceNeighbourSide.bTurb.val()[i] =
									boundaryConditionValueCalc.boundaryConditionValueCell(
											surfaceOwnerSide.bTurb.cval()[i],
											surfaceOwnerSide.bTurb.boundCond()[i],
											i, i);

							surfaceNeighbourSide.rhobTurb.val()[i] =
									surfaceNeighbourSide.density[0].cval()[i]
											* surfaceNeighbourSide.bTurb.cval()[i];
						}
					}
				}

				std::valarray<scalar> concentrations_i(
						surfaceNeighbourSide.concentration.v.size());
				for (std::size_t k = 0; k < concentrations_i.size(); ++k)
					concentrations_i[k] =
							surfaceNeighbourSide.concentration.v[k].cval()[i];

				surfaceNeighbourSide.internalEnergy.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->UvFromp(
								concentrations_i,
								surfaceNeighbourSide.pressure.cval()[i]);

				surfaceNeighbourSide.temperature.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->TFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy.cval()[i]);

				const scalar v2 { surfaceNeighbourSide.velocity.cval()[i]
						& surfaceNeighbourSide.velocity.cval()[i] };

				surfaceNeighbourSide.totalEnergy.val()[i] =
						surfaceNeighbourSide.internalEnergy.cval()[i]
								+ surfaceNeighbourSide.density[0].cval()[i] * v2
										* 0.5
								+ surfaceNeighbourSide.rhokTurb.cval()[i];

				surfaceNeighbourSide.HelmholtzEnergy.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Fv(
								concentrations_i,
								surfaceNeighbourSide.temperature.cval()[i]);

				surfaceNeighbourSide.entropy.val()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Sv(
								concentrations_i,
								surfaceNeighbourSide.temperature.cval()[i]);
			}
	}
}
