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
				divRhoVVHancock.r()[i] +=
						(((surfaceOwnerSide.momentum()[surfaceIndex]
								* surfaceOwnerSide.velocity()[surfaceIndex])
								& mesh_.surfaces()[surfaceIndex].N())
								+ (mesh_.surfaces()[surfaceIndex].N()
										* (surfaceOwnerSide.pressure()[surfaceIndex]
												+ twothirds
														* surfaceOwnerSide.rhokTurb()[surfaceIndex])))
								* mesh_.surfaces()[surfaceIndex].S();
			else if (i == mesh_.surfaceNeighbour()[surfaceIndex])
				divRhoVVHancock.r()[i] +=
						(((surfaceNeighbourSide.momentum()[surfaceIndex]
								* surfaceNeighbourSide.velocity()[surfaceIndex])
								& (mesh_.surfaces()[surfaceIndex].N() * (-1)))
								+ ((mesh_.surfaces()[surfaceIndex].N() * (-1))
										* (surfaceNeighbourSide.pressure()[surfaceIndex]
												+ twothirds
														* surfaceNeighbourSide.rhokTurb()[surfaceIndex])))
								* mesh_.surfaces()[surfaceIndex].S();
			else
				throw exception(
						"Cell is neither owner, nor neighbour to surface.",
						errors::systemError);
		}
		divRhoVVHancock.r()[i] /= mesh_.cells()[i].V();
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
				divRhoEVHancock.r()[i] +=
						(((surfaceOwnerSide.velocity()[surfaceIndex]
								* surfaceOwnerSide.totalEnergy()[surfaceIndex])
								& mesh_.surfaces()[surfaceIndex].N())
								+ ((surfaceOwnerSide.velocity()[surfaceIndex]
										* (surfaceOwnerSide.pressure()[surfaceIndex]
												+ twothirds
														* surfaceOwnerSide.rhokTurb()[surfaceIndex]))
										& mesh_.surfaces()[surfaceIndex].N()))
								* mesh_.surfaces()[surfaceIndex].S();
			else if (i == mesh_.surfaceNeighbour()[surfaceIndex])
				divRhoEVHancock.r()[i] +=
						(((surfaceNeighbourSide.velocity()[surfaceIndex]
								* surfaceNeighbourSide.totalEnergy()[surfaceIndex])
								& (mesh_.surfaces()[surfaceIndex].N() * (-1)))
								+ ((surfaceNeighbourSide.velocity()[surfaceIndex]
										* (surfaceNeighbourSide.pressure()[surfaceIndex]
												+ twothirds
														* surfaceNeighbourSide.rhokTurb()[surfaceIndex]))
										& (mesh_.surfaces()[surfaceIndex].N()
												* (-1))))
								* mesh_.surfaces()[surfaceIndex].S();
			else
				throw exception(
						"Cell is neither owner, nor neighbour to surface.",
						errors::systemError);
		}
		divRhoEVHancock.r()[i] /= mesh_.cells()[i].V();
	}
	if (surfaceOwnerSide.turbulenceSources->turbulence)
	{
		/*Rhok.*/
		divRhokVHancock = HancockDivergence(surfaceOwnerSide.velocity,
				surfaceNeighbourSide.velocity, surfaceOwnerSide.rhokTurb,
				surfaceNeighbourSide.rhokTurb);
		/*RhoEpsilon.*/
		divRhoEpsVHancock = HancockDivergence(surfaceOwnerSide.velocity,
				surfaceNeighbourSide.velocity, surfaceOwnerSide.rhoepsTurb,
				surfaceNeighbourSide.rhoepsTurb);
		if (surfaceOwnerSide.turbulenceSources->aField)
		{
			/*Rhoa.*/
			divRhoaVHancock = HancockDivergence(surfaceOwnerSide.velocity,
					surfaceNeighbourSide.velocity, surfaceOwnerSide.rhoaTurb,
					surfaceNeighbourSide.rhoaTurb);
			if (surfaceOwnerSide.turbulenceSources->bField)
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
	if (surfaceOwnerSide.turbulenceSources->turbulence)
	{
		/*Rhok.*/
		HancockTimeIntegration(divRhokVHancock, surfaceOwnerSide.rhokTurb,
				surfaceNeighbourSide.rhokTurb, halfTimestep);
		/*RhoEpsilon.*/
		HancockTimeIntegration(divRhoEpsVHancock, surfaceOwnerSide.rhoepsTurb,
				surfaceNeighbourSide.rhoepsTurb, halfTimestep);
		if (surfaceOwnerSide.turbulenceSources->aField)
		{
			/*Rhoa.*/
			HancockTimeIntegration(divRhoaVHancock, surfaceOwnerSide.rhoaTurb,
					surfaceNeighbourSide.rhoaTurb, halfTimestep);
			if (surfaceOwnerSide.turbulenceSources->bField)
				/*Rhob.*/
				HancockTimeIntegration(divRhobVHancock,
						surfaceOwnerSide.rhobTurb,
						surfaceNeighbourSide.rhobTurb, halfTimestep);
		}
	}

	/*Recalculation of other fields.*/
	{
		/*Owner side.*/
		surfaceOwnerSide.concentration.v[0].r() = 0;
		for (std::size_t k = 1; k < surfaceOwnerSide.concentration.v.size();
				++k)
		{
			surfaceOwnerSide.concentration.v[k].r() =
					surfaceOwnerSide.density[k]()
							/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			surfaceOwnerSide.concentration.v[0].r() +=
					surfaceOwnerSide.concentration.v[k]();
		}

		if (surfaceOwnerSide.turbulenceSources->turbulence)
		{
			surfaceOwnerSide.kTurb.r() = surfaceOwnerSide.rhokTurb()
					/ surfaceOwnerSide.density[0]();

			surfaceOwnerSide.epsTurb.r() = surfaceOwnerSide.rhoepsTurb()
					/ surfaceOwnerSide.density[0]();

			if (surfaceOwnerSide.turbulenceSources->aField)
			{
				surfaceOwnerSide.aTurb.r() = division(surfaceOwnerSide.rhoaTurb,
						surfaceOwnerSide.density[0])();

				if (surfaceOwnerSide.turbulenceSources->bField)
					surfaceOwnerSide.bTurb.r() = surfaceOwnerSide.rhobTurb()
							/ surfaceOwnerSide.density[0]();
			}
		}

		surfaceOwnerSide.velocity.r() = division(surfaceOwnerSide.momentum,
				surfaceOwnerSide.density[0])();

		{
			const auto v2 = ampProduct(surfaceOwnerSide.velocity,
					surfaceOwnerSide.velocity);

			surfaceOwnerSide.internalEnergy.r() = surfaceOwnerSide.totalEnergy()
					- surfaceOwnerSide.density[0]() * v2() * 0.5
					- surfaceOwnerSide.rhokTurb();
		}

		surfaceOwnerSide.pressure.r() =
				surfaceOwnerSide.phaseThermodynamics->pFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy());

		surfaceOwnerSide.temperature.r() =
				surfaceOwnerSide.phaseThermodynamics->TFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy());

		surfaceOwnerSide.HelmholtzEnergy.r() =
				surfaceOwnerSide.phaseThermodynamics->Fv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature());

		surfaceOwnerSide.entropy.r() = surfaceOwnerSide.phaseThermodynamics->Sv(
				surfaceOwnerSide.concentration.p,
				surfaceOwnerSide.temperature());

		/*Update parallel boundary conditions*/
		parallelism.correctBoundaryValues(surfaceOwnerSide);

		/*Neighbour side.*/
		for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
			if (mesh_.surfaceNeighbour()[i] != mesh_.nonexistCell())
			{
				surfaceNeighbourSide.concentration.v[0].r()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[k].r()[i] =
							surfaceNeighbourSide.density[k]()[i]
									/ surfaceNeighbourSide.phaseThermodynamics->Mv()[k
											- 1];

					surfaceNeighbourSide.concentration.v[0].r()[i] +=
							surfaceNeighbourSide.concentration.v[k]()[i];
				}

				surfaceNeighbourSide.velocity.r()[i] =
						surfaceNeighbourSide.momentum()[i]
								/ surfaceNeighbourSide.density[0]()[i];

				if (surfaceNeighbourSide.turbulenceSources->turbulence)
				{
					surfaceNeighbourSide.kTurb.r()[i] =
							surfaceNeighbourSide.rhokTurb()[i]
									/ surfaceNeighbourSide.density[0]()[i];

					surfaceNeighbourSide.epsTurb.r()[i] =
							surfaceNeighbourSide.rhoepsTurb()[i]
									/ surfaceNeighbourSide.density[0]()[i];

					if (surfaceNeighbourSide.turbulenceSources->aField)
					{
						surfaceNeighbourSide.aTurb.r()[i] =
								surfaceNeighbourSide.rhoaTurb()[i]
										/ surfaceNeighbourSide.density[0]()[i];

						if (surfaceNeighbourSide.turbulenceSources->bField)
							surfaceNeighbourSide.bTurb.r()[i] =
									surfaceNeighbourSide.rhobTurb()[i]
											/ surfaceNeighbourSide.density[0]()[i];
					}
				}

				const scalar v2 { surfaceNeighbourSide.velocity()[i]
						& surfaceNeighbourSide.velocity()[i] };

				surfaceNeighbourSide.internalEnergy.r()[i] =
						surfaceNeighbourSide.totalEnergy()[i]
								- surfaceNeighbourSide.density[0]()[i] * v2
										* 0.5
								- surfaceNeighbourSide.rhokTurb()[i];

				std::valarray<scalar> concentrations_i(
						surfaceNeighbourSide.concentration.v.size());
				for (std::size_t k = 0; k < concentrations_i.size(); ++k)
					concentrations_i[k] =
							surfaceNeighbourSide.concentration.v[k]()[i];

				surfaceNeighbourSide.pressure.r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->pFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy()[i]);

				surfaceNeighbourSide.temperature.r() =
						surfaceNeighbourSide.phaseThermodynamics->TFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy()[i]);

				surfaceNeighbourSide.HelmholtzEnergy.r() =
						surfaceNeighbourSide.phaseThermodynamics->Fv(
								concentrations_i,
								surfaceNeighbourSide.temperature()[i]);

				surfaceNeighbourSide.entropy.r() =
						surfaceNeighbourSide.phaseThermodynamics->Sv(
								concentrations_i,
								surfaceNeighbourSide.temperature()[i]);
			}
			else
			{
				/*Coping from surface owner side.*/
				surfaceNeighbourSide.concentration.v[0].r()[i] = 0;
				surfaceNeighbourSide.density[0].r()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[k].r()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.concentration.v[k]()[i],
									surfaceOwnerSide.concentration.v[k].boundCond()[i],
									i, i, k);

					surfaceNeighbourSide.density[k].r()[i] =
							surfaceNeighbourSide.concentration.v[k]()[i]
									* surfaceNeighbourSide.phaseThermodynamics->Mv()[k
											- 1];

					surfaceNeighbourSide.concentration.v[0].r()[i] +=
							surfaceNeighbourSide.concentration.v[k]()[i];

					surfaceNeighbourSide.density[0].r()[i] +=
							surfaceNeighbourSide.density[k]()[i];
				}

				surfaceNeighbourSide.velocity.r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.velocity()[i],
								surfaceOwnerSide.velocity.boundCond()[i], i, i);

				surfaceNeighbourSide.momentum.r()[i] =
						surfaceNeighbourSide.velocity.r()[i]
								* surfaceNeighbourSide.density[0]()[i];

				surfaceNeighbourSide.pressure.r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.pressure()[i],
								surfaceOwnerSide.pressure.boundCond()[i], i, i);

				if (surfaceNeighbourSide.turbulenceSources->turbulence)
				{
					surfaceNeighbourSide.kTurb.r()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.kTurb()[i],
									surfaceOwnerSide.kTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhokTurb.r()[i] =
							surfaceNeighbourSide.density[0].r()[i]
									* surfaceNeighbourSide.kTurb()[i];

					surfaceNeighbourSide.epsTurb.r()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.epsTurb()[i],
									surfaceOwnerSide.epsTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhoepsTurb.r()[i] =
							surfaceNeighbourSide.density[0].r()[i]
									* surfaceNeighbourSide.epsTurb()[i];

					if (surfaceNeighbourSide.turbulenceSources->aField)
					{
						surfaceNeighbourSide.aTurb.r()[i] =
								boundaryConditionValueCalc.boundaryConditionValueCell(
										surfaceOwnerSide.aTurb()[i],
										surfaceOwnerSide.aTurb.boundCond()[i],
										i, i);

						surfaceNeighbourSide.rhoaTurb.r()[i] =
								surfaceNeighbourSide.aTurb.r()[i]
										* surfaceNeighbourSide.density[0]()[i];

						if (surfaceNeighbourSide.turbulenceSources->bField)
						{
							surfaceNeighbourSide.bTurb.r()[i] =
									boundaryConditionValueCalc.boundaryConditionValueCell(
											surfaceOwnerSide.bTurb()[i],
											surfaceOwnerSide.bTurb.boundCond()[i],
											i, i);

							surfaceNeighbourSide.rhobTurb.r()[i] =
									surfaceNeighbourSide.density[0].r()[i]
											* surfaceNeighbourSide.bTurb()[i];
						}
					}
				}

				std::valarray<scalar> concentrations_i(
						surfaceNeighbourSide.concentration.v.size());
				for (std::size_t k = 0; k < concentrations_i.size(); ++k)
					concentrations_i[k] =
							surfaceNeighbourSide.concentration.v[k]()[i];

				surfaceNeighbourSide.internalEnergy.r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->UvFromp(
								concentrations_i,
								surfaceNeighbourSide.pressure()[i]);

				surfaceNeighbourSide.temperature.r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->TFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy()[i]);

				const scalar v2 { surfaceNeighbourSide.velocity()[i]
						& surfaceNeighbourSide.velocity()[i] };

				surfaceNeighbourSide.totalEnergy.r()[i] =
						surfaceNeighbourSide.internalEnergy.r()[i]
								+ surfaceNeighbourSide.density[0]()[i] * v2
										* 0.5
								+ surfaceNeighbourSide.rhokTurb()[i];

				surfaceNeighbourSide.HelmholtzEnergy.r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Fv(
								concentrations_i,
								surfaceNeighbourSide.temperature()[i]);

				surfaceNeighbourSide.entropy.r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Sv(
								concentrations_i,
								surfaceNeighbourSide.temperature()[i]);
			}
	}
}
