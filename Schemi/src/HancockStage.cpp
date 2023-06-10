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
	auto & mesh { surfaceOwnerSide.pressure.meshRef() };

	/*Creating Hancock divergences.*/
	std::vector<volumeField<scalar>> divRhoVHancock {
			surfaceOwnerSide.phaseThermodynamics->Mv().size() + 1, volumeField<
					scalar> { mesh, 0 } };

	volumeField<vector> divRhoVVHancock(
			volumeField<vector> { mesh, vector(0) });
	volumeField<scalar> divRhoEVHancock { mesh, 0 };
	volumeField<scalar> divRhokVHancock { mesh, 0 };
	volumeField<scalar> divRhoEpsVHancock { mesh, 0 };
	volumeField<vector> divRhoaVHancock { mesh, vector(0) };
	volumeField<scalar> divRhobVHancock { mesh, 0 };

	/*Hancock divergence calculation.*/
	/*Rho.*/
	for (std::size_t k = 0; k < divRhoVHancock.size(); ++k)
		divRhoVHancock[k] = HancockDivergence(surfaceOwnerSide.velocity,
				surfaceNeighbourSide.velocity, surfaceOwnerSide.density[k],
				surfaceNeighbourSide.density[k]);
	/*Momentum.*/
	for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
	{
		const std::vector<std::size_t> & surfacesOfCell_i {
				mesh.surfacesOfCells()[i] };

		for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
		{
			const std::size_t surfaceIndex { surfacesOfCell_i[j] };

			if (i == mesh.surfaceOwner()[surfaceIndex])
				divRhoVVHancock.ref_r()[i] +=
						(((surfaceOwnerSide.momentum.ref()[surfaceIndex]
								* surfaceOwnerSide.velocity.ref()[surfaceIndex])
								& mesh.surfaces()[surfaceIndex].N())
								+ (mesh.surfaces()[surfaceIndex].N()
										* (surfaceOwnerSide.pressure.ref()[surfaceIndex]
												+ twothirds
														* surfaceOwnerSide.rhokTurb.ref()[surfaceIndex])))
								* mesh.surfaces()[surfaceIndex].S();
			else if (i == mesh.surfaceNeighbour()[surfaceIndex])
				divRhoVVHancock.ref_r()[i] +=
						(((surfaceNeighbourSide.momentum.ref()[surfaceIndex]
								* surfaceNeighbourSide.velocity.ref()[surfaceIndex])
								& (mesh.surfaces()[surfaceIndex].N() * (-1)))
								+ ((mesh.surfaces()[surfaceIndex].N() * (-1))
										* (surfaceNeighbourSide.pressure.ref()[surfaceIndex]
												+ twothirds
														* surfaceNeighbourSide.rhokTurb.ref()[surfaceIndex])))
								* mesh.surfaces()[surfaceIndex].S();
			else
				throw exception(
						"Cell is neither owner, nor neighbour to surface.",
						errorsEnum::systemError);
		}
		divRhoVVHancock.ref_r()[i] /= mesh.cells()[i].V();
	}
	/*Total energy.*/
	for (std::size_t i = 0; i < mesh.cellsSize(); ++i)
	{
		const std::vector<std::size_t> & surfacesOfCell_i {
				mesh.surfacesOfCells()[i] };

		for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
		{
			const std::size_t surfaceIndex { surfacesOfCell_i[j] };

			if (i == mesh.surfaceOwner()[surfaceIndex])
				divRhoEVHancock.ref_r()[i] +=
						(((surfaceOwnerSide.velocity.ref()[surfaceIndex]
								* surfaceOwnerSide.totalEnergy.ref()[surfaceIndex])
								& mesh.surfaces()[surfaceIndex].N())
								+ ((surfaceOwnerSide.velocity.ref()[surfaceIndex]
										* (surfaceOwnerSide.pressure.ref()[surfaceIndex]
												+ twothirds
														* surfaceOwnerSide.rhokTurb.ref()[surfaceIndex]))
										& mesh.surfaces()[surfaceIndex].N()))
								* mesh.surfaces()[surfaceIndex].S();
			else if (i == mesh.surfaceNeighbour()[surfaceIndex])
				divRhoEVHancock.ref_r()[i] +=
						(((surfaceNeighbourSide.velocity.ref()[surfaceIndex]
								* surfaceNeighbourSide.totalEnergy.ref()[surfaceIndex])
								& (mesh.surfaces()[surfaceIndex].N() * (-1)))
								+ ((surfaceNeighbourSide.velocity.ref()[surfaceIndex]
										* (surfaceNeighbourSide.pressure.ref()[surfaceIndex]
												+ twothirds
														* surfaceNeighbourSide.rhokTurb.ref()[surfaceIndex]))
										& (mesh.surfaces()[surfaceIndex].N()
												* (-1))))
								* mesh.surfaces()[surfaceIndex].S();
			else
				throw exception(
						"Cell is neither owner, nor neighbour to surface.",
						errorsEnum::systemError);
		}
		divRhoEVHancock.ref_r()[i] /= mesh.cells()[i].V();
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
		if ((surfaceOwnerSide.turbulenceSources->model
				== turbulenceModelEnum::BHRSource)
				|| (surfaceOwnerSide.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (surfaceOwnerSide.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
		{
			/*Rhoa.*/
			divRhoaVHancock = HancockDivergence(surfaceOwnerSide.velocity,
					surfaceNeighbourSide.velocity, surfaceOwnerSide.rhoaTurb,
					surfaceNeighbourSide.rhoaTurb);
			if ((surfaceOwnerSide.turbulenceSources->model
					== turbulenceModelEnum::BHRSource)
					|| (surfaceOwnerSide.turbulenceSources->model
							== turbulenceModelEnum::BHRKLSource))
				/*Rhob.*/
				divRhobVHancock = HancockDivergence(surfaceOwnerSide.velocity,
						surfaceNeighbourSide.velocity,
						surfaceOwnerSide.rhobTurb,
						surfaceNeighbourSide.rhobTurb);
		}
	}

	/*Hancock time integration.*/
	const scalar halfTimestep { 0.5 * mesh.timestep() };

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
		if ((surfaceOwnerSide.turbulenceSources->model
				== turbulenceModelEnum::BHRSource)
				|| (surfaceOwnerSide.turbulenceSources->model
						== turbulenceModelEnum::BHRKLSource)
				|| (surfaceOwnerSide.turbulenceSources->model
						== turbulenceModelEnum::kEpsASource))
		{
			/*Rhoa.*/
			HancockTimeIntegration(divRhoaVHancock, surfaceOwnerSide.rhoaTurb,
					surfaceNeighbourSide.rhoaTurb, halfTimestep);
			if ((surfaceOwnerSide.turbulenceSources->model
					== turbulenceModelEnum::BHRSource)
					|| (surfaceOwnerSide.turbulenceSources->model
							== turbulenceModelEnum::BHRKLSource))
				/*Rhob.*/
				HancockTimeIntegration(divRhobVHancock,
						surfaceOwnerSide.rhobTurb,
						surfaceNeighbourSide.rhobTurb, halfTimestep);
		}
	}

	/*Recalculation of other fields.*/
	{
		/*Owner side.*/
		surfaceOwnerSide.concentration.v[0].ref_r() = 0;
		for (std::size_t k = 1; k < surfaceOwnerSide.concentration.v.size();
				++k)
		{
			surfaceOwnerSide.concentration.v[k].ref_r() =
					surfaceOwnerSide.density[k].ref()
							/ surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			surfaceOwnerSide.concentration.v[0].ref_r() +=
					surfaceOwnerSide.concentration.v[k].ref();
		}

		if (surfaceOwnerSide.turbulenceSources->turbulence)
		{
			surfaceOwnerSide.kTurb.ref_r() = surfaceOwnerSide.rhokTurb.ref()
					/ surfaceOwnerSide.density[0].ref();

			surfaceOwnerSide.epsTurb.ref_r() = surfaceOwnerSide.rhoepsTurb.ref()
					/ surfaceOwnerSide.density[0].ref();

			if ((surfaceOwnerSide.turbulenceSources->model
					== turbulenceModelEnum::BHRSource)
					|| (surfaceOwnerSide.turbulenceSources->model
							== turbulenceModelEnum::BHRKLSource)
					|| (surfaceOwnerSide.turbulenceSources->model
							== turbulenceModelEnum::kEpsASource))
			{
				surfaceOwnerSide.aTurb.ref_r() =
						division(surfaceOwnerSide.rhoaTurb,
								surfaceOwnerSide.density[0]).ref();

				if ((surfaceOwnerSide.turbulenceSources->model
						== turbulenceModelEnum::BHRSource)
						|| (surfaceOwnerSide.turbulenceSources->model
								== turbulenceModelEnum::BHRKLSource))
					surfaceOwnerSide.bTurb.ref_r() =
							surfaceOwnerSide.rhobTurb.ref()
									/ surfaceOwnerSide.density[0].ref();
			}
		}

		surfaceOwnerSide.velocity.ref_r() = division(surfaceOwnerSide.momentum,
				surfaceOwnerSide.density[0]).ref();

		{
			const auto v2 = ampProduct(surfaceOwnerSide.velocity,
					surfaceOwnerSide.velocity);

			surfaceOwnerSide.internalEnergy.ref_r() =
					surfaceOwnerSide.totalEnergy.ref()
							- surfaceOwnerSide.density[0].ref() * v2.ref() * 0.5
							- surfaceOwnerSide.rhokTurb.ref();
		}

		surfaceOwnerSide.pressure.ref_r() =
				surfaceOwnerSide.phaseThermodynamics->pFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy.ref());

		surfaceOwnerSide.temperature.ref_r() =
				surfaceOwnerSide.phaseThermodynamics->TFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy.ref());

		surfaceOwnerSide.HelmholtzEnergy.ref_r() =
				surfaceOwnerSide.phaseThermodynamics->Fv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature.ref());

		surfaceOwnerSide.entropy.ref_r() =
				surfaceOwnerSide.phaseThermodynamics->Sv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature.ref());

		/*Update parallel boundary conditions*/
		parallelism.correctBoundaryValues(surfaceOwnerSide);

		/*Neighbour side.*/
		for (std::size_t i = 0; i < mesh.surfacesSize(); ++i)
			if (mesh.surfaceNeighbour()[i] != mesh.nonexistCell())
			{
				surfaceNeighbourSide.concentration.v[0].ref_r()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[k].ref_r()[i] =
							surfaceNeighbourSide.density[k].ref()[i]
									/ surfaceNeighbourSide.phaseThermodynamics->Mv()[k
											- 1];

					surfaceNeighbourSide.concentration.v[0].ref_r()[i] +=
							surfaceNeighbourSide.concentration.v[k].ref()[i];
				}

				surfaceNeighbourSide.velocity.ref_r()[i] =
						surfaceNeighbourSide.momentum.ref()[i]
								/ surfaceNeighbourSide.density[0].ref()[i];

				if (surfaceNeighbourSide.turbulenceSources->turbulence)
				{
					surfaceNeighbourSide.kTurb.ref_r()[i] =
							surfaceNeighbourSide.rhokTurb.ref()[i]
									/ surfaceNeighbourSide.density[0].ref()[i];

					surfaceNeighbourSide.epsTurb.ref_r()[i] =
							surfaceNeighbourSide.rhoepsTurb.ref()[i]
									/ surfaceNeighbourSide.density[0].ref()[i];

					if ((surfaceNeighbourSide.turbulenceSources->model
							== turbulenceModelEnum::BHRSource)
							|| (surfaceNeighbourSide.turbulenceSources->model
									== turbulenceModelEnum::BHRKLSource)
							|| (surfaceNeighbourSide.turbulenceSources->model
									== turbulenceModelEnum::kEpsASource))
					{
						surfaceNeighbourSide.aTurb.ref_r()[i] =
								surfaceNeighbourSide.rhoaTurb.ref()[i]
										/ surfaceNeighbourSide.density[0].ref()[i];

						if ((surfaceNeighbourSide.turbulenceSources->model
								== turbulenceModelEnum::BHRSource)
								|| (surfaceNeighbourSide.turbulenceSources->model
										== turbulenceModelEnum::BHRKLSource))
							surfaceNeighbourSide.bTurb.ref_r()[i] =
									surfaceNeighbourSide.rhobTurb.ref()[i]
											/ surfaceNeighbourSide.density[0].ref()[i];
					}
				}

				const scalar v2 { surfaceNeighbourSide.velocity.ref()[i]
						& surfaceNeighbourSide.velocity.ref()[i] };

				surfaceNeighbourSide.internalEnergy.ref_r()[i] =
						surfaceNeighbourSide.totalEnergy.ref()[i]
								- surfaceNeighbourSide.density[0].ref()[i] * v2
										* 0.5
								- surfaceNeighbourSide.rhokTurb.ref()[i];

				std::valarray<scalar> concentrations_i(
						surfaceNeighbourSide.concentration.v.size());
				for (std::size_t k = 0; k < concentrations_i.size(); ++k)
					concentrations_i[k] =
							surfaceNeighbourSide.concentration.v[k].ref()[i];

				surfaceNeighbourSide.pressure.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->pFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy.ref()[i]);

				surfaceNeighbourSide.temperature.ref_r() =
						surfaceNeighbourSide.phaseThermodynamics->TFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy.ref()[i]);

				surfaceNeighbourSide.HelmholtzEnergy.ref_r() =
						surfaceNeighbourSide.phaseThermodynamics->Fv(
								concentrations_i,
								surfaceNeighbourSide.temperature.ref()[i]);

				surfaceNeighbourSide.entropy.ref_r() =
						surfaceNeighbourSide.phaseThermodynamics->Sv(
								concentrations_i,
								surfaceNeighbourSide.temperature.ref()[i]);
			}
			else
			{
				/*Coping from surface owner side.*/
				surfaceNeighbourSide.concentration.v[0].ref_r()[i] = 0;
				surfaceNeighbourSide.density[0].ref_r()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[k].ref_r()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.concentration.v[k].ref()[i],
									surfaceOwnerSide.concentration.v[k].boundCond()[i],
									i, i, k);

					surfaceNeighbourSide.density[k].ref_r()[i] =
							surfaceNeighbourSide.concentration.v[k].ref()[i]
									* surfaceNeighbourSide.phaseThermodynamics->Mv()[k
											- 1];

					surfaceNeighbourSide.concentration.v[0].ref_r()[i] +=
							surfaceNeighbourSide.concentration.v[k].ref()[i];

					surfaceNeighbourSide.density[0].ref_r()[i] +=
							surfaceNeighbourSide.density[k].ref()[i];
				}

				surfaceNeighbourSide.velocity.ref_r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.velocity.ref()[i],
								surfaceOwnerSide.velocity.boundCond()[i], i, i);

				surfaceNeighbourSide.momentum.ref_r()[i] =
						surfaceNeighbourSide.velocity.ref_r()[i]
								* surfaceNeighbourSide.density[0].ref()[i];

				surfaceNeighbourSide.pressure.ref_r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.pressure.ref()[i],
								surfaceOwnerSide.pressure.boundCond()[i], i, i);

				if (surfaceNeighbourSide.turbulenceSources->turbulence)
				{
					surfaceNeighbourSide.kTurb.ref_r()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.kTurb.ref()[i],
									surfaceOwnerSide.kTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhokTurb.ref_r()[i] =
							surfaceNeighbourSide.density[0].ref_r()[i]
									* surfaceNeighbourSide.kTurb.ref()[i];

					surfaceNeighbourSide.epsTurb.ref_r()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.epsTurb.ref()[i],
									surfaceOwnerSide.epsTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhoepsTurb.ref_r()[i] =
							surfaceNeighbourSide.density[0].ref_r()[i]
									* surfaceNeighbourSide.epsTurb.ref()[i];

					if ((surfaceNeighbourSide.turbulenceSources->model
							== turbulenceModelEnum::BHRSource)
							|| (surfaceNeighbourSide.turbulenceSources->model
									== turbulenceModelEnum::BHRKLSource)
							|| (surfaceNeighbourSide.turbulenceSources->model
									== turbulenceModelEnum::kEpsASource))
					{
						surfaceNeighbourSide.aTurb.ref_r()[i] =
								boundaryConditionValueCalc.boundaryConditionValueCell(
										surfaceOwnerSide.aTurb.ref()[i],
										surfaceOwnerSide.aTurb.boundCond()[i],
										i, i);

						surfaceNeighbourSide.rhoaTurb.ref_r()[i] =
								surfaceNeighbourSide.aTurb.ref_r()[i]
										* surfaceNeighbourSide.density[0].ref()[i];

						if ((surfaceNeighbourSide.turbulenceSources->model
								== turbulenceModelEnum::BHRSource)
								|| (surfaceNeighbourSide.turbulenceSources->model
										== turbulenceModelEnum::BHRKLSource))
						{
							surfaceNeighbourSide.bTurb.ref_r()[i] =
									boundaryConditionValueCalc.boundaryConditionValueCell(
											surfaceOwnerSide.bTurb.ref()[i],
											surfaceOwnerSide.bTurb.boundCond()[i],
											i, i);

							surfaceNeighbourSide.rhobTurb.ref_r()[i] =
									surfaceNeighbourSide.density[0].ref_r()[i]
											* surfaceNeighbourSide.bTurb.ref()[i];
						}
					}
				}

				std::valarray<scalar> concentrations_i(
						surfaceNeighbourSide.concentration.v.size());
				for (std::size_t k = 0; k < concentrations_i.size(); ++k)
					concentrations_i[k] =
							surfaceNeighbourSide.concentration.v[k].ref()[i];

				surfaceNeighbourSide.internalEnergy.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->UvFromp(
								concentrations_i,
								surfaceNeighbourSide.pressure.ref()[i]);

				surfaceNeighbourSide.temperature.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->TFromUv(
								concentrations_i,
								surfaceNeighbourSide.internalEnergy.ref()[i]);

				const scalar v2 { surfaceNeighbourSide.velocity.ref()[i]
						& surfaceNeighbourSide.velocity.ref()[i] };

				surfaceNeighbourSide.totalEnergy.ref_r()[i] =
						surfaceNeighbourSide.internalEnergy.ref_r()[i]
								+ surfaceNeighbourSide.density[0].ref()[i] * v2
										* 0.5
								+ surfaceNeighbourSide.rhokTurb.ref()[i];

				surfaceNeighbourSide.HelmholtzEnergy.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Fv(
								concentrations_i,
								surfaceNeighbourSide.temperature.ref()[i]);

				surfaceNeighbourSide.entropy.ref_r()[i] =
						surfaceNeighbourSide.phaseThermodynamics->Sv(
								concentrations_i,
								surfaceNeighbourSide.temperature.ref()[i]);
			}
	}
}
