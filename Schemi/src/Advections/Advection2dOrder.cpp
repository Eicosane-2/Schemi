/*
 * Advection2dOrder.cpp
 *
 *  Created on: 2024/11/20
 *      Author: Maxim Boldyrev
 */

#include "Advection2dOrder.hpp"

#include <chrono>

#include "divergence.hpp"
#include "gradient.hpp"
#include "TVDLimiter.hpp"

schemi::starFields schemi::Advection2dOrder(
		homogeneousPhase<cubicCell> & gasPhase, const abstractLimiter & limiter,
		const abstractFlowSolver & fsolver, std::pair<bool, vector> gravitation,
		const boundaryConditionValue & boundaryConditionValueCalc,
		scalar & timeForTVD, scalar & timeForHancock,
		scalar & timeForFlowCalculation, scalar & timeForTimeIntegration,
		const MPIHandler & parallelism)
{
	auto & mesh_ { gasPhase.pressure.meshRef() };

	/*Creating surface fields.*/
	homogeneousPhase<quadraticSurface> surfaceOwnerSide { bunchOfFields<
			quadraticSurface>(gasPhase),
			transportCoefficients<quadraticSurface>(mesh_),
			gasPhase.phaseThermodynamics, gasPhase.turbulence,
			gasPhase.transportModel };

	homogeneousPhase<quadraticSurface> surfaceNeighbourSide { bunchOfFields<
			quadraticSurface>(gasPhase),
			transportCoefficients<quadraticSurface>(mesh_),
			gasPhase.phaseThermodynamics, gasPhase.turbulence,
			gasPhase.transportModel };

	/*Creating limited gradients.*/
	std::vector<volumeField<vector>> concentrationTVDGradient {
			gasPhase.phaseThermodynamics->Mv().size(), volumeField<vector>(
					mesh_, vector(0)) };
	volumeField<tensor> velocityTVDGradient { mesh_, tensor(0) };
	volumeField<vector> pressureTVDGradient { mesh_, vector(0) };
	volumeField<vector> kTVDGradient { mesh_, vector(0) };
	volumeField<vector> epsilonTVDGradient { mesh_, vector(0) };
	volumeField<tensor> aTVDGradient { mesh_, tensor(0) };
	volumeField<vector> bTVDGradient { mesh_, vector(0) };

	/*TVD Reconstruction.*/
	const auto TVDStartTime = std::chrono::high_resolution_clock::now();
	{
		/*TVD limiters.*/
		for (std::size_t k = 0; k < concentrationTVDGradient.size(); ++k)
		{
			const auto concentrationGradient_k = grad(
					gasPhase.concentration.v[k + 1], boundaryConditionValueCalc,
					k + 1);

			concentrationTVDGradient[k] = TVDLimiter(concentrationGradient_k,
					gasPhase.concentration.v[k + 1], limiter,
					boundaryConditionValueCalc, parallelism, k + 1);
		}

		const auto velocityGradient = grad(gasPhase.velocity,
				boundaryConditionValueCalc);

		velocityTVDGradient = TVDLimiter(velocityGradient, gasPhase.velocity,
				limiter, boundaryConditionValueCalc, parallelism);

		const auto pressureGradient = grad(gasPhase.pressure,
				boundaryConditionValueCalc);

		pressureTVDGradient = TVDLimiter(pressureGradient, gasPhase.pressure,
				limiter, boundaryConditionValueCalc, parallelism);

		if (gasPhase.turbulence->turbulence())
		{
			const auto kGradient = grad(gasPhase.kTurb,
					boundaryConditionValueCalc);

			kTVDGradient = TVDLimiter(kGradient, gasPhase.kTurb, limiter,
					boundaryConditionValueCalc, parallelism);

			const auto epsilonGradient = grad(gasPhase.epsTurb,
					boundaryConditionValueCalc);

			epsilonTVDGradient = TVDLimiter(epsilonGradient, gasPhase.epsTurb,
					limiter, boundaryConditionValueCalc, parallelism);

			if (gasPhase.turbulence->aField())
			{
				const auto aGradient = grad(gasPhase.aTurb,
						boundaryConditionValueCalc);

				aTVDGradient = TVDLimiter(aGradient, gasPhase.aTurb, limiter,
						boundaryConditionValueCalc, parallelism);

				if (gasPhase.turbulence->bField())
				{
					const auto bGradient = grad(gasPhase.bTurb,
							boundaryConditionValueCalc);

					bTVDGradient = TVDLimiter(bGradient, gasPhase.bTurb,
							limiter, boundaryConditionValueCalc, parallelism);
				}
			}
		}

		for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		{
			const std::vector<std::size_t> & surfacesOfCell_i {
					mesh_.surfacesOfCells()[i] };

			for (std::size_t j = 0; j < surfacesOfCell_i.size(); ++j)
			{
				const std::size_t surfaceIndex { surfacesOfCell_i[j] };

				const vector deltaVector = mesh_.surfaces()[surfaceIndex].rC()
						- mesh_.cells()[i].rC();

				std::valarray<scalar> reconstructedConcentrationsValues(
						gasPhase.phaseThermodynamics->Mv().size());

				for (std::size_t k = 0;
						k < reconstructedConcentrationsValues.size(); ++k)
					reconstructedConcentrationsValues[k] =
							(concentrationTVDGradient[k]()[i] & deltaVector)
									+ gasPhase.concentration.v[k + 1]()[i];

				const vector reconstructedVelocityValue =
						(velocityTVDGradient()[i] & deltaVector)
								+ gasPhase.velocity()[i];

				const scalar reconstructedPressureValue {
						(pressureTVDGradient()[i] & deltaVector)
								+ gasPhase.pressure()[i] };

				scalar reconstructedkValue { 0 };
				scalar reconstructedEpsilon { 0 };
				vector reconstructedaValue { 0 };
				scalar reconstructedbValue { 0 };
				if (gasPhase.turbulence->turbulence())
				{
					reconstructedkValue = (kTVDGradient()[i] & deltaVector)
							+ gasPhase.kTurb()[i];

					reconstructedEpsilon = (epsilonTVDGradient()[i]
							& deltaVector) + gasPhase.epsTurb()[i];

					if (gasPhase.turbulence->aField())
					{
						reconstructedaValue = (aTVDGradient()[i] & deltaVector)
								+ gasPhase.aTurb()[i];

						if (gasPhase.turbulence->bField())
							reconstructedbValue = (bTVDGradient()[i]
									& deltaVector) + gasPhase.bTurb()[i];
					}
				}

				if (mesh_.surfaceOwner()[surfaceIndex] == i)
				{
					for (std::size_t k = 0;
							k < reconstructedConcentrationsValues.size(); ++k)
						surfaceOwnerSide.concentration.v[k + 1].r()[surfaceIndex] =
								reconstructedConcentrationsValues[k];

					surfaceOwnerSide.velocity.r()[surfaceIndex] =
							reconstructedVelocityValue;

					surfaceOwnerSide.pressure.r()[surfaceIndex] =
							reconstructedPressureValue;

					if (gasPhase.turbulence->turbulence())
					{
						surfaceOwnerSide.kTurb.r()[surfaceIndex] =
								reconstructedkValue;

						surfaceOwnerSide.epsTurb.r()[surfaceIndex] =
								reconstructedEpsilon;

						if (gasPhase.turbulence->aField())
						{
							surfaceOwnerSide.aTurb.r()[surfaceIndex] =
									reconstructedaValue;

							if (gasPhase.turbulence->bField())
								surfaceOwnerSide.bTurb.r()[surfaceIndex] =
										reconstructedbValue;
						}
					}
				}
				else if (mesh_.surfaceNeighbour()[surfaceIndex] == i)
				{
					for (std::size_t k = 0;
							k < reconstructedConcentrationsValues.size(); ++k)
						surfaceNeighbourSide.concentration.v[k + 1].r()[surfaceIndex] =
								reconstructedConcentrationsValues[k];

					surfaceNeighbourSide.velocity.r()[surfaceIndex] =
							reconstructedVelocityValue;

					surfaceNeighbourSide.pressure.r()[surfaceIndex] =
							reconstructedPressureValue;

					if (gasPhase.turbulence->turbulence())
					{
						surfaceNeighbourSide.kTurb.r()[surfaceIndex] =
								reconstructedkValue;

						surfaceNeighbourSide.epsTurb.r()[surfaceIndex] =
								reconstructedEpsilon;

						if (gasPhase.turbulence->aField())
						{
							surfaceNeighbourSide.aTurb.r()[surfaceIndex] =
									reconstructedaValue;

							if (gasPhase.turbulence->bField())
								surfaceNeighbourSide.bTurb.r()[surfaceIndex] =
										reconstructedbValue;
						}
					}
				}
				else
					[[unlikely]]
					throw exception("Couldn't choose side to add",
							errors::systemError);
			}
		}

		/*Recalculation of other quantities for owner side*/
		surfaceOwnerSide.concentration.v[0].r() = 0;
		surfaceOwnerSide.density[0].r() = 0;
		for (std::size_t k = 1; k < surfaceOwnerSide.concentration.v.size();
				++k)
		{
			surfaceOwnerSide.concentration.v[0].r() +=
					surfaceOwnerSide.concentration.v[k]();

			surfaceOwnerSide.density[k].r() =
					surfaceOwnerSide.concentration.v[k]()
							* surfaceOwnerSide.phaseThermodynamics->Mv()[k - 1];

			surfaceOwnerSide.density[0].r() += surfaceOwnerSide.density[k]();
		}

		if (gasPhase.turbulence->turbulence())
		{
			surfaceOwnerSide.rhokTurb.r() = astProduct(surfaceOwnerSide.kTurb,
					surfaceOwnerSide.density[0])();

			surfaceOwnerSide.rhoepsTurb.r() = surfaceOwnerSide.epsTurb()
					* surfaceOwnerSide.density[0]();

			if (gasPhase.turbulence->aField())
			{
				surfaceOwnerSide.rhoaTurb.r() = astProduct(
						surfaceOwnerSide.aTurb, surfaceOwnerSide.density[0])();

				if (gasPhase.turbulence->bField())
					surfaceOwnerSide.rhobTurb.r() = surfaceOwnerSide.bTurb()
							* surfaceOwnerSide.density[0]();
			}
		}

		surfaceOwnerSide.internalEnergy.r() =
				surfaceOwnerSide.phaseThermodynamics->UvFromp(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.pressure());

		surfaceOwnerSide.temperature.r() =
				surfaceOwnerSide.phaseThermodynamics->TFromUv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.internalEnergy());

		surfaceOwnerSide.momentum.r() = astProduct(surfaceOwnerSide.velocity,
				surfaceOwnerSide.density[0])();

		{
			const auto v2 = ampProduct(surfaceOwnerSide.velocity,
					surfaceOwnerSide.velocity);

			surfaceOwnerSide.totalEnergy.r() = surfaceOwnerSide.internalEnergy()
					+ surfaceOwnerSide.density[0]() * v2() * 0.5
					+ surfaceOwnerSide.rhokTurb();
		}

		surfaceOwnerSide.HelmholtzEnergy.r() =
				surfaceOwnerSide.phaseThermodynamics->Fv(
						surfaceOwnerSide.concentration.p,
						surfaceOwnerSide.temperature());

		surfaceOwnerSide.entropy.r() = surfaceOwnerSide.phaseThermodynamics->Sv(
				surfaceOwnerSide.concentration.p,
				surfaceOwnerSide.temperature());

		/*Recalculation of other quantities for neighbour side*/
		const std::size_t nonExistentCell = mesh_.nonexistCell();
		for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
			if (mesh_.surfaceNeighbour()[i] != nonExistentCell)
			{
				surfaceNeighbourSide.concentration.v[0].r()[i] = 0;
				surfaceNeighbourSide.density[0].r()[i] = 0;
				for (std::size_t k = 1;
						k < surfaceNeighbourSide.concentration.v.size(); ++k)
				{
					surfaceNeighbourSide.concentration.v[0].r()[i] +=
							surfaceNeighbourSide.concentration.v[k]()[i];

					surfaceNeighbourSide.density[k].r()[i] =
							surfaceNeighbourSide.concentration.v[k]()[i]
									* gasPhase.phaseThermodynamics->Mv()[k - 1];

					surfaceNeighbourSide.density[0].r()[i] +=
							surfaceNeighbourSide.density[k]()[i];
				}

				surfaceNeighbourSide.momentum.r()[i] =
						surfaceNeighbourSide.velocity()[i]
								* surfaceNeighbourSide.density[0]()[i];

				if (gasPhase.turbulence->turbulence())
				{
					surfaceNeighbourSide.rhokTurb.r()[i] =
							surfaceNeighbourSide.kTurb()[i]
									* surfaceNeighbourSide.density[0]()[i];

					surfaceNeighbourSide.rhoepsTurb.r()[i] =
							surfaceNeighbourSide.epsTurb()[i]
									* surfaceNeighbourSide.density[0]()[i];

					if (gasPhase.turbulence->aField())
					{
						surfaceNeighbourSide.rhoaTurb.r()[i] =
								surfaceNeighbourSide.aTurb()[i]
										* surfaceNeighbourSide.density[0]()[i];

						if (gasPhase.turbulence->bField())
							surfaceNeighbourSide.rhobTurb.r()[i] =
									surfaceNeighbourSide.bTurb()[i]
											* surfaceNeighbourSide.density[0]()[i];
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
						surfaceNeighbourSide.internalEnergy()[i]
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
	const auto TVDEndTime = std::chrono::high_resolution_clock::now();
	timeForTVD += std::chrono::duration_cast<std::chrono::milliseconds>(
			TVDEndTime - TVDStartTime).count();

	const auto HancockStartTime = std::chrono::high_resolution_clock::now();
	const auto HancockEndTime = std::chrono::high_resolution_clock::now();
	timeForHancock += std::chrono::duration_cast<std::chrono::milliseconds>(
			HancockEndTime - HancockStartTime).count();

	/*Update parallel boundary conditions*/
	parallelism.correctBoundaryValues(surfaceOwnerSide);

	/*Outer neighbour side.*/
	const std::size_t nonExistentCell = mesh_.nonexistCell();
	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
		if (mesh_.surfaceNeighbour()[i] == nonExistentCell)
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

			if (surfaceNeighbourSide.turbulence->turbulence())
			{
				surfaceNeighbourSide.kTurb.r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.kTurb()[i],
								surfaceOwnerSide.kTurb.boundCond()[i], i, i);

				surfaceNeighbourSide.rhokTurb.r()[i] =
						surfaceNeighbourSide.density[0].r()[i]
								* surfaceNeighbourSide.kTurb()[i];

				surfaceNeighbourSide.epsTurb.r()[i] =
						boundaryConditionValueCalc.boundaryConditionValueCell(
								surfaceOwnerSide.epsTurb()[i],
								surfaceOwnerSide.epsTurb.boundCond()[i], i, i);

				surfaceNeighbourSide.rhoepsTurb.r()[i] =
						surfaceNeighbourSide.density[0]()[i]
								* surfaceNeighbourSide.epsTurb()[i];

				if (surfaceNeighbourSide.turbulence->aField())
				{
					surfaceNeighbourSide.aTurb.r()[i] =
							boundaryConditionValueCalc.boundaryConditionValueCell(
									surfaceOwnerSide.aTurb()[i],
									surfaceOwnerSide.aTurb.boundCond()[i], i,
									i);

					surfaceNeighbourSide.rhoaTurb.r()[i] =
							surfaceNeighbourSide.aTurb.r()[i]
									* surfaceNeighbourSide.density[0]()[i];

					if (surfaceNeighbourSide.turbulence->bField())
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
							+ surfaceNeighbourSide.density[0]()[i] * v2 * 0.5
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

	/*Numerical fluxes.*/
	const auto FlowCalculationStartTime =
			std::chrono::high_resolution_clock::now();

	/*Creating fields of flux solver.*/
	const auto [NumFluxFlows, star] = fsolver.calculateFlows(surfaceOwnerSide,
			surfaceNeighbourSide);

	const auto FlowCalculationEndTime =
			std::chrono::high_resolution_clock::now();
	timeForFlowCalculation += std::chrono::duration_cast<
			std::chrono::milliseconds>(
			FlowCalculationEndTime - FlowCalculationStartTime).count();

	/*Time integration.*/
	const auto TimeIntegrationStartTime =
			std::chrono::high_resolution_clock::now();
	{
		const scalar timestep = mesh_.timestep();

		gasPhase.totalEnergy.r() -= divergence(NumFluxFlows.totalEnergy)()
				* timestep;

		gasPhase.momentum.r() -= astProduct(divergence(NumFluxFlows.momentum),
				timestep)();

		if (gravitation.first)
		{
			gasPhase.totalEnergy.r() += gasPhase.density[0]()
					* ampProduct(gasPhase.velocity, gravitation.second)()
					* timestep;

			gasPhase.momentum.r() += astProduct(
					astProduct(gasPhase.density[0], gravitation.second),
					timestep)();
		}

		for (std::size_t k = 0; k < gasPhase.density.size(); ++k)
			gasPhase.density[k].r() -= divergence(NumFluxFlows.density[k])()
					* timestep;

		if (gasPhase.turbulence->turbulence())
		{
			gasPhase.rhokTurb.r() -= divergence(NumFluxFlows.rhokTurb)()
					* timestep;

			gasPhase.rhoepsTurb.r() -= divergence(NumFluxFlows.rhoepsTurb)()
					* timestep;

			if (gasPhase.turbulence->aField())
			{
				gasPhase.rhoaTurb.r() -= astProduct(
						divergence(NumFluxFlows.rhoaTurb), timestep)();

				if (gasPhase.turbulence->bField())
					gasPhase.rhobTurb.r() -= divergence(NumFluxFlows.rhobTurb)()
							* timestep;
			}
		}

		gasPhase.concentration.v[0].r() = 0;
		for (std::size_t k = 1; k < gasPhase.concentration.v.size(); ++k)
		{
			gasPhase.concentration.v[k].r() = gasPhase.density[k]()
					/ gasPhase.phaseThermodynamics->Mv()[k - 1];

			gasPhase.concentration.v[0].r() += gasPhase.concentration.v[k]();
		}

		gasPhase.velocity.r() =
				division(gasPhase.momentum, gasPhase.density[0])();

		{
			const auto v2 = ampProduct(gasPhase.velocity, gasPhase.velocity);

			gasPhase.internalEnergy.r() = gasPhase.totalEnergy()
					- gasPhase.density[0]() * v2() * 0.5 - gasPhase.rhokTurb();
		}

		gasPhase.pressure.r() = gasPhase.phaseThermodynamics->pFromUv(
				gasPhase.concentration.p, gasPhase.internalEnergy());

		gasPhase.temperature.r() = gasPhase.phaseThermodynamics->TFromUv(
				gasPhase.concentration.p, gasPhase.internalEnergy());

		gasPhase.HelmholtzEnergy.r() = gasPhase.phaseThermodynamics->Fv(
				gasPhase.concentration.p, gasPhase.temperature());

		gasPhase.entropy.r() = gasPhase.phaseThermodynamics->Sv(
				gasPhase.concentration.p, gasPhase.temperature());

		if (gasPhase.turbulence->turbulence())
		{
			gasPhase.kTurb.r() = gasPhase.rhokTurb() / gasPhase.density[0]();

			gasPhase.epsTurb.r() = gasPhase.rhoepsTurb()
					/ gasPhase.density[0]();

			if (gasPhase.turbulence->aField())
			{
				gasPhase.aTurb.r() = division(gasPhase.rhoaTurb,
						gasPhase.density[0])();

				if (gasPhase.turbulence->bField())
					gasPhase.bTurb.r() = gasPhase.rhobTurb()
							/ gasPhase.density[0]();
			}
		}
	}
	const auto TimeIntegrationEndTime =
			std::chrono::high_resolution_clock::now();
	timeForTimeIntegration += std::chrono::duration_cast<
			std::chrono::milliseconds>(
			TimeIntegrationEndTime - TimeIntegrationStartTime).count();

	return star;
}
