/*
 * Diffusion.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "Diffusion.hpp"

#include <chrono>

#include "arithmeticAModel.hpp"
#include "harmonicInterpolateScalar.hpp"
#include "divergence.hpp"
#include "diffusiveFields.hpp"
#include "gradient.hpp"
#include "kEpsAModel.hpp"
#include "linearInterpolate.hpp"

void schemi::Diffusion(homogeneousPhase<cubicCell> & gasPhase,
		const abstractMatrixSolver & msolver,
		const abstractMatrixSolver & msolverForEnthalpy,
		const std::pair<scalar, scalar> & timestepCoeffs,
		scalar & timeForDiffusion,
		const std::vector<boundaryConditionType> & commonConditions,
		const starFields & star, const enthalpyFlow enthalpySolverFlag,
		const bool linearRec, const boundaryConditionValue & bncCalc,
		const volumeField<scalar> & minimalLengthScale,
		[[maybe_unused]] const MPIHandler & parallelism,
		const timestep sourceTimeFlag, const bool molMassDiffusionFlag,
		const bool nonLinearityIteratonsFlag)
{
	constexpr bool distributionOn { true };

	const auto DiffusionStartTime { std::chrono::high_resolution_clock::now() };

	auto & mesh_ { gasPhase.pressure.meshRef() };

	const scalar timestep = mesh_.timestep();

	/*Fields for diffusion stage.*/
	volumeField<scalar> avMolMass { mesh_, 0 };
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedAverageMolarMass);
		avMolMass = volumeField<scalar>(mesh_, 0, subPatchData<scalar> {
				bndCon[0] }, subPatchData<scalar> { bndCon[1] },
				subPatchData<scalar> { bndCon[2] }, subPatchData<scalar> {
						bndCon[3] }, subPatchData<scalar> { bndCon[4] },
				subPatchData<scalar> { bndCon[5] });
	}

	surfaceField<scalar> surfaceRho { mesh_, 0 };

	volumeField<scalar> CVOld { mesh_, 0 };
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedCv);
		CVOld = volumeField<scalar>(mesh_, 0,
				subPatchData<scalar> { bndCon[0] }, subPatchData<scalar> {
						bndCon[1] }, subPatchData<scalar> { bndCon[2] },
				subPatchData<scalar> { bndCon[3] }, subPatchData<scalar> {
						bndCon[4] }, subPatchData<scalar> { bndCon[5] });
	}

	volumeField<scalar> CCVOld { mesh_, 0 };
	volumeField<scalar> CCVNew { mesh_, 0 };
	volumeField<scalar> nonIdealCorrectionOld { mesh_, 0 };
	volumeField<scalar> nonIdealCorrectionNew { mesh_, 0 };
	volumeField<scalar> NonIdRho { mesh_, 0 };
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedNonidealityCorrectionPerDensity);
		NonIdRho = volumeField<scalar>(mesh_, 0, subPatchData<scalar> {
				bndCon[0] }, subPatchData<scalar> { bndCon[1] },
				subPatchData<scalar> { bndCon[2] }, subPatchData<scalar> {
						bndCon[3] }, subPatchData<scalar> { bndCon[4] },
				subPatchData<scalar> { bndCon[5] });
	}

	volumeField<scalar> CvM { mesh_, 0 };
	{
		std::vector<boundaryConditionType> bndCon(commonConditions);
		std::replace(bndCon.begin(), bndCon.end(),
				boundaryConditionType::calculated,
				boundaryConditionType::calculatedCvM);
		CvM = volumeField<scalar>(mesh_, 0, subPatchData<scalar> { bndCon[0] },
				subPatchData<scalar> { bndCon[1] }, subPatchData<scalar> {
						bndCon[2] }, subPatchData<scalar> { bndCon[3] },
				subPatchData<scalar> { bndCon[4] }, subPatchData<scalar> {
						bndCon[5] });
	}

	diffusiveFields diffFieldsOld { mesh_, gasPhase, commonConditions,
			gasPhase.turbulence->turbulence(), gasPhase.turbulence->aField(),
			gasPhase.turbulence->bField() };
	diffusiveFields diffFieldsNew(diffFieldsOld), diffFieldsCur(diffFieldsOld);

	effectiveTransportCoefficients<quadraticSurface> effectiveCoeffs { mesh_,
			gasPhase.phaseThermodynamics->Mv().size() };

	concentrationsPack<cubicCell> calculatedConcentration { mesh_,
			gasPhase.phaseThermodynamics->Mv().size() };

	std::size_t nonLinearIteration(1);

	while (true)
	{
		const auto isNotFirstIter = nonLinearIteration > 1;

		/*Creating empty objects for SLE matrix.*/
		std::vector<SLEMatrix> massFractionMatrix {
				diffFieldsOld.massFraction.size(), SLEMatrix(
						std::string("Mass fraction")) };
		SLEMatrix velocityMatrix { std::string("Velocity") };
		SLEMatrix temperatureMatrix { std::string("Temperature") };
		SLEMatrix kMatrix { std::string("k") };
		SLEMatrix epsMatrix { std::string("epsilon") };
		SLEMatrix aMatrix { std::string("a") };
		SLEMatrix bMatrix { std::string("b") };
		/**/

		/*Turbulent coefficients recalculation.*/
		gasPhase.calculateCoefficients(*gasPhase.transportModel,
				gasPhase.phaseThermodynamics->Mv(), diffFieldsCur.temperature,
				gasPhase.pressure, gasPhase.concentration);
		if (gasPhase.turbulence->turbulence())
			gasPhase.calculateCoefficients(diffFieldsCur.k, diffFieldsCur.eps,
					*(gasPhase.turbulence));

		/*Calculation of fields for SLE matrix calculation*/
		nonIdealCorrectionOld.val() = gasPhase.phaseThermodynamics->nonIdeality(
				gasPhase.concentration.p, diffFieldsCur.temperature.cval());

		avMolMass.val() =
				(gasPhase.density[0] / gasPhase.concentration.v[0]).cval();

		CVOld.val() = gasPhase.phaseThermodynamics->Cv(
				gasPhase.concentration.p);

		CvM.val() = (CVOld / avMolMass).cval();

		const auto surfaceCv = linearInterpolate(CVOld, bncCalc);

		CCVOld.val() = (gasPhase.concentration.v[0] * CVOld).cval();

		concentrationsPack<quadraticSurface> surfaceConcentration { mesh_,
				gasPhase.phaseThermodynamics->Mv().size() };

		surfaceRho.val() = 0;

		for (std::size_t k = 1; k < surfaceConcentration.v.size(); ++k)
		{
			surfaceConcentration.v[k].val() = linearInterpolate(
					gasPhase.concentration.v[k], bncCalc).cval();

			surfaceConcentration.v[0] += surfaceConcentration.v[k];

			surfaceRho += surfaceConcentration.v[k]
					* gasPhase.phaseThermodynamics->Mv()[k - 1];
		}

		NonIdRho.val() = (nonIdealCorrectionOld / gasPhase.density[0]).cval();

		const auto gradNonIdRho = surfGrad(NonIdRho, bncCalc);

		const auto gradCvM = surfGrad(CvM, bncCalc);

		volumeField<tensor> gradV { mesh_ };
		volumeField<vector> gradRho { mesh_ };
		volumeField<vector> gradP { mesh_ };
		volumeField<scalar> divV { mesh_ };

		if (linearRec || isNotFirstIter)
		{
			gradV = grad(diffFieldsCur.velocity, bncCalc);
			gradRho = grad(gasPhase.density[0], bncCalc);
			gradP = grad(gasPhase.pressure, bncCalc);
			divV = divergence(diffFieldsCur.velocity, bncCalc);
		}
		else
		{
			gradV = grad(star.v);
			gradRho = grad(star.rho);
			gradP = grad(star.p);
			divV = divergence(star.v);
		}

		auto gradVTr = gradV;
		for (auto & gradVTr_i : gradVTr.val())
			gradVTr_i.transpose();

		const auto gradV_s = surfGrad(diffFieldsCur.velocity, bncCalc);

		auto gradVTr_s = gradV_s;
		for (auto & gradVTr_s_i : gradVTr_s.val())
			gradVTr_s_i.transpose();

		const auto divV_s = surfDivergence(diffFieldsCur.velocity, bncCalc);

		const auto surfaceTemperature = linearInterpolate(
				diffFieldsCur.temperature, bncCalc);

		auto strainVol = gradV + gradVTr;
		auto strainSurf = gradV_s + gradVTr_s;

		auto strainSurf1 = gradVTr_s;

		for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		{
			const scalar lambdadivV_i { twothirds * divV.cval()[i] };

			strainVol.val()[i].wr()[0] -= lambdadivV_i;
			strainVol.val()[i].wr()[4] -= lambdadivV_i;
			strainVol.val()[i].wr()[8] -= lambdadivV_i;
		}

		for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
		{
			const scalar lambdadivV_s_i { twothirds * divV_s.cval()[i] };

			strainSurf.val()[i].wr()[0] -= lambdadivV_s_i;
			strainSurf.val()[i].wr()[4] -= lambdadivV_s_i;
			strainSurf.val()[i].wr()[8] -= lambdadivV_s_i;

			strainSurf1.val()[i].wr()[0] -= lambdadivV_s_i;
			strainSurf1.val()[i].wr()[4] -= lambdadivV_s_i;
			strainSurf1.val()[i].wr()[8] -= lambdadivV_s_i;
		}

		/*Calculation of all diffusion coefficients.*/
		effectiveCoeffs.calculateCoefficients(*gasPhase.transportModel,
				gasPhase.phaseThermodynamics->Mv(), surfaceTemperature, star.p,
				surfaceConcentration);

		volumeField<scalar> thetaS_R { mesh_ };
		if (gasPhase.turbulence->turbulence())
		{
			const auto k_surf = linearInterpolate(diffFieldsCur.k, bncCalc);
			const auto eps_surf = linearInterpolate(diffFieldsCur.eps, bncCalc);

			effectiveCoeffs.calculateCoefficients(k_surf, eps_surf,
					*gasPhase.turbulence);

			const auto thetaS_s = gasPhase.turbulence->thetaS_D(divV_s, k_surf,
					eps_surf);

			effectiveCoeffs.k_D *= thetaS_s;
			effectiveCoeffs.eps_D *= thetaS_s;

			thetaS_R = gasPhase.turbulence->thetaS_R(divV, diffFieldsCur.k,
					diffFieldsCur.eps);
		}

		effectiveCoeffs.calculateEffectiveCoefficients(surfaceRho,
				*(gasPhase.turbulence), surfaceConcentration.v[0], surfaceCv);

		if (!isNotFirstIter)
			gasPhase.turbulence->checkTransitionToTurbulenceModel(
					gasPhase.physMu / gasPhase.density[0],
					effectiveCoeffs.physMu / surfaceRho, diffFieldsOld.k,
					diffFieldsOld.eps, diffFieldsOld.a, diffFieldsOld.b,
					gasPhase.concentration, bncCalc, timestep);

		if (gasPhase.turbulence->turbulence())
		{
			if (msolver.solverType != matrixSolver::explicitSolver)
				kMatrix.generateDTimeLaplacian(diffFieldsCur.k,
						diffFieldsOld.k * gasPhase.density[0],
						gasPhase.density[0], effectiveCoeffs.rhoDk, timestep,
						bncCalc);
			else
				kMatrix.generateDTimeExplicitLaplacian(diffFieldsCur.k,
						diffFieldsOld.k * gasPhase.density[0],
						gasPhase.density[0], effectiveCoeffs.rhoDk, bncCalc);

			if (msolver.solverType != matrixSolver::explicitSolver)
				epsMatrix.generateDTimeLaplacian(diffFieldsCur.eps,
						diffFieldsOld.eps * gasPhase.density[0],
						gasPhase.density[0], effectiveCoeffs.rhoDeps, timestep,
						bncCalc);
			else
				epsMatrix.generateDTimeExplicitLaplacian(diffFieldsCur.eps,
						diffFieldsOld.eps * gasPhase.density[0],
						gasPhase.density[0], effectiveCoeffs.rhoDeps, bncCalc);

			if (gasPhase.turbulence->aField())
			{
				if (msolver.solverType != matrixSolver::explicitSolver)
					aMatrix.generateDTimeLaplacian(diffFieldsCur.a,
							gasPhase.density[0],
							gasPhase.density[0] * diffFieldsOld.a,
							effectiveCoeffs.rhoDa, timestep, bncCalc);
				else
					aMatrix.generateDTimeExplicitLaplacian(diffFieldsCur.a,
							gasPhase.density[0],
							gasPhase.density[0] * diffFieldsOld.a,
							effectiveCoeffs.rhoDa, bncCalc);

				if (gasPhase.turbulence->bField())
				{
					if (msolver.solverType != matrixSolver::explicitSolver)
						bMatrix.generateDTimeLaplacian(diffFieldsCur.b,
								diffFieldsOld.b * gasPhase.density[0],
								gasPhase.density[0], effectiveCoeffs.rhoDb,
								timestep, bncCalc);
					else
						bMatrix.generateDTimeExplicitLaplacian(diffFieldsCur.b,
								diffFieldsOld.b * gasPhase.density[0],
								gasPhase.density[0], effectiveCoeffs.rhoDb,
								bncCalc);
				}
			}
		}

		/*Generate SLE matrix for each component and calculate explicit diffusion fields.*/
		for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
		{
			if (msolver.solverType != matrixSolver::explicitSolver)
				massFractionMatrix[k].generateDTimeLaplacian(
						diffFieldsCur.massFraction[k],
						diffFieldsOld.massFraction[k] * gasPhase.density[0],
						gasPhase.density[0], effectiveCoeffs.rhoD[k], timestep,
						bncCalc, k + 1);
			else
				massFractionMatrix[k].generateDTimeExplicitLaplacian(
						diffFieldsCur.massFraction[k],
						diffFieldsOld.massFraction[k] * gasPhase.density[0],
						gasPhase.density[0], effectiveCoeffs.rhoD[k], bncCalc,
						k + 1);
		}

		std::vector<volumeField<scalar>> cellMolFrac(massFractionMatrix.size(),
				volumeField<scalar>(mesh_, 0));
		for (std::size_t k = 1; k < gasPhase.concentration.v.size(); ++k)
		{
			std::vector<boundaryConditionType> bndCon(commonConditions);
			std::replace(bndCon.begin(), bndCon.end(),
					boundaryConditionType::calculated,
					boundaryConditionType::calculatedMolarFraction);
			cellMolFrac[k - 1] = volumeField<scalar>(mesh_, 0,
					subPatchData<scalar> { bndCon[0] }, subPatchData<scalar> {
							bndCon[1] }, subPatchData<scalar> { bndCon[2] },
					subPatchData<scalar> { bndCon[3] }, subPatchData<scalar> {
							bndCon[4] }, subPatchData<scalar> { bndCon[5] });

			cellMolFrac[k - 1].val() = (gasPhase.concentration.v[k]
					/ gasPhase.concentration.v[0]).cval();
		}

		std::vector<surfaceField<vector>> explGradW { massFractionMatrix.size(),
				surfaceField<vector>(mesh_, vector(0)) };
		std::vector<surfaceField<vector>> explGradX { massFractionMatrix.size(),
				surfaceField<vector>(mesh_, vector(0)) };

		surfaceField<vector> entExplDiffFlow { mesh_, vector(0) };

		if (molMassDiffusionFlag
				|| (enthalpySolverFlag != enthalpyFlow::noSolve))
		{
			for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
				explGradW[k] = surfGrad(diffFieldsCur.massFraction[k], bncCalc,
						k + 1);

			if (molMassDiffusionFlag)
			{
				for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
					explGradX[k] = surfGrad(cellMolFrac[k], bncCalc, k + 1);

				effectiveCoeffs.caclulateDFluxes(explGradX,
						gasPhase.phaseThermodynamics->Mv(),
						surfaceConcentration,
						gasPhase.turbulence->turbulence());
			}

			for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
			{
				const auto explDiffFlowW_k = -1 * effectiveCoeffs.rhoD[k]
						* explGradW[k];

				/*Calculating molar mass correction for diffusion flow and explicit enthalpy flow.*/
				surfaceField<vector> massFrCorr_k { mesh_, vector(0) };
				if (molMassDiffusionFlag)
				{
					surfaceField<vector> explDiffFlowX_k(mesh_, vector(0));

					for (std::size_t i = 0; i < explDiffFlowX_k.size(); ++i)
						explDiffFlowX_k.val()[i] =
								effectiveCoeffs.DFlux.cval()[i][k];

					massFrCorr_k = explDiffFlowX_k - explDiffFlowW_k;

					const auto divMolMassFlow = divergence(massFrCorr_k);

					if constexpr (distributionOn)
						massFractionMatrix[k].distributeSourceTerm(
								divMolMassFlow * (-1.0),
								diffFieldsCur.massFraction[k]);
					else
						massFractionMatrix[k].SLE[0].freeTerm -=
								divMolMassFlow.cval();
				}

				if (enthalpySolverFlag != enthalpyFlow::noSolve)
				{
					/*Calculate enthalpy per mole per Kelvin (isobaric molar heat capacity) for k component*/
					surfaceField<scalar> h_k { mesh_, 0 };
					h_k.val() = gasPhase.phaseThermodynamics->hkT(
							surfaceConcentration.v[k + 1].cval(),
							surfaceTemperature.cval(), k);

					entExplDiffFlow += (massFrCorr_k + explDiffFlowW_k) * h_k
							/ gasPhase.phaseThermodynamics->Mv()[k];
				}
			}
		}

		/*Calculation of resulting nonideality correction in laplacian*/
		const auto nonIdFlow = -1
				* (gradCvM * surfaceTemperature + gradNonIdRho);
		const auto nonIdealityCorrectionLaplacian = nonIdFlow * surfaceRho
				* effectiveCoeffs.tLambda;

		/*Calculation of new mass fraction*/
		{
			volumeField<scalar> sumMassFrac(mesh_, 0.);

			for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
			{
				if (msolver.solverType != matrixSolver::explicitSolver)
					diffFieldsNew.massFraction[k].val() = msolver.solve(
							diffFieldsCur.massFraction[k].cval(),
							massFractionMatrix[k]);
				else
					diffFieldsNew.massFraction[k].val() =
							(massFractionMatrix[k].SLE[0].explOldTime
									+ massFractionMatrix[k].SLE[0].freeTerm
											* timestep)
									/ massFractionMatrix[k].SLE[0].centralDiagonale;

				std::replace_if(std::begin(diffFieldsNew.massFraction[k].val()),
						std::end(diffFieldsNew.massFraction[k].val()),
						[](const scalar value) 
						{
							return value < 0.;
						}, 0.0);

				sumMassFrac += diffFieldsNew.massFraction[k];
			}

			for (auto & newMassFrac_k : diffFieldsNew.massFraction)
				newMassFrac_k /= sumMassFrac;
		}

		/*Calculation of new concentrations*/
		calculatedConcentration.v[0].val() = 0;
		for (std::size_t k = 1; k < calculatedConcentration.v.size(); ++k)
		{
			calculatedConcentration.v[k].boundCond_wr() =
					gasPhase.concentration.v[k].boundCond();

			calculatedConcentration.v[k].val() = (diffFieldsNew.massFraction[k
					- 1].cval() * gasPhase.density[0].cval()
					/ gasPhase.phaseThermodynamics->Mv()[k - 1]);

			calculatedConcentration.v[0] += calculatedConcentration.v[k];
		}

		CCVNew.val() = gasPhase.phaseThermodynamics->Cv(
				calculatedConcentration.p)
				* calculatedConcentration.v[0].cval();
		nonIdealCorrectionNew.val() = gasPhase.phaseThermodynamics->nonIdeality(
				calculatedConcentration.p, diffFieldsCur.temperature.cval());

		const auto oldEnergyField = CCVOld * diffFieldsOld.temperature
				+ nonIdealCorrectionOld - nonIdealCorrectionNew;

		/*Generate SLE matrix for components of the velocity vector.*/
		if (msolver.solverType != matrixSolver::explicitSolver)
			velocityMatrix.generateDTimeLaplacian(diffFieldsCur.velocity,
					gasPhase.density[0],
					gasPhase.density[0] * diffFieldsOld.velocity,
					effectiveCoeffs.mu, timestep, bncCalc);
		else
			velocityMatrix.generateDTimeExplicitLaplacian(
					diffFieldsCur.velocity, gasPhase.density[0],
					gasPhase.density[0] * diffFieldsOld.velocity,
					effectiveCoeffs.mu, bncCalc);

		if (msolver.solverType != matrixSolver::explicitSolver)
			temperatureMatrix.generateDTimeLaplacian(diffFieldsCur.temperature,
					oldEnergyField, CCVNew, effectiveCoeffs.kappa, timestep,
					bncCalc);
		else
			temperatureMatrix.generateDTimeExplicitLaplacian(
					diffFieldsCur.temperature, oldEnergyField, CCVNew,
					effectiveCoeffs.kappa, bncCalc);

		/*Explicit deviatoric part of molecular viscosity tensor for energy computation.*/
		const auto devPhysVisc = strainVol * gasPhase.physMu;

		/*Deviatoric part of effective viscosity tensor minus velocity gradient for velocity matrix.*/
		const auto devTotVisc1Surf = strainSurf1 * effectiveCoeffs.mu;

		volumeField<tensor> grada { mesh_, tensor(0) };
		volumeField<scalar> diva { mesh_, 0 };
		volumeField<vector> gradb { mesh_, vector(0) };
		surfaceField<tensor> devPhysViscSurf { mesh_, tensor(0) };

		if (gasPhase.turbulence->turbulence())
		{
			switch (gasPhase.turbulence->model())
			{
			case turbulenceModel::kEpsAModel:
			{
				const auto & kEpsA =
						dynamic_cast<kEpsAModel&>(*gasPhase.turbulence);

				diffFieldsCur.b.val() = kEpsA.calculate_b(mesh_,
						gasPhase.concentration.v, gasPhase.density).cval();
				diffFieldsNew.b.val() = diffFieldsCur.b.cval();

				if (linearRec || isNotFirstIter)
				{
					grada = grad(diffFieldsCur.a, bncCalc);
					diva = divergence(diffFieldsCur.a, bncCalc);
				}
				else
				{
					grada = grad(star.a);
					diva = divergence(star.a);
				}
			}
				break;
			case turbulenceModel::arithmeticA1Model:
			case turbulenceModel::arithmeticA2Model:
			case turbulenceModel::arithmeticA3Model:
			{
				const auto & arithmetic =
						dynamic_cast<arithmeticAModel&>(*gasPhase.turbulence);

				diffFieldsCur.a.val() = arithmetic.calculate_a(
						gasPhase.turbulence->model(), mesh_, gasPhase, gradRho,
						gradP).cval();

				grada = grad(diffFieldsCur.a, bncCalc);
				diva = divergence(diffFieldsCur.a, bncCalc);
			}
				break;
			case turbulenceModel::BHRModel:
			case turbulenceModel::BHRKLModel:
			{
				if (linearRec || isNotFirstIter)
				{
					grada = grad(diffFieldsCur.a, bncCalc);
					diva = divergence(diffFieldsCur.a, bncCalc);
					gradb = grad(diffFieldsCur.b, bncCalc);
				}
				else
				{
					grada = grad(star.a);
					diva = divergence(star.a);
					gradb = grad(star.b);
				}
			}
				break;
			default:
				break;
			}

			devPhysViscSurf = strainSurf * effectiveCoeffs.physMu;
		}

		/*Reynolds tensor's components for turbulence model.*/
		volumeField<tensor> spherTurbR { mesh_, tensor(0) };
		volumeField<tensor> devTurbR { mesh_, tensor(0) };
		/*Part of Reynolds tensor, which converts kinetic energy directly into internal.*/
		volumeField<tensor> philtTurbVisc { mesh_, tensor(0) };

		if (gasPhase.turbulence->turbulence())
		{
			const auto tNu1 = (1. - thetaS_R) * gasPhase.tNu;

			philtTurbVisc = strainVol * gasPhase.density[0] * tNu1;

			devTurbR = strainVol * gasPhase.density[0] * gasPhase.tNu
					* thetaS_R;

			/*Calculate turbulent pressure part*/
			spherTurbR = tensor(1, 0, 0, 0, 1, 0, 0, 0, 1) * -twothirds
					* gasPhase.rhokTurb;
		}

		/*Accounting of explicit part of effective viscosity tensor*/
		const auto divDevTotVisc1 = divergence(devTotVisc1Surf);
		for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
			for (std::size_t j = 0; j < vector::vsize; ++j)
				velocityMatrix.SLE[j].freeTerm[i] +=
						divDevTotVisc1.cval()[i]()[j];

		const auto divNonIdealityCorrectionLaplacian = divergence(
				nonIdealityCorrectionLaplacian);

		if constexpr (distributionOn)
			temperatureMatrix.distributeSourceTerm(
					divNonIdealityCorrectionLaplacian * (-1.0),
					diffFieldsCur.temperature);
		else
			temperatureMatrix.SLE[0].freeTerm -=
					divNonIdealityCorrectionLaplacian.cval();

		if ((msolver.solverType != matrixSolver::explicitSolver)
				&& (enthalpySolverFlag != enthalpyFlow::noSolve))
		{
			surfaceField<vector> entExplDiffFlowAfDif { mesh_, vector(0) };

			concentrationsPack<quadraticSurface> surfaceConcentrationAfDif {
					mesh_, gasPhase.phaseThermodynamics->Mv().size() };

			for (std::size_t k = 1; k < surfaceConcentrationAfDif.v.size(); ++k)
			{
				surfaceConcentrationAfDif.v[k] = linearInterpolate(
						calculatedConcentration.v[k], bncCalc);

				surfaceConcentrationAfDif.v[0] +=
						surfaceConcentrationAfDif.v[k];
			}

			for (std::size_t k = 1; k < calculatedConcentration.v.size(); ++k)
				cellMolFrac[k - 1].val() = (calculatedConcentration.v[k]
						/ calculatedConcentration.v[0]).cval();

			for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
				explGradW[k] = surfGrad(diffFieldsCur.massFraction[k], bncCalc,
						k + 1);

			if (molMassDiffusionFlag)
			{
				for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
					explGradX[k] = surfGrad(cellMolFrac[k], bncCalc, k + 1);

				effectiveCoeffs.caclulateDFluxes(explGradX,
						gasPhase.phaseThermodynamics->Mv(),
						surfaceConcentrationAfDif,
						gasPhase.turbulence->turbulence());
			}

			for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
			{
				const auto explDiffFlowW_k = -1 * effectiveCoeffs.rhoD[k]
						* explGradW[k];

				/*Calculating molar mass correction for diffusion flow and explicit enthalpy flow.*/
				surfaceField<vector> massFrCorr_k { mesh_, vector(0) };
				if (molMassDiffusionFlag)
				{
					surfaceField<vector> explDiffFlowX_k(mesh_, vector(0));

					for (std::size_t i = 0; i < explDiffFlowX_k.size(); ++i)
						explDiffFlowX_k.val()[i] =
								effectiveCoeffs.DFlux.cval()[i][k];

					massFrCorr_k = explDiffFlowX_k - explDiffFlowW_k;
				}

				/*Calculate enthalpy per mole per Kelvin (isobaric molar heat capacity) for k component*/
				surfaceField<scalar> h_k { mesh_, 0 };
				h_k.val() = gasPhase.phaseThermodynamics->hkT(
						surfaceConcentrationAfDif.v[k + 1].cval(),
						surfaceTemperature.cval(), k);

				entExplDiffFlowAfDif += (massFrCorr_k + explDiffFlowW_k) * h_k
						/ gasPhase.phaseThermodynamics->Mv()[k];
			}

			constexpr scalar we = 0.5;

			entExplDiffFlow = (entExplDiffFlowAfDif * we)
					+ (entExplDiffFlow * (1 - we));
		}

		if (enthalpySolverFlag == enthalpyFlow::explicitSolve)
		{
			entExplDiffFlow *= surfaceTemperature;

			const auto divEnthFlow = divergence(entExplDiffFlow);

			auto & nonConstMesh = const_cast<mesh&>(mesh_);

			nonConstMesh.timestepSourceRef() =
					timestepCoeffs.first
							* (CCVOld.cval() * diffFieldsOld.temperature.cval()
									/ std::abs(
											divEnthFlow.cval() + stabilizator)).min();

			if constexpr (distributionOn)
				temperatureMatrix.distributeSourceTerm(divEnthFlow * -1,
						diffFieldsCur.temperature);
			else
				temperatureMatrix.SLE[0].freeTerm -= divEnthFlow.cval();
		}

		temperatureMatrix.SLE[0].freeTerm += (devPhysVisc && gradV).cval();

		/*Calculating turbulent sources.*/
		volumeField<scalar> sigmaSourcek { mesh_, scalar { 0 } };
		volumeField<scalar> sigmaSourceeps { mesh_, 0 };
		volumeField<vector> sigmaSourcea { mesh_, vector(0) };
		volumeField<scalar> sigmaSourceb { mesh_, 0 };

		volumeField<vector> gradMav_Mav { mesh_, vector { 0 } };

		if (gasPhase.turbulence->turbulence())
		{
			if (gasPhase.turbulence->model() != turbulenceModel::zeroModel)
				temperatureMatrix.SLE[0].freeTerm +=
						gasPhase.turbulence->rhoepsilon(gasPhase,
								*gasPhase.phaseThermodynamics, diffFieldsCur.k,
								diffFieldsCur.eps);

			if (gasPhase.turbulence->aField())
				gradMav_Mav = grad(surfaceRho / surfaceConcentration.v[0])
						/ avMolMass;

			temperatureMatrix.SLE[0].freeTerm +=
					(philtTurbVisc && gradV).cval();

			auto & nonConstMesh = const_cast<mesh&>(mesh_);

			const auto [SourcekSuSp, SourceepsSuSp, SourceaSuSp, SourcebSuSp,
					gravEnSink] = gasPhase.turbulence->calculate(
					nonConstMesh.timestepSourceRef(), timestepCoeffs.first,
					gasPhase, diffFieldsCur, gradV, divergence(devPhysViscSurf),
					gradP, gradRho, grada, diva, gradb, spherTurbR, devTurbR,
					gradMav_Mav, *(gasPhase.phaseThermodynamics), gasPhase.tNu);

			if (msolver.solverType == matrixSolver::explicitSolver)
			{
				sigmaSourcek = SourcekSuSp.first
						+ SourcekSuSp.second * diffFieldsCur.k;
				sigmaSourceeps = SourceepsSuSp.first
						+ SourceepsSuSp.second * diffFieldsCur.eps;

				for (std::size_t i = 0; i < sigmaSourcea.size(); ++i)
				{
					const vector vec = vector(
							SourceaSuSp.second.cval()[i]()[0]
									* diffFieldsCur.a.cval()[i]()[0],
							SourceaSuSp.second.cval()[i]()[1]
									* diffFieldsCur.a.cval()[i]()[1],
							SourceaSuSp.second.cval()[i]()[2]
									* diffFieldsCur.a.cval()[i]()[2]);

					sigmaSourcea.val()[i] = SourceaSuSp.first.cval()[i] + vec;
				}
			}
			else
			{
				kMatrix.freeSourceTerm(SourcekSuSp.first);
				kMatrix.diagonaleSourceTerm(SourcekSuSp.second);

				epsMatrix.freeSourceTerm(SourceepsSuSp.first);
				epsMatrix.diagonaleSourceTerm(SourceepsSuSp.second);

				if (gasPhase.turbulence->aField())
				{
					aMatrix.freeSourceTerm(SourceaSuSp.first);
					aMatrix.diagonaleSourceTerm(SourceaSuSp.second);

					if (gasPhase.turbulence->bField())
					{
						bMatrix.freeSourceTerm(SourcebSuSp.first);
						bMatrix.diagonaleSourceTerm(SourcebSuSp.second);
					}
				}
			}

			switch (gasPhase.turbulence->model())
			{
			case turbulenceModel::kEpsAModel:
			case turbulenceModel::BHRModel:
			case turbulenceModel::BHRKLModel:
			case turbulenceModel::arithmeticA1Model:
			case turbulenceModel::arithmeticA2Model:
			case turbulenceModel::arithmeticA3Model:
			{
				if constexpr (distributionOn)
					temperatureMatrix.distributeSourceTerm(gravEnSink * (-1.0),
							diffFieldsCur.temperature);
				else
					temperatureMatrix.SLE[0].freeTerm -= gravEnSink.cval();
			}
				break;
			default:
				break;
			}
		}

		/*Calculation of new velocity by components of vector*/
		if (msolver.solverType != matrixSolver::explicitSolver)
			diffFieldsNew.velocity.val() = msolver.solve(
					diffFieldsCur.velocity.cval(), velocityMatrix);
		else
			for (std::size_t j = 0; j < vector::vsize; ++j)
				for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
					diffFieldsNew.velocity.val()[i].wr()[j] =
							(velocityMatrix.SLE[j].explOldTime[i]
									+ velocityMatrix.SLE[j].freeTerm[i]
											* timestep)
									/ velocityMatrix.SLE[j].centralDiagonale[i];

		/*Calculation of new temperature*/
		if (enthalpySolverFlag == enthalpyFlow::implicitSolve)
		{
			temperatureMatrix.addNabla(diffFieldsCur.temperature,
					entExplDiffFlow, bncCalc);

			if (msolver.solverType != matrixSolver::explicitSolver)
				diffFieldsNew.temperature.val() = msolver.solve(
						diffFieldsCur.temperature.cval(), temperatureMatrix);
			else
				diffFieldsNew.temperature.val() = msolverForEnthalpy.solve(
						diffFieldsCur.temperature.cval(), temperatureMatrix);
		}
		else
		{
			if (msolver.solverType != matrixSolver::explicitSolver)
				diffFieldsNew.temperature.val() = msolver.solve(
						diffFieldsCur.temperature.cval(), temperatureMatrix);
			else
				diffFieldsNew.temperature.val() =
						(temperatureMatrix.SLE[0].explOldTime
								+ temperatureMatrix.SLE[0].freeTerm * timestep)
								/ temperatureMatrix.SLE[0].centralDiagonale;
		}

		if (gasPhase.turbulence->turbulence())
		{
			/*Explicit source integration*/
			/*Calculation of new k*/
			if (msolver.solverType != matrixSolver::explicitSolver)
				diffFieldsNew.k.val() = msolver.solve(diffFieldsCur.k.cval(),
						kMatrix);
			else
				diffFieldsNew.k.val() = (kMatrix.SLE[0].explOldTime
						+ kMatrix.SLE[0].freeTerm * timestep)
						/ kMatrix.SLE[0].centralDiagonale
						+ sigmaSourcek.cval() / gasPhase.density[0].cval()
								* timestep;

			/*Calculation of new epsilon*/
			if (msolver.solverType != matrixSolver::explicitSolver)
				diffFieldsNew.eps.val() = msolver.solve(
						diffFieldsCur.eps.cval(), epsMatrix);
			else
				diffFieldsNew.eps.val() = (epsMatrix.SLE[0].explOldTime
						+ epsMatrix.SLE[0].freeTerm * timestep)
						/ epsMatrix.SLE[0].centralDiagonale
						+ sigmaSourceeps.cval() / gasPhase.density[0].cval()
								* timestep;

			if (gasPhase.turbulence->aField())
			{
				/*Calculation of new a by components of vector*/
				if (msolver.solverType != matrixSolver::explicitSolver)
					diffFieldsNew.a.val() = msolver.solve(
							diffFieldsCur.a.cval(), aMatrix);
				else
				{
					for (std::size_t j = 0; j < vector::vsize; ++j)
						for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
							diffFieldsNew.a.val()[i].wr()[j] =
									(aMatrix.SLE[j].explOldTime[i]
											+ aMatrix.SLE[j].freeTerm[i]
													* timestep)
											/ aMatrix.SLE[j].centralDiagonale[i];

					diffFieldsNew.a.val() += (sigmaSourcea / gasPhase.density[0]
							* timestep).cval();
				}

				/*Calculation of new b*/
				if (gasPhase.turbulence->bField())
				{
					if (msolver.solverType != matrixSolver::explicitSolver)
						diffFieldsNew.b.val() = msolver.solve(
								diffFieldsCur.b.cval(), bMatrix);
					else
						diffFieldsNew.b.val() = (bMatrix.SLE[0].explOldTime
								+ bMatrix.SLE[0].freeTerm * timestep)
								/ bMatrix.SLE[0].centralDiagonale
								+ sigmaSourceb.cval()
										/ gasPhase.density[0].cval() * timestep;
				}
			}

			/*Check for non-negativity turbulent quantities.*/
			std::replace_if(std::begin(diffFieldsNew.k.val()),
					std::end(diffFieldsNew.k.val()),
					[&gasPhase](const scalar value) 
					{
						return value < gasPhase.turbulence->mink();
					}, gasPhase.turbulence->mink());

			std::replace_if(std::begin(diffFieldsNew.eps.val()),
					std::end(diffFieldsNew.eps.val()),
					[&gasPhase](const scalar value) 
					{
						return value < gasPhase.turbulence->minepsilon();
					}, gasPhase.turbulence->minepsilon());

			if (gasPhase.turbulence->bField())
				std::replace_if(std::begin(diffFieldsNew.b.val()),
						std::end(diffFieldsNew.b.val()),
						[&gasPhase](const scalar value) 
						{
							return value < gasPhase.turbulence->minb();
						}, gasPhase.turbulence->minb());
		}

		auto deltaEpsMax = std::abs(
				(diffFieldsNew.eps.cval() - diffFieldsCur.eps.cval())
						/ (diffFieldsOld.eps.cval() + stabilizator)).max();

#ifdef MPI_VERSION
		/*Calculate parallel timestep*/
		{
			const MPIHandler::mpi_scalar deltaEpsMax_arr[1] { deltaEpsMax };

			std::vector<MPIHandler::mpi_scalar> nodesDelta(0);
			if (parallelism.isRoot())
				nodesDelta.resize(parallelism.mpi_size);

			MPI_Gather(deltaEpsMax_arr, 1, schemi_MPI_SCALAR, nodesDelta.data(),
					1,
					schemi_MPI_SCALAR, parallelism.root,
					MPI_COMM_WORLD);

			MPIHandler::mpi_scalar maxOfAll[1] { -1.23 };
			if (parallelism.isRoot())
			{
				std::valarray<scalar> deltaEps_max(parallelism.mpi_size);

				for (std::size_t i = 0; i < deltaEps_max.size(); ++i)
					deltaEps_max[i] = nodesDelta[i];

				maxOfAll[0] = deltaEps_max.max();
			}
			MPI_Bcast(maxOfAll, 1, schemi_MPI_SCALAR, parallelism.root,
			MPI_COMM_WORLD);

			deltaEpsMax = maxOfAll[0];
		}
#endif

		if (deltaEpsMax < convergenceToleranceGlobal * 1E4
				|| nonLinearIteration > 100 * nonLinearityIteratonsFlag)
			break;
		else
		{
			diffFieldsCur = diffFieldsNew;

			nonLinearIteration++;
		}
	}

	/*Recalculation of values in cells*/
	gasPhase.concentration.v[0].val() = calculatedConcentration.v[0].cval();
	for (std::size_t k = 1; k < gasPhase.density.size(); ++k)
	{
		gasPhase.concentration.v[k].val() = calculatedConcentration.v[k].cval();

		gasPhase.density[k].val() = (gasPhase.concentration.v[k]
				* gasPhase.phaseThermodynamics->Mv()[k - 1]).cval();
	}

	gasPhase.velocity.val() = diffFieldsNew.velocity.cval();

	gasPhase.momentum.val() = (gasPhase.velocity * gasPhase.density[0]).cval();

	gasPhase.temperature.val() = diffFieldsNew.temperature.cval();

	gasPhase.internalEnergy.val() = gasPhase.phaseThermodynamics->UvcFromT(
			gasPhase.concentration.p, gasPhase.temperature.cval())
			* gasPhase.concentration.v[0].cval();

	gasPhase.pressure.val() = gasPhase.phaseThermodynamics->pFromUv(
			gasPhase.concentration.p, gasPhase.internalEnergy.cval());

	if (gasPhase.turbulence->turbulence())
	{
		gasPhase.kTurb.val() = diffFieldsNew.k.cval();
		gasPhase.epsTurb.val() = diffFieldsNew.eps.cval();
		if (gasPhase.turbulence->aField())
		{
			gasPhase.aTurb.val() = diffFieldsNew.a.cval();
			gasPhase.bTurb.val() = diffFieldsNew.b.cval();
		}

		gasPhase.rhokTurb.val() = (gasPhase.kTurb * gasPhase.density[0]).cval();
		gasPhase.rhoepsTurb.val() = gasPhase.epsTurb.cval()
				* gasPhase.density[0].cval();
		if (gasPhase.turbulence->aField())
		{
			gasPhase.rhoaTurb.val() =
					(gasPhase.aTurb * gasPhase.density[0]).cval();

			if (gasPhase.turbulence->bField())
				gasPhase.rhobTurb.val() =
						(gasPhase.bTurb * gasPhase.density[0]).cval();
		}
	}

	{
		const auto v2 = gasPhase.velocity & gasPhase.velocity;

		gasPhase.totalEnergy.val() = (gasPhase.internalEnergy
				+ gasPhase.density[0] * v2 * 0.5 + gasPhase.rhokTurb).cval();
	}

	gasPhase.HelmholtzEnergy.val() = gasPhase.phaseThermodynamics->Fv(
			gasPhase.concentration.p, gasPhase.temperature.cval());

	gasPhase.entropy.val() = gasPhase.phaseThermodynamics->Sv(
			gasPhase.concentration.p, gasPhase.temperature.cval());

	/*Diffusion time-step.*/
	if (sourceTimeFlag == timestep::CourantAndSourceAndDiffusionTimeStep)
	{
		const auto maxVal = effectiveCoeffs.maxValue(
				gasPhase.turbulence->turbulence(),
				gasPhase.turbulence->aField(), gasPhase.turbulence->bField(),
				surfaceRho);

		const scalar diffuse_dt = 0.5 * timestepCoeffs.second
				/ (linearInterpolate(maxVal).cval() * minimalLengthScale.cval()
						* minimalLengthScale.cval()).min();

		auto & nonConstMesh = const_cast<mesh&>(mesh_);

		nonConstMesh.timestepSourceRef() = std::min(
				nonConstMesh.timestepSource(), diffuse_dt);
	}

	if (gasPhase.temperature.cval().min() < 0.)
		throw exception("Negative temperature after diffusion.",
				errors::negativeTemperatureError);

	const auto DiffusionEndTime { std::chrono::high_resolution_clock::now() };
	timeForDiffusion += std::chrono::duration_cast<std::chrono::milliseconds>(
			DiffusionEndTime - DiffusionStartTime).count();
}
