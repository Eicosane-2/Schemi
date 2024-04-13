/*
 * Diffusion.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "Diffusion.hpp"

#include <chrono>

#include "harmonicInterpolateScalar.hpp"
#include "divergence.hpp"
#include "fieldProducts.hpp"
#include "gradient.hpp"
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
		const MPIHandler & parallelism, const timestep sourceTimeFlag,
		const bool molMassDiffusionFlag)
{
	const auto DiffusionStartTime { std::chrono::high_resolution_clock::now() };

	auto & mesh_ { gasPhase.pressure.meshRef() };

	const scalar timestep = mesh_.timestep();

	/*Fields for diffusion stage.*/
	surfaceField<vector> entExplDiffFlow { mesh_, vector(0) };
	volumeField<scalar> intermediateTemperature(gasPhase.temperature);

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
			gasPhase.turbulenceSources->turbulence,
			gasPhase.turbulenceSources->model }, diffFieldsNew(diffFieldsOld);

	effectiveTransportCoefficients<quadraticSurface> effectiveCoeffs { mesh_,
			gasPhase.phaseThermodynamics->Mv().size() };

	/*Creating empty objects for SLE matrix.*/
	std::vector<SLEMatrix> massFractionMatrix {
			diffFieldsOld.massFraction.size(), SLEMatrix(
					std::string("Mass fraction")) };
	SLEMatrix velocityMatrix { std::string("Velocity") };
	SLEMatrix temperatureMatrix { std::string("Temperature") };
	SLEMatrix temperatureEnthMatrix { std::string("Temperature-enthalpy") };
	SLEMatrix kMatrix { std::string("k") };
	SLEMatrix epsMatrix { std::string("epsilon") };
	SLEMatrix aMatrix { std::string("a") };
	SLEMatrix bMatrix { std::string("b") };
	/**/

	/*Turbulent coefficients recalculation.*/
	gasPhase.calculateCoefficients(*gasPhase.transportModel,
			gasPhase.phaseThermodynamics->Mv(), gasPhase.temperature,
			gasPhase.pressure, gasPhase.concentration);
	if (gasPhase.turbulenceSources->turbulence)
		gasPhase.calculateCoefficients(gasPhase.kTurb, gasPhase.epsTurb,
				*(gasPhase.turbulenceSources->turbPar));

	/*Calculation of fields for SLE matrix calculation*/
	nonIdealCorrectionOld.r() = gasPhase.phaseThermodynamics->nonIdeality(
			gasPhase.concentration.p, gasPhase.temperature());

	CVOld.r() = gasPhase.phaseThermodynamics->Cv(gasPhase.concentration.p);

	CvM.r() = CVOld() * gasPhase.concentration.v[0]() / gasPhase.density[0]();

	const auto surfaceCv = linearInterpolate(CVOld, bncCalc);

	CCVOld.r() = gasPhase.concentration.v[0]() * CVOld();

	const auto surfaceTemperature = linearInterpolate(gasPhase.temperature,
			bncCalc);

	concentrationsPack<quadraticSurface> surfaceConcentration { mesh_,
			gasPhase.phaseThermodynamics->Mv().size() };

	for (std::size_t k = 1; k < surfaceConcentration.v.size(); ++k)
	{
		surfaceConcentration.v[k] = linearInterpolate(
				gasPhase.concentration.v[k], bncCalc);

		surfaceConcentration.v[0].r() += surfaceConcentration.v[k]();

		surfaceRho.r() += surfaceConcentration.v[k]()
				* gasPhase.phaseThermodynamics->Mv()[k - 1];
	}

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
	avMolMass.r() = gasPhase.density[0]() / gasPhase.concentration.v[0]();

	const auto gradMav = surfGrad(avMolMass, bncCalc);

	NonIdRho.r() = nonIdealCorrectionOld() / gasPhase.density[0]();

	const auto gradNonIdRho = surfGrad(NonIdRho, bncCalc);

	const auto gradCvM = surfGrad(CvM, bncCalc);

	volumeField<tensor> gradV { mesh_ };
	volumeField<vector> gradRho { mesh_ };
	volumeField<vector> gradP { mesh_ };
	volumeField<scalar> divV { mesh_ };
	surfaceField<scalar> divV_s { mesh_ };
	volumeField<scalar> thetaS_R { mesh_ };

	if (linearRec)
	{
		gradV = grad(gasPhase.velocity, bncCalc);
		gradRho = grad(gasPhase.density[0], bncCalc);
		gradP = grad(gasPhase.pressure, bncCalc);
		divV = divergence(gasPhase.velocity, bncCalc);
	}
	else
	{
		gradV = grad(star.v);
		gradRho = grad(star.rho);
		gradP = grad(star.p);
		divV = divergence(star.v);
	}

	divV_s = surfDivergence(diffFieldsOld.velocity, bncCalc);

	/*Calculation of all diffusion coefficients.*/
	effectiveCoeffs.calculateCoefficients(*gasPhase.transportModel,
			gasPhase.phaseThermodynamics->Mv(), surfaceTemperature, star.p,
			surfaceConcentration);

	if (gasPhase.turbulenceSources->turbulence)
	{
		effectiveCoeffs.tNu = harmonicInterpolate(gasPhase.tNu, bncCalc);
		effectiveCoeffs.calculateCoefficients(
				*gasPhase.turbulenceSources->turbPar);

		const auto k_surf = linearInterpolate(gasPhase.kTurb, bncCalc);
		const auto eps_surf = linearInterpolate(gasPhase.epsTurb, bncCalc);

		const auto thetaS_s = gasPhase.turbulenceSources->turbPar->thetaS_D(
				divV_s, k_surf, eps_surf);

		effectiveCoeffs.k_D.r() *= thetaS_s();
		effectiveCoeffs.eps_D.r() *= thetaS_s();

		thetaS_R = gasPhase.turbulenceSources->turbPar->thetaS_R(divV,
				diffFieldsOld.k, diffFieldsOld.eps);
	}

	effectiveCoeffs.calculateEffectiveCoefficients(surfaceRho,
			*(gasPhase.turbulenceSources), surfaceConcentration.v[0],
			surfaceCv);

	if (gasPhase.turbulenceSources->turbulence)
	{
		if (msolver.solverType != matrixSolver::explicitSolver)
			kMatrix.generateDTimeLaplacian(diffFieldsOld.k,
					astProduct(diffFieldsOld.k, gasPhase.density[0]),
					gasPhase.density[0], effectiveCoeffs.rhoDk, timestep,
					bncCalc);
		else
			kMatrix.generateDTimeExplicitLaplacian(diffFieldsOld.k,
					astProduct(diffFieldsOld.k, gasPhase.density[0]),
					gasPhase.density[0], effectiveCoeffs.rhoDk, bncCalc);

		if (msolver.solverType != matrixSolver::explicitSolver)
			epsMatrix.generateDTimeLaplacian(diffFieldsOld.eps,
					astProduct(diffFieldsOld.eps, gasPhase.density[0]),
					gasPhase.density[0], effectiveCoeffs.rhoDeps, timestep,
					bncCalc);
		else
			epsMatrix.generateDTimeExplicitLaplacian(diffFieldsOld.eps,
					astProduct(diffFieldsOld.eps, gasPhase.density[0]),
					gasPhase.density[0], effectiveCoeffs.rhoDeps, bncCalc);

		if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::kEpsASource))
		{
			if (msolver.solverType != matrixSolver::explicitSolver)
				aMatrix.generateDTimeLaplacian(diffFieldsOld.a,
						gasPhase.density[0], effectiveCoeffs.rhoDa, timestep,
						bncCalc);
			else
				aMatrix.generateDTimeExplicitLaplacian(diffFieldsOld.a,
						gasPhase.density[0], effectiveCoeffs.rhoDa, bncCalc);

			if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::BHRKLSource))
			{
				if (msolver.solverType != matrixSolver::explicitSolver)
					bMatrix.generateDTimeLaplacian(diffFieldsOld.b,
							astProduct(diffFieldsOld.b, gasPhase.density[0]),
							gasPhase.density[0], effectiveCoeffs.rhoDb,
							timestep, bncCalc);
				else
					bMatrix.generateDTimeExplicitLaplacian(diffFieldsOld.b,
							astProduct(diffFieldsOld.b, gasPhase.density[0]),
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
					diffFieldsOld.massFraction[k],
					astProduct(diffFieldsOld.massFraction[k],
							gasPhase.density[0]), gasPhase.density[0],
					effectiveCoeffs.rhoD[k], timestep, bncCalc, k + 1);
		else
			massFractionMatrix[k].generateDTimeExplicitLaplacian(
					diffFieldsOld.massFraction[k],
					astProduct(diffFieldsOld.massFraction[k],
							gasPhase.density[0]), gasPhase.density[0],
					effectiveCoeffs.rhoD[k], bncCalc, k + 1);
	}

	std::vector<volumeField<vector>> molFracGrad(massFractionMatrix.size(),
			volumeField<vector>(mesh_, vector(0)));
	for (std::size_t k = 1; k < surfaceConcentration.v.size(); ++k)
	{
		surfaceField<scalar> surfMolFrac_k { mesh_, 0 };

		surfMolFrac_k.r() = surfaceConcentration.v[k]()
				/ surfaceConcentration.v[0]();
		molFracGrad[k - 1] = grad(surfMolFrac_k);
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

		cellMolFrac[k - 1].r() = gasPhase.concentration.v[k]()
				/ gasPhase.concentration.v[0]();
	}

	std::vector<surfaceField<vector>> explGradW { massFractionMatrix.size(),
			surfaceField<vector>(mesh_, vector(0)) };
	std::vector<surfaceField<vector>> explGradX { massFractionMatrix.size(),
			surfaceField<vector>(mesh_, vector(0)) };

	for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
	{
		explGradW[k] = surfGrad(diffFieldsOld.massFraction[k], bncCalc, k + 1);

		explGradX[k] = surfGrad(cellMolFrac[k], bncCalc, k + 1);
	}

	effectiveCoeffs.caclulateDFluxes(explGradX,
			gasPhase.phaseThermodynamics->Mv(), surfaceConcentration,
			gasPhase.turbulenceSources->turbulence);

	for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
	{
		auto explDiffFlowW_k = astProduct(effectiveCoeffs.rhoD[k],
				explGradW[k]);
		astProductSelf(explDiffFlowW_k, -1);

		/*Calculate enthalpy per mole per Kelvin (isobaric molar heat capacity) for k component*/
		surfaceField<scalar> h_k { mesh_, 0 };
		h_k.r() = gasPhase.phaseThermodynamics->hkT(
				surfaceConcentration.v[k + 1](), surfaceTemperature(), k);

		/*Calculating molar mass correction for diffusion flow and explicit enthalpy flow.*/
		surfaceField<vector> massFrCorr_k { mesh_, vector(0) };
		if (molMassDiffusionFlag)
		{
			surfaceField<vector> explDiffFlowX_k(mesh_, vector(0));

			for (std::size_t i = 0; i < explDiffFlowX_k.size(); ++i)
				explDiffFlowX_k.r()[i] = effectiveCoeffs.DFlux()[i][k];

			massFrCorr_k.r() = explDiffFlowX_k() - explDiffFlowW_k();
		}

		for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
		{
			const vector entExplDiffFlow_k_i = (massFrCorr_k()[i]
					+ explDiffFlowW_k()[i]) * h_k()[i]
					/ gasPhase.phaseThermodynamics->Mv()[k];

			entExplDiffFlow.r()[i] += entExplDiffFlow_k_i;
		}

		const auto divMolMassFlow = divergence(massFrCorr_k);

		massFractionMatrix[k].SLE[0].freeTerm -= divMolMassFlow();

		//massFractionMatrix[k].distributeSourceTerm(
		//		astProduct(divMolMassFlow, -1.0), diffFieldsOld.massFraction[k],
		//		timestep);
	}

	/*Calculation of resulting nonideality correction in laplacian*/
	surfaceField<vector> nonIdFlow { mesh_, vector(0) };

	nonIdFlow.r() = astProduct(gradCvM, surfaceTemperature)() + gradNonIdRho();
	astProductSelf(nonIdFlow, -1.0);
	const auto nonIdealityCorrectionLaplacian = astProduct(nonIdFlow,
			astProduct(surfaceRho, effectiveCoeffs.tLambda));

	/*Calculation of new mass fraction*/
	{
		volumeField<scalar> sumMassFrac(mesh_, 0.);

		for (std::size_t k = 0; k < massFractionMatrix.size(); ++k)
		{
			if (msolver.solverType != matrixSolver::explicitSolver)
				diffFieldsNew.massFraction[k].r() = msolver.solve(
						diffFieldsOld.massFraction[k](), massFractionMatrix[k]);
			else
				diffFieldsNew.massFraction[k].r() =
						(massFractionMatrix[k].SLE[0].explOldTime
								+ massFractionMatrix[k].SLE[0].freeTerm
										* timestep)
								/ massFractionMatrix[k].SLE[0].centralDiagonale;

			std::replace_if(std::begin(diffFieldsNew.massFraction[k].r()),
					std::end(diffFieldsNew.massFraction[k].r()),
					[](const scalar value) 
					{
						return value < 0.;
					}, 0.0);

			sumMassFrac.r() += diffFieldsNew.massFraction[k]();
		}

		for (auto & newMassFrac_k : diffFieldsNew.massFraction)
			newMassFrac_k.r() /= sumMassFrac();
	}

	/*Calculation of new concentrations*/
	concentrationsPack<cubicCell> intermediateConcentration { mesh_,
			gasPhase.phaseThermodynamics->Mv().size() };
	for (std::size_t k = 1; k < intermediateConcentration.v.size(); ++k)
	{
		intermediateConcentration.v[k].r() = diffFieldsNew.massFraction[k - 1]()
				* gasPhase.density[0]()
				/ gasPhase.phaseThermodynamics->Mv()[k - 1];

		intermediateConcentration.v[0].r() += intermediateConcentration.v[k]();
	}

	CCVNew.r() = gasPhase.phaseThermodynamics->Cv(intermediateConcentration.p)
			* intermediateConcentration.v[0]();
	nonIdealCorrectionNew.r() = gasPhase.phaseThermodynamics->nonIdeality(
			intermediateConcentration.p, gasPhase.temperature());

	/*Generate SLE matrix for components of the velocity vector.*/
	surfaceField<tensor> devPhysViscSurf { mesh_, tensor(0) };

	if (msolver.solverType != matrixSolver::explicitSolver)
		velocityMatrix.generateDTimeLaplacian(diffFieldsOld.velocity,
				gasPhase.density[0], effectiveCoeffs.mu, timestep, bncCalc);
	else
		velocityMatrix.generateDTimeExplicitLaplacian(diffFieldsOld.velocity,
				gasPhase.density[0], effectiveCoeffs.mu, bncCalc);

	volumeField<scalar> oldEnergyField(mesh_);
	oldEnergyField.r() = CCVOld() * diffFieldsOld.temperature()
			+ nonIdealCorrectionOld() - nonIdealCorrectionNew();

	if (msolver.solverType != matrixSolver::explicitSolver)
		temperatureMatrix.generateDTimeLaplacian(diffFieldsOld.temperature,
				oldEnergyField, CCVNew, effectiveCoeffs.kappa, timestep,
				bncCalc);
	else
		temperatureMatrix.generateDTimeExplicitLaplacian(
				diffFieldsOld.temperature, oldEnergyField, CCVNew,
				effectiveCoeffs.kappa, bncCalc);

	/*Deviatoric part of physical viscosity tensor, begins from velocity gradient.*/
	const auto gradV_s = surfGrad(gasPhase.velocity, bncCalc);
	surfaceField<tensor> devTotVisc1 = gradV_s;

	volumeField<tensor> grada { mesh_, tensor(0) };
	volumeField<scalar> diva { mesh_, 0 };
	volumeField<vector> gradb { mesh_, vector(0) };

	//volumeField<vector> divDeva1 { mesh_, vector(0) };

	if (gasPhase.turbulenceSources->turbulence)
	{
		switch (gasPhase.turbulenceSources->model)
		{
		case turbulenceModel::kEpsASource:
		{
			volumeField<scalar> bCalculated { mesh_, 0 };
			volumeField<scalar> numerator { mesh_, 0 };
			volumeField<scalar> denomenator { mesh_, 0 };
			for (std::size_t k = 1; k < gasPhase.concentration.v.size(); ++k)
			{
				std::valarray<scalar> x_k = gasPhase.concentration.v[k]()
						/ gasPhase.concentration.v[0]();

				numerator.r() += x_k
						/ (gasPhase.density[k]() + stabilizator
								+ gasPhase.turbulenceSources->turbPar->Cb1()
										* gasPhase.density[0]());

				denomenator.r() += x_k * gasPhase.density[k]()
						/ (gasPhase.density[k]() + stabilizator
								+ gasPhase.turbulenceSources->turbPar->Cb1()
										* gasPhase.density[0]());
			}

			bCalculated.r() = gasPhase.density[0]() * numerator()
					/ denomenator() - 1.;

			std::replace_if(std::begin(bCalculated.r()),
					std::end(bCalculated.r()), [](const scalar value) 
					{
						return value < 0.;
					}, 0.0);

			diffFieldsOld.b.r() = bCalculated();
			diffFieldsNew.b.r() = bCalculated();

			if (linearRec)
			{
				grada = grad(diffFieldsOld.a, bncCalc);
				diva = divergence(diffFieldsOld.a, bncCalc);
			}
			else
			{
				grada = grad(star.a);
				diva = divergence(star.a);
			}

			//auto aDiffSurface = surfGrad(diffFieldsOld.a, bncCalc);

			//for (auto & a_i : aDiffSurface.r())
			//	a_i.transpose();

			//aDiffSurface.r() = astProduct(aDiffSurface,
			//		effectiveCoeffs.rhoDa)();

			//divDeva1 = divergence(aDiffSurface);
		}
			break;
		case turbulenceModel::BHRSource:
		case turbulenceModel::BHRKLSource:
		{
			if (linearRec)
			{
				grada = grad(diffFieldsOld.a, bncCalc);
				diva = divergence(diffFieldsOld.a, bncCalc);
				gradb = grad(diffFieldsOld.b, bncCalc);
			}
			else
			{
				grada = grad(star.a);
				diva = divergence(star.a);
				gradb = grad(star.b);
			}

			//auto aDiffSurface = surfGrad(diffFieldsOld.a, bncCalc);

			//for (auto & a_i : aDiffSurface.r())
			//	a_i.transpose();

			//aDiffSurface.r() = astProduct(aDiffSurface,
			//		effectiveCoeffs.rhoDa)();

			//divDeva1 = divergence(aDiffSurface);
		}
			break;
		case turbulenceModel::arithmeticA1Source:
		{
			volumeField<scalar> sonicSpeed2 { mesh_, 0 };
			sonicSpeed2.r() = gasPhase.phaseThermodynamics->sqSonicSpeed(
					gasPhase.concentration.p, gasPhase.density[0](),
					gasPhase.internalEnergy(), gasPhase.pressure());

			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsOld.a.r()[i] = -gasPhase.tNu()[i]
						/ gasPhase.turbulenceSources->turbPar->Ca1()
						* (gradRho()[i] / gasPhase.density[0]()[i]
								- gradP()[i]
										/ (gasPhase.density[0]()[i]
												* sonicSpeed2()[i]));

			grada = grad(diffFieldsOld.a, bncCalc);
			diva = divergence(diffFieldsOld.a, bncCalc);
		}
			break;
		case turbulenceModel::arithmeticA2Source:
		{
			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsOld.a.r()[i] = -gasPhase.tNu()[i]
						/ gasPhase.turbulenceSources->turbPar->Ca1()
						* (gradRho()[i] / gasPhase.density[0]()[i]
								- gradP()[i] / gasPhase.pressure()[i]);

			grada = grad(diffFieldsOld.a, bncCalc);
			diva = divergence(diffFieldsOld.a, bncCalc);
		}
			break;
		case turbulenceModel::arithmeticA3Source:
		{
			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsOld.a.r()[i] = -gasPhase.tNu()[i]
						/ gasPhase.turbulenceSources->turbPar->Ca1()
						* gradRho()[i] / gasPhase.density[0]()[i];

			grada = grad(diffFieldsOld.a, bncCalc);
			diva = divergence(diffFieldsOld.a, bncCalc);
		}
			break;
		default:
			break;
		}

		devPhysViscSurf = devTotVisc1;
	}

	/*Transposition of velocity gradient.*/
	for (auto & devTotVisc1_i : devTotVisc1.r())
		devTotVisc1_i.transpose();

	if (gasPhase.turbulenceSources->turbulence)
		devPhysViscSurf.r() += devTotVisc1();

	/*Explicit deviatoric part of molecular viscosity tensor for energy computation.*/
	auto devPhysVisc = gradV;

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		devPhysVisc.r()[i].transpose();

	devPhysVisc.r() += gradV();

	/*Add buoyant viscosity component to diagonal of all three tensors.*/
	for (std::size_t i = 0; i < mesh_.surfacesSize(); ++i)
	{
		const scalar lambdadivV_s_i { twothirds * divV_s()[i] };

		devTotVisc1.r()[i].r()[0] -= lambdadivV_s_i;
		devTotVisc1.r()[i].r()[4] -= lambdadivV_s_i;
		devTotVisc1.r()[i].r()[8] -= lambdadivV_s_i;

		if (gasPhase.turbulenceSources->turbulence)
		{
			devPhysViscSurf.r()[i].r()[0] -= lambdadivV_s_i;
			devPhysViscSurf.r()[i].r()[4] -= lambdadivV_s_i;
			devPhysViscSurf.r()[i].r()[8] -= lambdadivV_s_i;
		}
	}

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar lambdadivV_i { twothirds * divV()[i] };

		devPhysVisc.r()[i].r()[0] -= lambdadivV_i;
		devPhysVisc.r()[i].r()[4] -= lambdadivV_i;
		devPhysVisc.r()[i].r()[8] -= lambdadivV_i;
	}

	volumeField<tensor> spherTurbR { mesh_, tensor(0) };
	volumeField<tensor> devTurbR { mesh_, tensor(0) };
	volumeField<tensor> philtTurbVisc { mesh_, tensor(0) };

	if (gasPhase.turbulenceSources->turbulence)
	{
		devTurbR = devPhysVisc;
		philtTurbVisc = devPhysVisc;
	}

	/*Multiplying on physical, effective or turbulent viscosity correspondingly.*/
	devPhysVisc.r() = astProduct(devPhysVisc, gasPhase.physMu)();

	if (gasPhase.turbulenceSources->turbulence)
	{
		auto turbViscCoeff1 = gasPhase.tNu;
		turbViscCoeff1.r() = (1. - thetaS_R()) * turbViscCoeff1();

		philtTurbVisc = astProduct(
				astProduct(turbViscCoeff1, gasPhase.density[0]), philtTurbVisc);

		devTurbR = astProduct(astProduct(devTurbR, gasPhase.density[0]),
				astProduct(gasPhase.tNu, thetaS_R));

		for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		{
			/*Calculate turbulent pressure part*/
			const scalar turbPressure { -twothirds * gasPhase.rhokTurb()[i] };
			spherTurbR.r()[i] = tensor(turbPressure, 0, 0, 0, turbPressure, 0,
					0, 0, turbPressure);
		}
	}

	astProductSelf(devTotVisc1, effectiveCoeffs.mu);

	astProductSelf(devPhysViscSurf, effectiveCoeffs.physMu);

	const auto divDevPhysVisc = divergence(devPhysViscSurf);

	/*Accounting of explicit part of effective viscosity tensor*/
	const auto divDevTotVisc1 = divergence(devTotVisc1);
	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		for (std::size_t j = 0; j < vector::vsize; ++j)
			velocityMatrix.SLE[j].freeTerm[i] += divDevTotVisc1()[i]()[j];

	const auto divNonIdealityCorrectionLaplacian = divergence(
			nonIdealityCorrectionLaplacian);

	temperatureMatrix.SLE[0].freeTerm -= divNonIdealityCorrectionLaplacian();

	//temperatureMatrix.distributeSourceTerm(
	//		astProduct(divNonIdealityCorrectionLaplacian, -1.0),
	//		diffFieldsOld.temperature, timestep);

	if (enthalpySolverFlag == enthalpyFlow::explicitSolve)
	{
		astProductSelf(entExplDiffFlow, surfaceTemperature);

		const auto divEnthFlowLaplacian = divergence(entExplDiffFlow);

		temperatureMatrix.SLE[0].freeTerm -= divEnthFlowLaplacian();

		//temperatureMatrix.distributeSourceTerm(
		//		astProduct(divEnthFlowLaplacian, -1), diffFieldsOld.temperature,
		//		timestep);
	}

	volumeField<vector> gradMav_Mav { mesh_, vector { 0 } };

	if (gasPhase.turbulenceSources->turbulence)
	{
		if (gasPhase.turbulenceSources->model != turbulenceModel::zeroSource)
			temperatureMatrix.SLE[0].freeTerm +=
					gasPhase.turbulenceSources->turbPar->rhoepsilon(gasPhase);

		if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::kEpsASource))
			gradMav_Mav.r() = division(
					grad(division(surfaceRho, surfaceConcentration.v[0])),
					avMolMass)();

		switch (gasPhase.turbulenceSources->model)
		{
		case turbulenceModel::kEpsASource:
		case turbulenceModel::BHRSource:
		case turbulenceModel::BHRKLSource:
		case turbulenceModel::arithmeticA1Source:
		case turbulenceModel::arithmeticA2Source:
		case turbulenceModel::arithmeticA3Source:
		{
			volumeField<scalar> pDiva { mesh_ };
			volumeField<vector> totalPdivergence { mesh_ };

			totalPdivergence.r() = gradP() - divDevPhysVisc();

			pDiva.r() = ampProduct(diffFieldsOld.a, totalPdivergence)();
			astProductSelf(pDiva, (-1.0));

			temperatureMatrix.SLE[0].freeTerm += pDiva();

			//temperatureMatrix.distributeSourceTerm(pDiva,
			//		diffFieldsOld.temperature, timestep);
		}
			break;
		default:
			break;
		}

		//FIXME Do not work correctly with MPI.
		//if ((source.model == BHRSource) || (source.model == BHRKLSource)
		//		|| (source.model == kEpsASource))
		//	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		//		for (std::size_t j = 0; j < vector::vsize; ++j)
		//			aMatrix.SLE[j].freeTerm[i] += divDeva1()[i].v()[j];
	}

	temperatureMatrix.SLE[0].freeTerm += dampProduct(devPhysVisc, gradV)();

	if (gasPhase.turbulenceSources->turbulence)
		temperatureMatrix.SLE[0].freeTerm +=
				dampProduct(philtTurbVisc, gradV)();

	/*Calculating turbulent sources.*/
	volumeField<scalar> sigmaSourcek { mesh_, scalar { 0 } };
	volumeField<scalar> sigmaSourceeps { mesh_, 0 };
	volumeField<vector> sigmaSourcea { mesh_, vector(0) };
	volumeField<scalar> sigmaSourceb { mesh_, 0 };

	if (gasPhase.turbulenceSources->turbulence)
	{
		auto & nonConstMesh = const_cast<mesh&>(mesh_);

		const auto [SourcekSuSp, SourceepsSuSp, SourceaSuSp, SourcebSuSp] =
				gasPhase.turbulenceSources->calculate(
						nonConstMesh.timestepSourceRef(), timestepCoeffs.first,
						gasPhase, diffFieldsOld, gradV, divDevPhysVisc, gradP,
						gradRho, grada, diva, gradb, spherTurbR, devTurbR,
						gradMav_Mav, *(gasPhase.phaseThermodynamics),
						gasPhase.tNu);

		if (msolver.solverType == matrixSolver::explicitSolver)
		{
			sigmaSourcek.r() = SourcekSuSp.first()
					+ SourcekSuSp.second() * diffFieldsOld.k();
			sigmaSourceeps.r() = SourceepsSuSp.first()
					+ SourceepsSuSp.second() * diffFieldsOld.eps();

			volumeField<vector> SpA(mesh_, vector(0));
			for (std::size_t i = 0; i < SpA.size(); ++i)
			{
				const vector vec = vector(
						SourceaSuSp.second()[i]()[0]
								* diffFieldsOld.a()[i]()[0],
						SourceaSuSp.second()[i]()[1]
								* diffFieldsOld.a()[i]()[1],
						SourceaSuSp.second()[i]()[2]
								* diffFieldsOld.a()[i]()[2]);

				sigmaSourcea.r()[i] = SourceaSuSp.first()[i] + vec;
			}
		}
		else
		{
			kMatrix.freeSourceTerm(SourcekSuSp.first);
			kMatrix.diagonaleSourceTerm(SourcekSuSp.second);

			epsMatrix.freeSourceTerm(SourceepsSuSp.first);
			epsMatrix.diagonaleSourceTerm(SourceepsSuSp.second);

			aMatrix.freeSourceTerm(SourceaSuSp.first);
			aMatrix.diagonaleSourceTerm(SourceaSuSp.second);

			bMatrix.freeSourceTerm(SourcebSuSp.first);
			bMatrix.diagonaleSourceTerm(SourcebSuSp.second);
		}
	}

	/*Calculation of new velocity by components of vector*/
	if (msolver.solverType != matrixSolver::explicitSolver)
		diffFieldsNew.velocity.r() = msolver.solve(diffFieldsOld.velocity(),
				velocityMatrix);
	else
		for (std::size_t j = 0; j < vector::vsize; ++j)
			for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
				diffFieldsNew.velocity.r()[i].r()[j] =
						(velocityMatrix.SLE[j].explOldTime[i]
								+ velocityMatrix.SLE[j].freeTerm[i] * timestep)
								/ velocityMatrix.SLE[j].centralDiagonale[i];

	/*Calculation of new temperature*/
	if (enthalpySolverFlag == enthalpyFlow::implicitSolve)
	{
		if (msolver.solverType != matrixSolver::explicitSolver)
			intermediateTemperature.r() = msolver.solve(
					diffFieldsOld.temperature(), temperatureMatrix);
		else
			intermediateTemperature.r() = (temperatureMatrix.SLE[0].explOldTime
					+ temperatureMatrix.SLE[0].freeTerm * timestep)
					/ temperatureMatrix.SLE[0].centralDiagonale;

		parallelism.correctBoundaryConditions(intermediateTemperature);

		nonIdealCorrectionOld.r() = gasPhase.phaseThermodynamics->nonIdeality(
				gasPhase.concentration.p, intermediateTemperature());

		oldEnergyField.r() = CCVOld() * intermediateTemperature()
				+ nonIdealCorrectionOld() - nonIdealCorrectionNew();

		temperatureEnthMatrix.generateDTimeNabla(intermediateTemperature,
				oldEnergyField, CCVNew, entExplDiffFlow, timestep, bncCalc);

		diffFieldsNew.temperature.r() = msolverForEnthalpy.solve(
				intermediateTemperature(), temperatureEnthMatrix);
	}
	else
	{
		if (msolver.solverType != matrixSolver::explicitSolver)
			diffFieldsNew.temperature.r() = msolver.solve(
					diffFieldsOld.temperature(), temperatureMatrix);
		else
			diffFieldsNew.temperature.r() =
					(temperatureMatrix.SLE[0].explOldTime
							+ temperatureMatrix.SLE[0].freeTerm * timestep)
							/ temperatureMatrix.SLE[0].centralDiagonale;
	}

	if (gasPhase.turbulenceSources->turbulence)
	{
		/*Explicit source integration*/
		/*Calculation of new k*/
		if (msolver.solverType != matrixSolver::explicitSolver)
			diffFieldsNew.k.r() = msolver.solve(diffFieldsOld.k(), kMatrix);
		else
			diffFieldsNew.k.r() = (kMatrix.SLE[0].explOldTime
					+ kMatrix.SLE[0].freeTerm * timestep)
					/ kMatrix.SLE[0].centralDiagonale
					+ sigmaSourcek() / gasPhase.density[0]() * timestep;

		/*Calculation of new epsilon*/
		if (msolver.solverType != matrixSolver::explicitSolver)
			diffFieldsNew.eps.r() = msolver.solve(diffFieldsOld.eps(),
					epsMatrix);
		else
			diffFieldsNew.eps.r() = (epsMatrix.SLE[0].explOldTime
					+ epsMatrix.SLE[0].freeTerm * timestep)
					/ epsMatrix.SLE[0].centralDiagonale
					+ sigmaSourceeps() / gasPhase.density[0]() * timestep;

		if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::kEpsASource))
		{
			/*Calculation of new a by components of vector*/
			if (msolver.solverType != matrixSolver::explicitSolver)
				diffFieldsNew.a.r() = msolver.solve(diffFieldsOld.a(), aMatrix);
			else
			{
				for (std::size_t j = 0; j < vector::vsize; ++j)
					for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
						diffFieldsNew.a.r()[i].r()[j] =
								(aMatrix.SLE[j].explOldTime[i]
										+ aMatrix.SLE[j].freeTerm[i] * timestep)
										/ aMatrix.SLE[j].centralDiagonale[i];

				diffFieldsNew.a.r() +=
						astProduct(division(sigmaSourcea, gasPhase.density[0]),
								timestep)();
			}

			if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::BHRKLSource))
			/*Calculation of new b*/
			{
				if (msolver.solverType != matrixSolver::explicitSolver)
					diffFieldsNew.b.r() = msolver.solve(diffFieldsOld.b(),
							bMatrix);
				else
					diffFieldsNew.b.r() = (bMatrix.SLE[0].explOldTime
							+ bMatrix.SLE[0].freeTerm * timestep)
							/ bMatrix.SLE[0].centralDiagonale
							+ sigmaSourceb() / gasPhase.density[0]() * timestep;
			}
		}

		/*Check for non-negativity turbulent quantities.*/
		std::replace_if(std::begin(diffFieldsNew.k.r()),
				std::end(diffFieldsNew.k.r()), [&gasPhase](const scalar value) 
				{
					return value < gasPhase.turbulenceSources->turbPar->mink();
				}, gasPhase.turbulenceSources->turbPar->mink());

		std::replace_if(std::begin(diffFieldsNew.eps.r()),
				std::end(diffFieldsNew.eps.r()),
				[&gasPhase](
						const scalar value) 
						{
							return value < gasPhase.turbulenceSources->turbPar->mineps();
						}, gasPhase.turbulenceSources->turbPar->mineps());

		//XXX Non-physical bounding of vector a.
		//if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
		//		|| (gasPhase.turbulenceSources->model
		//				== turbulenceModel::BHRKLSource)
		//		|| (gasPhase.turbulenceSources->model
		//				== turbulenceModel::kEpsASource))
		//	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
		//		diffFieldsNew.a.r()[i] =
		//				[&diffFieldsNew](const std::size_t i_l)
		//				{
		//					const auto module_a = diffFieldsNew.a()[i_l].mag();
		//					const auto module_turb_V = std::sqrt(
		//							2 * diffFieldsNew.k()[i_l]);
		//					if (module_a < module_turb_V)
		//						return diffFieldsNew.a()[i_l];
		//					else
		//						return diffFieldsNew.a()[i_l] / module_a
		//								* module_turb_V;
		//				}(i);

		if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::BHRKLSource))
			std::replace_if(std::begin(diffFieldsNew.b.r()),
					std::end(diffFieldsNew.b.r()),
					[&gasPhase](
							const scalar value) 
							{
								return value < gasPhase.turbulenceSources->turbPar->minb_value;
							}, gasPhase.turbulenceSources->turbPar->minb_value);
	}

	/*Recalculation of values in cells*/
	gasPhase.concentration.v[0].r() = intermediateConcentration.v[0]();
	for (std::size_t k = 1; k < gasPhase.density.size(); ++k)
	{
		gasPhase.concentration.v[k].r() = intermediateConcentration.v[k]();

		gasPhase.density[k].r() = gasPhase.concentration.v[k]()
				* gasPhase.phaseThermodynamics->Mv()[k - 1];
	}

	gasPhase.velocity.r() = diffFieldsNew.velocity();

	gasPhase.momentum.r() =
			astProduct(gasPhase.velocity, gasPhase.density[0])();

	gasPhase.temperature.r() = diffFieldsNew.temperature();

	gasPhase.internalEnergy.r() = gasPhase.phaseThermodynamics->UvcFromT(
			gasPhase.concentration.p, gasPhase.temperature())
			* gasPhase.concentration.v[0]();

	gasPhase.pressure.r() = gasPhase.phaseThermodynamics->pFromUv(
			gasPhase.concentration.p, gasPhase.internalEnergy());

	if (gasPhase.turbulenceSources->turbulence)
	{
		gasPhase.kTurb.r() = diffFieldsNew.k();
		gasPhase.epsTurb.r() = diffFieldsNew.eps();
		if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::kEpsASource))
		{
			gasPhase.aTurb.r() = diffFieldsNew.a();
			gasPhase.bTurb.r() = diffFieldsNew.b();
		}

		gasPhase.rhokTurb.r() =
				astProduct(gasPhase.kTurb, gasPhase.density[0])();
		gasPhase.rhoepsTurb.r() = gasPhase.epsTurb() * gasPhase.density[0]();
		if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::BHRKLSource)
				|| (gasPhase.turbulenceSources->model
						== turbulenceModel::kEpsASource))
		{
			gasPhase.rhoaTurb.r() = astProduct(gasPhase.aTurb,
					gasPhase.density[0])();

			if ((gasPhase.turbulenceSources->model == turbulenceModel::BHRSource)
					|| (gasPhase.turbulenceSources->model
							== turbulenceModel::BHRKLSource))
				gasPhase.rhobTurb.r() = gasPhase.bTurb()
						* gasPhase.density[0]();
		}
	}

	{
		const auto v2 = ampProduct(gasPhase.velocity, gasPhase.velocity);

		gasPhase.totalEnergy.r() = gasPhase.internalEnergy()
				+ gasPhase.density[0]() * v2() * 0.5 + gasPhase.rhokTurb();
	}

	gasPhase.HelmholtzEnergy.r() = gasPhase.phaseThermodynamics->Fv(
			gasPhase.concentration.p, gasPhase.temperature.r());

	gasPhase.entropy.r() = gasPhase.phaseThermodynamics->Sv(
			gasPhase.concentration.p, gasPhase.temperature());

	/*Diffusion time-step.*/
	if (sourceTimeFlag == timestep::CourantAndSourceAndDiffusionTimeStep)
	{
		const auto maxVal = effectiveCoeffs.maxValue(
				gasPhase.turbulenceSources->turbulence,
				gasPhase.turbulenceSources->model);

		const scalar diffuse_dt = 0.5 * timestepCoeffs.second
				/ (linearInterpolate(maxVal)() * minimalLengthScale()
						* minimalLengthScale()).min();

		auto & nonConstMesh = const_cast<mesh&>(mesh_);

		nonConstMesh.timestepSourceRef() = std::min(
				nonConstMesh.timestepSource(), diffuse_dt);
	}

	if (gasPhase.temperature().min() < 0.)
		throw exception("Negative temperature after diffusion.",
				errors::negativeTemperatureError);

	const auto DiffusionEndTime { std::chrono::high_resolution_clock::now() };
	timeForDiffusion += std::chrono::duration_cast<std::chrono::milliseconds>(
			DiffusionEndTime - DiffusionStartTime).count();
}
