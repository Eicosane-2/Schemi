/*
 * BHR2Gen.cpp
 *
 *  Created on: 2024/05/01
 *      Author: Maxim Boldyrev
 */

#include "BHR2Gen.hpp"

#include <algorithm>

#include "doubleDotProduct.hpp"
#include "turbulentParametersKEPS.hpp"

schemi::BHR2Gen::BHR2Gen(const mesh & meshIn, const bool turb_in,
		const turbulenceModel tm_in) noexcept :
		abstractTurbulenceGen(turb_in, tm_in)
{
	turbPar = std::make_unique<turbulentParametersKEPS>(meshIn);
}

std::tuple<
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		std::pair<schemi::volumeField<schemi::vector>,
				schemi::volumeField<schemi::vector>>,
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		schemi::volumeField<schemi::scalar>> schemi::BHR2Gen::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld,
		const volumeField<tensor> & gradV,
		const volumeField<vector> & divDevPhysVisc,
		const volumeField<vector> & gradP, const volumeField<vector> & gradRho,
		const volumeField<tensor> & grada, const volumeField<scalar>&,
		const volumeField<vector> & gradb, const volumeField<tensor> & spherR,
		const volumeField<tensor> & devR, const volumeField<vector> & gradMav_n,
		const abstractMixtureThermodynamics & mixture,
		const volumeField<scalar> & nu_t) const noexcept
{
	auto & mesh_ { cellFields.pressure.meshRef() };

	std::pair<volumeField<scalar>, volumeField<scalar>> Sourcek { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	std::pair<volumeField<scalar>, volumeField<scalar>> Sourceeps { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	std::pair<volumeField<vector>, volumeField<vector>> Sourcea { volumeField<
			vector>(mesh_, vector(0)), volumeField<vector>(mesh_, vector(0)) };
	std::pair<volumeField<scalar>, volumeField<scalar>> Sourceb { volumeField<
			scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	volumeField<scalar> gravGenField(mesh_, 0);

	std::valarray<scalar> modeps(diffFieldsOld.eps());
	const scalar maxeps { modeps.max() };
	std::replace_if(std::begin(modeps), std::end(modeps),
			[maxeps](const scalar value) 
			{
				return value < 1E-3 * maxeps;
			}, veryBig);

	const auto a_s2 = mixture.sqSonicSpeed(cellFields.concentration.p,
			cellFields.density[0](), cellFields.internalEnergy(),
			cellFields.pressure());

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps()[i] / diffFieldsOld.k()[i] };

		const auto thetaS_i = turbPar->thetaS_R(gradV()[i].trace(),
				diffFieldsOld.k()[i], diffFieldsOld.eps()[i]);

		const auto thetaB_i = turbPar->thetaB(diffFieldsOld.a()[i],
				diffFieldsOld.k()[i], diffFieldsOld.eps()[i], gradMav_n()[i],
				a_s2[i],
				std::pair<scalar, vector>(cellFields.density[0]()[i],
						gradRho()[i]),
				std::pair<scalar, vector>(cellFields.pressure()[i], gradP()[i]),
				nu_t()[i]);

		const auto thetaA_i = turbPar->thetaA(diffFieldsOld.a()[i],
				diffFieldsOld.k()[i], diffFieldsOld.b()[i]);

		const scalar rhoSpherRGen = spherR()[i] && gradV()[i];

		const scalar rhoDevRGen = thetaS_i * devR()[i] && gradV()[i];

		const scalar gravGen { (diffFieldsOld.a()[i]
				& (gradP()[i] - divDevPhysVisc()[i])) };

		const scalar dissip(
				-cellFields.rhoepsTurb()[i]
						* (1 + 2 * diffFieldsOld.k()[i] / a_s2[i]));

		Sourcek.first.r()[i] = rhoSpherRGen + rhoDevRGen + gravGen;
		Sourcek.second.r()[i] = dissip / cellFields.kTurb()[i];

		Sourceeps.first.r()[i] = turbPar->C1() * ek * rhoDevRGen
				+ turbPar->C3() * ek * rhoSpherRGen
				+ std::min(turbPar->C0() / (thetaB_i + stabilizator),
						turbPar->C0max()) * ek * gravGen;
		Sourceeps.second.r()[i] = turbPar->C2() * dissip
				/ cellFields.kTurb()[i];

		/*Time-step calculation*/
		modeps[i] = std::abs(
				sourceTimestepCoeff * modeps[i]
						/ ((Sourceeps.first()[i]
								+ Sourceeps.second()[i]
										* cellFields.epsTurb()[i])
								/ cellFields.density[0]()[i] + stabilizator));

		const vector bGradP(
				(gradP()[i] - divDevPhysVisc()[i]) * diffFieldsOld.b()[i]);

		const vector tauGradRho(
				(devR()[i] * thetaS_i + spherR()[i])
						/ cellFields.density[0]()[i] & gradRho()[i]);

		const vector rhoAgradV(
				cellFields.rhoaTurb()[i] & (grada()[i] - gradV()[i]));

		Sourcea.first.r()[i] = bGradP + tauGradRho + rhoAgradV;
		Sourcea.second.r()[i] = -cellFields.density[0]()[i] * ek
				* turbPar->Ca1() * thetaA_i;

		const scalar bagradRho { -2 * (diffFieldsOld.b()[i] + 1.)
				* (diffFieldsOld.a()[i] & gradRho()[i]) };

		const scalar redistribution_b = 2
				* (cellFields.rhoaTurb()[i] & gradb()[i]);

		Sourceb.first.r()[i] = bagradRho + redistribution_b;
		Sourceb.second.r()[i] = -cellFields.density[0]()[i] * ek
				* turbPar->Cb1() * thetaA_i;

		gravGenField.r()[i] = gravGen;
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}
