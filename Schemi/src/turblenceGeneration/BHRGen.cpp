/*
 * BHRGen.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "BHRGen.hpp"

#include <algorithm>

#include "turbulentParametersKEPS.hpp"
#include "doubleDotProduct.hpp"

schemi::BHRGen::BHRGen(const mesh & meshIn, const bool turb_in,
		const turbulenceModel tm_in) noexcept :
		abstractTurbulenceGen(turb_in, tm_in)
{
	turbPar = std::make_unique<turbulentParametersKEPS>(meshIn);
}

std::tuple<schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::vector>, schemi::volumeField<schemi::scalar>> schemi::BHRGen::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld,
		const volumeField<tensor> & gradV,
		const volumeField<vector> & divDevPhysVisc,
		const volumeField<vector> & gradP, const volumeField<vector> & gradRho,
		const volumeField<tensor> & grada, const volumeField<scalar> & diva,
		const volumeField<vector>&, const volumeField<tensor> & spherR,
		const volumeField<tensor> & devR, const volumeField<vector> & gradMav_n,
		const abstractMixtureThermodynamics & mixture,
		const volumeField<scalar> & nu_t) const noexcept
{
	auto & mesh_ { cellFields.pressure.meshRef() };

	volumeField<scalar> sigmaSourcek(mesh_, 0);
	volumeField<scalar> sigmaSourceeps(mesh_, 0);
	volumeField<vector> sigmaSourcea(mesh_, vector(0));
	volumeField<scalar> sigmaSourceb(mesh_, 0);

	std::valarray<scalar> modeps(diffFieldsOld.eps());
	const scalar maxeps { modeps.max() };
	std::replace_if(std::begin(modeps), std::end(modeps),
			[maxeps](const scalar value) 
			{
				return value < 1E-3 * maxeps;
			}, veryBig);

	volumeField<vector> divaa(mesh_, vector(0));
	divaa.r() = astProduct(diffFieldsOld.a, diva)()
			+ ampProduct(diffFieldsOld.a, grada)();

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

		const scalar rhoSpherRGen = spherR()[i] && gradV()[i];

		const scalar rhoDevRGen = thetaS_i * devR()[i] && gradV()[i];

		const scalar gravGen { diffFieldsOld.a()[i]
				& (gradP()[i] - divDevPhysVisc()[i]) };

		const scalar dissip(-cellFields.rhoepsTurb()[i]);

		sigmaSourcek.r()[i] = rhoSpherRGen + rhoDevRGen + gravGen + dissip;

		sigmaSourceeps.r()[i] = turbPar->C1() * ek * rhoDevRGen
				+ turbPar->C3() * ek * rhoSpherRGen
				+ turbPar->C0() / thetaB_i * ek * gravGen
				+ turbPar->C2() * ek * dissip;

		/*Time-step calculation*/
		modeps[i] = std::abs(
				sourceTimestepCoeff * modeps[i]
						/ (sigmaSourceeps()[i] / cellFields.density[0]()[i]
								+ stabilizator));

		const vector bGradP(
				(gradP()[i] - divDevPhysVisc()[i]) * diffFieldsOld.b()[i]);

		const vector tauGradRho(
				(devR()[i] * thetaS_i + spherR()[i])
						/ cellFields.density[0]()[i] & gradRho()[i]);

		const vector rhoAgradV(
				cellFields.rhoaTurb()[i] & (grada()[i] - gradV()[i]));

		//const vector redistribution_a = cellFields.density[0]()[i]
		//		* divaa()[i];

		sigmaSourcea.r()[i] = bGradP + tauGradRho + rhoAgradV
		//+ redistribution_a
				- cellFields.rhoaTurb()[i] * ek * turbPar->Ca1();

		const scalar bagradRho { -2. * (diffFieldsOld.b()[i] + 1.)
				* (diffFieldsOld.a()[i] & gradRho()[i]) };

		//const scalar redistribution_b = 2. * cellFields.rhoaTurb()[i]
		//		& gradb()[i];

		sigmaSourceb.r()[i] = bagradRho		//+ redistribution_b
		- cellFields.rhobTurb()[i] * ek * turbPar->Cb1();
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(sigmaSourcek, sigmaSourceeps, sigmaSourcea,
			sigmaSourceb);
}
