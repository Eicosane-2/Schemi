/*
 * BHRKLGen.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "BHRKLGen.hpp"

#include <algorithm>

#include "turbulentParametersKL.hpp"
#include "doubleDotProduct.hpp"

schemi::BHRKLGen::BHRKLGen(const bool turb_in,
		const turbulenceModelEnum tm_in) noexcept :
		abstractTurbulenceGen(turb_in, tm_in)
{
	turbPar = std::make_unique<turbulentParametersKL>();
}

std::tuple<schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::vector>, schemi::volumeField<schemi::scalar>> schemi::BHRKLGen::calculate(
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

	std::valarray<scalar> modeps(diffFieldsOld.eps.ref());
	const scalar maxeps { modeps.max() };
	std::replace_if(std::begin(modeps), std::end(modeps),
			[maxeps](const scalar value) 
			{
				return value < 1E-3 * maxeps;
			}, veryBig);

	volumeField<vector> divaa(mesh_, vector(0));
	divaa.ref_r() = astProduct(diffFieldsOld.a, diva).ref()
			+ ampProduct(diffFieldsOld.a, grada).ref();

	const auto a_s2 = mixture.sqSonicSpeed(cellFields.concentration.p,
			cellFields.density[0].ref(), cellFields.internalEnergy.ref(),
			cellFields.pressure.ref());

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps.ref()[i] / diffFieldsOld.k.ref()[i] };
		const scalar sqrke { std::sqrt(diffFieldsOld.k.ref()[i])
				/ diffFieldsOld.eps.ref()[i] };

		const auto thetaS_i = turbPar->thetaS_R(gradV.ref()[i].trace(),
				diffFieldsOld.k.ref()[i], diffFieldsOld.eps.ref()[i]);

		const auto thetaB_i = turbPar->thetaB(diffFieldsOld.a.ref()[i],
				diffFieldsOld.k.ref()[i], diffFieldsOld.eps.ref()[i],
				gradMav_n.ref()[i], a_s2[i],
				std::pair<scalar, vector>(cellFields.density[0].ref()[i],
						gradRho.ref()[i]),
				std::pair<scalar, vector>(cellFields.pressure.ref()[i],
						gradP.ref()[i]), nu_t.ref()[i]);

		const scalar rhoSpherRGen = spherR.ref()[i] && gradV.ref()[i];

		const scalar rhoDevRGen = thetaS_i * devR.ref()[i] && gradV.ref()[i];

		const scalar gravGen { diffFieldsOld.a.ref()[i]
				& (gradP.ref()[i] - divDevPhysVisc.ref()[i]) };

		const scalar dissip(
				-cellFields.density[0].ref()[i] * diffFieldsOld.k.ref()[i]
						* sqrke);

		sigmaSourcek.ref_r()[i] = rhoSpherRGen + rhoDevRGen + gravGen + dissip;

		sigmaSourceeps.ref_r()[i] = turbPar->C1() * ek * rhoDevRGen
				+ turbPar->C3() * ek * rhoSpherRGen
				+ turbPar->C0() / thetaB_i * ek * gravGen
				+ turbPar->C2() * ek * dissip;

		/*Time-step calculation*/
		modeps[i] =
				std::abs(
						sourceTimestepCoeff * modeps[i]
								/ (sigmaSourceeps.ref()[i]
										/ cellFields.density[0].ref()[i]
										+ stabilizator));

		const vector bGradP(
				(gradP.ref()[i] - divDevPhysVisc.ref()[i])
						* diffFieldsOld.b.ref()[i]);

		const vector tauGradRho(
				(devR.ref()[i] * thetaS_i + spherR.ref()[i])
						/ cellFields.density[0].ref()[i] & gradRho.ref()[i]);

		const vector rhoAgradV(
				cellFields.rhoaTurb.ref()[i]
						& (grada.ref()[i] - gradV.ref()[i]));

		//const vector redistribution_a = cellFields.density[0].ref()[i]
		//		* divaa.ref()[i];

		sigmaSourcea.ref_r()[i] = bGradP + tauGradRho + rhoAgradV
		//+ redistribution_a
				- cellFields.rhoaTurb.ref()[i] * sqrke * turbPar->Ca1();

		const scalar bagradRho { -2. * (diffFieldsOld.b.ref()[i] + 1.)
				* (diffFieldsOld.a.ref()[i] & gradRho.ref()[i]) };

		//const scalar redistribution_b = 2. * cellFields.rhoaTurb.ref()[i]
		//		& gradb.ref()[i];

		sigmaSourceb.ref_r()[i] = bagradRho		//+ redistribution_b
		- cellFields.rhobTurb.ref()[i] * sqrke * turbPar->Cb1();
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(sigmaSourcek, sigmaSourceeps, sigmaSourcea,
			sigmaSourceb);
}
