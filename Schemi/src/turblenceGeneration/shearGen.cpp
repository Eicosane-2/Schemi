/*
 * shearGen.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "shearGen.hpp"

#include <algorithm>

#include "turbulentParametersKEPS.hpp"
#include "doubleDotProduct.hpp"

schemi::shearGen::shearGen(const bool turb_in,
		const turbulenceModelEnum tm_in) noexcept :
		abstractTurbulenceGen(turb_in, tm_in)
{
	turbPar = std::make_unique<turbulentParametersKEPS>();
}

std::tuple<schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::vector>, schemi::volumeField<schemi::scalar>> schemi::shearGen::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld,
		const volumeField<tensor> & gradV, const volumeField<vector>&,
		const volumeField<vector>&, const volumeField<vector>&,
		const volumeField<tensor>&, const volumeField<scalar>&,
		const volumeField<vector>&, const volumeField<tensor> & spherR,
		const volumeField<tensor> & devR, const volumeField<vector>&,
		const abstractMixtureThermodynamics&,
		const volumeField<scalar>&) const noexcept
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

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps.ref()[i] / diffFieldsOld.k.ref()[i] };

		const scalar rhoSpherRGen = spherR.ref()[i] && gradV.ref()[i];

		const scalar rhoDevRGen = devR.ref()[i] && gradV.ref()[i];

		const scalar dissip(-cellFields.rhoepsTurb.ref()[i]);

		sigmaSourcek.ref_r()[i] = rhoSpherRGen + rhoDevRGen + dissip;

		sigmaSourceeps.ref_r()[i] = turbPar->C1() * ek * rhoDevRGen
				+ turbPar->C3() * ek * rhoSpherRGen
				+ turbPar->C2() * ek * dissip;

		/*Time-step calculation*/
		modeps[i] =
				std::abs(
						sourceTimestepCoeff * modeps[i]
								/ (sigmaSourceeps.ref()[i]
										/ cellFields.density[0].ref()[i]
										+ stabilizator));
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(sigmaSourcek, sigmaSourceeps, sigmaSourcea,
			sigmaSourceb);
}
