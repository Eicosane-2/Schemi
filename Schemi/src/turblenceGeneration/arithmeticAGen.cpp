/*
 * arithmeticAGen.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "arithmeticAGen.hpp"

#include <algorithm>

#include "turbulentParametersKEPS.hpp"
#include "vectorVectorDotProduct.hpp"
#include "doubleDotProduct.hpp"
#include "fieldProducts.hpp"

schemi::arithmeticAGen::arithmeticAGen(const mesh & meshIn, const bool turb_in,
		const turbulenceModel tm_in) noexcept :
		abstractTurbulenceGen(turb_in, tm_in)
{
	turbPar = std::make_unique<turbulentParametersKEPS>(meshIn);
}

std::tuple<schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::vector>, schemi::volumeField<schemi::scalar>> schemi::arithmeticAGen::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld,
		const volumeField<tensor> & gradV,
		const volumeField<vector> & divDevPhysVisc,
		const volumeField<vector> & gradP, const volumeField<vector>&,
		const volumeField<tensor> & grada, const volumeField<scalar> & diva,
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

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps()[i] / diffFieldsOld.k()[i] };

		const scalar rhoSpherRGen = spherR()[i] && gradV()[i];

		const scalar rhoDevRGen = devR()[i] && gradV()[i];

		const scalar gravGen { diffFieldsOld.a()[i]
				& (gradP()[i] - divDevPhysVisc()[i]) };

		const scalar dissip(-cellFields.rhoepsTurb()[i]);

		sigmaSourcek.r()[i] = rhoSpherRGen + rhoDevRGen + gravGen + dissip;

		sigmaSourceeps.r()[i] = turbPar->C1() * ek * rhoDevRGen
				+ turbPar->C3() * ek * rhoSpherRGen
				+ turbPar->C0() * ek * gravGen + turbPar->C2() * ek * dissip;

		/*Time-step calculation*/
		modeps[i] = std::abs(
				sourceTimestepCoeff * modeps[i]
						/ (sigmaSourceeps()[i] / cellFields.density[0]()[i]
								+ stabilizator));
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(sigmaSourcek, sigmaSourceeps, sigmaSourcea,
			sigmaSourceb);
}
