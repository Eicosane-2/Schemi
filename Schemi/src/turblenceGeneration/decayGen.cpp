/*
 * decayGen.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "decayGen.hpp"

#include <algorithm>

#include "turbulentParametersKEPS.hpp"

schemi::decayGen::decayGen(const mesh & meshIn, const bool turb_in,
		const turbulenceModel tm_in) noexcept :
		abstractTurbulenceGen(turb_in, tm_in)
{
	turbPar = std::make_unique<turbulentParametersKEPS>(meshIn);
}

std::tuple<schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::vector>, schemi::volumeField<schemi::scalar>> schemi::decayGen::calculate(
		scalar & sourceTimestep, const scalar sourceTimestepCoeff,
		const bunchOfFields<cubicCell> & cellFields,
		const diffusiveFields & diffFieldsOld, const volumeField<tensor>&,
		const volumeField<vector>&, const volumeField<vector>&,
		const volumeField<vector>&, const volumeField<tensor>&,
		const volumeField<scalar>&, const volumeField<vector>&,
		const volumeField<tensor>&, const volumeField<tensor>&,
		const volumeField<vector>&, const abstractMixtureThermodynamics&,
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

	for (std::size_t i = 0; i < mesh_.cellsSize(); ++i)
	{
		const scalar ek { diffFieldsOld.eps()[i] / diffFieldsOld.k()[i] };

		const scalar dissip(-cellFields.rhoepsTurb()[i]);

		sigmaSourcek.r()[i] = dissip;

		sigmaSourceeps.r()[i] = turbPar->C2() * ek * dissip;

		/*Time-step calculation*/
		modeps[i] = std::abs(
				sourceTimestepCoeff * modeps[i]
						/ (sigmaSourceeps()[i] / cellFields.density[0]()[i]
								+ stabilizator));

		sigmaSourcea.r()[i] = (-1.) * cellFields.rhoaTurb()[i] * ek
				* turbPar->Ca1();

		sigmaSourceb.r()[i] = -cellFields.rhobTurb()[i] * ek * turbPar->Cb1();
	}

	sourceTimestep = std::min(mesh_.timestepSource(), modeps.min());

	return std::make_tuple(sigmaSourcek, sigmaSourceeps, sigmaSourcea,
			sigmaSourceb);
}
