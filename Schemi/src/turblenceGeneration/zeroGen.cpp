/*
 * zeroGen.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "zeroGen.hpp"

#include "turbulentParametersKEPS.hpp"

schemi::zeroGen::zeroGen(const mesh & meshIn, const bool turb_in,
		const turbulenceModel tm_in) noexcept :
		abstractTurbulenceGen(turb_in, tm_in)
{
	turbPar = std::make_unique<turbulentParametersKEPS>(meshIn);
}

std::tuple<schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::scalar>,
		schemi::volumeField<schemi::vector>, schemi::volumeField<schemi::scalar>> schemi::zeroGen::calculate(
		scalar & sourceTimestep, const scalar,
		const bunchOfFields<cubicCell> & cellFields, const diffusiveFields&,
		const volumeField<tensor>&, const volumeField<vector>&,
		const volumeField<vector>&, const volumeField<vector>&,
		const volumeField<tensor>&, const volumeField<scalar>&,
		const volumeField<vector>&, const volumeField<tensor>&,
		const volumeField<tensor>&, const volumeField<vector>&,
		const abstractMixtureThermodynamics&,
		const volumeField<scalar>&) const noexcept
{
	auto & mesh_ { cellFields.pressure.meshRef() };

	volumeField<scalar> sigmaSourcek(mesh_, 0);
	volumeField<scalar> sigmaSourceeps(mesh_, 0);
	volumeField<vector> sigmaSourcea(mesh_, vector(0));
	volumeField<scalar> sigmaSourceb(mesh_, 0);

	sourceTimestep = std::min(mesh_.timestepSource(), veryBig);

	return std::make_tuple(sigmaSourcek, sigmaSourceeps, sigmaSourcea,
			sigmaSourceb);
}
