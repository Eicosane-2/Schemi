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

std::tuple<
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		std::pair<schemi::volumeField<schemi::vector>,
				schemi::volumeField<schemi::vector>>,
		std::pair<schemi::volumeField<schemi::scalar>,
				schemi::volumeField<schemi::scalar>>,
		schemi::volumeField<schemi::scalar>> schemi::zeroGen::calculate(
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

	const std::pair<volumeField<scalar>, volumeField<scalar>> Sourcek {
			volumeField<scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	const std::pair<volumeField<scalar>, volumeField<scalar>> Sourceeps {
			volumeField<scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	const std::pair<volumeField<vector>, volumeField<vector>> Sourcea {
			volumeField<vector>(mesh_, vector(0)), volumeField<vector>(mesh_,
					vector(0)) };
	const std::pair<volumeField<scalar>, volumeField<scalar>> Sourceb {
			volumeField<scalar>(mesh_, 0), volumeField<scalar>(mesh_, 0) };
	const volumeField<scalar> gravGenField(mesh_, 0);

	sourceTimestep = std::min(mesh_.timestepSource(), veryBig);

	return std::make_tuple(Sourcek, Sourceeps, Sourcea, Sourceb, gravGenField);
}
