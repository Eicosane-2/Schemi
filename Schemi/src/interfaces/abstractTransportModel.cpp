/*
 * abstractTransportModel.cpp
 *
 *  Created on: 2023/08/07
 *      Author: Maxim Boldyrev
 */

#include "abstractTransportModel.hpp"

schemi::scalar schemi::abstractTransportModel::muConst() const noexcept
{
	return mu;
}

schemi::scalar schemi::abstractTransportModel::DConst() const noexcept
{
	return D;
}

schemi::scalar schemi::abstractTransportModel::kappaConst() const noexcept
{
	return kappa;
}

schemi::abstractTransportModel::abstractTransportModel(const scalar mu_in,
		const scalar D_in, const scalar kappa_in) noexcept :
		mu(mu_in), D(D_in), kappa(kappa_in)
{
}

schemi::abstractTransportModel::~abstractTransportModel() noexcept
{
}

schemi::volumeField<schemi::scalar> schemi::abstractTransportModel::calculateMu(
		const std::valarray<scalar>&, const volumeField<scalar> & temperature,
		const concentrationsPack<cubicCell>&) const noexcept
{
	return volumeField<scalar>(temperature.meshRef(), mu);
}

schemi::volumeField<std::valarray<std::valarray<schemi::scalar>>> schemi::abstractTransportModel::calculateDm(
		const std::valarray<scalar>&, const volumeField<scalar> & temperature,
		const volumeField<scalar>&,
		const concentrationsPack<cubicCell> & concentrations) const noexcept
{
	volumeField<std::valarray<std::valarray<scalar>>> Dmatrix(
			temperature.meshRef());

	const std::size_t compNumber { concentrations.v.size() - 1 };

	Dmatrix.ref_r() = std::valarray<std::valarray<scalar>>(
			std::valarray<scalar>(0., compNumber), compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				if (k1 == k2)
					continue;

				Dmatrix.ref_r()[i][k1][k2] = D;
			}
	}

	return Dmatrix;
}

schemi::volumeField<std::valarray<schemi::scalar>> schemi::abstractTransportModel::calculateDs(
		const std::valarray<scalar>&, const volumeField<scalar> & temperature,
		const volumeField<scalar>&,
		const concentrationsPack<cubicCell> & concentrations) const noexcept
{
	return volumeField<std::valarray<scalar>>(temperature.meshRef(),
			std::valarray<scalar>(D, concentrations.v.size() - 1));
}

schemi::volumeField<schemi::scalar> schemi::abstractTransportModel::calculateKappa(
		const std::valarray<scalar>&, const volumeField<scalar> & temperature,
		const concentrationsPack<cubicCell>&) const noexcept
{
	return volumeField<scalar>(temperature.meshRef(), kappa);
}

schemi::surfaceField<schemi::scalar> schemi::abstractTransportModel::calculateMu(
		const std::valarray<scalar>&, const surfaceField<scalar> & temperature,
		const concentrationsPack<quadraticSurface>&) const noexcept
{
	return surfaceField<scalar>(temperature.meshRef(), mu);
}

schemi::surfaceField<std::valarray<std::valarray<schemi::scalar>>> schemi::abstractTransportModel::calculateDm(
		const std::valarray<scalar>&, const surfaceField<scalar> & temperature,
		const surfaceField<scalar>&,
		const concentrationsPack<quadraticSurface> & concentrations) const noexcept
{
	surfaceField<std::valarray<std::valarray<scalar>>> Dmatrix(
			temperature.meshRef());

	const std::size_t compNumber { concentrations.v.size() - 1 };

	Dmatrix.ref_r() = std::valarray<std::valarray<scalar>>(
			std::valarray<scalar>(0., compNumber), compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				if (k1 == k2)
					continue;

				Dmatrix.ref_r()[i][k1][k2] = D;
			}
	}

	return Dmatrix;
}

schemi::surfaceField<std::valarray<schemi::scalar>> schemi::abstractTransportModel::calculateDs(
		const std::valarray<scalar>&, const surfaceField<scalar> & temperature,
		const surfaceField<scalar>&,
		const concentrationsPack<quadraticSurface> & concentrations) const noexcept
{
	return surfaceField<std::valarray<scalar>>(temperature.meshRef(),
			std::valarray<scalar>(D, concentrations.v.size() - 1));
}

schemi::surfaceField<schemi::scalar> schemi::abstractTransportModel::calculateKappa(
		const std::valarray<scalar>&, const surfaceField<scalar> & temperature,
		const concentrationsPack<quadraticSurface>&) const noexcept
{
	return surfaceField<scalar>(temperature.meshRef(), kappa);
}
