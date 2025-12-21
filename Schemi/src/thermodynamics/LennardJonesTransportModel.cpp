/*
 * LennardJonesTransportModel.cpp
 *
 *  Created on: 2025/11/27
 *      Author: Maxim Boldyrev
 */

#include "LennardJonesTransportModel.hpp"

#include "intExpPow.hpp"

schemi::scalar schemi::LennardJonesTransportModel::OmegaCalcMuKappa(
		const scalar Ts) const noexcept
{
	constexpr scalar Ac { 1.16145 }, Bc { 0.14874 }, Cc { 0.52487 }, Dc {
			0.77320 }, Ec { 2.16178 }, Fc { 2.4378 };

	return Ac * std::pow(Ts, -Bc) + Cc * std::exp(-Dc * Ts)
			+ Ec * std::exp(-Fc * Ts);
}
schemi::scalar schemi::LennardJonesTransportModel::OmegaCalcD(
		const scalar Ts) const noexcept
{
	constexpr scalar Ac { 1.06036 }, Bc { 0.1561 }, Cc { 0.19300 },
			Dc { 0.47635 }, Ec { 1.03587 }, Fc { 1.52996 }, Gc { 1.76474 }, Hc {
					3.89411 };

	return Ac * std::pow(Ts, -Bc) + Cc * std::exp(-Dc * Ts)
			+ Ec * std::exp(-Fc * Ts) + Gc * std::exp(-Hc * Ts);
}

schemi::scalar schemi::LennardJonesTransportModel::muLJ(const scalar M,
		const scalar T, const scalar Omega, const scalar sigma) const noexcept
{
	const auto calcVal = 2.6693E-6 * std::sqrt(1E3) * std::sqrt(M * T) / Omega
			/ pow<scalar, 2>(sigma);

	return muConst() + calcVal;
}

schemi::scalar schemi::LennardJonesTransportModel::DLJ(const scalar M,
		const scalar T, const scalar p, const scalar Omega,
		const scalar sigma) const noexcept
{
	const auto calcVal = 2.6280E-7 / std::sqrt(1E3) / Omega
			* std::sqrt(pow<scalar, 3>(T) / M)
			/ (p / 101325. * pow<scalar, 2>(sigma));

	return DConst() + calcVal;
}

schemi::scalar schemi::LennardJonesTransportModel::kappaLJ(const scalar M,
		const scalar T, const scalar Omega, const scalar sigma) const noexcept
{
	const auto calcVal = 4.184 * 1.9891E-2 / std::sqrt(1E3) * std::sqrt(T / M)
			/ Omega / pow<scalar, 2>(sigma);

	return kappaConst() + calcVal;
}

schemi::LennardJonesTransportModel::LennardJonesTransportModel(
		const scalar mu_in, const scalar D_in, const scalar kappa_in,
		const std::valarray<scalar> & epsilonLJk_in,
		const std::valarray<scalar> & sigmaLJ_in) noexcept :
		abstractTransportModel(mu_in, D_in, kappa_in), epsilonLJk(
				epsilonLJk_in), sigmaLJ(sigmaLJ_in)
{
}

schemi::volumeField<schemi::scalar> schemi::LennardJonesTransportModel::calculateMu(
		const std::valarray<scalar> & M,
		const volumeField<scalar> & temperature,
		const concentrationsPack<cubicCell> & concentrations) const noexcept
{
	const std::size_t compNumber = concentrations.v.size() - 1;

	volumeField<scalar> muField(temperature.meshRef());

	std::valarray<scalar> xArray(0., compNumber);

	std::valarray<scalar> muArray(0., compNumber);

	auto Phi = std::valarray<std::valarray<scalar>>(
			std::valarray<scalar>(0., compNumber), compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k = 0; k < compNumber; ++k)
			xArray[k] = { concentrations.v[k + 1].cval()[i]
					/ concentrations.v[0].cval()[i] };

		for (std::size_t k = 0; k < compNumber; ++k)
		{
			const auto Omega = OmegaCalcMuKappa(
					temperature.cval()[i] / epsilonLJk[k]);

			muArray[k] = muLJ(M[k], temperature.cval()[i], Omega, sigmaLJ[k]);
		}

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				const scalar num { 1
						+ std::sqrt(muArray[k1] / muArray[k2])
								* std::sqrt(std::sqrt(M[k2] / M[k1])) };
				const scalar dem { 8 * (1 + M[k1] / M[k2]) };

				Phi[k1][k2] = num * num / std::sqrt(dem);
			}

		scalar muMix { 0 };

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
		{
			scalar Phi_k1 { 0 };

			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
				Phi_k1 += xArray[k2] * Phi[k1][k2];

			muMix += xArray[k1] * muArray[k1] / Phi_k1;
		}

		muField.val()[i] = muMix;
	}

	return muField;
}

schemi::volumeField<std::valarray<std::valarray<schemi::scalar>>> schemi::LennardJonesTransportModel::calculateDm(
		const std::valarray<scalar> & M,
		const volumeField<scalar> & temperature,
		const volumeField<scalar> & pressure,
		const concentrationsPack<cubicCell> & concentrations) const noexcept
{
	const std::size_t compNumber = concentrations.v.size() - 1;

	volumeField<std::valarray<std::valarray<scalar>> > Dmatrix(
			temperature.meshRef());

	Dmatrix.val() = std::valarray<std::valarray<scalar>>(
			std::valarray<scalar>(0., compNumber), compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				if (k1 == k2)
					continue;

				const scalar molMass12 = 2 * M[k1] * M[k2] / (M[k1] + M[k2]);
				const scalar Omega = OmegaCalcD(
						temperature.cval()[i]
								/ std::sqrt(epsilonLJk[k1] * epsilonLJk[k2]));
				const scalar sigma12 = 0.5 * (sigmaLJ[k1] + sigmaLJ[k2]);

				Dmatrix.val()[i][k1][k2] = DLJ(molMass12, temperature.cval()[i],
						pressure.cval()[i], Omega, sigma12);
			}
	}

	return Dmatrix;
}

schemi::volumeField<std::valarray<schemi::scalar>> schemi::LennardJonesTransportModel::calculateDs(
		const std::valarray<scalar> & M,
		const volumeField<scalar> & temperature,
		const volumeField<scalar> & pressure,
		const concentrationsPack<cubicCell> & concentrations) const noexcept
{
	const std::size_t compNumber = concentrations.v.size() - 1;

	volumeField<std::valarray<scalar>> Dsingle(temperature.meshRef(),
			std::valarray<scalar>(0., compNumber));

	std::valarray<scalar> xArray(0., compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k = 0; k < compNumber; ++k)
			xArray[k] = { concentrations.v[k + 1].cval()[i]
					/ concentrations.v[0].cval()[i] };

		std::valarray<std::valarray<scalar>> Dmatrix { std::valarray<scalar>(0.,
				compNumber), compNumber };

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				if (k1 == k2)
					continue;

				const scalar molMass12 = 2 * M[k1] * M[k2] / (M[k1] + M[k2]);
				const scalar Omega = OmegaCalcD(
						temperature.cval()[i]
								/ std::sqrt(epsilonLJk[k1] * epsilonLJk[k2]));
				const scalar sigma12 = 0.5 * (sigmaLJ[k1] + sigmaLJ[k2]);

				Dmatrix[k1][k2] = DLJ(molMass12, temperature.cval()[i],
						pressure.cval()[i], Omega, sigma12);
			}

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
		{
			scalar rDk { 0 };

			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				if (k1 == k2)
					continue;

				rDk += xArray[k2] / Dmatrix[k1][k2];
			}

			rDk *= 1.0 / (1.0 - xArray[k1]);

			Dsingle.val()[i][k1] = 1 / rDk;
		}
	}

	return Dsingle;
}

schemi::volumeField<schemi::scalar> schemi::LennardJonesTransportModel::calculateKappa(
		const std::valarray<scalar> & M,
		const volumeField<scalar> & temperature,
		const concentrationsPack<cubicCell> & concentrations) const noexcept
{
	const std::size_t compNumber = concentrations.v.size() - 1;

	volumeField<scalar> kappaField(temperature.meshRef());

	std::valarray<scalar> xArray(0., compNumber);

	std::valarray<scalar> kappaArray(0., compNumber);

	auto Phi = std::valarray<std::valarray<scalar>>(
			std::valarray<scalar>(0., compNumber), compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k = 0; k < compNumber; ++k)
			xArray[k] = { concentrations.v[k + 1].cval()[i]
					/ concentrations.v[0].cval()[i] };

		for (std::size_t k = 0; k < compNumber; ++k)
		{
			const auto Omega = OmegaCalcMuKappa(
					temperature.cval()[i] / epsilonLJk[k]);

			kappaArray[k] = kappaLJ(M[k], temperature.cval()[i], Omega,
					sigmaLJ[k]);
		}

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				const scalar num { 1
						+ std::sqrt(kappaArray[k1] / kappaArray[k2])
								* std::sqrt(std::sqrt(M[k2] / M[k1])) };
				const scalar dem { 8 * (1 + M[k1] / M[k2]) };

				Phi[k1][k2] = 0.85 * num * num / std::sqrt(dem);
			}

		scalar kappaMix { 0 };

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
		{
			scalar Phi_k1 { 0 };

			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
				Phi_k1 += xArray[k2] * Phi[k1][k2];

			kappaMix += xArray[k1] * kappaArray[k1] / Phi_k1;
		}

		kappaField.val()[i] = kappaMix;
	}

	return kappaField;
}

schemi::surfaceField<schemi::scalar> schemi::LennardJonesTransportModel::calculateMu(
		const std::valarray<scalar> & M,
		const surfaceField<scalar> & temperature,
		const concentrationsPack<quadraticSurface> & concentrations) const noexcept
{
	const std::size_t compNumber = concentrations.v.size() - 1;

	surfaceField<scalar> muField(temperature.meshRef());

	std::valarray<scalar> xArray(0., compNumber);

	std::valarray<scalar> muArray(0., compNumber);

	auto Phi = std::valarray<std::valarray<scalar>>(
			std::valarray<scalar>(0., compNumber), compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k = 0; k < compNumber; ++k)
			xArray[k] = { concentrations.v[k + 1].cval()[i]
					/ concentrations.v[0].cval()[i] };

		for (std::size_t k = 0; k < compNumber; ++k)
		{
			const auto Omega = OmegaCalcMuKappa(
					temperature.cval()[i] / epsilonLJk[k]);

			muArray[k] = muLJ(M[k], temperature.cval()[i], Omega, sigmaLJ[k]);
		}

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				const scalar num { 1
						+ std::sqrt(muArray[k1] / muArray[k2])
								* std::sqrt(std::sqrt(M[k2] / M[k1])) };
				const scalar dem { 8 * (1 + M[k1] / M[k2]) };

				Phi[k1][k2] = num * num / std::sqrt(dem);
			}

		scalar muMix { 0 };

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
		{
			scalar Phi_k1 { 0 };

			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
				Phi_k1 += xArray[k2] * Phi[k1][k2];

			muMix += xArray[k1] * muArray[k1] / Phi_k1;
		}

		muField.val()[i] = muMix;
	}

	return muField;
}

schemi::surfaceField<std::valarray<std::valarray<schemi::scalar>>> schemi::LennardJonesTransportModel::calculateDm(
		const std::valarray<scalar> & M,
		const surfaceField<scalar> & temperature,
		const surfaceField<scalar> & pressure,
		const concentrationsPack<quadraticSurface> & concentrations) const noexcept
{
	const std::size_t compNumber = concentrations.v.size() - 1;

	surfaceField<std::valarray<std::valarray<scalar>> > Dmatrix(
			temperature.meshRef());

	Dmatrix.val() = std::valarray<std::valarray<scalar>>(
			std::valarray<scalar>(0., compNumber), compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				if (k1 == k2)
					continue;

				const scalar molMass12 = 2 * M[k1] * M[k2] / (M[k1] + M[k2]);
				const scalar Omega = OmegaCalcD(
						temperature.cval()[i]
								/ std::sqrt(epsilonLJk[k1] * epsilonLJk[k2]));
				const scalar sigma12 = 0.5 * (sigmaLJ[k1] + sigmaLJ[k2]);

				Dmatrix.val()[i][k1][k2] = DLJ(molMass12, temperature.cval()[i],
						pressure.cval()[i], Omega, sigma12);
			}
	}

	return Dmatrix;
}

schemi::surfaceField<std::valarray<schemi::scalar>> schemi::LennardJonesTransportModel::calculateDs(
		const std::valarray<scalar> & M,
		const surfaceField<scalar> & temperature,
		const surfaceField<scalar> & pressure,
		const concentrationsPack<quadraticSurface> & concentrations) const noexcept
{
	const std::size_t compNumber = concentrations.v.size() - 1;

	surfaceField<std::valarray<scalar>> Dsingle(temperature.meshRef(),
			std::valarray<scalar>(0., compNumber));

	std::valarray<scalar> xArray(0., compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k = 0; k < compNumber; ++k)
			xArray[k] = { concentrations.v[k + 1].cval()[i]
					/ concentrations.v[0].cval()[i] };

		std::valarray<std::valarray<scalar>> Dmatrix { std::valarray<scalar>(0.,
				compNumber), compNumber };

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				if (k1 == k2)
					continue;

				const scalar molMass12 = 2 * M[k1] * M[k2] / (M[k1] + M[k2]);
				const scalar Omega = OmegaCalcD(
						temperature.cval()[i]
								/ std::sqrt(epsilonLJk[k1] * epsilonLJk[k2]));
				const scalar sigma12 = 0.5 * (sigmaLJ[k1] + sigmaLJ[k2]);

				Dmatrix[k1][k2] = DLJ(molMass12, temperature.cval()[i],
						pressure.cval()[i], Omega, sigma12);
			}

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
		{
			scalar rDk { 0 };

			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				if (k1 == k2)
					continue;

				rDk += xArray[k2] / Dmatrix[k1][k2];
			}

			if (rDk >= stabilizator)
			{
				rDk *= 1.0 / (1.0 - xArray[k1] + stabilizator);

				Dsingle.val()[i][k1] = 1 / rDk;
			}
			else
				Dsingle.val()[i][k1] = stabilizator;
		}
	}

	return Dsingle;
}

schemi::surfaceField<schemi::scalar> schemi::LennardJonesTransportModel::calculateKappa(
		const std::valarray<scalar> & M,
		const surfaceField<scalar> & temperature,
		const concentrationsPack<quadraticSurface> & concentrations) const noexcept
{
	const std::size_t compNumber = concentrations.v.size() - 1;

	surfaceField<scalar> kappaField(temperature.meshRef());

	std::valarray<scalar> xArray(0., compNumber);

	std::valarray<scalar> kappaArray(0., compNumber);

	auto Phi = std::valarray<std::valarray<scalar>>(
			std::valarray<scalar>(0., compNumber), compNumber);

	for (std::size_t i = 0; i < temperature.size(); ++i)
	{
		for (std::size_t k = 0; k < compNumber; ++k)
			xArray[k] = { concentrations.v[k + 1].cval()[i]
					/ concentrations.v[0].cval()[i] };

		for (std::size_t k = 0; k < compNumber; ++k)
		{
			const auto Omega = OmegaCalcMuKappa(
					temperature.cval()[i] / epsilonLJk[k]);

			kappaArray[k] = kappaLJ(M[k], temperature.cval()[i], Omega,
					sigmaLJ[k]);
		}

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
			{
				const scalar num { 1
						+ std::sqrt(kappaArray[k1] / kappaArray[k2])
								* std::sqrt(std::sqrt(M[k2] / M[k1])) };
				const scalar dem { 8 * (1 + M[k1] / M[k2]) };

				Phi[k1][k2] = 0.85 * num * num / std::sqrt(dem);
			}

		scalar kappaMix { 0 };

		for (std::size_t k1 = 0; k1 < compNumber; ++k1)
		{
			scalar Phi_k1 { 0 };

			for (std::size_t k2 = 0; k2 < compNumber; ++k2)
				Phi_k1 += xArray[k2] * Phi[k1][k2];

			kappaMix += xArray[k1] * kappaArray[k1] / Phi_k1;
		}

		kappaField.val()[i] = kappaMix;
	}

	return kappaField;
}
