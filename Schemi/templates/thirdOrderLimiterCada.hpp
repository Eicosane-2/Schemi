/*
 * thirdOrderLimiterCada.hpp
 *
 *  Created on: 2025/05/24
 *      Author: Maxim Boldyrev
 */

#ifndef THIRDORDERLIMITERCADA_HPP_
#define THIRDORDERLIMITERCADA_HPP_

#include "abstractLimiter.hpp"
#include "elementsProduct.hpp"
#include "intExpPow.hpp"
#include "volumeField.hpp"
#include "surfaceField.hpp"

namespace schemi
{
template<typename T>
volumeField<std::vector<T>> thirdOrderLimiterCada(
		const surfaceField<T> & surfGrad,
		const abstractLimiter & limiterObjectP,
		const volumeField<scalar> & minimalLengthScale)
{
	constexpr scalar alpha(pow<scalar, 2>(8 * Pi_number)), eps { 1E-6 };

	auto & mesh_ { surfGrad.meshRef() };

	volumeField<std::vector<T>> limitedGradients(mesh_, std::vector<T>(0));

	for (std::size_t i = 0; i < limitedGradients.size(); ++i)
	{
		const auto & surfacesOfCell_i = mesh_.surfacesOfCells()[i];

		const auto dx2 = pow<scalar, 2>(1 / minimalLengthScale.cval()[i]);

		for (std::size_t s1 = 0; s1 < surfacesOfCell_i.size(); ++s1)
		{
			const auto surfaceIndex { surfacesOfCell_i[s1] };

			std::pair<std::size_t, scalar> farSurface(0, 0);

			for (std::size_t s2 = 0; s2 < surfacesOfCell_i.size(); ++s2)
			{
				if (s1 == s2)
					continue;

				const auto secSurfaceIndex { surfacesOfCell_i[s2] };

				const scalar surfDistance =
						(mesh_.surfaces()[secSurfaceIndex].rC()
								- mesh_.surfaces()[surfaceIndex].rC()).mag();

				if (farSurface.second < surfDistance)
				{
					farSurface.second = surfDistance;
					farSurface.first = secSurfaceIndex;
				}
			}

			const auto & gradNear = surfGrad.cval()[surfaceIndex];
			const auto & gradOppose = surfGrad.cval()[farSurface.first];

			scalar delta2Near, delta2Oppose;

			if constexpr (std::is_same<T, vector>::value)
			{
				delta2Near = dx2 * (gradNear & gradNear);
				delta2Oppose = dx2 * (gradOppose & gradOppose);
			}
			else
			{
				delta2Near = dx2 * (gradNear && gradNear);
				delta2Oppose = dx2 * (gradOppose && gradOppose);
			}

			const scalar eta = (delta2Near + delta2Oppose)
					/ (std::sqrt(5.0 / 2.0) * alpha * dx2);

			if (eta < (1 - eps))
			{
				limitedGradients.val()[i].emplace_back(
						0.5
								* ((1 - onethirds) * gradOppose
										+ (1 + onethirds) * gradNear));
			}
			else if (eta > (1 + eps))
			{
				const auto r = elementsDivision(gradNear, gradOppose);

				const auto retr = elementsDivision(decltype(r)(1), r);

				const auto limGradSurf = limiterObjectP.calculateNoRSLimit(retr,
						gradNear);

				const auto limOppGradSurf = limiterObjectP.calculateNoRSLimit(r,
						gradOppose);

				limitedGradients.val()[i].emplace_back(
						0.5
								* ((1 - onethirds) * limOppGradSurf
										+ (1 + onethirds) * limGradSurf));
			}
			else
			{
				const auto r = elementsDivision(gradNear, gradOppose);

				const auto retr = elementsDivision(decltype(r)(1), r);

				const auto limGradSurf = limiterObjectP.calculateNoRSLimit(retr,
						gradNear);

				const auto limOppGradSurf = limiterObjectP.calculateNoRSLimit(r,
						gradOppose);

				const auto limitedGrad = 0.5
						* ((1 - onethirds) * limOppGradSurf
								+ (1 + onethirds) * limGradSurf);

				const auto unlimitedGrad = 0.5
						* ((1 - onethirds) * gradOppose
								+ (1 + onethirds) * gradNear);

				const auto eta1eps = (eta - 1) / eps;

				limitedGradients.val()[i].emplace_back(
						0.5
								* ((1 - eta1eps) * unlimitedGrad
										+ (1 + eta1eps) * limitedGrad));
			}
		}
	}

	return limitedGradients;
}
}  // namespace schemi

#endif /* THIRDORDERLIMITERCADA_HPP_ */
