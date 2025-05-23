/*
 * thirdOrderLimiter.hpp
 *
 *  Created on: 2023/11/04
 *      Author: Maxim Boldyrev
 */

#ifndef THIRDORDERLIMITER_HPP_
#define THIRDORDERLIMITER_HPP_

#include "abstractLimiter.hpp"
#include "elementsProduct.hpp"
#include "volumeField.hpp"
#include "surfaceField.hpp"

namespace schemi
{
template<typename T>
volumeField<std::vector<T>> thirdOrderLimiter(const surfaceField<T> & surfGrad,
		const abstractLimiter & limiterObjectP)
{
	auto & mesh_ { surfGrad.meshRef() };

	volumeField<std::vector<T>> limitedGradients(mesh_, std::vector<T>(0));

	for (std::size_t i = 0; i < limitedGradients.size(); ++i)
	{
		const auto & surfacesOfCell_i = mesh_.surfacesOfCells()[i];

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

			const auto r = elementsDivision(surfGrad()[surfaceIndex],
					surfGrad()[farSurface.first]);

			const auto retr = elementsDivision(decltype(r)(1), r);

			const auto limGradSurf = limiterObjectP.calculateNoRSLimit(retr,
					surfGrad()[surfaceIndex]);

			const auto limOppGradSurf = limiterObjectP.calculateNoRSLimit(r,
					surfGrad()[farSurface.first]);

			limitedGradients.r()[i].emplace_back(
					0.5
							* ((1 - onethirds) * limOppGradSurf
									+ (1 + onethirds) * limGradSurf));
		}
	}

	return limitedGradients;
}
}  // namespace schemi

#endif /* THIRDORDERLIMITER_HPP_ */
