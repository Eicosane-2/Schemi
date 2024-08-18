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

	volumeField<std::vector<T>> cellLimiters(mesh_, std::vector<T>(0));

	for (std::size_t i = 0; i < cellLimiters.size(); ++i)
	{
		const auto & surfacesOfCell_i = mesh_.surfacesOfCells()[i];

		for (std::size_t s1 = 0; s1 < surfacesOfCell_i.size(); ++s1)
		{
			const auto surfaceIndex { surfacesOfCell_i[s1] };

			std::size_t farSurfaceIndex { 0 };

			scalar farSurfaceDistance { 0 };

			for (std::size_t s2 = 0; s2 < surfacesOfCell_i.size(); ++s2)
			{
				if (s1 == s2)
					continue;

				const auto counterSurfaceIndex { surfacesOfCell_i[s2] };

				const scalar surfDistance =
						(mesh_.surfaces()[counterSurfaceIndex].rC()
								- mesh_.surfaces()[surfaceIndex].rC()).mag();

				if (farSurfaceDistance < surfDistance)
				{
					farSurfaceDistance = surfDistance;
					farSurfaceIndex = counterSurfaceIndex;
				}
			}

			const auto r = elementsDivision(surfGrad()[surfaceIndex],
					surfGrad()[farSurfaceIndex]);

			cellLimiters.r()[i].emplace_back(
					limiterObjectP.calculateNoRSLimit(r,
							surfGrad()[farSurfaceIndex]));
		}
	}

	return cellLimiters;
}
}  // namespace schemi

#endif /* THIRDORDERLIMITER_HPP_ */
