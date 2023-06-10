/*
 * thirdOrderAveraging.hpp
 *
 *  Created on: 2023/01/08
 *      Author: Maxim Boldyrev
 */

#ifndef THIRDORDERAVERAGING_HPP_
#define THIRDORDERAVERAGING_HPP_

#include "elementsProduct.hpp"

namespace schemi
{
template<typename T>
volumeField<std::vector<T>> thirdOrderAveraging(const surfaceField<T> & inField,
		const abstractLimiter & limiterObjectP)
{
	auto & mesh { inField.meshRef() };

	volumeField < std::vector < T >> returnField(mesh, std::vector < T > (0));

	for (std::size_t i = 0; i < returnField.size(); ++i)
	{
		const auto nSurfs = mesh.surfacesOfCells()[i].size();

		returnField.ref_r()[i].resize(nSurfs);

		for (std::size_t j = 0; j < nSurfs; ++j)
		{
			const auto surfIndex_j = mesh.surfacesOfCells()[i][j];

			std::pair < std::size_t, scalar > maxDistData { 0, 0. };

			for (std::size_t k = 0; k < nSurfs; ++k)
			{
				const auto surfIndex_k = mesh.surfacesOfCells()[i][k];

				if (surfIndex_k == surfIndex_j)
					continue;

				const scalar surfDelta = (mesh.surfaces()[surfIndex_j].rC()
						- mesh.surfaces()[surfIndex_k].rC()).mag();

				if (surfDelta > maxDistData.second)
					maxDistData = std::make_pair(surfIndex_k, surfDelta);
			}

			const auto r = elementsDivision(inField.ref()[surfIndex_j],
					inField.ref()[maxDistData.first]);

			const auto r1 = elementsDivision(T(1.), r);

			const auto surfTVDGrad = limiterObjectP.calculateNoRightLimit(r1,
					inField.ref()[surfIndex_j]);
			const auto farSurfTVDGrad = limiterObjectP.calculateNoRightLimit(r,
					inField.ref()[maxDistData.first]);

			returnField.ref_r()[i][j] = (1 - kappaPPM) / 2 * farSurfTVDGrad
					+ (1 + kappaPPM) / 2 * surfTVDGrad;
		}
	}

	return returnField;
}
}  // namespace schemi

#endif /* THIRDORDERAVERAGING_HPP_ */
