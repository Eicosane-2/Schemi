/*
 * SLEMatrixDotProduct.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "SLEMatrixDotProduct.hpp"

std::valarray<schemi::scalar> schemi::operator&(
		const SLEMatrix::SLEMatrixStorage & M,
		const std::valarray<scalar> & v) noexcept
{
	std::valarray<scalar> result(v.size());

	for (std::size_t i = 0; i < v.size(); ++i)
	{
		auto i_res = M.centralDiagonale[i] * v[i];

		const auto & lowTri = M.lowerTriangle[i];

		for (std::size_t j = 0; j < lowTri.size(); ++j)
		{
			const auto index = lowTri[j].second;

			i_res += lowTri[j].first * v[index];
		}

		const auto & upTri = M.upperTriangle[i];

		for (std::size_t j = 0; j < upTri.size(); ++j)
		{
			const auto index = upTri[j].second;

			i_res += upTri[j].first * v[index];
		}

		result[i] = i_res;
	}

	return result;
}
