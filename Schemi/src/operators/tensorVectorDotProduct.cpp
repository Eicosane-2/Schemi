/*
 * tensorVectorDotProduct.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensorVectorDotProduct.hpp"

schemi::vector schemi::operator&(const tensor & inTensor,
		const vector & inVector) noexcept
{
	return vector { std::get<0>(inTensor()) * std::get<0>(inVector())
			+ std::get<1>(inTensor()) * std::get<1>(inVector())
			+ std::get<2>(inTensor()) * std::get<2>(inVector()), std::get<3>(
			inTensor()) * std::get<0>(inVector())
			+ std::get<4>(inTensor()) * std::get<1>(inVector())
			+ std::get<5>(inTensor()) * std::get<2>(inVector()), std::get<6>(
			inTensor()) * std::get<0>(inVector())
			+ std::get<7>(inTensor()) * std::get<1>(inVector())
			+ std::get<8>(inTensor()) * std::get<2>(inVector()) };
}

schemi::vector schemi::operator&(const vector & inVector,
		const tensor & inTensor) noexcept
{
	tensor bufTensor(inTensor);
	bufTensor.transpose();

	return bufTensor & inVector;
}

schemi::tensor schemi::operator&(const tensor3 & inTensor,
		const vector & inVector) noexcept
{
	return tensor { std::get<0>(inTensor()) * std::get<0>(inVector())
			+ std::get<1>(inTensor()) * std::get<1>(inVector())
			+ std::get<2>(inTensor()) * std::get<2>(inVector()), std::get<3>(
			inTensor()) * std::get<0>(inVector())
			+ std::get<4>(inTensor()) * std::get<1>(inVector())
			+ std::get<5>(inTensor()) * std::get<2>(inVector()), std::get<6>(
			inTensor()) * std::get<0>(inVector())
			+ std::get<7>(inTensor()) * std::get<1>(inVector())
			+ std::get<8>(inTensor()) * std::get<2>(inVector()),

	std::get<9 + 0>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 1>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 2>(inTensor()) * std::get<2>(inVector()), std::get<
			9 + 3>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 4>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 5>(inTensor()) * std::get<2>(inVector()), std::get<
			9 + 6>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 7>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 8>(inTensor()) * std::get<2>(inVector()),

	std::get<18 + 0>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 1>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 2>(inTensor()) * std::get<2>(inVector()), std::get<
			18 + 3>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 4>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 5>(inTensor()) * std::get<2>(inVector()), std::get<
			18 + 6>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 7>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 8>(inTensor()) * std::get<2>(inVector()) };
}

schemi::tensor schemi::operator&(const vector & inVector,
		const tensor3 & inTensor) noexcept
{
	return tensor { std::get<0>(inTensor()) * std::get<0>(inVector())
			+ std::get<3>(inTensor()) * std::get<1>(inVector())
			+ std::get<6>(inTensor()) * std::get<2>(inVector()), std::get<1>(
			inTensor()) * std::get<0>(inVector())
			+ std::get<4>(inTensor()) * std::get<1>(inVector())
			+ std::get<7>(inTensor()) * std::get<2>(inVector()), std::get<2>(
			inTensor()) * std::get<0>(inVector())
			+ std::get<5>(inTensor()) * std::get<1>(inVector())
			+ std::get<8>(inTensor()) * std::get<2>(inVector()),

	std::get<9 + 0>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 3>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 6>(inTensor()) * std::get<2>(inVector()), std::get<
			9 + 1>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 4>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 7>(inTensor()) * std::get<2>(inVector()), std::get<
			9 + 2>(inTensor()) * std::get<0>(inVector())
			+ std::get<9 + 5>(inTensor()) * std::get<1>(inVector())
			+ std::get<9 + 8>(inTensor()) * std::get<2>(inVector()),

	std::get<18 + 0>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 3>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 6>(inTensor()) * std::get<2>(inVector()), std::get<
			18 + 1>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 4>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 7>(inTensor()) * std::get<2>(inVector()), std::get<
			18 + 2>(inTensor()) * std::get<0>(inVector())
			+ std::get<18 + 5>(inTensor()) * std::get<1>(inVector())
			+ std::get<18 + 8>(inTensor()) * std::get<2>(inVector()) };
}
