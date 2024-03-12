/*
 * dyadicProduct.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "dyadicProduct.hpp"

schemi::tensor schemi::operator*(const vector & vector1,
		const vector & vector2) noexcept
{
	return tensor { std::get<0>(vector1()) * std::get<0>(vector2()),
			std::get<0>(vector1()) * std::get<1>(vector2()), std::get<0>(
					vector1()) * std::get<2>(vector2()), std::get<1>(vector1())
					* std::get<0>(vector2()), std::get<1>(vector1())
					* std::get<1>(vector2()), std::get<1>(vector1())
					* std::get<2>(vector2()), std::get<2>(vector1())
					* std::get<0>(vector2()), std::get<2>(vector1())
					* std::get<1>(vector2()), std::get<2>(vector1())
					* std::get<2>(vector2()) };
}

schemi::tensor3 schemi::operator*(const tensor & inTensor,
		const vector & inVector) noexcept
{
	return tensor3 {

	std::get<0>(inTensor()) * std::get<0>(inVector()), std::get<0>(inTensor())
			* std::get<1>(inVector()), std::get<0>(inTensor())
			* std::get<2>(inVector()), std::get<1>(inTensor())
			* std::get<0>(inVector()), std::get<1>(inTensor())
			* std::get<1>(inVector()), std::get<1>(inTensor())
			* std::get<2>(inVector()), std::get<2>(inTensor())
			* std::get<0>(inVector()), std::get<2>(inTensor())
			* std::get<1>(inVector()), std::get<2>(inTensor())
			* std::get<2>(inVector()),

	std::get<3>(inTensor()) * std::get<0>(inVector()), std::get<3>(inTensor())
			* std::get<1>(inVector()), std::get<3>(inTensor())
			* std::get<2>(inVector()), std::get<4>(inTensor())
			* std::get<0>(inVector()), std::get<4>(inTensor())
			* std::get<1>(inVector()), std::get<4>(inTensor())
			* std::get<2>(inVector()), std::get<5>(inTensor())
			* std::get<0>(inVector()), std::get<5>(inTensor())
			* std::get<1>(inVector()), std::get<5>(inTensor())
			* std::get<2>(inVector()),

	std::get<6>(inTensor()) * std::get<0>(inVector()), std::get<6>(inTensor())
			* std::get<1>(inVector()), std::get<6>(inTensor())
			* std::get<2>(inVector()), std::get<7>(inTensor())
			* std::get<0>(inVector()), std::get<7>(inTensor())
			* std::get<1>(inVector()), std::get<7>(inTensor())
			* std::get<2>(inVector()), std::get<8>(inTensor())
			* std::get<0>(inVector()), std::get<8>(inTensor())
			* std::get<1>(inVector()), std::get<8>(inTensor())
			* std::get<2>(inVector()) };
}

schemi::tensor3 schemi::operator*(const vector & inVector,
		const tensor & inTensor) noexcept
{
	return tensor3 {

	std::get<0>(inVector()) * std::get<0>(inTensor()), std::get<0>(inVector())
			* std::get<1>(inTensor()), std::get<0>(inVector())
			* std::get<2>(inTensor()), std::get<0>(inVector())
			* std::get<3>(inTensor()), std::get<0>(inVector())
			* std::get<4>(inTensor()), std::get<0>(inVector())
			* std::get<5>(inTensor()), std::get<0>(inVector())
			* std::get<6>(inTensor()), std::get<0>(inVector())
			* std::get<7>(inTensor()), std::get<0>(inVector())
			* std::get<8>(inTensor()),

	std::get<1>(inVector()) * std::get<0>(inTensor()), std::get<1>(inVector())
			* std::get<1>(inTensor()), std::get<1>(inVector())
			* std::get<2>(inTensor()), std::get<1>(inVector())
			* std::get<3>(inTensor()), std::get<1>(inVector())
			* std::get<4>(inTensor()), std::get<1>(inVector())
			* std::get<5>(inTensor()), std::get<1>(inVector())
			* std::get<6>(inTensor()), std::get<1>(inVector())
			* std::get<7>(inTensor()), std::get<1>(inVector())
			* std::get<8>(inTensor()),

	std::get<2>(inVector()) * std::get<0>(inTensor()), std::get<2>(inVector())
			* std::get<1>(inTensor()), std::get<2>(inVector())
			* std::get<2>(inTensor()), std::get<2>(inVector())
			* std::get<3>(inTensor()), std::get<2>(inVector())
			* std::get<4>(inTensor()), std::get<2>(inVector())
			* std::get<5>(inTensor()), std::get<2>(inVector())
			* std::get<6>(inTensor()), std::get<2>(inVector())
			* std::get<7>(inTensor()), std::get<2>(inVector())
			* std::get<8>(inTensor()) };
}
