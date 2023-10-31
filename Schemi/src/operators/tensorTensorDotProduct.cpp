/*
 * tensorTensorDotProduct.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "tensorTensorDotProduct.hpp"

schemi::tensor schemi::operator&(const tensor & inTensor1,
		const tensor & inTensor2) noexcept
{
	return tensor { std::get<0>(inTensor1()) * std::get<0>(inTensor2())
			+ std::get<1>(inTensor1()) * std::get<3>(inTensor2())
			+ std::get<2>(inTensor1()) * std::get<6>(inTensor2()), std::get<0>(
			inTensor1()) * std::get<1>(inTensor2())
			+ std::get<1>(inTensor1()) * std::get<4>(inTensor2())
			+ std::get<2>(inTensor1()) * std::get<7>(inTensor2()), std::get<0>(
			inTensor1()) * std::get<2>(inTensor2())
			+ std::get<1>(inTensor1()) * std::get<5>(inTensor2())
			+ std::get<2>(inTensor1()) * std::get<8>(inTensor2()), std::get<3>(
			inTensor1()) * std::get<0>(inTensor2())
			+ std::get<4>(inTensor1()) * std::get<3>(inTensor2())
			+ std::get<5>(inTensor1()) * std::get<6>(inTensor2()), std::get<3>(
			inTensor1()) * std::get<1>(inTensor2())
			+ std::get<4>(inTensor1()) * std::get<4>(inTensor2())
			+ std::get<5>(inTensor1()) * std::get<7>(inTensor2()), std::get<3>(
			inTensor1()) * std::get<2>(inTensor2())
			+ std::get<4>(inTensor1()) * std::get<5>(inTensor2())
			+ std::get<5>(inTensor1()) * std::get<8>(inTensor2()), std::get<6>(
			inTensor1()) * std::get<0>(inTensor2())
			+ std::get<7>(inTensor1()) * std::get<3>(inTensor2())
			+ std::get<8>(inTensor1()) * std::get<6>(inTensor2()), std::get<6>(
			inTensor1()) * std::get<1>(inTensor2())
			+ std::get<7>(inTensor1()) * std::get<4>(inTensor2())
			+ std::get<8>(inTensor1()) * std::get<7>(inTensor2()), std::get<6>(
			inTensor1()) * std::get<2>(inTensor2())
			+ std::get<7>(inTensor1()) * std::get<5>(inTensor2())
			+ std::get<8>(inTensor1()) * std::get<8>(inTensor2()) };
}
