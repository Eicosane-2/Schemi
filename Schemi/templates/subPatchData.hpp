/*
 * subPatchData.hpp
 *
 *  Created on: 2023/06/17
 *      Author: Maxim Boldyrev
 */

#ifndef SUBPATCHDATA_HPP_
#define SUBPATCHDATA_HPP_

#include "boundaryConditionTypesEnum.hpp"
#include "vector.hpp"
#include "globalConstants.hpp"

namespace schemi
{
template<typename T>
struct subPatchData
{
	subPatchData() :
			bType(boundaryConditionType::blank), fixVal(0)
	{
	}

	explicit subPatchData(boundaryConditionType boundType, const T & value = T {
			0 }) :
			bType(boundType), fixVal(value)
	{
	}

	subPatchData(boundaryConditionType boundType, const vector & b,
			const vector & e, T value = T { 0 }) :
			bType(boundType), fixVal(value), patchBeg(b), patchEnd(e)
	{
	}

	boundaryConditionType bType { boundaryConditionType::blank };

	T fixVal { 0 };

	vector patchBeg { 0.0, 0.0, 0.0 }, patchEnd { veryBig, veryBig, veryBig };
};
}  // namespace schemi

#endif /* SUBPATCHDATA_HPP_ */
