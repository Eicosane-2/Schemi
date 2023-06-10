/*
 * boundaryConditionFromString.hpp
 *
 *  Created on: 2019/12/01
 *      Author: Maxim Boldyrev
 *
 *      Functions converting strings to boundary condition type.
 */

#ifndef BOUNDARYCONDITIONFROMSTRING_HPP_
#define BOUNDARYCONDITIONFROMSTRING_HPP_

#include <string>
#include <vector>

#include "boundaryConditionTypesEnum.hpp"

#include "scalar.hpp"
#include "tensor.hpp"
#include "tensor3.hpp"
#include "vector.hpp"

namespace schemi
{
void boundaryConditionFromString(const std::string & boundaryConditionString,
		const std::string & boundaryConditionValueString,
		boundaryConditionType & boundaryCondition,
		std::vector<scalar> & boundaryConditionValue);

void boundaryConditionFromString(const std::string & boundaryConditionString,
		const std::string & boundaryConditionValueStringX,
		const std::string & boundaryConditionValueStringY,
		const std::string & boundaryConditionValueStringZ,
		boundaryConditionType & boundaryCondition,
		std::vector<vector> & boundaryConditionValue);

void boundaryConditionFromString(const std::string & boundaryConditionString,
		const std::string & boundaryConditionValueStringXX,
		const std::string & boundaryConditionValueStringXY,
		const std::string & boundaryConditionValueStringXZ,
		const std::string & boundaryConditionValueStringYX,
		const std::string & boundaryConditionValueStringYY,
		const std::string & boundaryConditionValueStringYZ,
		const std::string & boundaryConditionValueStringZX,
		const std::string & boundaryConditionValueStringZY,
		const std::string & boundaryConditionValueStringZZ,
		boundaryConditionType & boundaryCondition,
		std::vector<tensor> & boundaryConditionValue);

void boundaryConditionFromString(const std::string & boundaryConditionString,

const std::string & boundaryConditionValueStringXXX,
		const std::string & boundaryConditionValueStringXXY,
		const std::string & boundaryConditionValueStringXXZ,
		const std::string & boundaryConditionValueStringXYX,
		const std::string & boundaryConditionValueStringXYY,
		const std::string & boundaryConditionValueStringXYZ,
		const std::string & boundaryConditionValueStringXZX,
		const std::string & boundaryConditionValueStringXZY,
		const std::string & boundaryConditionValueStringXZZ,

		const std::string & boundaryConditionValueStringYXX,
		const std::string & boundaryConditionValueStringYXY,
		const std::string & boundaryConditionValueStringYXZ,
		const std::string & boundaryConditionValueStringYYX,
		const std::string & boundaryConditionValueStringYYY,
		const std::string & boundaryConditionValueStringYYZ,
		const std::string & boundaryConditionValueStringYZX,
		const std::string & boundaryConditionValueStringYZY,
		const std::string & boundaryConditionValueStringYZZ,

		const std::string & boundaryConditionValueStringZXX,
		const std::string & boundaryConditionValueStringZXY,
		const std::string & boundaryConditionValueStringZXZ,
		const std::string & boundaryConditionValueStringZYX,
		const std::string & boundaryConditionValueStringZYY,
		const std::string & boundaryConditionValueStringZYZ,
		const std::string & boundaryConditionValueStringZZX,
		const std::string & boundaryConditionValueStringZZY,
		const std::string & boundaryConditionValueStringZZZ,

		boundaryConditionType & boundaryCondition,
		std::vector<tensor3> & boundaryConditionValue);
}  // namespace schemi

#endif /* BOUNDARYCONDITIONFROMSTRING_HPP_ */
