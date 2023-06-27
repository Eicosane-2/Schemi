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

#include "boundaryConditionTypesEnum.hpp"

namespace schemi
{
boundaryConditionType boundaryConditionFromString(
		const std::string & boundaryConditionString);
}  // namespace schemi

#endif /* BOUNDARYCONDITIONFROMSTRING_HPP_ */
