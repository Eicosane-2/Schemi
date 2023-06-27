/*
 * boundaryConditionFromString.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "boundaryConditionFromString.hpp"

#include "exception.hpp"

schemi::boundaryConditionType schemi::boundaryConditionFromString(
		const std::string & boundaryConditionString)
{
	if (boundaryConditionString == "blank")
		return boundaryConditionType::blank;
	else if (boundaryConditionString == "freeBoundary")
		return boundaryConditionType::freeBoundary;
	else if (boundaryConditionString == "slip")
		return boundaryConditionType::slip;
	else if (boundaryConditionString == "fixedValue")
		return boundaryConditionType::fixedValue;
	else if (boundaryConditionString == "innerSurface")
		throw exception("<<innerSurface>> can not be boundary surface type.",
				errorsEnum::boundaryConditionError);
	else
		throw exception("Unknown type of boundary condition",
				errorsEnum::boundaryConditionError);
}
