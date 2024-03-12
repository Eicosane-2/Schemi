/*
 * boundaryConditionFromString.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "boundaryConditionFromString.hpp"

#include "exception.hpp"

schemi::boundaryConditionType schemi::boundaryConditionFromString(
		const std::string_view boundaryConditionString)
{
	if (boundaryConditionString == "blank")
		return boundaryConditionType::blank;
	else if (boundaryConditionString == "freeBoundary")
		return boundaryConditionType::freeBoundary;
	else if (boundaryConditionString == "slip")
		return boundaryConditionType::slip;
	else if (boundaryConditionString == "fixedValueCell")
		return boundaryConditionType::fixedValueCell;
	else if (boundaryConditionString == "fixedValueSurface")
		return boundaryConditionType::fixedValueSurface;
	else if (boundaryConditionString == "innerSurface")
		throw exception("<<innerSurface>> can not be boundary surface type.",
				errors::boundaryConditionError);
	else
		throw exception("Unknown type of boundary condition",
				errors::boundaryConditionError);
}
