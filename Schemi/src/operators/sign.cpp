/*
 * sign.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "sign.hpp"

#include "exception.hpp"
#include "globalConstants.hpp"

int schemi::sign(const scalar inValue)
{
	if (inValue > zeroLevel)
		return 1;
	else if (inValue < -zeroLevel)
		return -1;
	else if ((inValue <= zeroLevel) && (inValue >= -zeroLevel))
		return 0;
	else
		throw exception(
				std::string("Signum. Unknown value: ")
						+ std::to_string(inValue), errorsEnum::NaNError);
}
