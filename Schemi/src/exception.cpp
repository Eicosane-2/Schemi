/*
 * exception.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "exception.hpp"

schemi::exception::exception(const std::string & in_string,
		const errorsEnum in_errType) noexcept :
		ex_txt(in_string), errType(in_errType)
{
}

const char* schemi::exception::what() const noexcept
{
	return ex_txt.c_str();
}
