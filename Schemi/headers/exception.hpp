/*
 * exception.hpp
 *
 *  Created on: 2020/07/07
 *      Author: Maxim Boldyrev
 *
 *      Project's exception class.
 */

#ifndef EXCEPTION_HPP_
#define EXCEPTION_HPP_

#include <exception>
#include <string>

#include "errorsEnum.hpp"

namespace schemi
{
class exception: public std::exception
{
	const std::string ex_txt;

public:
	const errorsEnum errType;

	exception(const std::string & in_string,
			const errorsEnum in_errType) noexcept;

	const char* what() const noexcept override;
};
}  // namespace schemi

#endif /* EXCEPTION_HPP_ */
