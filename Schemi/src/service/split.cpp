/*
 * split.cpp
 *
 *  Created on: 2026/06/25
 *      Author: Maxim Boldyrev
 */

#include "split.hpp"

#include <algorithm>

#include "exception.hpp"

std::array<std::string, 2> schemi::split(const std::string & str,
		const char delim)
{
	const auto delimPos = std::find(str.cbegin(), str.cend(), delim);

	if (delimPos == str.cend())
		throw exception("Could not find delimeter.",
				errors::initialisationError);

	auto next = delimPos;
	next++;

	return
	{	std::string(str.cbegin(), delimPos), std::string(next, str.cend())};
}
