/*
 * line.hpp
 *
 *  Created on: 2024/10/10
 *      Author: Maxim Boldyrev
 */

#ifndef LINE_HPP_
#define LINE_HPP_

#include "vector.hpp"

namespace schemi
{
class line
{
	vector tangent { 0, 0, 0 };
	vector point { 0, 0, 0 };
public:
	line() noexcept = default;

	line(const vector & t, const vector & p) noexcept;

	scalar distanceFromPoint(const vector & p) const noexcept;
};
}

#endif /* LINE_HPP_ */
