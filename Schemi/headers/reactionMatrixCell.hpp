/*
 * reactionMatrixCell.hpp
 *
 *  Created on: 2026/07/07
 *      Author: Maxim Boldyrev
 */

#ifndef REACTIONMATRIXCELL_HPP_
#define REACTIONMATRIXCELL_HPP_

#include <array>
#include <vector>

#include "scalar.hpp"

namespace schemi
{
namespace chemicalKinetics
{

struct reactionMatrixCell
{
	std::array<std::size_t, 2> indexes { 0, 0 };
	bool nullCell { true };
	bool diagonalCell { false };
	std::pair<scalar, std::vector<scalar>> molMass {0.0, 0};
	std::vector<std::array<scalar, 3>> reactParams { };
	std::vector<scalar> reactWeight { };
	std::vector<std::vector<std::pair<std::size_t, int>>> comp { };
	std::vector<int> coeffs { };
};

}
}

#endif /* REACTIONMATRIXCELL_HPP_ */
