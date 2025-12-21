/*
 * output.hpp
 *
 *  Created on: 2019/12/07
 *      Author: Maxim Boldyrev
 *
 *      Functions for output to text files.
 */

#ifndef OUTPUT_HPP_
#define OUTPUT_HPP_

#include <cstddef>

#include "abstractTurbulenceModel.hpp"
#include "scalar.hpp"
#include "structForOutput.hpp"

namespace schemi
{
namespace output
{
void dataOutput(const structForOutput & outputData, const std::size_t nOutput,
		const scalar Time, const abstractTurbulenceModel & turb);

void mixedZoneWidth1D(const structForOutput & outputData, const scalar Time);
}  // namespace output
}  // namespace schemi

#endif /* OUTPUT_HPP_ */
