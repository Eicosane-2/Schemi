/*
 * turbulenceModelEnum.hpp
 *
 *  Created on: 2020/03/12
 *      Author: Maxim Boldyrev
 *
 *      Types of turbulence model's turbulence generation.
 */

#ifndef TURBULENCEMODELENUM_HPP_
#define TURBULENCEMODELENUM_HPP_

namespace schemi
{
enum class turbulenceModel
{
	BHRSource,
	zeroSource,
	decaySource,
	shearSource,
	arithmeticA1Source,
	arithmeticA2Source,
	kEpsASource,
	arithmeticA3Source,
	BHRKLSource
};
}  // namespace schemi

#endif /* TURBULENCEMODELENUM_HPP_ */
