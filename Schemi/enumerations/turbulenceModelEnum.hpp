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
	BHRModel,
	zeroModel,
	decayModel,
	shearModel,
	arithmeticA1Model,
	arithmeticA2Model,
	arithmeticA3Model,
	kEpsAModel,
	BHRKLModel,
	BHR2Model,
	BHR3Model,
	unknownModel
};
}  // namespace schemi

#endif /* TURBULENCEMODELENUM_HPP_ */
