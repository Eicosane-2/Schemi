/*
 * diffusiveFields.hpp
 *
 *  Created on: 2020/01/25
 *      Author: Maxim Boldyrev
 *
 *      Structure storing fields for which diffusion equations are solved.
 */

#ifndef DIFFUSIVEFIELDS_HPP_
#define DIFFUSIVEFIELDS_HPP_

#include <vector>

#include "turbulenceModelEnum.hpp"

#include "bunchOfFields.hpp"
#include "volumeField.hpp"

namespace schemi
{
struct diffusiveFields
{
	std::vector<volumeField<scalar>> massFraction;
	volumeField<vector> velocity;
	volumeField<scalar> temperature;
	volumeField<scalar> k;
	volumeField<scalar> eps;
	volumeField<vector> a;
	volumeField<scalar> b;

	diffusiveFields(const mesh & meshRef,
			const bunchOfFields<cubicCell> & cellFields,
			const std::vector<boundaryConditionType> & commBoundCond,
			const bool turbulenceFlag, const turbulenceModelEnum sourceFlag);
};
}  // namespace schemi

#endif /* DIFFUSIVEFIELDS_HPP_ */
