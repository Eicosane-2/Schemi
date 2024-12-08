/*
 * boundaryConditionValue.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "boundaryConditionValue.hpp"

schemi::boundaryConditionValue::boundaryConditionValue(
		const abstractTurbulenceModel & turb_in,
		const bunchOfFields<cubicCell> & cellFields_in,
		const abstractMixtureThermodynamics & mix_in,
		const MPIHandler & parallelism_in) noexcept :
		turb(turb_in), cellFields(cellFields_in), mix(mix_in), meshReference(
				cellFields.pressure.meshRef()), parallelism(parallelism_in)
{
}

