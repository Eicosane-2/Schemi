/*
 * boundaryConditionValue.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "boundaryConditionValue.hpp"

schemi::boundaryConditionValue::boundaryConditionValue(
		const abstractTurbulentParameters & turbPar_in,
		const bunchOfFields<cubicCell> & cellFields_in,
		const abstractMixtureThermodynamics & mix_in,
		const MPIHandler & parallelism_in) noexcept :
		turbPar(turbPar_in), cellFields(cellFields_in), mix(mix_in), meshReference(
				cellFields.pressure.meshRef()), parallelism(parallelism_in)
{
}

