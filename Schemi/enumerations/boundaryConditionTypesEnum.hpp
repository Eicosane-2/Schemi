/*
 * boundaryConditionTypesEnum.hpp
 *
 *  Created on: 2019/11/14
 *      Author: Maxim Boldyrev
 *
 *      Types of boundary conditions.
 */

#ifndef BOUNDARYCONDITIONTYPESENUM_HPP_
#define BOUNDARYCONDITIONTYPESENUM_HPP_

namespace schemi
{
enum class boundaryConditionType
{
	blank,
	freeBoundary,
	slip,
	fixedValue,
	innerSurface,
	calculated,
	calculatedTurbulentViscosity,
	calculatedTemperature,
	calculatedMassFraction,
	calculatedAverageMolarMass,
	calculatedNonidealityCorrectionPerDensity,
	calculatedCv,
	calculatedCvM,
	calculatedDensity,
	calculatedParallelBoundary,
	calculatedMolarFraction
};
}  // namespace schemi

#endif /* BOUNDARYCONDITIONTYPESENUM_HPP_ */
