/*
 * errorsEnum.hpp
 *
 *  Created on: 2020/06/12
 *      Author: Maxim Boldyrev
 *
 *      Types of errors.
 */

#ifndef ERRORSENUM_HPP_
#define ERRORSENUM_HPP_

namespace schemi
{
enum class errors
{
	noErrors,
	RiemannSolverError,
	negativeTemperatureError,
	tooBigOutputNumberError,
	systemError,
	boundaryConditionError,
	cubicEquationError,
	positivnessError,
	MPIError,
	fieldInitializationError,
	meshGenerationError,
	initializationError,
	NaNError
};
}  // namespace schemi

#endif /* ERRORSENUM_HPP_ */
