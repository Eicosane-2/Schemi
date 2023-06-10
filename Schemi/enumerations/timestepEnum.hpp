/*
 * timestepEnum.hpp
 *
 *  Created on: 2022/01/06
 *      Author: Maxim Boldyrev
 *
 *      Enum for types of timestep calculation.
 */

#ifndef TIMESTEPENUM_HPP_
#define TIMESTEPENUM_HPP_

enum class timestepEnum
{
	CourantTimeStep,
	CourantAndSourceTimeStep,
	CourantAndSourceAndDiffusionTimeStep
};

#endif /* TIMESTEPENUM_HPP_ */
