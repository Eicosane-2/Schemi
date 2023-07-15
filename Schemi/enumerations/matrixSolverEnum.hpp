/*
 * matrixSolverEnum.hpp
 *
 *  Created on: 2019/12/06
 *      Author: Maxim Boldyrev
 *
 *      Type of matrix diagonalization solvers.
 */

#ifndef MATRIXSOLVERENUM_HPP_
#define MATRIXSOLVERENUM_HPP_

namespace schemi
{
enum class matrixSolverEnum
{
	GaussSeidel,
	Thomas,
	noSovler,
	explicitSolver,
	conjugateGradient,
	JacobiConjugateGradient
};
}  // namespace schemi

#endif /* MATRIXSOLVERENUM_HPP_ */
