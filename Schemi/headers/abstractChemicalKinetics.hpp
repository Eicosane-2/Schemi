/*
 * abstractChemicalKinetics.hpp
 *
 *  Created on: 2023/05/08
 *      Author: Maxim Boldyrev
 */

#ifndef ABSTRACTCHEMICALKINETICS_HPP_
#define ABSTRACTCHEMICALKINETICS_HPP_

#include <cstddef>
#include <vector>
#include <utility>

#include "cubicCell.hpp"
#include "scalar.hpp"
#include "homogeneousPhase.hpp"

namespace schemi
{
class abstractChemicalKinetics
{
protected:
	typedef std::vector<std::pair<scalar, std::size_t>> triangleList;

	constexpr static scalar convergenceTolerance { 1E-10 };
	constexpr static scalar massFracTolerance { 1E-2 };
	std::size_t maxIterationNumber { 0 };
	const scalar minTimestep { 0 };

	enum class iterativeSolver
	{
		noSolver,
		GaussSeidel,
		ConjugateGradient,
		JacobiConjugateGradient,
		Jacobi,
		GaussElimination
	};

	struct cellReactingFields
	{
		scalar internalEnergy;

		scalar temperature;

		std::valarray<scalar> concentration;

		std::valarray<scalar> density;

		cellReactingFields(scalar internalEnergy_in, scalar temperature_in,
				std::valarray<scalar> concentration_in,
				std::valarray<scalar> density_in) noexcept :
				internalEnergy(internalEnergy_in), temperature(temperature_in), concentration(
						concentration_in), density(density_in)
		{
		}
	};

public:
	const bool chemicalReaction;

	abstractChemicalKinetics(const bool flag, const scalar mt) noexcept;

	virtual ~abstractChemicalKinetics() noexcept =0;

	virtual void solveChemicalKinetics(homogeneousPhase<cubicCell>&) const =0;
};
}  // namespace schemi

#endif /* ABSTRACTCHEMICALKINETICS_HPP_ */
