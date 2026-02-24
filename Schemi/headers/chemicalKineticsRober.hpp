/*
 * chemicalKineticsRober.hpp
 *
 *  Created on: 2024/06/01
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSROBER_HPP_
#define CHEMICALKINETICSROBER_HPP_

#include "abstractChemicalKinetics.hpp"

#include <functional>

namespace schemi
{
namespace chemicalKinetics
{
class Rober: public abstractChemicalKinetics
{
	static constexpr std::size_t N { 3 };

	/*(/ 1E6) --- cm^3/mole/s converted to m^3/mole/s*/
	/*(/ 1E12) --- cm^6/mole^2/s converted to m^6/mole^2/s*/
	scalar A_1 { 0.04 };
	scalar n_1 { 0 };
	scalar E_1 { 0 };

	scalar A_2 { 3E7 };
	scalar n_2 { 0 };
	scalar E_2 { 0 };

	scalar A_3 { 1E4 };
	scalar n_3 { 0 };
	scalar E_3 { 0 };

	iterativeSolver itSolv;

	class cellReactionMatrix
	{
		iterativeSolver solverFlag;

		struct reactionMatrix
		{
			std::array<scalar, N> Diagonale { 0.0, 0.0, 0.0 };

			std::array<triangleList, N> LeftTriangle {

			triangleList(0),

			triangleList(1, std::make_pair(0.0, 0)),

			triangleList(1, std::make_pair(0.0, 1))

			};

			std::array<triangleList, N> RightTriangle {

			triangleList( { std::make_pair(0.0, 1), std::make_pair(0.0, 2) }),

			triangleList(1, std::make_pair(0.0, 2)),

			triangleList(0)

			};

			std::array<scalar, N> FreeTerm { 0.0, 0.0, 0.0 };

			void transpose() noexcept;
		} matrix;

		std::array<scalar, N> (*my_solveJ)(const reactionMatrix&,
				const std::array<scalar, N>&,
				const std::size_t) = solveJ<reactionMatrix, N>;
		std::array<scalar, N> (*my_solveGS)(const reactionMatrix&,
				const std::array<scalar, N>&,
				const std::size_t) = solveGS<reactionMatrix, N>;
		std::array<scalar, N> (*my_solveCG)(const reactionMatrix&,
				const std::array<scalar, N>&,
				const std::size_t) = solveCG<reactionMatrix, N>;
		std::array<scalar, N> (*my_solveJCG)(const reactionMatrix&,
				const std::array<scalar, N>&,
				const std::size_t) = solveJCG<reactionMatrix, N>;
		std::array<scalar, N> (*my_solveGE)(const reactionMatrix&,
				const std::array<scalar, N>&,
				const std::size_t) = solveGE<reactionMatrix, N>;

		std::function<
				std::array<scalar, N>(const reactionMatrix&,
						const std::array<scalar, N>&, const std::size_t)> solverF;

		void setMatrix(const scalar timeStep, const scalar k_1,
				const scalar k_2, const scalar k_3, const scalar C_1_0,
				const scalar C_2_0, const scalar C_3_0, const scalar rho_0,
				const std::array<scalar, N> & molMass) noexcept;
	public:
		cellReactionMatrix() noexcept;

		explicit cellReactionMatrix(const iterativeSolver solverType);

		cellReactionMatrix(const scalar timeStep, const scalar k_1,
				const scalar k_2, const scalar k_3, const scalar C_1_0,
				const scalar C_2_0, const scalar C_3_0, const scalar rho_0,
				const std::array<scalar, N> & molMass,
				const iterativeSolver solverType);

		auto solve(const std::array<scalar, N> & oldField,
				const std::size_t maxIterationNumber) -> std::array<scalar, N>;

		void velocityCalculation(const scalar timestep, const scalar T,
				const std::array<scalar, N + 1> & concentrations,
				const std::array<scalar, N> & molarMasses, const scalar rho,
				const scalar R, const kineticParams & eq1,
				const kineticParams & eq2, const kineticParams & eq3) noexcept;
	};

	cellReactionMatrix cellReactionVel;

	void timeStepIntegration(homogeneousPhase<cubicCell> & phaseN);
public:
	Rober(const homogeneousPhase<cubicCell> & phaseIn, const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) override;
};
}
}  // namespace schemi

#endif /* CHEMICALKINETICSROBER_HPP_ */
