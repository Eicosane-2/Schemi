/*
 * chemicalKineticsRober.hpp
 *
 *  Created on: 2024/06/01
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSROBER_HPP_
#define CHEMICALKINETICSROBER_HPP_

#include "abstractChemicalKinetics.hpp"

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

		std::valarray<scalar> matrixDotProduct(const reactionMatrix & m,
				const std::valarray<scalar> & v) const noexcept;

		auto solveJ(const std::array<scalar, N> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, N>;

		auto solveGS(const std::array<scalar, N> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, N>;

		auto solveCG(const std::array<scalar, N> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, N>;

		auto solveJCG(const std::array<scalar, N> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, N>;

		auto solveGE() const -> std::array<scalar, N>;
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const scalar timeStep, const scalar k_1,
				const scalar k_2, const scalar k_3, const scalar C_1_0,
				const scalar C_2_0, const scalar C_3_0, const scalar M_0,
				const scalar rho_0, const std::array<scalar, N> & molMass,
				const iterativeSolver solverType);

		auto solve(const std::array<scalar, N> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, N>;
	};

	cellReactionMatrix velocityCalculation(const scalar timestep,
			const scalar T, const std::array<scalar, N + 1> & concentrations,
			const std::array<scalar, N> & molarMasses, const scalar rho,
			const scalar R) const noexcept;

	void timeStepIntegration(homogeneousPhase<cubicCell> & phaseN) const;
public:
	Rober(const homogeneousPhase<cubicCell> & phaseIn, const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) const
			override;
};
}
}  // namespace schemi

#endif /* CHEMICALKINETICSROBER_HPP_ */
