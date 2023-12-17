/*
 * chemicalKineticsСhlorumDissociation.hpp
 *
 *  Created on: 2023/05/09
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSCHLORUMDISSOCIATION_
#define CHEMICALKINETICSCHLORUMDISSOCIATION_

#include "abstractChemicalKinetics.hpp"

namespace schemi
{
class chemicalKineticsChlorumDissociation: public abstractChemicalKinetics
{
	// (/ 1E6) --- cm^3/mole/s converted to m^3/mole/s
	scalar A_forw { 5.74 * 1E15 / 1E6 };
	scalar n_forw { 0 };
	scalar E_forw { 55.596 * 1000 * 4.184 };

	scalar A_backw { 2.3 * 1E19 / 1E12 };
	scalar n_backw { -1.5 };
	scalar E_backw { 0 };

	iterativeSolver itSolv;

	class cellReactionMatrix
	{
		iterativeSolver solverFlag;

		struct reactionMatrix
		{
			std::array<scalar, 2> Diagonale { 0.0, 0.0 };

			std::array<triangleList, 2> LeftTriangle {

			triangleList(0),

			triangleList(1, std::make_pair(0.0, 0))

			};

			std::array<triangleList, 2> RightTriangle {

			triangleList(1, std::make_pair(0.0, 1)),

			triangleList(0)

			};

			std::array<scalar, 2> FreeTerm { 0.0, 0.0 };

			void transpose() noexcept;
		} matrix;

		std::valarray<scalar> matrixDotProduct(const reactionMatrix & m,
				const std::valarray<scalar> & v) const noexcept;

		auto solveJ(const std::array<scalar, 2> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 2>;

		auto solveGS(const std::array<scalar, 2> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 2>;

		auto solveCG(const std::array<scalar, 2> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 2>;

		auto solveJCG(const std::array<scalar, 2> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 2>;

		auto solveGE() const -> std::array<scalar, 2>;
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const scalar timeStep, const scalar k_diss,
				const scalar k_recomb, const scalar C_Cl2_0,
				const scalar C_Cl_0, const scalar M_0, const scalar rho_0,
				const std::array<scalar, 2> & molMass,
				const iterativeSolver solverType);

		auto solve(const std::array<scalar, 2> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 2>;
	};

	cellReactionMatrix velocityCalculation(const scalar timestep,
			const scalar T, const std::array<scalar, 3> & concentrations,
			const std::array<scalar, 2> & molarMasses, const scalar rho,
			const scalar R) const noexcept;

	void timeStepIntegration(homogeneousPhase<cubicCell> & phaseN) const;
public:
	chemicalKineticsChlorumDissociation(
			const homogeneousPhase<cubicCell> & phaseIn, const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) const
			override;
};
}  // namespace schemi

#endif /* CHEMICALKINETICSCHLORUMDISSOCIATION_ */
