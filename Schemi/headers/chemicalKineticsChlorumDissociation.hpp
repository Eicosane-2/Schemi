/*
 * chemicalKineticsСhlorumDissociation.hpp
 *
 *  Created on: 2023/05/09
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSCHLORUMDISSOCIATION_
#define CHEMICALKINETICSCHLORUMDISSOCIATION_

#include "abstractChemicalKinetics.hpp"

#include <functional>

namespace schemi
{
namespace chemicalKinetics
{
class ChlorumDissociation: public abstractChemicalKinetics
{
	static constexpr std::size_t N { 2 };

	/*(/ 1E6) --- cm^3/mole/s converted to m^3/mole/s*/
	/*(/ 1E12) --- cm^6/mole^2/s converted to m^6/mole^2/s*/
	scalar A_forw { 5.74E15 / 1E6 };
	scalar n_forw { 0 };
	scalar E_forw { 55.596 * 1000 * 4.184 };

	scalar A_backw { 2.3E19 / 1E12 };
	scalar n_backw { -1.5 };
	scalar E_backw { 0 };

	iterativeSolver itSolv;

	class cellReactionMatrix
	{
		iterativeSolver solverFlag;

		struct reactionMatrix
		{
			std::array<scalar, N> Diagonale { 0.0, 0.0 };

			std::array<triangleList, N> LeftTriangle {

			triangleList(0),

			triangleList(1, std::make_pair(0.0, 0))

			};

			std::array<triangleList, N> RightTriangle {

			triangleList(1, std::make_pair(0.0, 1)),

			triangleList(0)

			};

			std::array<scalar, N> FreeTerm { 0.0, 0.0 };

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

		void setMatrix(const scalar timeStep, const scalar k_diss,
				const scalar k_recomb, const scalar C_Cl2_0,
				const scalar C_Cl_0, const scalar M_0, const scalar rho_0,
				const std::array<scalar, N> & molMas) noexcept;
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const iterativeSolver solverType);

		cellReactionMatrix(const scalar timeStep, const scalar k_diss,
				const scalar k_recomb, const scalar C_Cl2_0,
				const scalar C_Cl_0, const scalar M_0, const scalar rho_0,
				const std::array<scalar, N> & molMass,
				const iterativeSolver solverType);

		auto solve(const std::array<scalar, N> & oldField,
				const std::size_t maxIterationNumber) -> std::array<scalar, N>;

		void velocityCalculation(const scalar timestep, const scalar T,
				const std::array<scalar, N + 1> & concentrations,
				const std::array<scalar, N> & molarMasses, const scalar rho,
				const scalar R, const kineticParams & forw,
				const kineticParams & backw) noexcept;
	};

	cellReactionMatrix cellReactionVel;

	void timeStepIntegration(homogeneousPhase<cubicCell> & phaseN);
public:
	ChlorumDissociation(const homogeneousPhase<cubicCell> & phaseIn,
			const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) override;
};
}
}  // namespace schemi

#endif /* CHEMICALKINETICSCHLORUMDISSOCIATION_ */
