/*
 * chemicalKineticsNO2Disproportionation.hpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSNO2DISPROPORTIONATION_HPP_
#define CHEMICALKINETICSNO2DISPROPORTIONATION_HPP_

#include "abstractChemicalKinetics.hpp"

#include <functional>

#include "globalConstants.hpp"

namespace schemi
{
namespace chemicalKinetics
{
class NO2Disproportionation: public abstractChemicalKinetics
{
	static constexpr std::size_t N { 4 };

	/*(/ 1E6) --- cm^3/mole/s converted to m^3/mole/s*/
	scalar A_forward { 5.36E-50 * NAvogardro * NAvogardro / 1E12 };
	scalar n_forward { 3.95 };
	scalar E_forward { -1825 * R_SI };

	scalar A_backward { 3.31E-19 / 1E6 * NAvogardro };
	scalar n_backward { 2.478 };
	scalar E_backward { 3199 * R_SI };

	/*ΔH*/
	scalar deltaH_298 = (-76.73 - 134.31 + 241.826 - 2 * 33.1) * 1000;
	scalar deltaU_298 = deltaH_298;
	constexpr static scalar deltan = -1;

	iterativeSolver itSolv;

	class cellReactionMatrix
	{
		iterativeSolver solverFlag;

		struct reactionMatrix
		{
			std::array<scalar, N> Diagonale { 0.0, 0.0, 0.0, 0.0 };

			std::array<triangleList, N> LeftTriangle {

			triangleList(0),

			triangleList(1, std::make_pair(0.0, 0)),

			triangleList { std::make_pair(0.0, 0), std::make_pair(0.0, 1) },

			triangleList { std::make_pair(0.0, 0), std::make_pair(0.0, 1),
					std::make_pair(0.0, 2) }

			};

			std::array<triangleList, N> RightTriangle {

			triangleList { std::make_pair(0.0, 1), std::make_pair(0.0, 2),
					std::make_pair(0.0, 3) },

			triangleList { std::make_pair(0.0, 2), std::make_pair(0.0, 3) },

			triangleList { std::make_pair(0.0, 3) },

			triangleList(0)

			};

			std::array<scalar, N> FreeTerm { 0.0, 0.0, 0.0, 0.0 };

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

		void setMatrix(const scalar timeStep, const scalar k_forw,
				const scalar k_backw, const scalar C_NO2_0,
				const scalar C_H2O_0, const scalar C_HNO2_0,
				const scalar C_HNO3_0, const scalar rho_0,
				const std::array<scalar, N> & molMass) noexcept;
	public:
		cellReactionMatrix() noexcept;

		explicit cellReactionMatrix(const iterativeSolver solverType);

		cellReactionMatrix(const scalar timeStep, const scalar k_forw,
				const scalar k_backw, const scalar C_NO2_0,
				const scalar C_H2O_0, const scalar C_HNO2_0,
				const scalar C_HNO3_0, const scalar rho_0,
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
	NO2Disproportionation(const homogeneousPhase<cubicCell> & phaseIn,
			const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) override;
};
}
}  // namespace schemi

#endif /* CHEMICALKINETICSNO2DISPROPORTIONATION_HPP_ */
