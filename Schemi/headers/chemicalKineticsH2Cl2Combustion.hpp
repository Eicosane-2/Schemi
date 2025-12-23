/*
 * chemicalKineticsH2Cl2Combustion.hpp
 *
 *  Created on: 2023/05/18
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSH2CL2COMBUSTION_HPP_
#define CHEMICALKINETICSH2CL2COMBUSTION_HPP_

#include "abstractChemicalKinetics.hpp"

#include <functional>

namespace schemi
{
namespace chemicalKinetics
{
class H2Cl2Combustion: public abstractChemicalKinetics
{
	static constexpr std::size_t N { 5 };

	/*(/ 1E6) --- cm^3/mole/s converted to m^3/mole/s*/
	/*(/ 1E12) --- cm^6/mole^2/s converted to m^6/mole^2/s*/
	scalar A_Cl2_diss { 5.74E15 / 1E6 };
	scalar n_Cl2_diss { 0 };
	scalar E_Cl2_diss { 55.596 * 1000 * 4.184 };

	scalar A_Cl_recomb { 2.3E19 / 1E12 };
	scalar n_Cl_recomb { -1.5 };
	scalar E_Cl_recomb { 0 };

	scalar A_H2_diss { 1.8E17 / 1E6 };
	scalar n_H2_diss { -0.5 };
	scalar E_H2_diss { 454.0 * 1000 };

	scalar A_H_recomb { 6.53E17 / 1E12 };
	scalar n_H_recomb { -1 };
	scalar E_H_recomb { 0 };

	scalar A_Cl_H2_prop { 2.88E8 / 1E6 };
	scalar n_Cl_H2_prop { 1.58 };
	scalar E_Cl_H2_prop { 3.199 * 1000 * 4.184 };

	scalar A_H_Cl2_prop { 9.64E13 / 1E6 };
	scalar n_H_Cl2_prop { 0 };
	scalar E_H_Cl2_prop { 1.68 * 1000 * 4.184 };

	scalar A_H_Cl_recomb { 2E23 / 1E12 };
	scalar n_H_Cl_recomb { -2.45 };
	scalar E_H_Cl_recomb { 0 };

	scalar A_HCl_diss { 4.4E13 / 1E6 };
	scalar n_HCl_diss { 0 };
	scalar E_HCl_diss { 81.753 * 1000 * 4.184 };

	/*ΔH*/
	scalar deltaH_298 = -92300;
	scalar deltaU_298 = deltaH_298;
	constexpr static scalar deltan = 0;

	iterativeSolver itSolv;

	class cellReactionMatrix
	{
		iterativeSolver solverFlag;

		struct reactionMatrix
		{
			std::array<scalar, N> Diagonale { 0.0, 0.0, 0.0, 0.0, 0.0 };

			std::array<triangleList, N> LeftTriangle {

			triangleList(0), // A1

			triangleList(1, std::make_pair(0.0, 0)), // A2

			triangleList(1, std::make_pair(0.0, 1)), // A3

			triangleList { std::make_pair(0.0, 0), std::make_pair(0.0, 1),
					std::make_pair(0.0, 2) }, // A4

					triangleList { std::make_pair(0.0, 0), std::make_pair(0.0,
							1), std::make_pair(0.0, 2), std::make_pair(0.0, 3) } // A5

			};

			std::array<triangleList, N> RightTriangle {

			triangleList { std::make_pair(0.0, 1), std::make_pair(0.0, 3) }, // A1

					triangleList { std::make_pair(0.0, 2), std::make_pair(0.0,
							3), std::make_pair(0.0, 4) }, // A2

					triangleList(1, std::make_pair(0.0, 3)), // A3

					triangleList(1, std::make_pair(0.0, 4)), // A4

					triangleList(0) // A5

			};

			std::array<scalar, N> FreeTerm { 0.0, 0.0, 0.0, 0.0, 0.0 };

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
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const iterativeSolver solverType);

		cellReactionMatrix(const scalar timeStep, const scalar k_diss_Cl2,
				const scalar k_recomb_Cl, const scalar k_diss_H2,
				const scalar k_recomb_H, const scalar k_prop_Cl_H2,
				const scalar k_prop_H_Cl2, const scalar k_recomb_H_Cl,
				const scalar k_diss_HCl, const scalar C_Cl2_0,
				const scalar C_Cl_0, const scalar C_H2_0, const scalar C_H_0,
				const scalar C_HCl_0, const scalar M_0, const scalar rho_0,
				const std::array<scalar, N> & molMass,
				const iterativeSolver solverType);

		auto solve(const std::array<scalar, N> & oldField,
				const std::size_t maxIterationNumber) -> std::array<scalar, N>;

		const reactionMatrix& getMaxtrix() const noexcept
		{
			return matrix;
		}

		void extractMatrix(const cellReactionMatrix & inReactionMatrix)
		{
			matrix = inReactionMatrix.getMaxtrix();
		}
	};

	cellReactionMatrix cellReactionVel;

	cellReactionMatrix velocityCalculation(const scalar timestep,
			const scalar T, const std::array<scalar, N + 1> & concentrations,
			const std::array<scalar, N> & molarMasses, const scalar rho,
			const scalar R) noexcept;

	void timeStepIntegration(homogeneousPhase<cubicCell> & phaseN);
public:
	H2Cl2Combustion(const homogeneousPhase<cubicCell> & phaseIn,
			const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) override;
};
}
}  // namespace schemi

#endif /* CHEMICALKINETICSH2CL2COMBUSTION_HPP_ */
