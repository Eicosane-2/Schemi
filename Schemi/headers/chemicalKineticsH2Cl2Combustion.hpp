/*
 * chemicalKineticsH2Cl2Combustion.hpp
 *
 *  Created on: 2023/05/18
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSH2CL2COMBUSTION_HPP_
#define CHEMICALKINETICSH2CL2COMBUSTION_HPP_

#include "abstractChemicalKinetics.hpp"

namespace schemi
{
class chemicalKineticsH2Cl2Combustion: public abstractChemicalKinetics
{
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
	scalar E_H2_diss { 454 * 1000 };

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

	scalar ΔН_298 = -92300;
	scalar ΔU_298 = ΔН_298;
	constexpr static scalar Δn = 0;

	iterativeSolver itSolv;

	class cellReactionMatrix
	{
		iterativeSolver solverFlag;

		struct reactionMatrix
		{
			std::array<scalar, 5> Diagonale { 0.0, 0.0, 0.0, 0.0, 0.0 };

			std::array<triangleList, 5> LeftTriangle {

			triangleList(0), // A1

			triangleList(1, std::make_pair(0.0, 0)), // A2

			triangleList(1, std::make_pair(0.0, 1)), // A3

			triangleList { std::make_pair(0.0, 0), std::make_pair(0.0, 1),
					std::make_pair(0.0, 2) }, // A4

					triangleList { std::make_pair(0.0, 0), std::make_pair(0.0,
							1), std::make_pair(0.0, 2), std::make_pair(0.0, 3) } // A5

			};

			std::array<triangleList, 5> RightTriangle {

			triangleList { std::make_pair(0.0, 1), std::make_pair(0.0, 3) }, // A1

					triangleList { std::make_pair(0.0, 2), std::make_pair(0.0,
							3), std::make_pair(0.0, 4) }, // A2

					triangleList(1, std::make_pair(0.0, 3)), // A3

					triangleList(1, std::make_pair(0.0, 4)), // A4

					triangleList(0) // A5

			};

			std::array<scalar, 5> FreeTerm { 0.0, 0.0, 0.0, 0.0, 0.0 };

			void transpose() noexcept;
		} matrix;

		std::valarray<scalar> matrixDotProduct(const reactionMatrix & m,
				const std::valarray<scalar> & v) const noexcept;

		auto solveJ(const std::array<scalar, 5> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 5>;

		auto solveGS(const std::array<scalar, 5> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 5>;

		auto solveCG(const std::array<scalar, 5> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 5>;

		auto solveJCG(const std::array<scalar, 5> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 5>;

		auto solveGE() const -> std::array<scalar, 5>;
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const scalar timeStep, const scalar k_diss_Cl2,
				const scalar k_recomb_Cl, const scalar k_diss_H2,
				const scalar k_recomb_H, const scalar k_prop_Cl_H2,
				const scalar k_prop_H_Cl2, const scalar k_recomb_H_Cl,
				const scalar k_diss_HCl, const scalar C_Cl2_0,
				const scalar C_Cl_0, const scalar C_H2_0, const scalar C_H_0,
				const scalar C_HCl_0, const scalar M_0, const scalar rho_0,
				const std::array<scalar, 5> & molMass,
				const iterativeSolver solverType);

		auto solve(const std::array<scalar, 5> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 5>;
	};

	cellReactionMatrix velocityCalculation(const scalar timestep,
			const scalar T, const std::array<scalar, 6> & concentrations,
			const std::array<scalar, 5> & molarMasses, const scalar rho,
			const scalar R) const noexcept;

	void timeStepIntegration(homogeneousPhase<cubicCell> & phaseN) const;
public:
	chemicalKineticsH2Cl2Combustion(const homogeneousPhase<cubicCell> & phaseIn,
			const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) const
			override;
};
}  // namespace schemi

#endif /* CHEMICALKINETICSH2CL2COMBUSTION_HPP_ */
