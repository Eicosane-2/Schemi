/*
 * chemicalKineticsH2O2Combustion.hpp
 *
 *  Created on: 2024/02/11
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSH2O2COMBUSTION_HPP_
#define CHEMICALKINETICSH2O2COMBUSTION_HPP_

#include "abstractChemicalKinetics.hpp"

namespace schemi
{
class chemicalKineticsH2O2Combustion: public abstractChemicalKinetics
{
	/*(/ 1E6) --- cm^3/mole/s converted to m^3/mole/s*/
	/*(/ 1E12) --- cm^6/mole^2/s converted to m^6/mole^2/s*/
	scalar A_R1 { 2.2E13 / 1E6 };
	scalar n_R1 { 0.5 };
	scalar E_R1 { 21500 * R_SI };

	scalar A_R2 { 7.6E13 / 1E6 };
	scalar n_R2 { 0.5 };
	scalar E_R2 { 29500 * R_SI };

	scalar A_R3 { 2.2E12 / 1E6 };
	scalar n_R3 { 0.5 };
	scalar E_R3 { 46600 * R_SI };

	scalar A_R4 { 1.8E11 / 1E6 };
	scalar n_R4 { 0.5 };
	scalar E_R4 { 48100 * R_SI };

	scalar A_R5 { 1.81E10 / 1E6 };
	scalar n_R5 { 1 };
	scalar E_R5 { 4480 * R_SI };

	scalar A_R6 { 1.4E14 / 1E6 };
	scalar n_R6 { 0 };
	scalar E_R6 { 8250 * R_SI };

	scalar A_R7 { 6.02E11 / 1E6 };
	scalar n_R7 { 0 };
	scalar E_R7 { 9440 * R_SI };

	scalar A_R8 { 7.25E15 / 1E12 };
	scalar n_R8 { 0 };
	scalar E_R8 { 500 * R_SI };

	scalar A_R9 { 1.14E9 / 1E6 };
	scalar n_R9 { 1.3 };
	scalar E_R9 { 1825 * R_SI };

	scalar A_R10 { 6.53E5 };
	scalar n_R10 { -1 };
	scalar E_R10 { 0 };

	scalar A_R11 { 7.8E7 };
	scalar n_R11 { -1.04 };
	scalar E_R11 { -1.75 * R_SI };

	scalar A_R12 { 2.54E6 };
	scalar n_R12 { -1.5 };
	scalar E_R12 { 0 };

	scalar ΔН_298 = -242000;
	scalar ΔU_298 = ΔН_298;
	constexpr static scalar Δn = -1;

	iterativeSolver itSolv;

	class cellReactionMatrix
	{
		iterativeSolver solverFlag;

		struct reactionMatrix
		{
			std::array<scalar, 7> Diagonale { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

			std::array<triangleList, 7> LeftTriangle {

			triangleList(0), // A1

					triangleList { std::make_pair(0.0, 0) }, // A2

					triangleList { std::make_pair(0.0, 0), std::make_pair(0.0,
							1) }, // A3

					triangleList { std::make_pair(0.0, 0), std::make_pair(0.0,
							1), std::make_pair(0.0, 2) }, // A4

					triangleList { std::make_pair(0.0, 0), std::make_pair(0.0,
							1), std::make_pair(0.0, 2), std::make_pair(0.0, 3) }, // A5

					triangleList { std::make_pair(0.0, 0), std::make_pair(0.0,
							2), std::make_pair(0.0, 3) }, // A6

					triangleList { std::make_pair(0.0, 2), std::make_pair(0.0,
							3), std::make_pair(0.0, 4), std::make_pair(0.0, 5) } // A7
			};

			std::array<triangleList, 7> RightTriangle {

			triangleList { std::make_pair(0.0, 1), std::make_pair(0.0, 2),
					std::make_pair(0.0, 3) }, // A1

					triangleList { std::make_pair(0.0, 2), std::make_pair(0.0,
							3) }, // A2

					triangleList { std::make_pair(0.0, 3), std::make_pair(0.0,
							4), std::make_pair(0.0, 5) }, // A3

					triangleList { std::make_pair(0.0, 4) }, // A4

					triangleList(0), // A5

					triangleList(0), // A6

					triangleList(0) // A7
			};

			std::array<scalar, 7> FreeTerm { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

			void transpose() noexcept;
		} matrix;

		std::valarray<scalar> matrixDotProduct(const reactionMatrix & m,
				const std::valarray<scalar> & v) const noexcept;

		auto solveJ(const std::array<scalar, 7> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 7>;

		auto solveGS(const std::array<scalar, 7> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 7>;

		auto solveCG(const std::array<scalar, 7> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 7>;

		auto solveJCG(const std::array<scalar, 7> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 7>;

		auto solveGE() const -> std::array<scalar, 7>;
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const scalar timeStep, const scalar k_R1,
				const scalar k_R2, const scalar k_R3, const scalar k_R4,
				const scalar k_R5, const scalar k_R6, const scalar k_R7,
				const scalar k_R8, const scalar k_R9, const scalar k_R10,
				const scalar k_R11, const scalar k_R12, const scalar C_O2_0,
				const scalar C_O_0, const scalar C_H2_0, const scalar C_H_0,
				const scalar C_OH_0, const scalar C_HO2_0, const scalar C_H2O_0,
				const scalar M_0, const scalar rho_0,
				const std::array<scalar, 7> & molMass,
				const iterativeSolver solverType);

		auto solve(const std::array<scalar, 7> & oldField,
				const std::size_t maxIterationNumber) const ->
						std::array<scalar, 7>;
	};

	cellReactionMatrix velocityCalculation(const scalar timestep,
			const scalar T, const std::array<scalar, 8> & concentrations,
			const std::array<scalar, 7> & molarMasses, const scalar rho,
			const scalar R) const noexcept;

	void timeStepIntegration(homogeneousPhase<cubicCell> & phaseN) const;
public:
	chemicalKineticsH2O2Combustion(const homogeneousPhase<cubicCell> & phaseIn,
			const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) const
			override;
};
}  // namespace schemi

#endif /* CHEMICALKINETICSH2O2COMBUSTION_HPP_ */
