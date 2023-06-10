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
	// (/ 1E6) --- cm^3/mole/s converted to m^3/mole/s
	scalar A_Cl2_diss { 5.74 * 1E15 / 1E6 };
	scalar n_Cl2_diss { 0 };
	scalar E_Cl2_diss { 55.596 * 1000 * 4.184 };

	scalar A_Cl_recomb { 2.3 * 1E19 / 1E12 };
	scalar n_Cl_recomb { -1.5 };
	scalar E_Cl_recomb { 0 };

	scalar A_H2_diss { 1.8 * 1E17 / 1E6 };
	scalar n_H2_diss { -0.5 };
	scalar E_H2_diss { 454 * 1000 };

	scalar A_H_recomb { 6.53 * 1E17 / 1E12 };
	scalar n_H_recomb { -1 };
	scalar E_H_recomb { 0 };

	scalar A_Cl_H2_prop { 2.88 * 1E8 / 1E6 };
	scalar n_Cl_H2_prop { 1.58 };
	scalar E_Cl_H2_prop { 3.199 * 1000 * 4.184 };

	scalar A_H_Cl2_prop { 9.64 * 1E13 / 1E6 };
	scalar n_H_Cl2_prop { 0 };
	scalar E_H_Cl2_prop { 1.68 * 1000 * 4.184 };

	scalar A_H_Cl_recomb { 2 * 1E23 / 1E12 };
	scalar n_H_Cl_recomb { -2.45 };
	scalar E_H_Cl_recomb { 0 };

	scalar A_HCl_diss { 4.4 * 1E13 / 1E6 };
	scalar n_HCl_diss { 0 };
	scalar E_HCl_diss { 81.753 * 1000 * 4.184 };

	scalar ΔН_298 = -92300;
	scalar ΔU_298 = ΔН_298;
	constexpr static scalar Δn = 0;

	class cellReactionMatrix
	{
		std::array<scalar, 5> matrixDiagonale { 0.0, 0.0, 0.0, 0.0, 0.0 };

		std::array<triangleList, 5> matrixLeftTriangle {

		triangleList(0),

		triangleList(1, std::make_pair(0.0, 0)),

		triangleList(1, std::make_pair(0.0, 1)),

		triangleList { std::make_pair(0.0, 0), std::make_pair(0.0, 1),
				std::make_pair(0.0, 2) },

		triangleList { std::make_pair(0.0, 0), std::make_pair(0.0, 1),
				std::make_pair(0.0, 2), std::make_pair(0.0, 3) }

		};

		std::array<triangleList, 5> matrixRightTriangle {

		triangleList { std::make_pair(0.0, 1), std::make_pair(0.0, 3) },

		triangleList { std::make_pair(0.0, 2), std::make_pair(0.0, 3),
				std::make_pair(0.0, 4) },

		triangleList(1, std::make_pair(0.0, 3)),

		triangleList(1, std::make_pair(0.0, 4)),

		triangleList(0)

		};

		std::array<scalar, 5> matrixFreeTerm { 0.0, 0.0, 0.0, 0.0, 0.0 };

		void normalize(std::valarray<scalar> & res) const noexcept;
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const scalar timeStep, const scalar k_diss_Cl2,
				const scalar k_recomb_Cl, const scalar k_diss_H2,
				const scalar k_recomb_H, const scalar k_prop_Cl_H2,
				const scalar k_prop_H_Cl2, const scalar k_recomb_H_Cl,
				const scalar k_diss_HCl, const scalar C_Cl2_0,
				const scalar C_Cl_0, const scalar C_H2_0, const scalar C_H_0,
				const scalar C_HCl_0, const scalar M_0, const scalar rho_0,
				const std::valarray<scalar> & molMass) noexcept;

		std::array<scalar, 5> solve(const std::array<scalar, 5> & oldField,
				const std::size_t maxIterationNumber) const noexcept;
	};

	std::vector<cellReactionMatrix> velocityCalculation(const scalar timestep,
			const homogeneousPhase<cubicCell> & phase) const noexcept;

	void timeStepIntegration(
			homogeneousPhase<cubicCell> & phaseN) const noexcept;
public:
	chemicalKineticsH2Cl2Combustion(const homogeneousPhase<cubicCell> & phaseIn,
			const std::size_t itNumber);

	void solveChemicalKinetics(
			homogeneousPhase<cubicCell> & phaseIn) const noexcept override;
};
}  // namespace schemi

#endif /* CHEMICALKINETICSH2CL2COMBUSTION_HPP_ */
