/*
 * chemicalKineticsNO2Disproportionation.hpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSNO2DISPROPORTIONATION_HPP_
#define CHEMICALKINETICSNO2DISPROPORTIONATION_HPP_

#include "abstractChemicalKinetics.hpp"
#include "globalConstants.hpp"

namespace schemi
{
class chemicalKineticsNO2Disproportionation: public abstractChemicalKinetics
{
	// (/ 1E6) --- cm^3/mole/s converted to m^3/mole/s
	scalar A_forward { 5.36 * 1E-50 * NAvogardro * NAvogardro / 1E12 };
	scalar n_forward { 3.95 };
	scalar E_forward { -1825 * R_SI };

	scalar A_backward { 3.31 * 1E-19 / 1E6 * NAvogardro };
	scalar n_backward { 2.478 };
	scalar E_backward { 3199 * R_SI };

	scalar ΔН_298 = (-76.73 - 134.31 + 241.826 - 2 * 33.1) * 1000;
	scalar ΔU_298 = ΔН_298;
	constexpr static scalar Δn = -1;

	class cellReactionMatrix
	{
		std::array<scalar, 4> matrixDiagonale { 0.0, 0.0, 0.0, 0.0 };

		std::array<triangleList, 4> matrixLeftTriangle {

		triangleList(0),

		triangleList(1, std::make_pair(0.0, 0)),

		triangleList { std::make_pair(0.0, 0), std::make_pair(0.0, 1) },

		triangleList { std::make_pair(0.0, 0), std::make_pair(0.0, 1),
				std::make_pair(0.0, 2) }

		};

		std::array<triangleList, 4> matrixRightTriangle {

		triangleList { std::make_pair(0.0, 1), std::make_pair(0.0, 2),
				std::make_pair(0.0, 3) },

		triangleList { std::make_pair(0.0, 2), std::make_pair(0.0, 3) },

		triangleList { std::make_pair(0.0, 3) },

		triangleList(0)

		};

		std::array<scalar, 4> matrixFreeTerm { 0.0, 0.0, 0.0, 0.0 };

		void normalize(std::valarray<scalar> & res) const noexcept;
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const scalar timeStep, const scalar k_f,
				const scalar k_b, const scalar C_NO2_0, const scalar C_H2O_0,
				const scalar C_HNO2_0, const scalar C_HNO3_0,
				const scalar rho_0,
				const std::valarray<scalar> & molMass) noexcept;

		std::array<scalar, 4> solve(const std::array<scalar, 4> & oldField,
				const std::size_t maxIterationNumber) const noexcept;
	};

	std::vector<cellReactionMatrix> velocityCalculation(const scalar timestep,
			const homogeneousPhase<cubicCell> & phase) const noexcept;

	void timeStepIntegration(
			homogeneousPhase<cubicCell> & phaseN) const noexcept;
public:
	chemicalKineticsNO2Disproportionation(
			const homogeneousPhase<cubicCell> & phaseIn,
			const std::size_t itNumber);

	void solveChemicalKinetics(
			homogeneousPhase<cubicCell> & phaseIn) const noexcept override;
};
}  // namespace schemi

#endif /* CHEMICALKINETICSNO2DISPROPORTIONATION_HPP_ */
