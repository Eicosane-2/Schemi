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
namespace chemicalKinetics
{
class NO2Disproportionation: public abstractChemicalKinetics
{
	static constexpr std::size_t N { 4 };

	// (/ 1E6) --- cm^3/mole/s converted to m^3/mole/s
	scalar A_forward { 5.36E-50 * NAvogardro * NAvogardro / 1E12 };
	scalar n_forward { 3.95 };
	scalar E_forward { -1825 * R_SI };

	scalar A_backward { 3.31E-19 / 1E6 * NAvogardro };
	scalar n_backward { 2.478 };
	scalar E_backward { 3199 * R_SI };

	scalar ΔН_298 = (-76.73 - 134.31 + 241.826 - 2 * 33.1) * 1000;
	scalar ΔU_298 = ΔН_298;
	constexpr static scalar Δn = -1;

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
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const scalar timeStep, const scalar k_f,
				const scalar k_b, const scalar C_NO2_0, const scalar C_H2O_0,
				const scalar C_HNO2_0, const scalar C_HNO3_0,
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
	NO2Disproportionation(const homogeneousPhase<cubicCell> & phaseIn,
			const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) const
			override;
};
}
}  // namespace schemi

#endif /* CHEMICALKINETICSNO2DISPROPORTIONATION_HPP_ */
