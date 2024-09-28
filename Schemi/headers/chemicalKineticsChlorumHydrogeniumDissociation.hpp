/*
 * chemicalKineticsChlorumHydrogeniumDissociation.hpp
 *
 *  Created on: 2023/05/10
 *      Author: Maxim Boldyrev
 */

#ifndef CHEMICALKINETICSCHLORUMHYDROGENIUMDISSOCIATION_HPP_
#define CHEMICALKINETICSCHLORUMHYDROGENIUMDISSOCIATION_HPP_

#include "abstractChemicalKinetics.hpp"

namespace schemi
{
namespace chemicalKinetics
{
class ChlorumHydrogeniumDissociation: public abstractChemicalKinetics
{
	static constexpr std::size_t N { 4 };

	/*(/ 1E6) --- cm^3/mole/s converted to m^3/mole/s*/
	/*(/ 1E12) --- cm^6/mole^2/s converted to m^6/mole^2/s*/
	scalar A_Cl2_forw { 5.74E15 / 1E6 };
	scalar n_Cl2_forw { 0 };
	scalar E_Cl2_forw { 55.596 * 1000 * 4.184 };

	scalar A_Cl_backw { 2.3E19 / 1E12 };
	scalar n_Cl_backw { -1.5 };
	scalar E_Cl_backw { 0 };

	scalar A_H2_forw { 1.8E17 / 1E6 };
	scalar n_H2_forw { -0.5 };
	scalar E_H2_forw { 454.0 * 1000 };

	scalar A_H2_backw { 6.53E17 / 1E12 };
	scalar n_H2_backw { -1 };
	scalar E_H2_backw { 0 };

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

			triangleList(0),

			triangleList(1, std::make_pair(0.0, 2))

			};

			std::array<triangleList, N> RightTriangle {

			triangleList(1, std::make_pair(0.0, 1)),

			triangleList(0),

			triangleList(1, std::make_pair(0.0, 3)),

			triangleList(0)

			};

			std::array<scalar, N> FreeTerm { 0.0, 0.0, 0.0, 0.0 };

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

		cellReactionMatrix(const scalar timeStep, const scalar k_diss_Cl2,
				const scalar k_recomb_Cl2, const scalar k_diss_H2,
				const scalar k_recomb_H2, const scalar C_Cl2_0,
				const scalar C_Cl_0, const scalar C_H2_0, const scalar C_H_0,
				const scalar M_0, const scalar rho_0,
				const std::array<scalar, N> & molMass,
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
	ChlorumHydrogeniumDissociation(const homogeneousPhase<cubicCell> & phaseIn,
			const scalar mt);

	void solveChemicalKinetics(homogeneousPhase<cubicCell> & phaseIn) const
			override;
};
}
}  // namespace schemi

#endif /* CHEMICALKINETICSCHLORUMHYDROGENIUMDISSOCIATION_HPP_ */
