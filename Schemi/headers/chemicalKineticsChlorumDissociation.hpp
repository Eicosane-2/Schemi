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

	class cellReactionMatrix
	{
		std::array<scalar, 2> matrixDiagonale { 0.0, 0.0 };

		std::array<triangleList, 2> matrixLeftTriangle {

		triangleList(0),

		triangleList(1, std::make_pair(0.0, 0))

		};

		std::array<triangleList, 2> matrixRightTriangle {

		triangleList(1, std::make_pair(0.0, 1)),

		triangleList(0)

		};

		std::array<scalar, 2> matrixFreeTerm { 0.0, 0.0 };

		void normalize(std::valarray<scalar> & res) const noexcept;
	public:
		cellReactionMatrix() noexcept;

		cellReactionMatrix(const scalar timeStep, const scalar k_diss,
				const scalar k_recomb, const scalar C_Cl2_0,
				const scalar C_Cl_0, const scalar M_0, const scalar rho_0,
				const std::valarray<scalar> & molMass) noexcept;

		std::array<scalar, 2> solve(const std::array<scalar, 2> & oldField,
				const std::size_t maxIterationNumber) const noexcept;
	};

	std::vector<cellReactionMatrix> velocityCalculation(const scalar timestep,
			const homogeneousPhase<cubicCell> & phase) const noexcept;

	void timeStepIntegration(
			homogeneousPhase<cubicCell> & phaseN) const noexcept;
public:
	chemicalKineticsChlorumDissociation(
			const homogeneousPhase<cubicCell> & phaseIn,
			const std::size_t itNumber);

	void solveChemicalKinetics(
			homogeneousPhase<cubicCell> & phaseIn) const noexcept override;
};
}  // namespace schemi

#endif /* CHEMICALKINETICSCHLORUMDISSOCIATION_ */
