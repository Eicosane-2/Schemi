/*
 * SLEMatrix.hpp
 *
 *  Created on: 2020/01/16
 *      Author: Maxim Boldyrev
 *
 *      Class generating and storing matrix for system of linear equations
 */

#ifndef SLEMATRIX_HPP_
#define SLEMATRIX_HPP_

#include <string>
#include <valarray>
#include <vector>

#include "boundaryConditionValue.hpp"
#include "mesh.hpp"
#include "scalar.hpp"
#include "tensor.hpp"
#include "surfaceField.hpp"
#include "volumeField.hpp"

namespace schemi
{
class SLEMatrix
{
public:
	explicit SLEMatrix(const std::string & stringIn) noexcept;

	struct SLEMatrixStorage
	{
		std::valarray<scalar> centralDiagonale;
		std::valarray<scalar> freeTerm;
		std::vector<std::vector<std::pair<scalar, std::size_t>>> lowerTriangle,
				upperTriangle;

		SLEMatrixStorage();

		explicit SLEMatrixStorage(const mesh & meshRef) noexcept;

		void transpose() noexcept;
	};

	const std::string name;

	std::vector<SLEMatrixStorage> SLE { };

	void generateNabla(const volumeField<scalar> & vField,
			const volumeField<scalar> & rFieldOld,
			const volumeField<scalar> & rFieldNew,
			const surfaceField<vector> & additionalField, const scalar timestep,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	void generateLaplacian(const volumeField<scalar> & vField,
			const volumeField<scalar> & rFieldOld,
			const volumeField<scalar> & rFieldNew,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateLaplacian2TO(const volumeField<scalar> & vField,
			const volumeField<scalar> & rFieldOld,
			const volumeField<scalar> & rFieldNew,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateLaplacian(const volumeField<vector> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateLaplacian2TO(const volumeField<vector> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateLaplacian(const volumeField<tensor> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateLaplacian2TO(const volumeField<tensor> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateExplicitLaplacian(const volumeField<scalar> & vField,
			const volumeField<scalar> & rFieldOld,
			const volumeField<scalar> & rFieldNew,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const surfaceField<scalar> & surfaceVField, const scalar timestep,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	void generateExplicitLaplacian(const volumeField<vector> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const surfaceField<vector> & surfaceVField, const scalar timestep,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	void generateExplicitLaplacian(const volumeField<tensor> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const surfaceField<tensor> & surfaceVField, const scalar timestep,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	/*Such distribution is only possible for non-negative scalars, so no overloads for vector and tensor.*/
	void distributeMinusSourceTerm(const volumeField<scalar> & source,
			const volumeField<scalar> & basicField) noexcept;

	void distributePlusSourceTerm(const volumeField<scalar> & source,
			const volumeField<scalar> & basicField) noexcept;
};
}  // namespace schemi

#endif /* SLEMATRIX_HPP_ */
