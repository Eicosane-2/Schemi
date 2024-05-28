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
class boundaryConditionValue;

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

		std::valarray<scalar> explOldTime = std::valarray<scalar>(0);

		SLEMatrixStorage();

		explicit SLEMatrixStorage(const mesh & meshRef) noexcept;

		void transpose() noexcept;
	};

	const std::string name;

	std::vector<SLEMatrixStorage> SLE { };

	void generateLaplacianSurfaceBoundary(const volumeField<scalar> & vField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	void generateDTimeNabla(const volumeField<scalar> & vField,
			const volumeField<scalar> & rFieldOld,
			const volumeField<scalar> & rFieldNew,
			const surfaceField<vector> & additionalField, const scalar timestep,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	void generateDTimeLaplacian(const volumeField<scalar> & vField,
			const volumeField<scalar> & rFieldOld,
			const volumeField<scalar> & rFieldNew,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateDTimeLaplacian2TO(const volumeField<scalar> & vField,
			const volumeField<scalar> & rFieldOld,
			const volumeField<scalar> & rFieldNew,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateDTimeLaplacian(const volumeField<vector> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateDTimeLaplacian2TO(const volumeField<vector> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateDTimeLaplacian(const volumeField<tensor> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateDTimeLaplacian2TO(const volumeField<tensor> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const scalar timestep, const boundaryConditionValue & bncCalc,
			const std::size_t compt = componentPlaceholder);

	void generateDTimeExplicitLaplacian(const volumeField<scalar> & vField,
			const volumeField<scalar> & rFieldOld,
			const volumeField<scalar> & rFieldNew,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	void generateDTimeExplicitLaplacian(const volumeField<vector> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	void generateDTimeExplicitLaplacian(const volumeField<tensor> & vField,
			const volumeField<scalar> & rField,
			const surfaceField<scalar> & effectiveDiffusionCoefficient,
			const boundaryConditionValue & bncCalc, const std::size_t compt =
					componentPlaceholder);

	void freeSourceTerm(const volumeField<scalar> & source) noexcept;

	void freeSourceTerm(const volumeField<vector> & source) noexcept;

	void freeSourceTerm(const volumeField<tensor> & source) noexcept;

	void distributeSourceTerm(const volumeField<scalar> & source,
			const volumeField<scalar> & basicField) noexcept;

	void distributeSourceTerm(const volumeField<vector> & source,
			const volumeField<vector> & basicField) noexcept;

	void distributeSourceTerm(const volumeField<tensor> & source,
			const volumeField<tensor> & basicField) noexcept;

	void diagonaleSourceTerm(const volumeField<scalar> & source) noexcept;

	void diagonaleSourceTerm(const volumeField<vector> & source) noexcept;

	void diagonaleSourceTerm(const volumeField<tensor> & source) noexcept;
};
}  // namespace schemi

#endif /* SLEMATRIX_HPP_ */
