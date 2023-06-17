/*
 * field.hpp
 *
 *  Created on: 2019/11/13
 *      Author: Maxim Boldyrev
 *
 *      Structure storing physical fields.
 */
#ifndef FIELD_HPP_
#define FIELD_HPP_

#include <valarray>

#include "exception.hpp"
#include "globalConstants.hpp"
#include "mesh.hpp"

namespace schemi
{
template<typename typeOfValue, typename typeOfEntity>
struct field
{
	field(const field<typeOfValue, typeOfEntity>&) = default;

	auto& operator=(const field<typeOfValue, typeOfEntity> & f)
	{
		if ((f.boundCond().size() != boundaryConditionInfo.size())
				&& (boundaryConditionInfo.size() != 0))
			throw exception("Fields have different boundary condition size.",
					errorsEnum::fieldInitializationError);

		valueField = f.ref();
		boundaryConditionInfo = f.boundCond();

		fieldSize = f.size();

		return *this;
	}

	explicit field(const mesh & meshIn,
			const typeOfValue & value = typeOfValue(0)) :
			valueField(0), boundaryConditionInfo(0), meshReference { meshIn }, fieldSize {
					0 }
	{
		is_mesh_initialised();

		if constexpr (std::is_same_v<typeOfEntity, cubicCell>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.cellsSize());

			boundaryConditionInfo.resize(meshReference.surfacesSize());
			for (std::size_t i = 0; i < boundaryConditionInfo.size(); ++i)
			{
				boundaryConditionInfo[i].first = meshReference.bndType()[i];
				boundaryConditionInfo[i].second = value;
			}

			fieldSize = meshReference.cellsSize();
		}
		else if constexpr (std::is_same_v<typeOfEntity, quadraticSurface>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.surfacesSize());

			fieldSize = meshReference.surfacesSize();
		}
		else
			throw exception("Unknown type of field.",
					errorsEnum::fieldInitializationError);
	}

	field(const mesh & meshIn, const typeOfValue & value,
			const boundaryConditionType & b1, const typeOfValue & bv1,
			const boundaryConditionType & b2, const typeOfValue & bv2,
			const boundaryConditionType & b3, const typeOfValue & bv3,
			const boundaryConditionType & b4, const typeOfValue & bv4,
			const boundaryConditionType & b5, const typeOfValue & bv5,
			const boundaryConditionType & b6, const typeOfValue & bv6) :
			valueField(0), boundaryConditionInfo(0), meshReference { meshIn }, fieldSize {
					0 }
	{
		is_mesh_initialised();

		if constexpr (std::is_same_v<typeOfEntity, cubicCell>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.cellsSize());

			const std::size_t surfSize = meshReference.surfacesSize();
			boundaryConditionInfo.resize(surfSize,
					std::pair<boundaryConditionType, typeOfValue>(
							boundaryConditionType::innerSurface,
							typeOfValue(0)));

			std::size_t prev { 0 };

			for (std::size_t j = prev; j < (prev + meshReference.tailNumber());
					++j)
			{
				boundaryConditionInfo[j].first = b1;
				boundaryConditionInfo[j].second = bv1;
			}
			prev += meshReference.tailNumber() + meshReference.innerNumber();

			for (std::size_t j = prev; j < (prev + meshReference.pointNumber());
					++j)
			{
				boundaryConditionInfo[j].first = b2;
				boundaryConditionInfo[j].second = bv2;
			}
			prev += meshReference.pointNumber();

			for (std::size_t j = prev;
					j < (prev + meshReference.bottomNumber()); ++j)
			{
				boundaryConditionInfo[j].first = b3;
				boundaryConditionInfo[j].second = bv3;
			}
			prev += meshReference.bottomNumber();

			for (std::size_t j = prev; j < (prev + meshReference.rightNumber());
					++j)
			{
				boundaryConditionInfo[j].first = b4;
				boundaryConditionInfo[j].second = bv4;
			}
			prev += meshReference.rightNumber();

			for (std::size_t j = prev; j < (prev + meshReference.leftNumber());
					++j)
			{
				boundaryConditionInfo[j].first = b5;
				boundaryConditionInfo[j].second = bv5;
			}
			prev += meshReference.leftNumber();

			for (std::size_t j = prev; j < (prev + meshReference.topNumber());
					++j)
			{
				boundaryConditionInfo[j].first = b6;
				boundaryConditionInfo[j].second = bv6;
			}

			fieldSize = meshReference.cellsSize();
		}
		else
			throw exception("Field is not volumeField.",
					errorsEnum::fieldInitializationError);
	}

	field(const mesh & meshIn, const typeOfValue & value,
			const boundaryConditionType & b1,
			const std::vector<typeOfValue> & bv1,
			const boundaryConditionType & b2,
			const std::vector<typeOfValue> & bv2,
			const boundaryConditionType & b3,
			const std::vector<typeOfValue> & bv3,
			const boundaryConditionType & b4,
			const std::vector<typeOfValue> & bv4,
			const boundaryConditionType & b5,
			const std::vector<typeOfValue> & bv5,
			const boundaryConditionType & b6,
			const std::vector<typeOfValue> & bv6, const std::size_t mpiSize) :
			valueField(0), boundaryConditionInfo(0), meshReference { meshIn }, fieldSize {
					0 }
	{
		is_mesh_initialised();

		if constexpr (std::is_same_v<typeOfEntity, cubicCell>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.cellsSize());

			const std::size_t surfSize = meshReference.surfacesSize();
			boundaryConditionInfo.resize(surfSize,
					std::pair<boundaryConditionType, typeOfValue>(
							boundaryConditionType::innerSurface,
							typeOfValue(0)));

			const auto right0Index = meshReference.tailNumber()
					+ meshReference.innerNumber() + meshReference.pointNumber()
					+ meshReference.bottomNumber();
			const auto left0Index = right0Index + meshReference.rightNumber();

			const std::size_t tailOIndex = 0;
			const auto point0Index = meshReference.tailNumber()
					+ meshReference.innerNumber();

			scalar YLength { 1. / stabilizator };
			scalar XLengthLocal { 1. / stabilizator };

			if (meshReference.taskDimension() != dimensionsEnum::task1D)
			{
				YLength = (meshReference.surfaces()[left0Index].rC()
						- meshReference.surfaces()[right0Index].rC()).mag();
				XLengthLocal = (meshReference.surfaces()[point0Index].rC()
						- meshReference.surfaces()[tailOIndex].rC()).mag();
			}

			const auto XLength = XLengthLocal * mpiSize;

			const auto YLengthPartitionTail = YLength / bv1.size();
			const auto YLengthPartitionPoint = YLength / bv2.size();

			const auto XLengthPartitionBottom = XLength / bv3.size();
			const auto XLengthPartitionRight = XLength / bv4.size();
			const auto XLengthPartitionLeft = XLength / bv5.size();
			const auto XLengthPartitionTop = XLength / bv6.size();

			std::size_t prev { 0 };

			for (std::size_t j = prev; j < (prev + meshReference.tailNumber());
					++j)
			{
				const std::size_t boundaryNumber =
						static_cast<std::size_t>(meshReference.surfaces()[j].rC().v()[1]
								/ YLengthPartitionTail);

				boundaryConditionInfo[j].first = b1;
				boundaryConditionInfo[j].second = bv1[boundaryNumber];
			}
			prev += meshReference.tailNumber() + meshReference.innerNumber();

			for (std::size_t j = prev; j < (prev + meshReference.pointNumber());
					++j)
			{
				const std::size_t boundaryNumber =
						static_cast<std::size_t>(meshReference.surfaces()[j].rC().v()[1]
								/ YLengthPartitionPoint);

				boundaryConditionInfo[j].first = b2;
				boundaryConditionInfo[j].second = bv2[boundaryNumber];
			}
			prev += meshReference.pointNumber();

			for (std::size_t j = prev;
					j < (prev + meshReference.bottomNumber()); ++j)
			{
				const std::size_t boundaryNumber =
						static_cast<std::size_t>(meshReference.surfaces()[j].rC().v()[0]
								/ XLengthPartitionBottom);

				boundaryConditionInfo[j].first = b3;
				boundaryConditionInfo[j].second = bv3[boundaryNumber];
			}
			prev += meshReference.bottomNumber();

			for (std::size_t j = prev; j < (prev + meshReference.rightNumber());
					++j)
			{
				const std::size_t boundaryNumber =
						static_cast<std::size_t>(meshReference.surfaces()[j].rC().v()[0]
								/ XLengthPartitionRight);

				boundaryConditionInfo[j].first = b4;
				boundaryConditionInfo[j].second = bv4[boundaryNumber];
			}
			prev += meshReference.rightNumber();

			for (std::size_t j = prev; j < (prev + meshReference.leftNumber());
					++j)
			{
				const std::size_t boundaryNumber =
						static_cast<std::size_t>(meshReference.surfaces()[j].rC().v()[0]
								/ XLengthPartitionLeft);

				boundaryConditionInfo[j].first = b5;
				boundaryConditionInfo[j].second = bv5[boundaryNumber];
			}
			prev += meshReference.leftNumber();

			for (std::size_t j = prev; j < (prev + meshReference.topNumber());
					++j)
			{
				const std::size_t boundaryNumber =
						static_cast<std::size_t>(meshReference.surfaces()[j].rC().v()[0]
								/ XLengthPartitionTop);

				boundaryConditionInfo[j].first = b6;
				boundaryConditionInfo[j].second = bv6[boundaryNumber];
			}

			fieldSize = meshReference.cellsSize();
		}
		else
			throw exception("Field is not volumeField.",
					errorsEnum::fieldInitializationError);
	}

	field(const mesh & meshIn, const typeOfValue & value,
			const std::vector<std::pair<boundaryConditionType, typeOfValue>> & boundCondtIn) :
			valueField(0), boundaryConditionInfo(0), meshReference { meshIn }, fieldSize {
					0 }
	{
		is_mesh_initialised();

		if constexpr (std::is_same_v<typeOfEntity, cubicCell>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.cellsSize());

			fieldSize = meshReference.cellsSize();
		}
		else if constexpr (std::is_same_v<typeOfEntity, quadraticSurface>)
		{
			valueField = std::valarray<typeOfValue>(value,
					meshReference.surfacesSize());

			fieldSize = meshReference.surfacesSize();
		}
		else
			throw exception("Unknown type of field.",
					errorsEnum::fieldInitializationError);

		boundaryConditionInfo = boundCondtIn;
	}

	const std::valarray<typeOfValue>& ref() const noexcept
	{
		return valueField;
	}

	std::valarray<typeOfValue>& ref_r() noexcept
	{
		return valueField;
	}

	const std::vector<std::pair<boundaryConditionType, typeOfValue>>& boundCond() const noexcept
	{
		return boundaryConditionInfo;
	}

	std::vector<std::pair<boundaryConditionType, typeOfValue>>& boundCond_r() noexcept
	{
		return boundaryConditionInfo;
	}

	const mesh& meshRef() const noexcept
	{
		return meshReference;
	}

	std::size_t size() const noexcept
	{
		return fieldSize;
	}
private:
	std::valarray<typeOfValue> valueField;
	std::vector<std::pair<boundaryConditionType, typeOfValue>> boundaryConditionInfo;
	const mesh & meshReference;
	std::size_t fieldSize;

	void is_mesh_initialised() const
	{
		if (!meshReference.is_initialized())
			throw exception("Mesh has not yet been initialized.",
					errorsEnum::fieldInitializationError);
	}
};
}  // namespace schemi

#endif /* FIELD_HPP_ */
