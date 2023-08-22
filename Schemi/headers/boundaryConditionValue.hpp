/*
 * boundaryConditionValue.hpp
 *
 *  Created on: 2019/12/20
 *      Author: Maxim Boldyrev
 *
 *      Class calculating boundary conditions.
 */

#ifndef BOUNDARYCONDITIONVALUE_HPP_
#define BOUNDARYCONDITIONVALUE_HPP_

#include "abstractMixtureThermodynamics.hpp"
#include "abstractTurbulentParameters.hpp"
#include "slipFunction.hpp"

namespace schemi
{
class boundaryConditionValue
{
	const abstractTurbulentParameters & turbPar;
	const bunchOfFields<cubicCell> & cellFields;
	const abstractMixtureThermodynamics & mix;
	const mesh & meshReference;
public:
	boundaryConditionValue(const abstractTurbulentParameters & turbPar_in,
			const bunchOfFields<cubicCell> & cellFields_in,
			const abstractMixtureThermodynamics & mix_in) noexcept;

	template<typename typeOfValue>
	typeOfValue boundaryConditionValueSurface(const typeOfValue & inValue,
			const std::pair<boundaryConditionType, typeOfValue> & bndCondInfo,
			const std::size_t cellIndex, const std::size_t surfaceIndex,
			const std::size_t componentIndex = componentPlaceholder) const
	{
		typeOfValue retValue { 0 };

		switch (bndCondInfo.first)
		{
		case boundaryConditionType::freeBoundary:
			retValue = inValue;
			break;
		case boundaryConditionType::slip:
			retValue = slipFunction(inValue,
					meshReference.surfaces()[surfaceIndex].N());
			break;
		case boundaryConditionType::calculatedParallelBoundary:
		case boundaryConditionType::fixedValue:
			retValue = (inValue + bndCondInfo.second) * 0.5;
			break;
		case boundaryConditionType::innerSurface:
			throw exception(
					"<<innerSurface>> is incorrect variant for returning boundary condition value.",
					errors::boundaryConditionError);
			break;
		case boundaryConditionType::calculated:
			throw exception(
					"<<calculated>> is incorrect variant for returning boundary condition value.",
					errors::boundaryConditionError);
			break;
		case boundaryConditionType::calculatedTurbulentViscosity:
		{
			const scalar kBoundary(
					boundaryConditionValueSurface(
							cellFields.kTurb.ref()[cellIndex],
							cellFields.kTurb.boundCond()[surfaceIndex],
							cellIndex, surfaceIndex));

			const scalar epsilonBoundary { boundaryConditionValueSurface(
					cellFields.epsTurb.ref()[cellIndex],
					cellFields.epsTurb.boundCond()[surfaceIndex], cellIndex,
					surfaceIndex) };

			retValue = turbPar.calculateNut(kBoundary, epsilonBoundary);
		}
			break;
		case boundaryConditionType::calculatedTemperature:
		{
			std::valarray<scalar> concentrationsBoundary(0.,
					cellFields.concentration.v.size());

			for (std::size_t k = 1; k < concentrationsBoundary.size(); ++k)
			{
				concentrationsBoundary[k] = boundaryConditionValueSurface(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

				concentrationsBoundary[0] += concentrationsBoundary[k];
			}

			const scalar pressureBoundary { boundaryConditionValueSurface(
					cellFields.pressure.ref()[cellIndex],
					cellFields.pressure.boundCond()[surfaceIndex], cellIndex,
					surfaceIndex) };

			const scalar internalEnergyBoundary { mix.UvFromp(
					concentrationsBoundary, pressureBoundary) };

			retValue = mix.TFromUv(concentrationsBoundary,
					internalEnergyBoundary);
		}
			break;
		case boundaryConditionType::calculatedMassFraction:
		{
			const scalar density {
					boundaryConditionValueSurface(
							cellFields.concentration.v[componentIndex].ref()[cellIndex],
							cellFields.concentration.v[componentIndex].boundCond()[surfaceIndex],
							cellIndex, surfaceIndex)
							* mix.Mv()[componentIndex - 1] };

			scalar sumDensity { 0 };

			for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
				sumDensity += boundaryConditionValueSurface(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex) * mix.Mv()[k - 1];

			retValue = density / sumDensity;
		}
			break;
		case boundaryConditionType::calculatedMolarFraction:
		{
			const scalar concentration {
					boundaryConditionValueSurface(
							cellFields.concentration.v[componentIndex].ref()[cellIndex],
							cellFields.concentration.v[componentIndex].boundCond()[surfaceIndex],
							cellIndex, surfaceIndex) };

			scalar sumConcentration { 0 };

			for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
				sumConcentration += boundaryConditionValueSurface(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

			retValue = concentration / sumConcentration;
		}
			break;
		case boundaryConditionType::calculatedAverageMolarMass:
		{
			scalar sumDensity { 0 }, sumConcentration { 0 };

			for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
			{
				const scalar boundaryConcentration_k {
						boundaryConditionValueSurface(
								cellFields.concentration.v[k].ref()[cellIndex],
								cellFields.concentration.v[k].boundCond()[surfaceIndex],
								cellIndex, surfaceIndex) };

				sumConcentration += boundaryConcentration_k;

				sumDensity += boundaryConcentration_k * mix.Mv()[k - 1];
			}

			retValue = sumDensity / sumConcentration;
		}
			break;
		case boundaryConditionType::calculatedNonidealityCorrectionPerDensity:
		{
			std::valarray<scalar> concentrationsBoundary(0.,
					cellFields.concentration.v.size());

			scalar density { 0 };

			for (std::size_t k = 1; k < concentrationsBoundary.size(); ++k)
			{
				concentrationsBoundary[k] = boundaryConditionValueSurface(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

				concentrationsBoundary[0] += concentrationsBoundary[k];

				density += concentrationsBoundary[k] * mix.Mv()[k - 1];
			}

			const scalar temperature { boundaryConditionValueSurface(
					cellFields.temperature.ref()[cellIndex],
					cellFields.temperature.boundCond()[surfaceIndex], cellIndex,
					surfaceIndex) };

			const scalar nonIdealityCorrection { mix.nonIdeality(
					concentrationsBoundary, temperature) };

			retValue = nonIdealityCorrection / density;
		}
			break;
		case boundaryConditionType::calculatedCv:
		{
			std::valarray<scalar> concentrationsBoundary(0.,
					cellFields.concentration.v.size());

			for (std::size_t k = 1; k < concentrationsBoundary.size(); ++k)
			{
				concentrationsBoundary[k] = boundaryConditionValueSurface(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

				concentrationsBoundary[0] += concentrationsBoundary[k];
			}

			retValue = mix.Cv(concentrationsBoundary);
		}
			break;
		case boundaryConditionType::calculatedCvM:
		{
			std::valarray<scalar> concentrationsBoundary(0.,
					cellFields.concentration.v.size());

			scalar density { 0 };

			for (std::size_t k = 1; k < concentrationsBoundary.size(); ++k)
			{
				concentrationsBoundary[k] = boundaryConditionValueSurface(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

				concentrationsBoundary[0] += concentrationsBoundary[k];

				density += concentrationsBoundary[k] * mix.Mv()[k - 1];
			}

			retValue = mix.Cv(concentrationsBoundary)
					* concentrationsBoundary[0] / density;
		}
			break;
		case boundaryConditionType::calculatedDensity:
		{
			scalar density { 0 };

			for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
			{
				const scalar concentrationsBoundary_k {
						boundaryConditionValueSurface(
								cellFields.concentration.v[k].ref()[cellIndex],
								cellFields.concentration.v[k].boundCond()[surfaceIndex],
								cellIndex, surfaceIndex) };

				density += concentrationsBoundary_k * mix.Mv()[k - 1];
			}

			retValue = density;
		}
			break;
		default:
			throw exception("Unknown type of boundary condition.",
					errors::boundaryConditionError);
			break;
		}

		return retValue;
	}

	template<typename typeOfValue>
	typeOfValue boundaryConditionValueCell(const typeOfValue & inValue,
			const std::pair<boundaryConditionType, typeOfValue> & bndCondInfo,
			const std::size_t cellIndex, const std::size_t surfaceIndex,
			const std::size_t componentIndex = componentPlaceholder) const
	{
		typeOfValue retValue { 0 };

		switch (bndCondInfo.first)
		{
		case boundaryConditionType::freeBoundary:
			retValue = inValue;
			break;
		case boundaryConditionType::slip:
			retValue = slipFunction(inValue,
					meshReference.surfaces()[surfaceIndex].N()) * 2 - inValue;
			break;
		case boundaryConditionType::calculatedParallelBoundary:
		case boundaryConditionType::fixedValue:
			retValue = bndCondInfo.second;
			break;
		case boundaryConditionType::innerSurface:
			throw exception(
					"<<innerSurface>> is incorrect variant for returning boundary condition value.",
					errors::boundaryConditionError);
			break;
		case boundaryConditionType::calculated:
			throw exception(
					"<<calculated>> is incorrect variant for returning boundary condition value.",
					errors::boundaryConditionError);
			break;
		case boundaryConditionType::calculatedTurbulentViscosity:
		{
			const scalar kBoundary(
					boundaryConditionValueCell(
							cellFields.kTurb.ref()[cellIndex],
							cellFields.kTurb.boundCond()[surfaceIndex],
							cellIndex, surfaceIndex));

			const scalar epsilonBoundary { boundaryConditionValueCell(
					cellFields.epsTurb.ref()[cellIndex],
					cellFields.epsTurb.boundCond()[surfaceIndex], cellIndex,
					surfaceIndex) };

			retValue = turbPar.calculateNut(kBoundary, epsilonBoundary);
		}
			break;
		case boundaryConditionType::calculatedTemperature:
		{
			std::valarray<scalar> concentrationsBoundary(0.,
					cellFields.concentration.v.size());

			for (std::size_t k = 1; k < concentrationsBoundary.size(); ++k)
			{
				concentrationsBoundary[k] = boundaryConditionValueCell(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

				concentrationsBoundary[0] += concentrationsBoundary[k];
			}

			const scalar pressureBoundary { boundaryConditionValueCell(
					cellFields.pressure.ref()[cellIndex],
					cellFields.pressure.boundCond()[surfaceIndex], cellIndex,
					surfaceIndex) };

			const scalar internalEnergyBoundary { mix.UvFromp(
					concentrationsBoundary, pressureBoundary) };

			retValue = mix.TFromUv(concentrationsBoundary,
					internalEnergyBoundary);
		}
			break;
		case boundaryConditionType::calculatedMassFraction:
		{
			const scalar density {
					boundaryConditionValueCell(
							cellFields.concentration.v[componentIndex].ref()[cellIndex],
							cellFields.concentration.v[componentIndex].boundCond()[surfaceIndex],
							cellIndex, surfaceIndex)
							* mix.Mv()[componentIndex - 1] };

			scalar sumDensity { 0 };

			for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
				sumDensity += boundaryConditionValueCell(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex) * mix.Mv()[k - 1];

			retValue = density / sumDensity;
		}
			break;
		case boundaryConditionType::calculatedMolarFraction:
		{
			const scalar concentration {
					boundaryConditionValueCell(
							cellFields.concentration.v[componentIndex].ref()[cellIndex],
							cellFields.concentration.v[componentIndex].boundCond()[surfaceIndex],
							cellIndex, surfaceIndex) };

			scalar sumConcentration { 0 };

			for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
				sumConcentration += boundaryConditionValueCell(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

			retValue = concentration / sumConcentration;
		}
			break;
		case boundaryConditionType::calculatedAverageMolarMass:
		{
			scalar sumDensity { 0 }, sumConcentration { 0 };

			for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
			{
				const scalar boundaryConcentration_k {
						boundaryConditionValueCell(
								cellFields.concentration.v[k].ref()[cellIndex],
								cellFields.concentration.v[k].boundCond()[surfaceIndex],
								cellIndex, surfaceIndex) };

				sumConcentration += boundaryConcentration_k;

				sumDensity += boundaryConcentration_k * mix.Mv()[k - 1];
			}

			retValue = sumDensity / sumConcentration;
		}
			break;
		case boundaryConditionType::calculatedNonidealityCorrectionPerDensity:
		{
			std::valarray<scalar> concentrationsBoundary(0.,
					cellFields.concentration.v.size());

			scalar density { 0 };

			for (std::size_t k = 1; k < concentrationsBoundary.size(); ++k)
			{
				concentrationsBoundary[k] = boundaryConditionValueCell(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

				concentrationsBoundary[0] += concentrationsBoundary[k];

				density += concentrationsBoundary[k] * mix.Mv()[k - 1];
			}

			const scalar temperature { boundaryConditionValueSurface(
					cellFields.temperature.ref()[cellIndex],
					cellFields.temperature.boundCond()[surfaceIndex], cellIndex,
					surfaceIndex) };

			const scalar nonIdealityCorrection { mix.nonIdeality(
					concentrationsBoundary, temperature) };

			retValue = nonIdealityCorrection / density;
		}
			break;
		case boundaryConditionType::calculatedCv:
		{
			std::valarray<scalar> concentrationsBoundary(0.,
					cellFields.concentration.v.size());

			for (std::size_t k = 1; k < concentrationsBoundary.size(); ++k)
			{
				concentrationsBoundary[k] = boundaryConditionValueCell(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

				concentrationsBoundary[0] += concentrationsBoundary[k];
			}

			retValue = mix.Cv(concentrationsBoundary);
		}
			break;
		case boundaryConditionType::calculatedCvM:
		{
			std::valarray<scalar> concentrationsBoundary(0.,
					cellFields.concentration.v.size());

			scalar density { 0 };

			for (std::size_t k = 1; k < concentrationsBoundary.size(); ++k)
			{
				concentrationsBoundary[k] = boundaryConditionValueCell(
						cellFields.concentration.v[k].ref()[cellIndex],
						cellFields.concentration.v[k].boundCond()[surfaceIndex],
						cellIndex, surfaceIndex);

				concentrationsBoundary[0] += concentrationsBoundary[k];

				density += concentrationsBoundary[k] * mix.Mv()[k - 1];
			}

			retValue = mix.Cv(concentrationsBoundary)
					* concentrationsBoundary[0] / density;
		}
			break;
		case boundaryConditionType::calculatedDensity:
		{
			scalar density { 0 };

			for (std::size_t k = 1; k < cellFields.concentration.v.size(); ++k)
			{
				const scalar concentrationsBoundary_k {
						boundaryConditionValueCell(
								cellFields.concentration.v[k].ref()[cellIndex],
								cellFields.concentration.v[k].boundCond()[surfaceIndex],
								cellIndex, surfaceIndex) };

				density += concentrationsBoundary_k * mix.Mv()[k - 1];
			}

			retValue = density;
		}
			break;
		default:
			throw exception("Unknown type of boundary condition.",
					errors::boundaryConditionError);
			break;
		}

		return retValue;
	}
};
}

#endif /* BOUNDARYCONDITIONVALUE_HPP_ */
