/*
 * GoncharovTracerModel.hpp
 *
 *  Created on: 2025/06/28
 *      Author: Maxim Boldyrev
 */

#ifndef GONCHAROVTRACERMODEL_HPP_
#define GONCHAROVTRACERMODEL_HPP_

#include <functional>

#include "boundaryConditionValue.hpp"
#include "concentrationsPack.hpp"
#include "interfaceStatusEnum.hpp"
#include "scalar.hpp"
#include "vector.hpp"
#include "tracerParticle.hpp"
#include "volumeField.hpp"

namespace schemi
{
class GoncharovTracerModel: public tracerParticle
{
	constexpr static scalar ReCriterion { 300 };
	constexpr static scalar relativeCriterion { 5.0 / 3.0 };

	constexpr static scalar Ck { 0.35 }, Ceps { 1.05 }, Cb { 1 };

	/*Initial conditions*/
	std::array<std::size_t, 2> s12 { componentPlaceholder, componentPlaceholder };
	int c { 0 };
	scalar eta_0 { 0 };
	scalar lambda { 1E15 }, k { 0 };
	scalar radiusOfIfluence { 1E15 };

	/*Calculated data*/
	scalar eta { 0 }, eta_1 { 0 }, eta_2 { 0 };
	scalar rho1 { 0 }, rho2 { 0 };
	scalar k0 { 0 }, eps0 { 0 }, b0 { 0 };
	vector a0 { 0, 0, 0 };
	scalar Re { 0 }, relDeltaEta { 0 };
	scalar timestep_1 { 0 };

	interfaceStatus status { interfaceStatus::notDeveloped };

	std::function<bool(const scalar, const scalar)> transitioniCriterion {
			nullptr };

	void initialiseTransitionCheck(const std::string & initCheckMethod);

protected:
	constexpr static scalar cellCoefficient { 6 };

	void setEta(const scalar etaIn) noexcept
	{
		eta = etaIn;
	}
	void setStatus(const interfaceStatus statusIn) noexcept
	{
		status = statusIn;
	}
	void setk0(const scalar kIn) noexcept
	{
		k0 = kIn;
	}
	void seteps0(const scalar epsIn) noexcept
	{
		eps0 = epsIn;
	}
	void seta0(const vector & aIn) noexcept
	{
		a0 = aIn;
	}
	void setb0(const scalar bIn) noexcept
	{
		b0 = bIn;
	}
public:
	constexpr static scalar C01 { 4 }, C02 { 0.5 };
	GoncharovTracerModel() noexcept = default;
	GoncharovTracerModel(const std::string & initCheckMethod,
			const vector & inPos, const vector & inVelocity,
			const std::size_t sub1, const std::size_t sub2, const int pertType,
			const scalar eta0In, const scalar lambdaIn,
			const scalar radiusOfIfluenceIn);
	GoncharovTracerModel(const std::string & initCheckMethod,
			const vector & inPos, const vector & inPos1,
			const std::array<vector, 4> & inVelocity, const std::size_t inStep,
			const std::size_t sub1, const std::size_t sub2, const int pertType,
			const scalar eta0In, const scalar lambdaIn,
			const scalar radiusOfIfluenceIn, const scalar etaCur,
			const scalar eta1Cur, const scalar eta2Cur, const scalar rho1Cur,
			const scalar rho2Cur, const scalar k0Cur, const scalar eps0Cur,
			const scalar b0Cur, const vector & a0Cur, const scalar timeStepOld,
			const interfaceStatus curStatus) noexcept;
	virtual ~GoncharovTracerModel() noexcept override = default;

	const std::array<std::size_t, 2>& substances() const noexcept
	{
		return s12;
	}
	scalar getRadiusOfInfluence() const noexcept
	{
		return radiusOfIfluence;
	}
	scalar getEta() const noexcept
	{
		return eta;
	}
	interfaceStatus getStatus() const noexcept
	{
		return status;
	}
	scalar getk0() const noexcept
	{
		return k0;
	}
	scalar geteps0() const noexcept
	{
		return eps0;
	}
	const vector& geta0() const noexcept
	{
		return a0;
	}
	scalar getb0() const noexcept
	{
		return b0;
	}
	scalar getCurRe() const noexcept
	{
		return Re;
	}
	scalar getCurRelDeltaEta() const noexcept
	{
		return relDeltaEta;
	}

	virtual void writeOutput(std::ofstream & output) const override;

	virtual void timeIntegration(const scalar density1, const scalar density2,
			const vector & gradRho, const vector & u, const vector & g,
			const scalar timestep, const vector & gradP = vector(0),
			const scalar rho = 0, const scalar divU = 0, const tensor & gradU =
					tensor(0));

	virtual interfaceStatus checkTransition(const scalar nu,
			const scalar timestep, const vector & cellRadius,
			const vector & surfaceRadius) noexcept;
};
}  // namespace schemi

#endif /* GONCHAROVTRACERMODEL_HPP_ */
