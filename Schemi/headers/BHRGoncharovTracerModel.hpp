/*
 * BHRGoncharovTracerModel.hpp
 *
 *  Created on: 2025/12/06
 *      Author: Maxim Boldyrev
 */

#ifndef BHRGONCHAROVTRACERMODEL_HPP_
#define BHRGONCHAROVTRACERMODEL_HPP_

#include "GoncharovTracerModel.hpp"

namespace schemi
{
class BHRGoncharovTracerModel: public GoncharovTracerModel
{
	scalar Cmu { 0 }, C0 { 0 }, C2 { 0 }, C3 { 0 }, Ca1 { 0 }, Cb1 { 0 };
	constexpr static scalar Ceta { 1.75 };
public:
	BHRGoncharovTracerModel() noexcept;
	BHRGoncharovTracerModel(const std::string & initCheckMethod,
			const vector & inPos, const vector & inVelocity,
			const std::size_t sub1, const std::size_t sub2, const int pertType,
			const scalar eta0In, const scalar lambdaIn,
			const scalar radiusOfIfluenceIn, const scalar CmuIn,
			const scalar C0In, const scalar C2In, const scalar C3In,
			const scalar Ca1In, const scalar Cb1In);
	BHRGoncharovTracerModel(const std::string & initCheckMethod,
			const vector & inPos, const vector & inPos1,
			const std::array<vector, 4> & inVelocity, const std::size_t inStep,
			const std::size_t sub1, const std::size_t sub2, const int pertType,
			const scalar eta0In, const scalar lambdaIn,
			const scalar radiusOfIfluenceIn, const scalar etaCur,
			const scalar eta1Cur, const scalar eta2Cur, const scalar rho1Cur,
			const scalar rho2Cur, const scalar k0Cur, const scalar eps0Cur,
			const scalar b0Cur, const vector & a0Cur,
			const interfaceStatus curStatus, const scalar CmuIn,
			const scalar C0In, const scalar C2In, const scalar C3In,
			const scalar Ca1In, const scalar Cb1In);
	virtual ~BHRGoncharovTracerModel() noexcept override;

	virtual void timeIntegration(const scalar density1, const scalar density2,
			const vector & gradRho, const vector & u, const scalar timestep,
			const vector & gradP, const scalar rho, const scalar divU,
			const tensor & gradU) override;

	virtual interfaceStatus checkTransition(const scalar nu,
			const scalar timestep, const vector & cellRadius,
			const vector & surfaceRadius) noexcept override;
};
}

#endif /* BHRGONCHAROVTRACERMODEL_HPP_ */
