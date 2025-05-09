/*
 * diffusiveFields.cpp
 *
 *  Created on: 2023/06/04
 *      Author: Maxim Boldyrev
 */

#include "diffusiveFields.hpp"

#include <algorithm>

schemi::diffusiveFields::diffusiveFields(const mesh & meshRef,
		const bunchOfFields<cubicCell> & cellFields,
		const std::vector<boundaryConditionType> & commBoundCond,
		const bool turbulenceFlag, const bool aField, const bool bField) :
		massFraction { cellFields.concentration.v.size() - 1,
				volumeField<scalar> { meshRef, 0 } },

		velocity { meshRef, vector(0), cellFields.velocity.boundCond() },

		temperature { meshRef, 0, cellFields.temperature.boundCond() },

		k { meshRef, scalar { 0 }, cellFields.kTurb.boundCond() },

		eps { meshRef, 0, cellFields.epsTurb.boundCond() },

		a { meshRef, vector(0), cellFields.aTurb.boundCond() },

		b { meshRef, 0, cellFields.bTurb.boundCond() }
{
	std::vector<boundaryConditionType> bndCon(commBoundCond);
	std::replace(bndCon.begin(), bndCon.end(),
			boundaryConditionType::calculated,
			boundaryConditionType::calculatedMassFraction);
	std::fill(massFraction.begin(), massFraction.end(),
			volumeField<scalar>(meshRef, 0, subPatchData<scalar> { bndCon[0] },
					subPatchData<scalar> { bndCon[1] }, subPatchData<scalar> {
							bndCon[2] }, subPatchData<scalar> { bndCon[3] },
					subPatchData<scalar> { bndCon[4] }, subPatchData<scalar> {
							bndCon[5] }));

	std::valarray<scalar> sumMassFrac(scalar(0), meshRef.cellsSize());
	for (std::size_t k_ind = 0; k_ind < massFraction.size(); ++k_ind)
	{
		massFraction[k_ind].r() = cellFields.density[k_ind + 1]()
				/ cellFields.density[0]();
		sumMassFrac += massFraction[k_ind]();
	}

	for (auto & diffFieldMassFrac_k : massFraction)
		diffFieldMassFrac_k.r() /= sumMassFrac;

	velocity.r() = cellFields.velocity();
	temperature.r() = cellFields.temperature();
	if (turbulenceFlag)
	{
		k.r() = cellFields.kTurb();
		eps.r() = cellFields.epsTurb();
		if (aField)
		{
			a.r() = cellFields.aTurb();
			if (bField)
				b.r() = cellFields.bTurb();
		}
	}
}
