/*
 * transportCoefficients.hpp
 *
 *  Created on: 2020/01/25
 *      Author: Maxim Boldyrev
 *
 *      Structure storing fields for diffusion coefficients.
 */

#ifndef TRANSPORTCOEFFICIENTS_HPP_
#define TRANSPORTCOEFFICIENTS_HPP_

#include <iostream>

#include "turbulenceModelEnum.hpp"
#include "abstractTurbulentParameters.hpp"

namespace schemi
{
template<typename typeOfEntity>
struct transportCoefficients
{
	field<scalar, typeOfEntity> tNu;
	field<scalar, typeOfEntity> tD;
	field<scalar, typeOfEntity> tKappa;
	field<scalar, typeOfEntity> tLambda;
	field<scalar, typeOfEntity> k_D;
	field<scalar, typeOfEntity> eps_D;
	field<scalar, typeOfEntity> a_D;
	field<scalar, typeOfEntity> b_D;

	field<scalar, typeOfEntity> physMu;
	field<std::valarray<std::valarray<scalar>>, typeOfEntity> physDm;
	field<std::valarray<scalar>, typeOfEntity> physDs;
	field<scalar, typeOfEntity> physKappa;

	explicit transportCoefficients(const mesh & meshRef) noexcept :
			tNu { meshRef, 0 },

			tD { meshRef, 0 },

			tKappa { meshRef, 0 },

			tLambda { meshRef, 0 },

			k_D { meshRef, 0 },

			eps_D { meshRef, 0 },

			a_D { meshRef, 0 },

			b_D { meshRef, 0 },

			physMu { meshRef, 0 },

			physDm { meshRef, std::valarray<std::valarray<scalar>>(0) },

			physDs { meshRef, std::valarray<scalar>(0) },

			physKappa { meshRef, 0 }
	{
	}

	void calculateCoefficients(const field<scalar, typeOfEntity> & kField,
			const field<scalar, typeOfEntity> & epsField,
			const abstractTurbulentParameters & tp) noexcept
	{
		tNu.ref_r() = tp.calculateNut(kField.ref(), epsField.ref());

		calculateCoefficients(tp);
	}

	void calculateCoefficients(const abstractTurbulentParameters & tp) noexcept
	{
		tD.ref_r() = tNu.ref() / tp.sigmaSc();
		tKappa.ref_r() = tNu.ref() / tp.sigmaT();
		tLambda.ref_r() = tNu.ref() / tp.sigmaE();
		k_D.ref_r() = tNu.ref() / tp.sigmaK();
		eps_D.ref_r() = tNu.ref() / tp.sigmaEps();
		a_D.ref_r() = tNu.ref() / tp.sigmaa();
		b_D.ref_r() = tNu.ref() / tp.sigmab();
	}

	void tAssign(const scalar val) noexcept
	{
		tNu.ref_r() = val;
		tD.ref_r() = val;
		tKappa.ref_r() = val;
		tLambda.ref_r() = val;
		k_D.ref_r() = val;
		eps_D.ref_r() = val;
		a_D.ref_r() = val;
		b_D.ref_r() = val;
	}

	auto& operator=(const scalar val) noexcept
	{
		tNu.ref_r() = val;
		tD.ref_r() = val;
		tKappa.ref_r() = val;
		tLambda.ref_r() = val;
		k_D.ref_r() = val;
		eps_D.ref_r() = val;
		a_D.ref_r() = val;
		b_D.ref_r() = val;

		physMu.ref_r() = val;
		for (auto & dm_i : physDm.ref_r())
			for (std::size_t k1 = 0; k1 < dm_i.size(); ++k1)
				for (std::size_t k2 = 0; k2 < dm_i[k1].size(); ++k2)
				{
					if (k1 == k2)
						continue;

					dm_i[k1][k2] = val;
				}
		physDs.ref_r() = std::valarray<scalar>(val, physDs.ref()[0].size);
		physKappa.ref_r() = val;

		return *this;
	}

	void calculateCoefficients(const abstractTransportModel & model,
			const std::valarray<scalar> & M,
			const field<scalar, typeOfEntity> & temperature,
			const field<scalar, typeOfEntity> & pressure,
			const concentrationsPack<typeOfEntity> & concentrations) noexcept
	{
		physMu = model.calculateMu(M, temperature, concentrations);
		physDm = model.calculateDm(M, temperature, pressure, concentrations);
		physDs = model.calculateDs(M, temperature, pressure, concentrations);
		physKappa = model.calculateKappa(M, temperature, concentrations);
	}
};

template<typename typeOfEntity>
struct effectiveTransportCoefficients: transportCoefficients<typeOfEntity>
{
	field<scalar, typeOfEntity> mu;
	field<std::valarray<vector>, typeOfEntity> DFlux;
	std::vector<field<scalar, typeOfEntity>> rhoD;
	field<scalar, typeOfEntity> kappa;
	field<scalar, typeOfEntity> rhoDk;
	field<scalar, typeOfEntity> rhoDeps;
	field<scalar, typeOfEntity> rhoDa;
	field<scalar, typeOfEntity> rhoDb;

	effectiveTransportCoefficients(const mesh & meshRef,
			const std::size_t compt) noexcept :
			transportCoefficients<typeOfEntity>(meshRef),

			mu { meshRef, 0 },

			DFlux(meshRef, std::valarray<vector>(vector(0.), compt)),

			rhoD(compt, field<scalar, typeOfEntity>(meshRef, 0.)),

			kappa { meshRef, 0 },

			rhoDk { meshRef, 0 },

			rhoDeps { meshRef, 0 },

			rhoDa { meshRef, 0 },

			rhoDb { meshRef, 0 }
	{
	}

	field<scalar, typeOfEntity> maxValue(const bool turbulenceFlag,
			const turbulenceModel sourceFlag) const noexcept
	{
		auto maxFieldVal = this->physKappa;

		for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
		{
			const auto D = this->physDs.ref()[i];

			maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i], D.max());
			maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
					this->physMu.ref()[i]);
		}

		if (turbulenceFlag)
		{
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						this->tNu.ref()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						this->tD.ref()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						2 * this->tLambda.ref()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						this->k_D.ref()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
						this->eps_D.ref()[i]);

			if ((sourceFlag == turbulenceModel::BHRSource)
					|| (sourceFlag == turbulenceModel::BHRKLSource)
					|| (sourceFlag == turbulenceModel::kEpsASource))
			{
				for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
					maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
							this->a_D.ref()[i]);

				if ((sourceFlag == turbulenceModel::BHRSource)
						|| (sourceFlag == turbulenceModel::BHRKLSource))
					for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
						maxFieldVal.ref_r()[i] = std::max(maxFieldVal.ref()[i],
								this->b_D.ref()[i]);
			}
		}

		return maxFieldVal;
	}

	void calculateEffectiveCoefficients(const field<scalar, typeOfEntity> & rho,
			const abstractTurbulenceGen & t,
			const field<scalar, typeOfEntity> & conc,
			const field<scalar, typeOfEntity> & Cv) noexcept
	{
		if (t.turbulence)
		{
			mu.ref_r() = this->physMu.ref() + rho.ref() * this->tNu.ref();

			for (std::size_t k = 0; k < rhoD.size(); ++k)
				for (std::size_t i = 0; i < rho.size(); ++i)
					rhoD[k].ref_r()[i] = rho.ref()[i]
							* (this->tD.ref()[i] + this->physDs.ref()[i][k]);

			kappa.ref_r() = this->physKappa.ref()
					+ rho.ref() * this->tLambda.ref()
					+ conc.ref() * Cv.ref() * this->tKappa.ref();

			rhoDk.ref_r() = this->physMu.ref() + rho.ref() * this->k_D.ref();
			rhoDeps.ref_r() = this->physMu.ref() + rho.ref() * this->k_D.ref();
			if ((t.model == turbulenceModel::BHRSource)
					|| (t.model == turbulenceModel::BHRKLSource)
					|| (t.model == turbulenceModel::kEpsASource))
			{
				rhoDa.ref_r() = this->physMu.ref()
						+ rho.ref() * this->k_D.ref();
				rhoDb.ref_r() = this->physMu.ref()
						+ rho.ref() * this->k_D.ref();
			}
		}
		else
		{
			mu.ref_r() = this->physMu.ref();

			for (std::size_t k = 0; k < rhoD.size(); ++k)
				for (std::size_t i = 0; i < rho.size(); ++i)
					rhoD[k].ref_r()[i] = rho.ref()[i]
							* this->physDs.ref()[i][k];

			kappa.ref_r() = this->physKappa.ref();
		}
	}

	void caclulateDFluxes(
			const std::vector<field<vector, typeOfEntity>> & gradX,
			const std::valarray<scalar> & M,
			const concentrationsPack<typeOfEntity> & N, const bool turb)
	{
		const std::size_t Ncomp = M.size(), Ncomp1 = Ncomp + 1;

		if (!turb)
			for (std::size_t i = 0; i < DFlux.size(); ++i)
			{
				std::valarray<scalar> c(Ncomp);
				for (std::size_t k = 1; k < Ncomp1; ++k)
					c[k - 1] = N.v[k].ref()[i];

				std::size_t replaceIndex { 0 };
				const scalar replaceC { c.max() };
				for (std::size_t k = 0; k < c.size(); ++k)
					if (replaceC == c[k])
					{
						replaceIndex = k + 1;
						break;
					}

				for (std::size_t f = 0; f < vector::vsize; ++f)
				{
					std::valarray<std::valarray<scalar>> cellDFluxesMatrix(
							std::valarray<scalar>(0., Ncomp), Ncomp);
					std::valarray<scalar> freeTerm(0., Ncomp);

					for (std::size_t k = 0; k < Ncomp; ++k)
						freeTerm[k] = -gradX[k].ref()[i].v()[f];

					for (std::size_t k1 = 1; k1 < Ncomp1; ++k1)
						if (k1 != replaceIndex)
							for (std::size_t k2 = 1; k2 < Ncomp1; ++k2)
							{
								if (k1 == k2)
									continue;

								const auto & N1 = N.v[k1].ref()[i];
								const auto & N2 = N.v[k2].ref()[i];
								const auto & N0 = N.v[0].ref()[i];

								const scalar Aij =
										N1 * N2
												/ (N0 * N0
														* this->physDm.ref()[i][k1
																- 1][k2 - 1]);

								cellDFluxesMatrix[k1 - 1][k2 - 1] = -Aij;
								cellDFluxesMatrix[k1 - 1][k1 - 1] += Aij;
							}
						else
						{
							freeTerm[k1 - 1] = 0;

							for (std::size_t k2 = 1; k2 < Ncomp1; ++k2)
								cellDFluxesMatrix[k1 - 1][k2 - 1] =
										N.v[k2].ref()[i] * M[k2 - 1];
						}

					const std::valarray<scalar> resFlows =
							GaussEliminationSolver(cellDFluxesMatrix, freeTerm)
									* M * c;

					for (std::size_t k = 0; k < Ncomp; ++k)
						DFlux.ref_r()[i][k].v_r()[f] = resFlows[k];
				}
			}
		else
			for (std::size_t i = 0; i < DFlux.size(); ++i)
			{
				std::valarray<scalar> c(Ncomp);
				for (std::size_t k = 1; k < Ncomp1; ++k)
					c[k - 1] = N.v[k].ref()[i];

				std::size_t replaceIndex { 0 };
				const scalar replaceC { c.max() };
				for (std::size_t k = 0; k < c.size(); ++k)
					if (replaceC == c[k])
					{
						replaceIndex = k + 1;
						break;
					}

				for (std::size_t f = 0; f < vector::vsize; ++f)
				{
					std::valarray<std::valarray<scalar>> cellDFluxesMatrix(
							std::valarray<scalar>(0., Ncomp), Ncomp);
					std::valarray<scalar> freeTerm(0., Ncomp);

					for (std::size_t k = 0; k < Ncomp; ++k)
						freeTerm[k] = -gradX[k].ref()[i].v()[f];

					for (std::size_t k1 = 1; k1 < Ncomp1; ++k1)
						if (k1 != replaceIndex)
							for (std::size_t k2 = 1; k2 < Ncomp1; ++k2)
							{
								if (k1 == k2)
									continue;

								const auto & N1 = N.v[k1].ref()[i];
								const auto & N2 = N.v[k2].ref()[i];
								const auto & N0 = N.v[0].ref()[i];

								const scalar Aij =
										N1 * N2
												/ (N0 * N0
														* (this->physDm.ref()[i][k1
																- 1][k2 - 1]
																+ this->tD.ref()[i]));

								cellDFluxesMatrix[k1 - 1][k2 - 1] = -Aij;
								cellDFluxesMatrix[k1 - 1][k1 - 1] += Aij;
							}
						else
						{
							freeTerm[k1 - 1] = 0;

							for (std::size_t k2 = 1; k2 < Ncomp1; ++k2)
								cellDFluxesMatrix[k1 - 1][k2 - 1] =
										N.v[k2].ref()[i] * M[k2 - 1];
						}

					const std::valarray<scalar> resFlows =
							GaussEliminationSolver(cellDFluxesMatrix, freeTerm)
									* M * c;

					for (std::size_t k = 0; k < Ncomp; ++k)
						DFlux.ref_r()[i][k].v_r()[f] = resFlows[k];
				}
			}
	}
private:
	std::valarray<scalar> GaussEliminationSolver(
			std::valarray<std::valarray<scalar>> A,
			std::valarray<scalar> b) const noexcept
	{
		const std::size_t N = b.size();

		for (std::size_t k = 0; k < N - 1; ++k)
			for (std::size_t i = k + 1; i < N; ++i)
			{
				const auto ratio = A[i][k] / (A[k][k] + stabilizator);

				for (std::size_t j = 0; j < N; ++j)
					A[i][j] = A[i][j] - ratio * A[k][j];

				b[i] = b[i] - ratio * b[k];
			}

		std::valarray<scalar> phi(b.size());

		phi[N - 1] = b[N - 1] / (A[N - 1][N - 1] + stabilizator);

		for (std::size_t i = N - 2;; --i)
		{
			scalar term = 0.0;

			for (std::size_t j = i + 1; j < N; ++j)
				term += A[i][j] * phi[j];

			phi[i] = (b[i] - term) / (A[i][i] + stabilizator);

			if (i == 0)
				break;
		}

		normalize(phi);

		return phi;
	}

	std::valarray<scalar> GaussSeidelSolver(
			const std::valarray<std::valarray<scalar>> & A,
			const std::valarray<scalar> & b) const
	{
		std::valarray<scalar> oldIteration(0., b.size());

		std::valarray<scalar> newIteration(oldIteration);

		std::size_t nIterations { 0 };

		while (true)
		{
			nIterations++;

			for (std::size_t i = 0; i < oldIteration.size(); ++i)
			{
				const scalar aii { 1. / (A[i][i] + stabilizator) };

				const scalar bi { b[i] };

				newIteration[i] = bi * aii;

				for (std::size_t j = 0; j < i; ++j)
					newIteration[i] -= A[i][j] * newIteration[j] * aii;

				for (std::size_t j = i + 1; j < oldIteration.size(); ++j)
					newIteration[i] -= A[i][j] * oldIteration[j] * aii;
			}

			for (std::size_t i = oldIteration.size() - 1;; --i)
			{
				const scalar aii { 1. / (A[i][i] + stabilizator) };

				const scalar bi { b[i] };

				newIteration[i] = bi * aii;

				if (i != 0)
					for (std::size_t j = i - 1;; --j)
					{
						newIteration[i] -= A[i][j] * oldIteration[j] * aii;

						if (j == 0)
							break;
					}

				if (i != oldIteration.size() - 1)
					for (std::size_t j = oldIteration.size() - 1;; --j)
					{
						newIteration[i] -= A[i][j] * newIteration[j] * aii;

						if (j == i + 1)
							break;
					}

				if (i == 0)
					break;
			}

			const scalar diff {
					100.
							* std::abs(
									(newIteration - oldIteration)
											/ (std::abs(newIteration)
													+ stabilizator)).max() };

			if (diff < 1E-12)
			{
				normalize(newIteration);
				return newIteration;
			}
			else if (nIterations >= 100)
			{
				std::clog
						<< "Gauss-Seidel algorithm for diffusive flows did not converged. Difference is: "
						<< diff << std::endl;

				normalize(newIteration);
				return newIteration;
			}
			else
				oldIteration = newIteration;
		}
	}

	void normalize(std::valarray<scalar> & res) const noexcept
	{
		for (auto & i : res)
			if (std::abs(i) < std::numeric_limits<scalar>::epsilon())
				i = 0;
	}
};
}  // namespace schemi

#endif /* TRANSPORTCOEFFICIENTS_HPP_ */
