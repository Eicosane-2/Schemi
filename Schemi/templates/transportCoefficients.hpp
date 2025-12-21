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

#include <algorithm>
#include <iostream>

#include "turbulenceModelEnum.hpp"
#include "abstractTurbulenceModel.hpp"

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
			const abstractTurbulenceModel & t) noexcept
	{
		tNu.val() = t.calculateNut(kField.cval(), epsField.cval());

		calculateCoefficients(t);
	}

	void calculateCoefficients(const abstractTurbulenceModel & t) noexcept
	{
		tD = tNu / t.sigmaSc();
		tKappa = tNu / t.sigmaT();
		tLambda = tNu / t.sigmaE();
		k_D = tNu / t.sigmaK();
		eps_D = tNu / t.sigmaEps();
		a_D = tNu / t.sigmaa();
		b_D = tNu / t.sigmab();
	}

	void tAssign(const scalar val) noexcept
	{
		tNu.val() = val;
		tD.val() = val;
		tKappa.val() = val;
		tLambda.val() = val;
		k_D.val() = val;
		eps_D.val() = val;
		a_D.val() = val;
		b_D.val() = val;
	}

	auto& operator=(const scalar val) noexcept
	{
		tNu.val() = val;
		tD.val() = val;
		tKappa.val() = val;
		tLambda.val() = val;
		k_D.val() = val;
		eps_D.val() = val;
		a_D.val() = val;
		b_D.val() = val;

		physMu.val() = val;
		for (auto & dm_i : physDm.val())
			for (std::size_t k1 = 0; k1 < dm_i.size(); ++k1)
				for (std::size_t k2 = 0; k2 < dm_i[k1].size(); ++k2)
				{
					if (k1 == k2)
						continue;

					dm_i[k1][k2] = val;
				}
		physDs.val() = std::valarray<scalar>(val, physDs()[0].size);
		physKappa.val() = val;

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
			const bool aField, const bool bField,
			const field<scalar, typeOfEntity> & rho) const noexcept
	{
		auto maxFieldVal = this->physKappa;

		for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
		{
			const auto & D = this->physDs.cval()[i];

			maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
					rho.cval()[i] * D.max());
			maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
					this->physMu.cval()[i]);
		}

		if (turbulenceFlag)
		{
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
						rho.cval()[i] * this->tNu.cval()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
						rho.cval()[i] * this->tD.cval()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
						rho.cval()[i] * this->tLambda.cval()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
						rho.cval()[i] * this->k_D.cval()[i]);
			for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
				maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
						rho.cval()[i] * this->eps_D.cval()[i]);

			if (aField)
			{
				for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
					maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
							rho.cval()[i] * this->a_D.cval()[i]);

				if (bField)
					for (std::size_t i = 0; i < maxFieldVal.size(); ++i)
						maxFieldVal.val()[i] = std::max(maxFieldVal.cval()[i],
								rho.cval()[i] * this->b_D.cval()[i]);
			}
		}

		return maxFieldVal;
	}

	void calculateEffectiveCoefficients(const field<scalar, typeOfEntity> & rho,
			const abstractTurbulenceModel & t,
			const field<scalar, typeOfEntity> & conc,
			const field<scalar, typeOfEntity> & Cv) noexcept
	{
		if (t.turbulence())
		{
			mu = this->physMu + rho * this->tNu;

			for (std::size_t k = 0; k < rhoD.size(); ++k)
				for (std::size_t i = 0; i < rho.size(); ++i)
					rhoD[k].val()[i] = rho.cval()[i]
							* (this->physDs.cval()[i][k] + this->tD.cval()[i]);

			kappa = this->physKappa + conc * Cv * this->tLambda
					+ conc * Cv * this->tKappa;

			rhoDk = this->physMu + rho * this->k_D;
			rhoDeps = this->physMu + rho * this->eps_D;
			if (t.aField())
			{
				rhoDa = this->physMu + rho * this->a_D;
				if (t.bField())
					rhoDb = this->physMu + rho * this->b_D;
			}
		}
		else
		{
			mu = this->physMu;

			for (std::size_t k = 0; k < rhoD.size(); ++k)
				for (std::size_t i = 0; i < rho.size(); ++i)
					rhoD[k].val()[i] = rho.cval()[i]
							* this->physDs.cval()[i][k];

			kappa = this->physKappa;
		}
	}

	void caclulateDFluxes(
			const std::vector<field<vector, typeOfEntity>> & gradX,
			const std::valarray<scalar> & M,
			const concentrationsPack<typeOfEntity> & N, const bool turb)
	{
		const std::size_t Ncomp = M.size(), Ncomp1 = Ncomp + 1;

		if (Ncomp > 1)
			if (!turb)
				for (std::size_t i = 0; i < DFlux.size(); ++i)
				{
					std::valarray<scalar> c(Ncomp);
					for (std::size_t k = 1; k < Ncomp1; ++k)
						c[k - 1] = N.v[k].cval()[i];

					std::size_t replaceIndex { 0 };
					const scalar replaceC { c.min() };
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
							freeTerm[k] = -gradX[k].cval()[i]()[f];

						for (std::size_t k1 = 1; k1 < Ncomp1; ++k1)
							if (k1 != replaceIndex)
								for (std::size_t k2 = 1; k2 < Ncomp1; ++k2)
								{
									if (k1 == k2)
										continue;

									const auto & N1 = N.v[k1].cval()[i];
									const auto & N2 = N.v[k2].cval()[i];
									const auto & N0 = N.v[0].cval()[i];

									const scalar Aij = N1 * N2
											/ (N0 * N0
													* this->physDm.cval()[i][k1
															- 1][k2 - 1]);

									cellDFluxesMatrix[k1 - 1][k2 - 1] = -Aij;
									cellDFluxesMatrix[k1 - 1][k1 - 1] += Aij;
								}
							else
							{
								freeTerm[k1 - 1] = 0;

								for (std::size_t k2 = 1; k2 < Ncomp1; ++k2)
									cellDFluxesMatrix[k1 - 1][k2 - 1] =
											N.v[k2].cval()[i] * M[k2 - 1];
							}

						const std::valarray<scalar> resFlows =
								GaussEliminationSolver(cellDFluxesMatrix,
										freeTerm) * M * c;

						for (std::size_t k = 0; k < Ncomp; ++k)
							DFlux.val()[i][k].wr()[f] = resFlows[k];
					}
				}
			else
				for (std::size_t i = 0; i < DFlux.size(); ++i)
				{
					std::valarray<scalar> c(Ncomp);
					for (std::size_t k = 1; k < Ncomp1; ++k)
						c[k - 1] = N.v[k].cval()[i];

					std::size_t replaceIndex { 0 };
					const scalar replaceC { c.min() };
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
							freeTerm[k] = -gradX[k].cval()[i]()[f];

						for (std::size_t k1 = 1; k1 < Ncomp1; ++k1)
							if (k1 != replaceIndex)
								for (std::size_t k2 = 1; k2 < Ncomp1; ++k2)
								{
									if (k1 == k2)
										continue;

									const auto & N1 = N.v[k1].cval()[i];
									const auto & N2 = N.v[k2].cval()[i];
									const auto & N0 = N.v[0].cval()[i];

									const scalar Aij =
											N1 * N2
													/ (N0 * N0
															* (this->physDm.cval()[i][k1
																	- 1][k2 - 1]
																	+ this->tD.cval()[i]));

									cellDFluxesMatrix[k1 - 1][k2 - 1] = -Aij;
									cellDFluxesMatrix[k1 - 1][k1 - 1] += Aij;
								}
							else
							{
								freeTerm[k1 - 1] = 0;

								for (std::size_t k2 = 1; k2 < Ncomp1; ++k2)
									cellDFluxesMatrix[k1 - 1][k2 - 1] =
											N.v[k2].cval()[i] * M[k2 - 1];
							}

						const std::valarray<scalar> resFlows =
								GaussEliminationSolver(cellDFluxesMatrix,
										freeTerm) * M * c;

						for (std::size_t k = 0; k < Ncomp; ++k)
							DFlux.val()[i][k].wr()[f] = resFlows[k];
					}
				}
		else
			for (std::size_t i = 0; i < DFlux.size(); ++i)
				DFlux.val()[i][0].wr() = { 0, 0, 0 };
	}
private:
	std::valarray<scalar> GaussEliminationSolver(
			std::valarray<std::valarray<scalar>> A,
			std::valarray<scalar> b) const noexcept
	{
		const std::size_t N = b.size();

		/*Pivoting*/
		for (std::size_t k = 0; k < (N - 1); ++k)
		{
			/*Find maximum value in column's next values.*/
			std::pair<scalar, std::size_t> maxVal(std::make_pair(A[k][k], k));
			for (std::size_t c = k + 1; c < N; ++c)
				if (A[c][k] > maxVal.first)
					maxVal = std::make_pair(A[c][k], c);

			const auto c = maxVal.second;

			if (c != k)
			{
				auto & origArr = A[c];
				auto & destArr = A[k];

				std::swap(origArr, destArr);
				std::swap(b[c], b[k]);
			}
		}
		/**/

		for (std::size_t k = 0; k < N - 1; ++k)
			for (std::size_t i = k + 1; i < N; ++i)
			{
				const auto ratio = A[i][k] / (A[k][k] + stabilizator);

				for (std::size_t j = k + 1; j < N; ++j)
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

			if (diff < convergenceToleranceGlobal)
			{
				normalize(newIteration);
				return newIteration;
			}
			else if (nIterations >= 100)
			//[[unlikely]]
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
		std::replace_if(std::begin(res), std::end(res), [](const auto & i) 
		{
			return std::abs(i) < std::numeric_limits<scalar>::epsilon();
		}, 0);
	}
};
}  // namespace schemi

#endif /* TRANSPORTCOEFFICIENTS_HPP_ */
