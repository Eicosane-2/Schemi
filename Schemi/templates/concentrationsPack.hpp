/*
 * concentrationsPack.hpp
 *
 *  Created on: 2020/02/19
 *      Author: Maxim Boldyrev
 *
 *      Structure for storing concentrations' fields.
 */

#ifndef CONCENTRATIONSPACK_HPP_
#define CONCENTRATIONSPACK_HPP_

#include "field.hpp"

namespace schemi
{
template<typename typeOfEntity>
struct concentrationsPack
{
	std::vector<field<scalar, typeOfEntity>> v;
	std::vector<const std::valarray<scalar>*> p;

	concentrationsPack(const mesh & meshRef,
			const std::size_t numberOfComponents) noexcept :
			v { numberOfComponents + 1, field<scalar, typeOfEntity>(meshRef, 0) },

			p { numberOfComponents + 1 }
	{
		for (std::size_t i = 0; i < p.size(); ++i)
			p[i] = &v[i]();
	}

	concentrationsPack(const concentrationsPack<typeOfEntity> & in) noexcept :
			v(in.v), p(in.p.size())
	{
		for (std::size_t i = 0; i < p.size(); ++i)
			p[i] = &v[i]();
	}

	auto& operator=(const concentrationsPack<typeOfEntity> & in) = delete;
};
}  // namespace schemi

#endif /* CONCENTRATIONSPACK_HPP_ */
