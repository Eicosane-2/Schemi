/*
 * starFields.cpp
 *
 *  Created on: 2023/06/03
 *      Author: Maxim Boldyrev
 */

#include "starFields.hpp"

schemi::starFields::starFields(const mesh & meshRef,
		const std::size_t numberOfcomponents) :
		c { meshRef, numberOfcomponents },

		rho { meshRef, 0 },

		v { meshRef, vector(0) },

		p { meshRef, 0 },

		a { meshRef, vector(0) },

		b { meshRef, 0 }
{
}

