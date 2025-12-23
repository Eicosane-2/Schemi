/*
 * intExpPow.hpp
 *
 *  Created on: 2023/01/19
 *      Author: Maxim Boldyrev
 */

#ifndef INTEXPPOW_HPP_
#define INTEXPPOW_HPP_

namespace schemi
{
#if __cplusplus > 201703L
template<typename T>
concept nonPowIntegral = std::same_as<T, bool> || std::same_as<T, char>
|| std::same_as<T, unsigned char> || std::same_as<T, signed char>
|| std::same_as<T, char8_t> || std::same_as<T, char16_t>
|| std::same_as<T, char32_t> || std::same_as<T, wchar_t>;

template<typename T>
concept powNumber = std::floating_point<T>
|| (std::integral<T> && !nonPowIntegral<T> );

template<powNumber T, int expon>
#else
		template<typename T, int expon>
#endif
		constexpr T pow(const T val) noexcept
		{
			if constexpr (expon > 0)
			{
				T retVal
				{	val};

				for (std::size_t p = 1; p < expon; ++p)
				retVal *= val;

				return retVal;
			}
			else if constexpr (expon < 0)
			{
				constexpr int posExpon = -expon;

				T retVal
				{	val};

				for (std::size_t p = 1; p < posExpon; ++p)
				retVal *= val;

				return T(1) / retVal;
			}
			else
			return T(1);
		}
	}
	// namespace schemi

#endif /* INTEXPPOW_HPP_ */
