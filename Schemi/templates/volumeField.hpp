/*
 * volumeField.hpp
 *
 *  Created on: 2019/11/23
 *      Author: Maxim Boldyrev
 *
 *      Alias template for volume field.
 */

#ifndef VOLUMEFIELD_HPP_
#define VOLUMEFIELD_HPP_

#include "cubicCell.hpp"
#include "field.hpp"

namespace schemi
{
template<typename typeOfValue> using volumeField = field<typeOfValue, cubicCell>;
}  // namespace schemi

#endif /* VOLUMEFIELD_HPP_ */
