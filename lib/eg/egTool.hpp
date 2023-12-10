/**
 * @file egTools.hpp
 * @author Daniel Cho
 * @date 2023.12.1
 * @version 0.0.1
 */
#ifndef __EGTYPES_H
#include "egTypes.hpp"
#endif
#define __EGTYPES_H

using namespace eg;

namespace eg::tool {
Dots merge(const Paths & a);

/**
 * @brief get borders aligned by cw from suzuki-ed matrix
 * @attention ans is int-valued matrix. not double!
 */
Paths Mat2dToBorders(const Eigen::Tensor<int, 2> & ans, int nbd, int h, int w);

/**
 * @brief do line clipping using Cohen-Sutherland algorithm
 * @see https://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland_algorithm
 * @attention rd is not included. [lu, rd)
 */
Segment clip(const Segment & s, const Dot & lu, const Dot & rd);
}
