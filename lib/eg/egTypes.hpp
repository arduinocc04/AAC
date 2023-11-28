/**
 * @file egTypes.hpp
 * @author Daniel Cho
 * @date 2023.11.24
 * @version 0.0.1
 */

#include <vector>

#ifndef __TENSOR
#include "unsupported/Eigen/CXX11/Tensor"
#endif
#define __TENSOR

namespace eg {
using Image = Eigen::Tensor<unsigned char, 3>;
using Mat2d = Eigen::Tensor<double, 2>;

using Dot = std::pair<int, int>;
using Dots = std::vector<Dot>;
using Path = Dots;
using Paths = std::vector<Path>;
}
