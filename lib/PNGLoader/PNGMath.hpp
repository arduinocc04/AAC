#ifndef __PNGEXCEPTIONS_H
#include "PNGExceptions.hpp"
#endif
#define __PNGEXCEPTIONS_H

#ifndef __Tensor
#include <unsupported/Eigen/CXX11/Tensor>
#endif
#define __Tensor

namespace eg::math {
/**
 * @todo Study other's implementation. I don't know how to behave at corners.
 */
Eigen::Tensor<double, 2> conv(
                    Eigen::Tensor<double, 2> & input,
                    Eigen::Tensor<double, 2> & kernel
                    );

}
