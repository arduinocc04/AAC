/**
 * @file egTrace.hpp
 * @author Daniel Cho
 * @date 2023.11.29
 * @version 0.0.1
 */
#include "egProcessing.hpp"

namespace eg::trace {
/**
 * @param a path must given in ccw or cw.
 */
Segments decomposePathToSegments(const Path & a, int method);
}
