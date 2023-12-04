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
Segments approxPathToSegments(const Path & a, int method);

Path approxPath(const Path & a);
}
