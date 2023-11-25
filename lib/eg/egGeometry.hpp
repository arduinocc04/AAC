/**
 * @file egGeometry.hpp
 * @author Daniel Cho
 * @date 2023.11.25
 * @version 0.0.1
 */
#include "egTypes.hpp"

using namespace eg;

namespace eg::geo {

int ccw(const Dot & a, const Dot & b, const Dot & c);

Dots getConvexHull(Dots & a);

bool isPointInsideHull(const Dot & p, const std::vector<Dot> & hull, bool includeBorder);
}
