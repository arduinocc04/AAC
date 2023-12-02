/**
 * @file egGeometry.hpp
 * @author Daniel Cho
 * @date 2023.11.25
 * @version 0.0.1
 */
#ifndef __EGTYPES_H
#include "egTypes.hpp"
#endif
#define __EGTYPES_H

using namespace eg;

namespace eg::geo {

int ccw(const Dot & a, const Dot & b, const Dot & c);

double dot(const Vec2 & a, const Vec2 & b);

/**
 * @param a 2d vector
 * @param b 2d vector
 */
double cross(const Vec2 & a, const Vec2 & b);

double euclideDist(const Dot & a, const Dot & b);

double logEuclideDist(const Dot & a, const Dot & b);

Dots getConvexHull(const Dots & a);

bool isPointInsideHull(const Dot & p, const std::vector<Dot> & hull, bool includeBorder);

double distSegDot(const Segment & a, const Dot & p);
}
