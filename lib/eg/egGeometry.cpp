/**
 * @file egGeometry.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include "egTypes.hpp"
#include "egGeometry.hpp"

using namespace eg;

int eg::geo::ccw(const Dot & a, const Dot & b, const Dot & c) {
    int x1 = a.first;
    int y1 = a.second;
    int x2 = b.first;
    int y2 = b.second;
    int x3 = c.first;
    int y3 = c.second;

    int  t = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
    if(t > 0) return 1;
    else if(t < 0) return -1;
    return 0;
}

bool eg::geo::isPointInsideHull(const Dot & p, const std::vector<Dot> & hull, bool includeBorder) {
    int prev = ccw(p, hull[0], hull[1]);
    for(int i = 1; i < hull.size() -1; i++) {
        int t = ccw(p, hull[i], hull[i + 1]);
        if(t*prev < 0)
            return false;
        else if(!includeBorder && t == 0)
            return false;
        prev = t;
    }
    return true;
}

/**
 * @brief calculate convex hull using monotone chain algorithm.
 * @param a x-sorted <int, int> pair vector.
 */
Dots eg::geo::getConvexHull(Dots & a) {
    Dots rHull;
    Dots lHull;

    for(int i = 0; i < a.size(); i++) {
        Dot p = a[i];
        while(rHull.size() >= 2 && ccw(rHull[rHull.size() - 2], rHull[rHull.size() - 1], p) <= 0)
            rHull.pop_back();
        rHull.push_back(p);
    }

    for(int i = a.size() - 1; i >= 0; i--) {
        int n = lHull.size();
        Dot p = a[i];
        while(lHull.size() >= 2 && ccw(lHull[lHull.size() - 2], lHull[lHull.size() - 1], p) <= 0)
            lHull.pop_back();
        lHull.push_back(p);
    }

    std::vector<Dot> hull;
    for(int i = 0; i < rHull.size() - 1; i++)
        hull.push_back(rHull[i]);
    for(int i = 0; i < lHull.size() - 1; i++)
        hull.push_back(lHull[i]);

    return hull;
}

