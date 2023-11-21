/**
 * @file PNGMath.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include "PNGMath.hpp"
#include <iostream>
using namespace Eigen;
using Mat2d = Tensor<double, 2>;
using Dot = std::pair<int, int>;

Mat2d eg::math::conv(Mat2d & input, Mat2d & kernel) {
    int kh = kernel.dimensions()[0];
    int kw = kernel.dimensions()[1];
    int h = input.dimensions()[0];
    int w = input.dimensions()[1];
    if((kh % 2 == 0) || (kw % 2 == 0)) // shape of kernel must odd * odd.
        throw eg::exceptions::InvalidParameter();

    Mat2d out(h, w);
    out.setConstant(0);

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            for(int dx = -kh/2; dx < kh/2 + 1; dx++) {
                if(i + dx < 0 || i + dx >= h)
                    continue;
                for(int dy = -kw/2; dy < kw/2 + 1; dy++) {
                    if(j + dy < 0 || j + dy >= w)
                        continue;
                    out(i, j) += kernel(dx + kh/2, dy + kw/2)*input(i + dx, j + dy);
                }
            }
        }
    }
    return out;
}

double eg::math::rmse(Mat2d & a, Mat2d & b) {
    if(a.dimensions()[0] != b.dimensions()[0] ||
        a.dimensions()[1] != b.dimensions()[1])
            throw exceptions::InvalidParameter();
    Tensor<double, 0> tmp = (a-b).square().sum();
    return std::sqrt(tmp(0));
}

Mat2d eg::math::inflate(Mat2d & a, int h, int w) {
    int ah = a.dimensions()[0], aw = a.dimensions()[1];
    if(ah > h || aw > w)
        throw exceptions::InvalidParameter();
    if(ah == h && aw == w)
        return a;

    Mat2d res(h, w);
    res.setConstant(0);
    for(int i = 0; i < ah; i++) {
        for(int j = 0; j < aw; j++) {
                res(i, j) = a(i, j);
        }
    }

    return res;
}

/**
 * @attention the border of mask must zero. If not, if will raise segfault.
 */
Mat2d eg::math::grassfire(Mat2d & a, Mat2d & mask) {
    if(a.dimensions() != mask.dimensions())
        throw exceptions::InvalidParameter();
    Mat2d ans(a.dimensions());
    ans.setConstant(0);
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(mask(i, j) && a(i, j) > 0)
                ans(i, j) = 1 + std::min(ans(i - 1, j), ans(i, j - 1)); // use taxi distance
        }
    }
    for(int i = h - 1; i >= 0; i--) {
        for(int j = w - 1; j >= 0; j--) {
            if(mask(i, j) && a(i, j) > 0)
                ans(i, j) = std::min(ans(i, j),
                                     1 + std::min(ans(i + 1, j), ans(i, j + 1)));
        }
    }

    return ans;
}

int ccw(const Dot & a, const Dot & b, const Dot & c) {
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

bool isPointInsideHull(const Dot & p, const std::vector<Dot> & hull, bool includeBorder) {
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
std::vector<Dot> eg::math::getConvexHull(std::vector<Dot> & a) {
    std::vector<Dot> rHull;
    std::vector<Dot> lHull;

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

Mat2d eg::math::getMask(Mat2d & a) {
    std::vector<Dot> candidates;
    for(int i = 0; i < a.dimensions()[0]; i++)
        for(int j = 0; j < a.dimensions()[1]; j++)
            if(a(i, j))
                candidates.push_back(std::make_pair(i, j));
    std::vector<Dot> hull = getConvexHull(candidates);

    for(int i = 0; i < hull.size(); i++) {
        std::cout << hull[i].first << " " << hull[i].second << std::endl;
    }
    Mat2d ans = Mat2d(a.dimensions());
    ans.setConstant(0);
    for(int i = 0; i < a.dimensions()[0]; i++) {
        for(int j = 0; j < a.dimensions()[1]; j++) {
            if(isPointInsideHull(std::make_pair(i, j), hull, true))
                ans(i, j) = 1;
        }
    }
    return ans;
}
