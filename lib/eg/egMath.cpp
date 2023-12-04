/**
 * @file egMath.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#include "egGeometry.hpp"
#include "egMath.hpp"
#include "egTool.hpp"

#define ZERO 1e-6

using namespace eg;
using namespace Eigen;

/**
 * @todo padding 붙여서 구현하는게 더 깔끔할듯
 */
Mat2d eg::math::conv(const Mat2d & input, const Mat2d & kernel) {
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

// I wanted to name this function as rmse, but stupid compiler thought it's ambiguous.
double calcrmse(const Mat2d & a, const Mat2d & b) {
    Tensor<double, 0> tmp = (a-b).square().sum();
    return std::sqrt(tmp(0));
}

double bhattacharyyaDist(const Mat2d & ha, const Mat2d & hb) {
    int h = ha.dimensions()[0];
    int w = ha.dimensions()[1];
    Eigen::Tensor<double, 0> sa = ha.sum();
    Eigen::Tensor<double, 0> sb = hb.sum();
    if(sa(0) == 0 || sb(0) == 0) return calcrmse(ha, hb); // todo change..
    double n = std::sqrt(sa(0))*std::sqrt(sb(0));
    double s = 0;
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            s += std::sqrt(ha(i, j))*std::sqrt(hb(i, j));
        }
    }
    return std::sqrt(std::abs(1 - 1/n*s));
}

double calcAE(const Mat2d & ha, const Mat2d & hb) { // absolute error
    int h = ha.dimensions()[0];
    int w = ha.dimensions()[1];
    Mat2d t = ha - hb;
    int res = 0;
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            res += abs(t(i, j));
    return res;
}

double eg::math::compareHistogram(const Mat2d & ha, const Mat2d & hb, int method) {
    if(ha.dimensions() != hb.dimensions())
        throw eg::exceptions::InvalidParameter();

    switch(method) {
        case eg::histCmpMethod::bhattacharyya:
            return bhattacharyyaDist(ha, hb);
        case eg::histCmpMethod::absolute:
            return calcAE(ha, hb);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

double eg::math::compareMat2d(const Mat2d & a, const Mat2d & b, int method) {
    if(a.dimensions() != b.dimensions())
        throw eg::exceptions::InvalidParameter();
    switch(method) {
        case eg::matCmpMethod::rmse:
            return calcrmse(a, b);
        case eg::matCmpMethod::shape: {
            std::pair<Paths, std::vector<int>> tmpA = eg::imgproc::getContours(a, eg::contourMethod::suzuki);
            std::pair<Paths, std::vector<int>> tmpB = eg::imgproc::getContours(b, eg::contourMethod::suzuki);
            Dots dotsConsistA = eg::tool::merge(tmpA.first);
            Dots dotsConsistB = eg::tool::merge(tmpB.first);
            Mat2d histogramA = eg::imgproc::logpolarAll(dotsConsistA);
            Mat2d histogramB = eg::imgproc::logpolarAll(dotsConsistB);
            return compareHistogram(histogramA, histogramB, eg::histCmpMethod::bhattacharyya);
        }
        case eg::matCmpMethod::logpolar: {
            int h = a.dimensions()[0];
            int w = a.dimensions()[1];
            Dots candidates;
            for(int i = 0; i < h; i += 2) // i += 2 is suggested in the paper structure-based ascii art
                for(int j = 0; j < w; j += 2) // j += 2 is suggested in the paper structure-based ascii art
                    candidates.push_back(std::make_pair(i, j));

            double res = 0;
            Eigen::Tensor<double, 0> n = a.sum();
            Eigen::Tensor<double, 0> np = b.sum();
            double M = n(0) + np(0);
            if(M == 0) return 0;

            for(int i = 0; i < candidates.size(); i++) {
                Mat2d histogramA = eg::imgproc::logpolarForMat2d(a, candidates[i]);
                Mat2d histogramB = eg::imgproc::logpolarForMat2d(b, candidates[i]);
                res += calcrmse(histogramA, histogramB);
            }

            return res/M;
        }
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

double calcDeformAngle(const Segment & a, const Segment & b) {
    const double lambda = 8/M_PI;
    const double aLength = eg::geo::euclideDist(a.first, a.second);
    const double bLength = eg::geo::euclideDist(b.first, b.second);
    const double t = eg::geo::dot(a.first - a.second, b.first - b.second)/(aLength*bLength);
    if(t != t) return 987;
    const double theta = std::acos(t);
    // std::cout << t << std::endl;
    return std::exp(lambda*theta);
}

double calcDeformLength(const Segment & a, const Segment & b) {
    const double lambda = (double)2/15;
    const double delta = 0.5;

    const double r = eg::geo::euclideDist(a.first, a.second);
    const double rr = eg::geo::euclideDist(b.first, b.second);
    const double lengthDelta = abs(r - rr);
    const double shortL = std::min(r, rr);
    // if(shortL < ZERO)  return 987;
    const double longL = std::max(r, rr);
    const double t = std::max(std::exp(lambda*lengthDelta),
                     std::exp(delta*longL/shortL));
    if(t != t) return 987;
    return t;
}

double eg::math::calcDeformLocal(const Segment & before, const Segment & after) {
    double t = calcDeformAngle(before, after);
    double tt = calcDeformLength(before, after);
    // std::cout << "ANGLE " << t << " LEN " << tt << std::endl;
    return std::max(t, tt);
}

/**
 * @attention this implementation isn't same in the paper.
 */
double eg::math::calcAccess(const Segment & before, const Segment & after, const Segments & ss) {
    Vec2 u = before.second - before.first;
    Dot mid = before.first + u/2;
    std::vector<std::pair<double, int>> calced;
    for(int i = 0; i < ss.size(); i++) {
        if(ss[i] == before) continue;
        calced.push_back(std::make_pair(eg::geo::distSegDot(ss[i], mid), i));
    }
    std::sort(calced.begin(), calced.end());
    std::vector<std::pair<double, double>> used;
    std::vector<int> candidates;
    for(int j = 0; j < calced.size(); j++) {
        int i = calced[j].second;
        double startAngle, endAngle;
        if(ss[i].first.first == mid.first)
            startAngle = 0;
        else
            startAngle = std::atan2(ss[i].first.second - mid.second, ss[i].first.first - mid.first) + M_PI;
        if(ss[i].second.first == mid.first)
            endAngle = 0;
        else
            endAngle = std::atan2(ss[i].second.second - mid.second, ss[i].second.first - mid.first) + M_PI;
        if(startAngle > endAngle) {
            double tmp = startAngle;
            startAngle = endAngle;
            endAngle = startAngle;
        }

        bool insert = true;
        for(int k = 0; k < used.size(); k++) {
            if(used[k].first <= startAngle && endAngle <= used[k].second) {
                insert = false;
                break;
            }
        }

        if(insert) {
            candidates.push_back(i);
            double l = startAngle;
            double r = endAngle;
            std::vector<int> toRemove;
            for(int k = 0; k < used.size(); k++) { // merge overlapped intervals
                double tl = used[k].first, tr = used[k].second;
                if((tl <= l && l <= tr) || (tl <= r && r <= tr) || (l <= tl && tl <= r) || (l <= tr && tr <= r)) {
                        toRemove.push_back(k);
                        l = std::min(tl, l);
                        r = std::max(tr, r);
                }
            }
            used.push_back(std::make_pair(l, r));
            for(int k = 0; k < toRemove.size(); k++) {
                used.erase(used.begin() + toRemove[k] - k);
            }
        }
    }

    u = after.second - after.first;
    Dot midAfter = after.first + u/2;

    double lengthSum = 0;
    for(int i = 0; i < candidates.size(); i++) {
        lengthSum += eg::geo::distSegDot(ss[i], mid);
    }
    // if(lengthSum < ZERO) return 987;
    double res = 0;
    for(int i = 0; i < candidates.size(); i++) {
        double w = eg::geo::distSegDot(ss[i], mid)/lengthSum;
        Dot mostCloseToMid;
        if(eg::geo::dot(ss[i].second - ss[i].first, mid - ss[i].first) < 0) {
            mostCloseToMid = ss[i].first;
        }
        else if(eg::geo::dot(ss[i].first - ss[i].second, mid - ss[i].second) < 0) {
            mostCloseToMid = ss[i].second;
        }
        else {
            u = ss[i].second - ss[i].first;
            Vec2 v = mid - ss[i].first;
            // if(eg::geo::norm(u) < ZERO) return 987;
            mostCloseToMid = ss[i].first + eg::geo::dot(u, v)/eg::geo::norm(u)*u;
        }
        double t = w*calcDeformLocal(std::make_pair(mid, mostCloseToMid), std::make_pair(midAfter, mostCloseToMid));
        res += t;
    }
    if(res != res) return 987;
    return res;
}

double eg::math::calcDeform(const Segment & before, const Segment & after, const Segments & ss) {
    if(eg::geo::euclideDist(after.first, after.second) < 1e-8) return  987;
    if(eg::geo::euclideDist(before.first, before.second) < 1e-8) return  987;
    double t = calcDeformLocal(before, after);
    double tt = calcAccess(before, after, ss);
    // std::cout << "LO " << t << " ACC " << tt << std::endl;
    double c = std::max(t, tt);
    if(c != c) return 987;
    return c;
}
