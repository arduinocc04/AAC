/**
 * @file egProcessing.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <queue>
#include <iostream>

#include "egProcessing.hpp"
#include "egTypes.hpp"
#include "egGeometry.hpp"
#include "egMath.hpp"

using namespace eg;
using namespace eg::imgproc;
using namespace eg::geo;

Mat2d eg::imgproc::binary(Mat2d & gray, double threshold) {
    int h = gray.dimensions()[0];
    int w = gray.dimensions()[1];
    Mat2d bin(h, w);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            bin(i, j) = (double)(gray(i, j) > threshold);

    return bin;
}

Mat2d eg::imgproc::getMask(Mat2d & gray) {
    int h = gray.dimensions()[0];
    int w = gray.dimensions()[1];
    Dots candidates;
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            if(gray(i, j))
                candidates.push_back(std::make_pair(i, j));
    Dots hull = getConvexHull(candidates);

    Mat2d ans = Mat2d(h, w);
    ans.setConstant(0);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(isPointInsideHull(std::make_pair(i, j), hull, true))
                ans(i, j) = 1;
        }
    }
    return ans;
}

Mat2d blurGaussian(Mat2d & gray) {
    Mat2d kernel(3, 3);
    kernel(0, 0) = (double)1/16, kernel(0, 1) = (double)2/16, kernel(0, 2) = (double)1/16;
    kernel(1, 0) = (double)2/16, kernel(1, 1) = (double)4/16, kernel(1, 2) = (double)2/16;
    kernel(2, 0) = (double)1/16, kernel(2, 1) = (double)2/16, kernel(2, 2) = (double)1/16;
    return eg::math::conv(gray, kernel);
}

Mat2d eg::imgproc::blur(Mat2d & gray, int method) {
    switch(method) {
        case eg::blurMethod::gaussian:
            return blurGaussian(gray);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d cvtGrayMean(Image & imageRGBA) {
    if(!imageRGBA.data())
        throw eg::exceptions::ImageNotOpened();

    int h = imageRGBA.dimensions()[0];
    int w = imageRGBA.dimensions()[1];
    Mat2d gray(h, w);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            double mean = 0;
            for(int k = 0; k < 3; k++)
                mean += imageRGBA(i, j, k);
            mean /= 3;
            gray(i, j) = mean;
        }
    }
    return gray;
}

Mat2d eg::imgproc::cvtGray(Image & imageRGBA, int method) {
    switch(method) {
        case eg::grayCvtMethod::mean:
            return cvtGrayMean(imageRGBA);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d getEdgeGrad(Mat2d & gray) {
    Eigen::Tensor<double, 2> kernel(3, 3);
    kernel(0, 0) = -1, kernel(0, 1) = -1, kernel(0, 2) = -1;
    kernel(1, 0) = -1, kernel(1, 1) =  8, kernel(1, 2) = -1;
    kernel(2, 0) = -1, kernel(2, 1) = -1, kernel(2, 2) = -1;

    return eg::math::conv(gray, kernel);
}

Mat2d eg::imgproc::getEdge(Mat2d & gray, int method) {
    switch(method) {
        case eg::edgeDetectMethod::gradient:
            return getEdgeGrad(gray);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d extCenterlineGrassfire(Mat2d & bin) {
    Mat2d mask = imgproc::getMask(bin);
    return eg::imgproc::grassfire(bin, mask);
}

Mat2d eg::imgproc::extractCenterline(Mat2d & bin, int method) {
    switch(method) {
        case eg::centerlineMethod::grassfire:
            return extCenterlineGrassfire(bin);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d addZeroPadding(Mat2d & a, int h, int j, int k, int l) {
    int ah = a.dimensions()[0];
    int aw = a.dimensions()[1];
    Mat2d res(ah + j + k, aw + h + l);
    res.setConstant(0);
    for(int r = 0; r < ah; r++)
        for(int c = h; c < aw; c++)
            res(r + k, c + h) = a(r, c);
    return res;
}

Mat2d eg::imgproc::addPadding(Mat2d & a, int h, int j, int k, int l, int method) {
    if(h < 0 || j < 0 || k < 0 || l < 0)
        throw eg::exceptions::InvalidParameter();
    switch(method) {
        case eg::paddingMethod::zero:
            return addZeroPadding(a, h, j, k, l);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d eg::imgproc::inflate(Mat2d & a, int h, int w) {
    int ah = a.dimensions()[0], aw = a.dimensions()[1];
    return imgproc::addPadding(a, 0, h - ah, 0, w - aw,
                     eg::paddingMethod::zero);
}

unsigned char saturate(double x) {
    if(x < 0) return 0;
    if(x > 255) return 255;
    return (unsigned char)x;
}

Image eg::imgproc::mat2dToImage(Mat2d & a) {
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];
    Image res(h, w, 4);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++) {
            for(int k = 0; k < 3; k++)
                res(i, j, k) = std::round(a(i, j));
            res(i, j, 3) = 255;
        }
    return res;
}

std::pair<Paths, std::vector<int>> getContourSuzuki(Mat2d & bin) {
    int h = bin.dimensions()[0];
    int w = bin.dimensions()[1];
    Mat2d ans(h, w);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            ans(i, j) = (double)(bin(i, j) > 0);

    int nbd = 1;
    int lnbd;
    int dir;
    int dx8ccw[] = {-1, -1, 0, 1, 1, 1, 0, -1};
    int dy8ccw[] = {0, -1, -1, -1, 0, 1, 1, 1};
    int dx8cw[] = {-1, -1, 0, 1, 1, 1, 0, -1};
    int dy8cw[] = {0, 1, 1, 1, 0, -1, -1, -1};
    std::vector<Path> borders;
    borders.push_back({});
    borders.push_back({});

    std::vector<int> p;
    p.push_back(1);
    p.push_back(1);
    std::vector<bool> is_hole;
    is_hole.push_back(true);
    is_hole.push_back(true);

    for(int i = 0; i < h; i++) {
        lnbd = 1; // 1 means frame
        for(int j = 0; j < w; j++) {
            if(ans(i, j) == 0) continue;
            bool flag = true;
            if((j == 0 || ans(i, j - 1) == 0) && ans(i, j) == 1) { // outer
                nbd++;
                if(nbd >= p.size()) {
                    if(is_hole[lnbd])
                        p.push_back(lnbd);
                    else
                        p.push_back(p[lnbd]);
                }
                if(nbd >= is_hole.size()) is_hole.push_back(false);

                dir = 6;
            }
            else if(ans(i, j) >= 1 && (j == w - 1 || ans(i, j + 1) == 0)) { // hole
                nbd++;
                if(ans(i, j) > 1) lnbd = ans(i, j);
                if(nbd >= p.size()) {
                    if(is_hole[lnbd])
                        p.push_back(p[lnbd]);
                    else
                        p.push_back(lnbd);
                }
                if(nbd >= is_hole.size()) is_hole.push_back(true);
                dir = 3;
            }
            else flag = false;

            if(flag) {
                bool found = false;
                int i1, j1;
                for(int k = 0; k < 8; k++) {
                    int tx = i + dx8cw[(dir + k) % 8];
                    int ty = j + dy8cw[(dir + k) % 8];
                    if(tx < 0 || tx >= h || ty < 0 || ty >= w)
                        continue;
                    if(ans(tx, ty) != 0) {
                        i1 = tx, j1 = ty;
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    ans(i, j) = -nbd;
                    break;
                }

                int i2, j2, i3, j3;
                i2 = i1, j2 = j1;
                i3 = i, j3 = j;

                while(true) {
                    for(int k = 0; k < 8; k++) {
                        if(dx8ccw[k] == i2 - i3 && dy8ccw[k] == j2 - j3) {
                            dir = k;
                            break;
                        }
                    }
                    int i4, j4;
                    bool examined = false;
                    for(int k = 1; k <= 8; k++) {
                        int tx = i3 + dx8ccw[(dir + k) % 8];
                        int ty = j3 + dy8ccw[(dir + k) % 8];
                        if((dir + k) % 8 == 6)
                            examined = true;
                        if(tx < 0 || tx >= h || ty < 0 || ty >= w)
                            continue;
                        if(ans(tx, ty)) {
                            i4 = tx, j4 = ty;
                            break;
                        }
                    }

                    bool rcond = (j3 == w - 1 || ans(i3, j3 + 1) == 0);

                    if(rcond && examined) {
                        ans(i3, j3) = -nbd;
                    }
                    else if(!(rcond && examined) && ans(i3, j3) == 1) {
                        ans(i3, j3) = nbd;
                    }

                    if(i4 == i && j4 == j && i3 == i1 && j3 == j1)
                        break;
                    i2 = i3, j2 = j3;
                    i3 = i4, j3 = j4;
                } // while(true)
            } // if(flag)
            if(ans(i, j) != 1) lnbd = abs(ans(i, j));
        } // for j
    } // for i
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(ans(i, j) != 0) {
                int t = abs(ans(i, j));
                if(t >= borders.size()) borders.push_back({std::make_pair(i, j)});
                else borders[t].push_back(std::make_pair(i, j));
            }
        }
    }
    return std::make_pair(borders, p);
}

std::pair<Paths, std::vector<int>> eg::imgproc::getContours(Mat2d & bin, int method) {
    switch(method) {
        case eg::contourMethod::suzuki:
            return getContourSuzuki(bin);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d eg::imgproc::saturate(Mat2d & gray) {
    int h = gray.dimensions()[0];
    int w = gray.dimensions()[1];
    Mat2d ans(h, w);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(gray(i, j) < 0)
                ans(i, j) = 0;
            else if(gray(i, j) > 255)
                ans(i, j) = 255;
            else
                ans(i, j) = gray(i, j);
        }
    }
    return ans;
}

Mat2d eg::imgproc::markOutlier(Mat2d & gray, double threshold) {
    int h = gray.dimensions()[0];
    int w = gray.dimensions()[1];
    Mat2d ans(h, w);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(abs(gray(i, j)) > threshold)
                ans(i, j) = 1;
            else
                ans(i, j) = 0;
        }
    }
    return ans;
}

Mat2d eg::imgproc::reverse(Mat2d & bin) {
    int h = bin.dimensions()[0];
    int w = bin.dimensions()[1];
    Mat2d ans(h, w);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(bin(i, j))
                ans(i, j) = 0;
            else
                ans(i, j) = 1;
        }
    }
    return ans;
}

double logEuclideDist(Dot & a, Dot & b) {
    int distSquared = (a.first - b.first)*(a.first - b.first) + (a.second - b.second)*(a.second - b.second);
    return std::log(std::sqrt(distSquared));
}

Mat2d eg::imgproc::logpolar(Dots & dots, Dot & p) {
    const int tbinCnt = 20;
    const int rbinCnt = 20;
    const double tbinSize = 2*M_PI/tbinCnt;
    const double rbinSize = 0.2;

    Mat2d histogram(tbinCnt, rbinCnt);
    histogram.setConstant(0);
    for(int i = 0; i < dots.size(); i++) {
        double rho = logEuclideDist(dots[i], p);
        if(dots[i].first == p.first) continue;
        double theta = std::atan2(dots[i].second - p.second, dots[i].first - p.first) + M_PI;
        int ri, ti;
        if(rho > rbinSize*rbinCnt)
            ri = rbinCnt - 1;
        else {
            for(int i = 1; i <= rbinCnt; i++) {
                if(rho <= i*rbinSize) {
                    ri = i - 1;
                    break;
                }
            }
        }
        if(theta > tbinSize*tbinCnt)
            ti = tbinCnt - 1;
        else {
            for(int i = 1; i <= tbinCnt; i++) {
                if(theta <= i*tbinSize) {
                    ti = i - 1;
                    break;
                }
            }
        }
        histogram(ri, ti)++;
    }
    return histogram;
}

Mat2d eg::imgproc::logpolarAll(Dots & dots) {
    const int tbinCnt = 20;
    const int rbinCnt = 20;
    Mat2d histogram(tbinCnt, rbinCnt);
    histogram.setConstant(0);
    for(int i = 0; i < dots.size(); i++) {
        histogram += logpolar(dots, dots[i]);
    }
    return histogram;
}

Mat2d eg::imgproc::grassfire(Mat2d & a, Mat2d & mask) {
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

Mat2d eg::imgproc::cycle(const Mat2d & a, int stride) {
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];
    Mat2d ans(h, w);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w - stride; j++)
            ans(i, j + stride) = a(i, j);
        for(int j = w - stride; j < w; j++)
            ans(i, j - (w - stride)) = a(i, j);
    }
    return ans;
}
