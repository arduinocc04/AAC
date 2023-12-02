/**
 * @file egProcessing.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <queue>
#include <algorithm>
#include <iostream>

#include "egProcessing.hpp"
#include "egGeometry.hpp"
#include "egMath.hpp"
#include "egTrace.hpp"
#include "egTool.hpp"

using namespace eg;
using namespace eg::imgproc;
using namespace eg::geo;

Mat2d eg::imgproc::binary(const Mat2d & gray, double threshold) {
    int h = gray.dimensions()[0];
    int w = gray.dimensions()[1];
    Mat2d bin(h, w);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            bin(i, j) = (double)(gray(i, j) > threshold);

    return bin;
}

Mat2d eg::imgproc::getMask(const Mat2d & gray) {
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

Mat2d blurGaussian(const Mat2d & gray) {
    Mat2d kernel(3, 3);
    kernel(0, 0) = (double)1/16, kernel(0, 1) = (double)2/16, kernel(0, 2) = (double)1/16;
    kernel(1, 0) = (double)2/16, kernel(1, 1) = (double)4/16, kernel(1, 2) = (double)2/16;
    kernel(2, 0) = (double)1/16, kernel(2, 1) = (double)2/16, kernel(2, 2) = (double)1/16;
    return eg::math::conv(gray, kernel);
}

Mat2d eg::imgproc::blur(const Mat2d & gray, int method) {
    switch(method) {
        case eg::blurMethod::gaussian:
            return blurGaussian(gray);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d cvtGrayMean(const Image & imageRGBA) {
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

Mat2d eg::imgproc::cvtGray(const Image & imageRGBA, int method) {
    switch(method) {
        case eg::grayCvtMethod::mean:
            return cvtGrayMean(imageRGBA);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d getEdgeGrad(const Mat2d & gray) {
    Eigen::Tensor<double, 2> kernel(3, 3);
    kernel(0, 0) = -1, kernel(0, 1) = -1, kernel(0, 2) = -1;
    kernel(1, 0) = -1, kernel(1, 1) =  8, kernel(1, 2) = -1;
    kernel(2, 0) = -1, kernel(2, 1) = -1, kernel(2, 2) = -1;

    return eg::math::conv(gray, kernel);
}

Mat2d eg::imgproc::getEdge(const Mat2d & gray, int method) {
    switch(method) {
        case eg::edgeDetectMethod::gradient:
            return getEdgeGrad(gray);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d extCenterlineGrassfire(const Mat2d & bin) {
    Mat2d mask = imgproc::getMask(bin);
    return eg::imgproc::grassfire(bin, mask);
}

Mat2d eg::imgproc::extractCenterline(const Mat2d & bin, int method) {
    switch(method) {
        case eg::centerlineMethod::grassfire:
            return extCenterlineGrassfire(bin);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d addZeroPadding(const Mat2d & a, int h, int j, int k, int l) {
    int ah = a.dimensions()[0];
    int aw = a.dimensions()[1];
    Mat2d res(ah + j + k, aw + h + l);
    res.setConstant(0);
    for(int r = 0; r < ah; r++)
        for(int c = h; c < aw; c++)
            res(r + k, c + h) = a(r, c);
    return res;
}

Mat2d eg::imgproc::addPadding(const Mat2d & a, int h, int j, int k, int l, int method) {
    if(h < 0 || j < 0 || k < 0 || l < 0)
        throw eg::exceptions::InvalidParameter();
    switch(method) {
        case eg::paddingMethod::zero:
            return addZeroPadding(a, h, j, k, l);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d eg::imgproc::inflate(const Mat2d & a, int h, int w) {
    int ah = a.dimensions()[0], aw = a.dimensions()[1];
    return imgproc::addPadding(a, 0, h - ah, 0, w - aw,
            eg::paddingMethod::zero);
}

unsigned char saturate(double x) {
    if(x < 0) return 0;
    if(x > 255) return 255;
    return (unsigned char)x;
}

Image eg::imgproc::mat2dToImage(const Mat2d & a) {
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

std::pair<Paths, std::vector<int>> getContourSuzuki(const Mat2d & bin) {
    int h = bin.dimensions()[0];
    int w = bin.dimensions()[1];
    Eigen::Tensor<int, 2> ans(h, w);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            ans(i, j) = bin(i, j) > 0;

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
    std::vector<bool> calced(nbd);
    std::vector<std::vector<bool>> used(h);
    for(int i = 0; i < h; i++)
        used[i].resize(w);
    calced[0] = calced[1] = true;
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            int cnbd = abs(ans(i, j));
            if(!calced[cnbd]) {
                for(int k = 0; k < h; k++)
                    for(int l = 0; l < w; l++)
                        used[k][l] = false;
                std::queue<Dot> q;

                used[i][j] = true;
                calced[cnbd] = true;
                int last;
                Path res;
                res.push_back(std::make_pair(i, j));
                for(int k = 3; k < 7; k++) { // go below or right
                    int tx = i + dx8ccw[k];
                    int ty = j + dy8ccw[k];
                    if(tx < 0 || tx >= h || ty < 0 || ty >= w)
                        continue;
                    if(abs(ans(tx, ty)) == cnbd) {
                        q.push(std::make_pair(tx, ty));
                        last = k;
                        break;
                    }
                }
                while(q.size()) {
                    Dot tmp = q.front();
                    res.push_back(tmp);
                    q.pop();
                    for(int k = 0; k < 8; k++) {
                        int tx = tmp.first + dx8ccw[k];
                        int ty = tmp.second + dy8ccw[k];
                        if(tx < 0 || tx >= h || ty < 0 || ty >= w || used[tx][ty])
                            continue;
                        if(abs(ans(tx, ty)) == cnbd) {
                            q.push(std::make_pair(tx, ty));
                            used[tx][ty] = true;
                            break;
                        }
                    }
                }
                while(q.size()) q.pop();
                bool flag = false;
                for(int k = last + 1; k < 7; k++) {
                    int tx = i + dx8ccw[k];
                    int ty = j + dy8ccw[k];
                    if(tx < 0 || tx >= h || ty < 0 || ty >= w)
                        continue;
                    if(abs(ans(tx, ty)) == cnbd && !used[tx][ty]) {
                        flag = true;
                        q.push(std::make_pair(tx, ty));
                        used[tx][ty] = true;
                        break;
                    }
                }
                if(flag) { // this will change order to clock wise.
                    std::reverse(res.begin(), res.end());
                    while(q.size()) {
                        Dot tmp = q.front();
                        res.push_back(tmp);
                        q.pop();
                        for(int k = 0; k < 8; k++) {
                            int tx = tmp.first + dx8ccw[k];
                            int ty = tmp.second + dy8ccw[k];
                            if(tx < 0 || tx >= h || ty < 0 || ty >= w || used[tx][ty])
                                continue;
                            if(abs(ans(tx, ty)) == cnbd) {
                                q.push(std::make_pair(tx, ty));
                                used[tx][ty] = true;
                                break;
                            }
                        }
                    }
                }
                borders.push_back(res);
            }
        }
    }
    return std::make_pair(borders, p);
}

std::pair<Paths, std::vector<int>> eg::imgproc::getContours(const Mat2d & bin, int method) {
    switch(method) {
        case eg::contourMethod::suzuki:
            return getContourSuzuki(bin);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

Mat2d eg::imgproc::saturate(const Mat2d & gray) {
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

Mat2d eg::imgproc::markOutlier(const Mat2d & gray, double threshold) {
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

Mat2d eg::imgproc::reverse(const Mat2d & bin) {
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

Mat2d eg::imgproc::logpolarForMat2d(const Mat2d & a, const Dot & p) {
    const int h = a.dimensions()[0];
    const int w = a.dimensions()[1];
    const int tbinCnt = 12;
    const int rbinCnt = 5;
    const double tbinSize = 2*M_PI/tbinCnt;
    const double rbinSize = std::log((double)(16/2)/rbinCnt);

    Mat2d histo(tbinCnt, rbinCnt);
    histo.setConstant(0);

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            const double rho = eg::geo::logEuclideDist(std::make_pair(i, j), p);
            if(i == p.first) continue;
            const double theta = std::atan2(j - p.second, i - p.first) + M_PI;
            int ri, ti;
            if(rho > rbinSize*rbinCnt)
                ri = rbinCnt - 1;
            else {
                for(int k = 1; k <= rbinCnt; k++) {
                    if(rho <= k*rbinSize) {
                        ri = k - 1;
                        break;
                    }
                }
            }
            if(theta > tbinSize*tbinCnt)
                ti = tbinCnt - 1;
            else {
                for(int k = 1; k <= tbinCnt; k++) {
                    if(theta <= k*tbinSize) {
                        ti = k - 1;
                        break;
                    }
                }
            }
            histo(ti, ri) += a(i, j);
        }
    }
    return histo;
}

Mat2d eg::imgproc::logpolar(const Dots & dots, const Dot & p) {
    const int tbinCnt = 12; // value in the paper structure-based ascii art
    const int rbinCnt = 5; // value in the paper structure-based ascii art
    const double tbinSize = 2*M_PI/tbinCnt;
    const double rbinSize = std::log((double)(16/2)/rbinCnt); // value in the paper structure-based ascii art

    Mat2d histo(tbinCnt, rbinCnt);
    histo.setConstant(0);
    for(int i = 0; i < dots.size(); i++) {
        double rho = eg::geo::logEuclideDist(dots[i], p);
        if(dots[i].first == p.first) continue;
        double theta = std::atan2(dots[i].second - p.second, dots[i].first - p.first) + M_PI;
        int ri, ti;
        if(rho > rbinSize*rbinCnt)
            ri = rbinCnt - 1;
        else {
            for(int j = 1; j <= rbinCnt; j++) {
                if(rho <= j*rbinSize) {
                    ri = j - 1;
                    break;
                }
            }
        }

        if(theta > tbinSize*tbinCnt)
            ti = tbinCnt - 1;
        else {
            for(int j = 1; j <= tbinCnt; j++) {
                if(theta <= j*tbinSize) {
                    ti = j - 1;
                    break;
                }
            }
        }
        ++histo(ti, ri);
    }
    return histo;
}

Mat2d eg::imgproc::logpolarAll(const Dots & dots) {
    if(dots.size() == 0)
        throw eg::exceptions::InvalidParameter();
    Mat2d histogram = logpolar(dots, dots[0]);
    for(int i = 1; i < dots.size(); i++) {
        histogram += logpolar(dots, dots[i]);
    }
    return histogram;
}

Mat2d eg::imgproc::grassfire(const Mat2d & a, const Mat2d & mask) {
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

void eg::imgproc::drawSegmentDirect(Mat2d & a, const Segment & s, int val) {
    int dx = abs((s.first - s.second).first);
    int sx = (s.first.first < s.second.first)?1:-1;
    int dy = -abs((s.first - s.second).second);
    int sy = (s.first.second < s.second.second)?1:-1;
    int err = dx + dy;
    Dot tmp = s.first;
    while(true) {
        a(tmp.first, tmp.second) = val;
        if(tmp == s.second) break;

        int e2 = 2*err;
        if(e2 >= dy) {
            if(tmp.first == s.second.first) break;
            err = err + dy;
            tmp.first += sx;
        }
        if(e2 <= dx) {
            if(tmp.second == s.second.second) break;
            err = err + dx;
            tmp.second += sy;
        }
    }

}

Mat2d eg::imgproc::drawSegment(const Mat2d & a, const Segment & s, int val) {
    Mat2d ans(a.dimensions());
    for(int i = 0; i < a.dimensions()[0]; i++)
        for(int j = 0; j < a.dimensions()[1]; j++)
            ans(i, j) = a(i, j);
    drawSegmentDirect(ans, s, val);
    return ans;
}

Mat2d eg::imgproc::drawSegments(const Mat2d & a, const Segments & ss, int val) {
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];
    Mat2d ans(h, w);
    for(int i = 0; i < h; i++)
        for(int j = 0; j < w; j++)
            ans(i, j) = a(i, j);

    for(int i = 0; i < ss.size(); i++) {
        drawSegmentDirect(ans, ss[i], val);
    }
    return ans;
}

Mat2d eg::imgproc::erode(const Mat2d & bin, int kh, int kw) {
    int h = bin.dimensions()[0];
    int w = bin.dimensions()[1];
    Mat2d ans(h, w);
    if((kh % 2 == 0) || (kw % 2 == 0))
        throw eg::exceptions::InvalidParameter();
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            bool flag = true;
            for(int dx = -kh/2; dx <= kh/2; dx++) {
                if(i + dx < 0 || i + dx >= h)
                    continue;
                for(int dy = -kw/2; dy <= kw/2; dy++) {
                    if(j + dy < 0 || j + dy >= w)
                        continue;
                    if(bin(i + dx, j + dy) == 0) {
                        flag = false;
                    }
                }
                if(!flag) break;
            }
            if(flag) ans(i, j) = 1;
            else ans(i, j) = 0;
        }
    }
    return ans;
}

Mat2d eg::imgproc::dilate(const Mat2d & bin, int kh, int kw) {
    int h = bin.dimensions()[0];
    int w = bin.dimensions()[1];
    Mat2d ans(h, w);
    if((kh % 2 == 0) || (kw % 2 == 0))
        throw eg::exceptions::InvalidParameter();
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            bool flag = false;
            for(int dx = -kh/2; dx <= kh/2; dx++) {
                if(i + dx < 0 || i + dx >= h)
                    continue;
                for(int dy = -kw/2; dy <= kw/2; dy++) {
                    if(j + dy < 0 || j + dy >= w)
                        continue;
                    if(bin(i + dx, j + dy) > 0) {
                        flag = true;
                        break;
                    }
                }
                if(flag) break;
            }
            if(flag) ans(i, j) = 1;
            else ans(i, j) = 0;
        }
    }
    return ans;
}

Mat2d eg::imgproc::approxUsingSegments(const Mat2d & a) {
    int h = a.dimensions()[0];
    int w = a.dimensions()[1];
    Mat2d ans(h, w);
    ans.setConstant(0);

    auto ttmp = getContours(a, eg::contourMethod::suzuki);
    for(int i = 0; i < ttmp.first.size(); i++) {
        if(!ttmp.first[i].size()) continue;
        Segments tmp = eg::trace::decomposePathToSegments(ttmp.first[i], eg::pathDecomMethod::greedy);
        for(int j = 0; j < tmp.size(); j++)
            drawSegmentDirect(ans, tmp[j], 1);
    }
    return ans;
}

Mat2d resizeVector(const Mat2d & bin, int targetH, int targetW) {
    int h = bin.dimensions()[0];
    int w = bin.dimensions()[1];

    Mat2d ans(targetH, targetW);
    ans.setConstant(0);

    double hRatio = (double)targetH/h;
    double wRatio = (double)targetW/w;

    auto ttmp = getContours(bin, eg::contourMethod::suzuki);
    for(int i = 0; i < ttmp.first.size(); i++) {
        if(!ttmp.first[i].size()) continue;
        Segments tmp = eg::trace::decomposePathToSegments(ttmp.first[i], eg::pathDecomMethod::greedy);
        for(int j = 0; j < tmp.size(); j++) {
            tmp[j].first.first = round(tmp[j].first.first*hRatio);
            tmp[j].second.first = round(tmp[j].second.first*hRatio);
            tmp[j].first.second = round(tmp[j].first.second*wRatio);
            tmp[j].second.second = round(tmp[j].second.second*wRatio);
            drawSegmentDirect(ans, tmp[j], 1);
        }
    }
    return ans;
}

Mat2d eg::imgproc::resizeImage(const Mat2d & bin, int method, int targetH, int targetW) {
    switch(method) {
        case eg::resizeMethod::vector: {
            return resizeVector(bin, targetH, targetW);
        }
        default:
            throw eg::exceptions::InvalidParameter();
    }
}
