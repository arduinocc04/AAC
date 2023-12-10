/**
 * @file access.cpp
 * @author Daniel Cho
 * @date 2023.12.10
 * @version 0.0.1
 */
#include <iostream>
#include <ctime>
#include <random>

#include "egLoader.hpp"
#include "egMath.hpp"
#include "egProcessing.hpp"
#include "egTrace.hpp"
#include "egTool.hpp"
#include "egGeometry.hpp"

#define INF 987654321

using namespace eg::imgproc;

template<typename T>
using ArrayMat = std::vector<std::vector<std::vector<T>>>;

using Index2d = std::pair<int, int>;

using SegWithOri = std::pair<Segment, Index2d>; // Segment with its origin location.

int main(int argc, char * argv[]) {
    if(argc != 2) {
        std::cout << "Use Program Properly! program INPUT_IMAGE_PATH" << std::endl;
        return -1;
    }
    std::string inputImagePath = argv[1];

    const int asciih = 22;
    const int asciiw = 16;

    time_t seed = time(nullptr);
    std::cout << "seed: " << seed << std::endl;
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> distribution(0, INF);
    auto mt = std::bind(distribution, engine);

    eg::PNG inputImage;
    std::cout << "Opening input image" << std::endl;
    inputImage.openImage(inputImagePath);
    Image i = inputImage.copy();
    const int h = i.dimensions()[0];
    const int w = i.dimensions()[1];
    std::cout << std::endl;
    std::cout << "Converting input image gray" << std::endl;
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    std::cout << "Getting Edge of input image" << std::endl;
    t = getEdge(t, eg::edgeDetectMethod::gradient);

    t = markOutlier(t, 50);
    t = dilate(t, 1, 3);
    t = erode(t, 1, 3);
    t = dilate(t, 3, 1);
    t = erode(t, 3, 1);
    i = mat2dToImage(t);

    std::cout << "getContours" << std::endl;
    auto contours = getContours(t, eg::contourMethod::suzuki);

    Segments segs;
    std::cout << "Merging contours" << std::endl;
    Paths paths;
    for(int i = 0; i < contours.first.size(); i++) {
        if(contours.first[i].size() == 0) continue;
        Path tmp = eg::trace::approxPath(contours.first[i]);
        paths.push_back(tmp);
        for(int j = 1; j < tmp.size(); j++)
            segs.push_back(std::make_pair(tmp[j - 1], tmp[j]));
    }
    Segments originalSegs = segs;
    Paths originalPaths = paths;

    inputImage.divideImageByLength(asciih, asciiw);
    int gcCnt = inputImage.getGridColCnt();
    int grCnt = inputImage.getGridRowCnt();

    std::cout << grCnt << "x" << gcCnt << std::endl;

    char devnull[10];

    std::cout << "segs size: " << segs.size() << std::endl;

    std::cout << "Decomposing lines into grids.." << std::endl;
    ArrayMat<SegWithOri> linesPerGrid(grCnt);
    for(int i = 0; i < grCnt; i++)
        linesPerGrid[i].resize(gcCnt);

    for(int k = 0; k < paths.size(); k++) {
        for(int l = 1; l < paths[k].size(); l++) {
            for(int i = 0; i < grCnt; i++) {
                for(int j = 0; j < gcCnt; j++) {
                    const Dot lu = std::make_pair(i*asciih, j*asciiw);
                    const Dot rd = std::make_pair((i + 1)*asciih, (j + 1)*asciiw);
                    const Segment original = std::make_pair(paths[k][l - 1], paths[k][l]);
                    const Segment clipped = eg::tool::clip(original, lu, rd);
                    if(clipped.first == clipped.second) continue;
                    if(clipped.first.first == -1) continue;

                    linesPerGrid[i][j].push_back(std::make_pair(clipped, std::make_pair(k, l)));
                }
            }
        }
    }

    ArrayMat<SegWithOri> originalLinesPerGrid = linesPerGrid;

    Mat2d tmpMat2d(h, w);
    tmpMat2d.setConstant(0);
    Image tmpImg(h, w, 4);
    tmpImg.setConstant(0);
    /// grid
    for(int i = 0; i < grCnt; i++) {
        if((i + 1)*asciih > h) break;
        for(int j = 0; j < w; j++) {
            for(int k = 0; k < 3; k++)
                tmpImg((i + 1)*asciih, j, k) = 127;
        }
    }
    for(int i = 0; i < gcCnt; i++) {
        if((i + 1)*asciiw > w) break;
        for(int j = 0; j < h; j++) {
            for(int k = 0; k < 3; k++)
                tmpImg(j, (i + 1)*asciiw, k) = 127;
        }
    }
    /// grid
    
    /// line
    tmpMat2d.setConstant(0);
    for(int i = 0; i < linesPerGrid.size(); i++) {
        for(int j = 0; j < linesPerGrid[i].size(); j++) {
            for(int k = 0; k < linesPerGrid[i][j].size(); k++)
                drawSegmentDirect(tmpMat2d, linesPerGrid[i][j][k].first, 127);
        }
    }
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            tmpImg(i, j, 3) = 255;
            for(int k = 0; k < 3; k++) {
                if(tmpMat2d(i, j) > 0) tmpImg(i, j, k) = tmpMat2d(i, j);
            }
        }
    }
    /// line
    Segment before;
    bool flag = true;
    do{
        flag = false;
        int i = mt() % linesPerGrid.size();
        if(linesPerGrid[i].size() == 0) {
            flag = true;
            continue;
        }
        int j = mt() % linesPerGrid[i].size();
        if(linesPerGrid[i][j].size() == 0) {
            flag = true;
            continue;
        }
        int k = mt() % linesPerGrid[i][j].size();
        before = linesPerGrid[i][j][k].first;
    }while(flag);

    Vec2 u = before.second - before.first;
    Dot mid = before.first + u/2;

    std::vector<std::pair<double, int>> calced;
    for(int i = 0; i < originalSegs.size(); i++) {
        double dist = eg::geo::distSegDot(originalSegs[i], mid);
        if(dist < 0.5) continue;
        calced.push_back(std::make_pair(dist, i));
    }

    std::sort(calced.begin(), calced.end());
    std::vector<std::pair<double, double>> used;
    std::vector<int> candidates;
    bool * usedBool = new bool[originalSegs.size()];
    for(int i = 0; i < originalSegs.size(); i++) usedBool[i] = false;

    for(double theta = 0; theta < 2*M_PI; theta += 0.01) {
        double minVal = INF;
        int minIndex = -1;
        const Vec2 direc = std::make_pair(5000*std::cos(theta), 5000*std::sin(theta));
        const Segment ray = std::make_pair(mid, mid + direc);
        std::cout << "theta" << theta << std::endl;
        for(int i = 0; i < originalSegs.size(); i++) {
            const auto & tmp = originalSegs[i];
            if(ray.first.first < std::min(tmp.first.first, tmp.second.first) && ray.second.first < std::min(tmp.first.first, tmp.second.first)) continue;
            if(ray.first.first > std::max(tmp.first.first, tmp.second.first) && ray.second.first > std::max(tmp.first.first, tmp.second.first)) continue;
            if(eg::geo::ccw(ray.first, ray.second, tmp.first)*eg::geo::ccw(ray.first, ray.second, tmp.second) <= 0 &&
               eg::geo::ccw(tmp.first, tmp.second, ray.first)*eg::geo::ccw(tmp.first, tmp.second, ray.second) <= 0) {
                double dist = eg::geo::distSegDot(tmp, mid);
                std::cout << "DIST" << dist << std::endl;
                if(dist < 1) continue;
                if(minVal > dist) {
                    minVal = dist;
                    minIndex = i;
                }
            }
        }
        if(minIndex != -1 && !usedBool[minIndex]) {
            candidates.push_back(minIndex);
            usedBool[minIndex] = true;
        }
    }
    delete[] usedBool;

    // for(int j = 0; j < calced.size(); j++) {
    //     int i = calced[j].second;
    //     double startAngle, endAngle;
    //     startAngle = std::atan2(originalSegs[i].first.second - mid.second, originalSegs[i].first.first - mid.first) + M_PI;
    //     endAngle = std::atan2(originalSegs[i].second.second - mid.second, originalSegs[i].second.first - mid.first) + M_PI;
    //     if(startAngle > endAngle) endAngle += 2*M_PI;

    //     bool insert = true;
    //     for(int k = 0; k < used.size(); k++) {
    //         if(used[i].first <= startAngle && endAngle <= used[k].second) {
    //             insert = false;
    //             break;
    //         }
    //     }

    //     std::vector<int> toRemove;
    //     if(insert) {
    //         candidates.push_back(i);

    //         double l = startAngle;
    //         double r = endAngle;
    //         for(int k = 0; k < used.size(); k++) { // merge overlapped intervals
    //             double tl = used[k].first, tr = used[k].second;
    //             if((tl <= l && l <= tr) || (tl <= r && r <= tr) || (l <= tl && tl <= r) || (l <= tr && tr <= r)) {
    //                     toRemove.push_back(k);
    //                     l = std::min(tl, l);
    //                     r = std::max(tr, r);
    //             }
    //         }
    //         used.push_back(std::make_pair(l, r));
    //         for(int k = 0; k < toRemove.size(); k++) {
    //             used.erase(used.begin() + toRemove[k] - k);
    //         }
    //     }
    // }

    /// calced -access -line
    tmpMat2d.setConstant(0);
    for(int i = 0; i < candidates.size(); i++) {
        drawSegmentDirect(tmpMat2d, originalSegs[candidates[i]], 255);
    }
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(tmpMat2d(i, j) > 0) {
                tmpImg(i, j, 0) = 0;
                tmpImg(i, j, 1) = 255;
                tmpImg(i, j, 2) = 0;
            }
        }
    }
    /// calced -access -line

    
    tmpMat2d.setConstant(0);
    for(int j = 0; j < candidates.size(); j++) {
        int i = candidates[j];
        std::cout << eg::geo::distSegDot(originalSegs[i], mid) << std::endl;
        Dot mostCloseToMid;
		Dot mostCloseToMidAfter;
        if(eg::geo::dot(originalSegs[i].second - originalSegs[i].first, mid - originalSegs[i].first) < 0) {
            // std::cout << "TYPE 0" << std::endl;
            mostCloseToMid = originalSegs[i].first;
        }
        else if(eg::geo::dot(originalSegs[i].first - originalSegs[i].second, mid - originalSegs[i].second) < 0) {
            // std::cout << "TYPE 1" << std::endl;
            mostCloseToMid = originalSegs[i].second;
        }
        else {
            // std::cout << "TYPE 2" << std::endl;
            u = originalSegs[i].second - originalSegs[i].first;
            Vec2 v = mid - originalSegs[i].first;
            /*
            std::cout << "u " << u.first << " " << u.second << std::endl;
            std::cout << "v " << v.first << " " << v.second << std::endl;
            std::cout << "dot " << eg::geo::dot(u, v) << std::endl;
            std::cout << "norm " << eg::geo::norm(u) << std::endl;
            std::cout << "originalifirst " << originalSegs[i].first.first << " " << originalSegs[i].first.second << std::endl;
            */
            // if(eg::geo::norm(u) < ZERO) return 987;
            mostCloseToMid = originalSegs[i].first + eg::geo::dot(u, v)/eg::geo::dot(u, u)*u;
        }
        drawSegmentDirect(tmpMat2d, std::make_pair(mid, mostCloseToMid), 255);
    }
    /// mid-line
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(tmpMat2d(i, j) > 0) {
                tmpImg(i, j, 0) = 0;
                tmpImg(i, j, 1) = 0;
                tmpImg(i, j, 2) = 255;
            }
        }
    }
    /// mid-line

    /// vertices
    for(int i = 0; i < paths.size(); i++) {
        for(int j = 0; j < paths[i].size(); j++) {
            int x = paths[i][j].first;
            int y = paths[i][j].second;
            for(int k = 0; k < 3; k++)
                tmpImg(x, y, k) = 255;
        }
    }
    /// vertices


    /// before(one segment)
    tmpMat2d.setConstant(0);
    drawSegmentDirect(tmpMat2d, before, 255);
    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            if(tmpMat2d(i, j) > 0) {
                tmpImg(i, j, 0) = 255;
                tmpImg(i, j, 1) = 0;
                tmpImg(i, j, 2) = 0;
            }
        }
    }
    /// before

    /// mid
    tmpImg(mid.first, mid.second, 0) = 255;
    tmpImg(mid.first, mid.second, 1) = 255;
    tmpImg(mid.first, mid.second, 2) = 255;
    /// mid

    PNG tmpPNG;
    tmpPNG.openImage(inputImagePath);
    tmpPNG.setImage(tmpImg);
    tmpPNG.saveImage("access-" + inputImagePath);
}