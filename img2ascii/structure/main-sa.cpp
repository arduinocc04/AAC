/**
 * @file main.cpp
 * @author Daniel Cho
 * @date 2023.12.3
 * @version 0.0.1
 */
#include <iostream>
#include <filesystem>
#include <ctime>
#include <random>
#include <thread>
#include <vector>
#include <set>

#include "egLoader.hpp"
#include "egProcessing.hpp"
#include "egMath.hpp"
#include "egGeometry.hpp"
#include "egTrace.hpp"
#include "egTool.hpp"

#define INF 987654321
#define THREAD_CNT 4

#define DEBUG_WITH_TMP_IMAGE
#define DEBUG_WITH_DECOM_IMAGE
#define DEBUG_WITH_NODECOM_IMAGE
#define DEBUG_DRAW_GRID

using namespace eg::imgproc;
namespace fs = std::filesystem;

std::string * names;

template<typename T>
using ArrayMat = std::vector<std::vector<std::vector<T>>>;

using Index2d = std::pair<int, int>;

using SegWithOri = std::pair<Segment, Index2d>; // Segment with its origin location.

int getFileCount(std::string path) {
    int fCnt = 0;
    for(const auto & entry : fs::directory_iterator(path))
        fCnt++;
    return fCnt;
}

void print(Eigen::Tensor<double, 2> & a) {
    for(int i = 0; i < a.dimensions()[0]; i++) {
        for(int j = 0; j < a.dimensions()[1]; j++) {
            std::cout << (a(i, j) > 0);
        }
        std::cout << std::endl;
    }
}

Mat2d * getAllImages(std::string path, int fCnt) {
    Mat2d * ans = new Mat2d[fCnt];

    int i = 0;
    eg::PNG png;
    for(const auto & entry : fs::directory_iterator(path)) {
        png.openImage(entry.path());
        names[i] = entry.path();
        Image t = png.copy();
        ans[i] = cvtGray(t, eg::grayCvtMethod::mean);
        ans[i] = markOutlier(ans[i], 70); // 70 is founded by experiment.
        ans[i] = eg::imgproc::grassfire(ans[i], ans[i]);
        ans[i] = binary(ans[i], 1);
        i += 1;
    }
    std::cout << std::endl;

    return ans;
}

std::string getAsciiFromPath(std::string path) {
    std::string ans = path.substr(path.find_last_of("/\\") + 1);
    std::string nameWithoutExt = ans.substr(0, ans.find_last_of('.'));
    if(nameWithoutExt[0] == '=') {
        char s[2];
        s[0] = std::stoi(nameWithoutExt.substr(1), nullptr, 16);
        s[1] = '\0';
        return std::string(s);
    }
    return nameWithoutExt;
}

void fillDist(double * distances, int iStart, int iEnd, Mat2d & tmp, Mat2d * asciiPNGs) {
    for(int i = iStart; i < iEnd; i++) {
        distances[i] = eg::math::compareMat2d(asciiPNGs[i], tmp, eg::matCmpMethod::logpolar);
    }
}

bool isDotInAABB(const Dot & p, const Dot & lu, const Dot & rd) { // Don't include border
    return lu.first < p.first && p.first < rd.first && lu.second < p.second && p.second < rd.second;
}

std::pair<int, int> getChanged(const std::pair<int, int> & coord, const std::pair<int, int> & lineIdx, const Paths & original, const Paths & paths) {
    Segment u = std::make_pair(original[lineIdx.first][lineIdx.second - 1], original[lineIdx.first][lineIdx.second]);
    Segment v = std::make_pair(paths[lineIdx.first][lineIdx.second - 1], paths[lineIdx.first][lineIdx.second]);
    Vec2 vecU = u.second - u.first;
    Vec2 vecV = v.second - v.first;
    Vec2 t = coord - u.first;
    if(eg::geo::norm(vecU) == 0) std::cout << "ZERO" << std::endl;
    double ratio = eg::geo::norm(t)/eg::geo::norm(vecU);
    double x = v.first.first + ratio*vecV.first;
    double y = v.first.second + ratio*vecV.second;
    return std::make_pair(round(x), round(y));
}

void updateLine(const std::pair<int, int> & toChange, int dx, int dy, Segments & segs, Paths & paths) {
    if(toChange.second == 0) {
        paths.at(toChange.first).at(0).first += dx;
        paths.at(toChange.first).at(0).second += dy;
        if(paths[toChange.first][0] == paths[toChange.first][paths[toChange.first].size() - 1]) {
            paths[toChange.first][paths[toChange.first].size() - 1].first += dx;
            paths[toChange.first][paths[toChange.first].size() - 1].second += dy;
        }
    }
    else if(toChange.second == paths[toChange.first].size() - 1) {
        paths[toChange.first][paths[toChange.first].size() - 1].first += dx;
        paths[toChange.first][paths[toChange.first].size() - 1].second += dy;
        if(paths[toChange.first][0] == paths[toChange.first][paths[toChange.first].size() - 1]) {
            paths[toChange.first][0].first += dx;
            paths[toChange.first][0].second += dy;
        }
    }
    else {
        paths.at(toChange.first).at(toChange.second).first += dx;
        paths.at(toChange.first).at(toChange.second).second += dy;
    }

    // std::cout << "HI" << std::endl;
    int ssss = 0;
    for(int i = 0; i < toChange.first - 1; i++)
        ssss += paths.at(i).size() - 1;
    if(toChange.second == 0) {
        segs.at(ssss + toChange.second).first.first += dx;
        segs.at(ssss + toChange.second).first.second += dy;
        if(paths[toChange.first][0] == paths[toChange.first][paths[toChange.first].size() - 1]) {
            segs[ssss + paths[toChange.first].size() - 1].second.first += dx;
            segs[ssss + paths[toChange.first].size() - 1].second.second += dy;
        }
    }
    else if(toChange.second == paths[toChange.first].size() - 1) {
        segs[ssss + toChange.second - 1].second.first += dx;
        segs[ssss + toChange.second - 1].second.second += dy;
        if(paths[toChange.first][0] == paths[toChange.first][paths[toChange.first].size() - 1]) {
            segs[ssss].second.first += dx;
            segs[ssss].second.second += dy;
        }
    }
    else {
        segs.at(ssss + toChange.second - 1).second.first += dx;
        segs.at(ssss + toChange.second - 1).second.second += dy;
        segs.at(ssss + toChange.second).first.first += dx;
        segs.at(ssss + toChange.second).first.second += dy;
    }
}

std::pair<Index2d, std::pair<int, int>> chooseRandomValue(const Paths & paths, std::_Bind<std::uniform_int_distribution<int> (std::mt19937)> & mt, int asciih, int asciiw, int h, int w) {
    Index2d toChange;
    int randI;
    // std::cout << paths.size() << std::endl;
    do{
        randI = mt() % paths.size();
    }while(paths.at(randI).size() < 2);

    // std::cout << randI << "SS" << paths.at(randI).size() << std::endl;
    int randJ = mt() % paths[randI].size();

    // std::cout << randI << " " << randJ << " selected." << std::endl;
    toChange = std::make_pair(randI, randJ);

    int dx, dy;
    bool lengthNotZero = false;

    do{
        lengthNotZero = true;
        dx = (mt() % asciih) - asciih/2;
        dy = (mt() % asciiw) - asciiw/2;
        if(paths[toChange.first][toChange.second].first + dx < 0 || paths[toChange.first][toChange.second].first + dx >= h) {
            lengthNotZero = false;
            continue;
        }
        if(paths[toChange.first][toChange.second].second + dy < 0 || paths[toChange.first][toChange.second].second + dy >= w) {
            lengthNotZero = false;
            continue;
        }

        if(toChange.second == 0) {
            const Dot a = paths[toChange.first][toChange.second];
            const Dot b = paths[toChange.first][toChange.second + 1];
            if(a.first + dx == b.first) {
                if(a.second + dy == b.second) {
                    lengthNotZero = false;
                }
            }
        }
        else if(toChange.second == paths[toChange.first].size() - 1) {
            const Dot a = paths[toChange.first][toChange.second - 1];
            const Dot b = paths[toChange.first][toChange.second];
            if(a.first == b.first + dx) {
                if(a.second == b.second + dy) {
                    lengthNotZero = false;
                }
            }
        }
        else {
            const Dot a = paths[toChange.first][toChange.second - 1];
            const Dot b = paths[toChange.first][toChange.second];
            const Dot c = paths[toChange.first][toChange.second + 1];
            if(a.first == b.first + dx) {
                if(a.second == b.second + dy) {
                    lengthNotZero = false;
                }
            }
            if(b.first + dx == c.first) {
                if(b.second + dy == c.second) {
                    lengthNotZero = false;
                }
            }
        }
    }while(dx*dx + dy*dy > std::max(asciih, asciiw)*std::max(asciih, asciiw) || !lengthNotZero);

    return std::make_pair(toChange, std::make_pair(dx, dy));
}

std::pair<double, int> calcAISS(int i, int j, const std::vector<SegWithOri> & linesInGrid, int asciih, int asciiw, int fCnt, double * distances, Mat2d * asciiPNGs) {
    Mat2d sample(asciih, asciiw);
    sample.setConstant(0);

    for(int k = 0; k < linesInGrid.size(); k++) {
        Segment lineInGrid = linesInGrid[k].first;
        lineInGrid.first.first -= i*asciih;
        lineInGrid.first.second -= j*asciiw;
        lineInGrid.second.first -= i*asciih;
        lineInGrid.second.second -= j*asciiw;
        drawSegmentDirect(sample, lineInGrid, 1);
    }

    Eigen::Tensor<double, 0> tmp = sample.sum();
    if(tmp(0) < 1) {
        return std::pair(-1, -1);
    }

    std::thread ts[THREAD_CNT];
    for(int k = 0; k < THREAD_CNT; k++) {
        ts[k] = std::thread(fillDist, std::ref(distances), fCnt*((double)k/THREAD_CNT), fCnt*((double)(k+1)/THREAD_CNT),
                            std::ref(sample), std::ref(asciiPNGs));
        ts[k].join();
    }

    int minIndex = 0;
    for(int k = 0; k < fCnt; k++) {
        if(distances[k] < distances[minIndex]) {
            minIndex = k;
        }
    }
    return std::make_pair(distances[minIndex], minIndex);
}

int main(int argc, char * argv[]) {
    if(argc != 3) {
        std::cout << "Use Program Properly! program ASCII_IMAGES_PATH INPUT_IMAGE_PATH" << std::endl;
        return -1;
    }
    std::string asciiImagesPath = argv[1];
    std::string inputImagePath = argv[2];

    int fCnt = getFileCount(asciiImagesPath);
    std::cout << "Opening ascii images" << std::endl;

    names = new std::string[fCnt];

    Mat2d * asciiPNGs = getAllImages(asciiImagesPath, fCnt);
    std::cout << "Opened ascii images" << std::endl;
    const int asciih = asciiPNGs[0].dimensions()[0];
    const int asciiw = asciiPNGs[0].dimensions()[1];

    eg::PNG inputImage;
    std::cout << "Opening input image" << std::endl;
    inputImage.openImage(inputImagePath);
    Image i = inputImage.copy();
    // for(int j = 0; j < i.dimensions()[0]; j++) {
    //     for(int k = 0; k < i.dimensions()[1]; k++) {
    //         std::cout << "/";
    //         for(int l = 0; l < i.dimensions()[2]; l++) {
    //             std::cout << (int)i(j, k, l) << " ";
    //         }
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << std::endl;
    std::cout << "Converting input image gray" << std::endl;
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    // print(t);
    std::cout << "Getting Edge of input image" << std::endl;
    t = getEdge(t, eg::edgeDetectMethod::gradient);

    t = markOutlier(t, 50);
    // print(t);
    t = dilate(t, 1, 3);
    t = erode(t, 1, 3);
    t = dilate(t, 3, 1);
    t = erode(t, 3, 1);
    i = mat2dToImage(t);

    std::cout << "getContours" << std::endl;
    auto contours = getContours(t, eg::contourMethod::suzuki);

    Segments segs;
    // for(int i = 0; i < contours.first[9].size(); i++) {
    //     std::cout << contours.first[9][i].first << " " << contours.first[9][i].second << std::endl;
    // }
    std::cout << "Merging contours" << std::endl;
    Paths paths;
    for(int i = 0; i < contours.first.size(); i++) {
        if(contours.first[i].size() == 0) continue;
        Path tmp = eg::trace::approxPath(contours.first[i]);
        paths.push_back(tmp);
        for(int j = 1; j < tmp.size(); j++)
            segs.push_back(std::make_pair(tmp[j - 1], tmp[j]));
    }
    // for(int i = 0; i < paths[9].size(); i++) {
    //     std::cout << paths[9][i].first << " " << paths[9][i].second << std::endl;
    // }
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
    std::cout << "Decompose finished" << std::endl;

#ifdef DEBUG_WITH_NODECOM_IMAGE
    Mat2d noDecomTestMat(inputImage.info.height, inputImage.info.width);
    noDecomTestMat.setConstant(0);
    for(int i = 0; i < segs.size(); i++) {
        drawSegmentDirect(noDecomTestMat, segs[i], 255);
    }
    Image noDecomImage = mat2dToImage(noDecomTestMat);
    PNG noDecom;
    noDecom.openImage(inputImage.getInputPath());
    noDecom.setImage(noDecomImage);
    noDecom.saveImage("noDecomed-approx.png");
#endif

#ifdef DEBUG_WITH_DECOM_IMAGE
    Mat2d decomTestMat(inputImage.info.height, inputImage.info.width);
    decomTestMat.setConstant(0);
    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            for(int k = 0; k < linesPerGrid[i][j].size(); k++) {
                drawSegmentDirect(decomTestMat, linesPerGrid[i][j][k].first, 255);
            }
        }
    }
    Image decomImage = mat2dToImage(decomTestMat);
    PNG decom;
    decom.openImage(inputImage.getInputPath());
    decom.setImage(decomImage);
    decom.saveImage("decomed-approx.png");
#endif

    std::cout << "Calculating initial E" << std::endl;

    int co = 0;
    int c = 1; // iteration index

    Mat2d sample(asciih, asciiw);

    double * distances = new double[fCnt];
    char ** buffer = new char * [grCnt];
    char ** best = new char * [grCnt];
    for(int i = 0; i < grCnt; i++) {
        buffer[i] = new char[gcCnt];
        best[i] = new char[gcCnt];
    }

    Mat2d errCell(grCnt, gcCnt);
    int K = grCnt*gcCnt;
    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            std::cout << "\rCalculating " << i << " " << j << std::flush;
            const std::pair<double, int> tmp = calcAISS(i, j, linesPerGrid[i][j], asciih, asciiw, fCnt, distances, asciiPNGs);
            const double minVal = tmp.first;
            const int minIndex = tmp.second;

            if(minIndex == -1) {
                buffer[i][j] = ' ';
                errCell(i, j) = 0;
                K--;
            }
            else {
                errCell(i, j) = minVal;
                buffer[i][j] = getAsciiFromPath(names[minIndex])[0];
            }
        }
    }

    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            best[i][j] = buffer[i][j];
            std::cout << buffer[i][j];
        }
        std::cout << "\n";
    }
    std::cout << std::endl;

    Eigen::Tensor<double, 0> tmpErr = errCell.sum();
    double prevE = tmpErr(0);
    double ta = prevE/(grCnt*gcCnt);
    if(K > 0) prevE /= K;
    double lastMinE = prevE;

    std::cout << "Finished calc. E: " << prevE << std::endl;

    time_t seed = time(nullptr);
    std::cout << "seed: " << seed << std::endl;
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> distribution(0, INF);
    auto mt = std::bind(distribution, engine);

    const int h = inputImage.info.height;
    const int w = inputImage.info.width;

#ifdef DEBUG_WITH_TMP_IMAGE
    PNG tmpPNG;
    tmpPNG.openImage(inputImagePath);
    Mat2d tmpMat2d(asciih*grCnt, asciiw*gcCnt);
    Image tmpImage(h, w, 4);
#endif

    bool possible = false;
    for(int i = 0; i < paths.size(); i++) {
        if(paths[i].size() >= 2) {
            possible = true;
            break;
        }
    }
    
    ArrayMat<std::vector<int>> candidates(grCnt);
    for(int i = 0; i < grCnt; i++) candidates[i].resize(gcCnt);
    for(int i = 0; i < grCnt; i++)
        for(int j = 0; j < gcCnt; j++)
            candidates[i][j].resize(linesPerGrid[i][j].size());
    bool * usedBool = new bool[originalSegs.size()];
    for(int k = 0; k < grCnt; k++) {
        for(int l = 0; l < gcCnt; l++) {
            for(int m = 0; m < linesPerGrid[k][l].size(); m++) {
                const Dot mid = std::make_pair((linesPerGrid[k][l][m].first.first.first + linesPerGrid[k][l][m].first.second.first)/2, (linesPerGrid[k][l][m].first.first.second + linesPerGrid[k][l][m].first.second.second)/2);
                for(int i = 0; i < originalSegs.size(); i++) usedBool[i] = false;
                for(double theta = 0; theta < 2*M_PI; theta += 0.01) {
                    double minVal = 987654321;
                    int minIndex = -1;
                    const Vec2 direc = std::make_pair(18000*std::cos(theta), 18000*std::sin(theta));
                    const Segment ray = std::make_pair(mid, mid + direc);
                    for(int i = 0; i < originalSegs.size(); i++) {
                        const auto & tmp = originalSegs[i];
                        if(ray.first.first < std::min(tmp.first.first, tmp.second.first) && ray.second.first < std::min(tmp.first.first, tmp.second.first)) continue;
                        if(ray.first.first > std::max(tmp.first.first, tmp.second.first) && ray.second.first > std::max(tmp.first.first, tmp.second.first)) continue;
                        if(eg::geo::ccw(ray.first, ray.second, tmp.first)*eg::geo::ccw(ray.first, ray.second, tmp.second) <= 0 &&
                        eg::geo::ccw(tmp.first, tmp.second, ray.first)*eg::geo::ccw(tmp.first, tmp.second, ray.second) <= 0) {
                            double dist = eg::geo::distSegDot(tmp, mid);
                            if(dist < 1) continue;
                            if(minVal > dist) {
                                minVal = dist;
                                minIndex = i;
                            }
                        }
                    }
                    if(minIndex != -1 && !usedBool[minIndex]) {
                        candidates[k][l][m].push_back(minIndex);
                        usedBool[minIndex] = true;
                    }
                }
            }
        }
    }
    delete[] usedBool;

    while(possible) {
        const auto tmp = chooseRandomValue(paths, mt, asciih, asciiw, h, w);
        const Index2d toChange = tmp.first;
        const int dx = tmp.second.first;
        const int dy = tmp.second.second;
        
        updateLine(toChange, dx, dy, segs, paths);

#ifdef DEBUG_WITH_TMP_IMAGE
        tmpMat2d.setConstant(0);
        tmpImage.setConstant(0);
#endif

        ArrayMat<SegWithOri> removedLines(grCnt);
        for(int i = 0; i < grCnt; i++) removedLines[i].resize(gcCnt);

        std::set<Index2d> needCalc;
        int KMinused = 0;
        for(int i = 0; i < grCnt; i++) { // remove lines in grid originated from updated line.
            for(int j = 0; j < gcCnt; j++) {
                std::vector<SegWithOri> & linesInGrid = linesPerGrid[i][j];
                for(int k = 0; k < linesInGrid.size(); k++) {
                    Index2d origin = linesInGrid[k].second;
                    if(origin.first != toChange.first) continue;
                    if((toChange.second == 0 && origin.second == 1) || (origin.second == toChange.second + 1) || (origin.second == toChange.second)) {
                        needCalc.insert(std::make_pair(i, j));
                        // std::cout << i << " " << j << " " << k << " removed" << std::endl;

#ifdef DEBUG_WITH_TMP_IMAGE
                        drawSegmentDirect(tmpMat2d, linesInGrid[k].first, 255);
#endif

                        removedLines[i][j].push_back(linesInGrid[k]);
                        linesInGrid.erase(linesInGrid.begin() + k--);

                        --K;
                        ++KMinused;
                    }
                }
            }
        }
#ifdef DEBUG_DRAW_GRID
        for(int i = 0; i < grCnt; i++) {
            if((i + 1)*asciih > h) break;
            for(int j = 0; j < w; j++) {
                for(int k = 0; k < 3; k++)
                    tmpImage((i + 1)*asciih, j, k) = 127;
            }
        }
        for(int i = 0; i < gcCnt; i++) {
            if((i + 1)*asciiw > w) break;
            for(int j = 0; j < h; j++) {
                for(int k = 0; k < 3; k++)
                    tmpImage(j, (i + 1)*asciiw, k) = 127;
            }
        }
#endif
#ifdef DEBUG_WITH_TMP_IMAGE
        for(const auto & p: needCalc) {
            int i = p.first, j = p.second;
            for(int k = 0; k < asciih; k++) {
                if(k + i*asciih < h) tmpImage(i*asciih + k, j*asciiw, 0) = 255;
                if((j + 1)*asciiw < w && k + i*asciih < h)tmpImage(k + i*asciih, (j + 1)*asciiw, 0) = 255;
            }
            for(int k = 0; k < asciiw; k++) {
                if(k + j*asciiw < w) tmpImage(i*asciih, k + j*asciiw, 0) = 255;
                if((i + 1)*asciih < h && k + j*asciiw < w)tmpImage((i + 1)*asciih, k + j*asciiw, 0) = 255;
            }
        }
        for(int i = 0; i < h; i++) {
            for(int j = 0; j < w; j++) {
                tmpImage(i, j, 3) = 255;
                if(tmpMat2d(i, j) > 0) tmpImage(i, j, 2) = 255; // blue: line removed
            }
        }
        tmpMat2d.setConstant(0);
        for(int i = 0; i < grCnt; i++) {
            for(int j = 0; j < gcCnt; j++) {
                for(int k = 0; k < linesPerGrid[i][j].size(); k++) {
                    drawSegmentDirect(tmpMat2d, linesPerGrid[i][j][k].first, 127);
                }
            }
        }
        for(int i = 0; i < h; i++) {
            for(int j = 0; j < w; j++) {
                if(tmpMat2d(i, j) > 0) tmpImage(i, j, 0) = 255; //red: line exist
            }
        }
        tmpMat2d.setConstant(0);
#endif
        // std::cout << "HI" << std::endl;
        std::vector<std::vector<int>> addedLinesCnt(grCnt);
        for(int i = 0; i < grCnt; i++) addedLinesCnt[i].resize(gcCnt);
        for(int i = 0; i < grCnt; i++)
            for(int j = 0; j < gcCnt; j++)
                addedLinesCnt[i][j] = 0;

        for(int i = 0; i < grCnt; i++) { // collect grids that updated line passes by.
            for(int j = 0; j < gcCnt; j++) {
                const Dot lu = std::make_pair(i*asciih, j*asciiw);
                const Dot rd = std::make_pair((i + 1)*asciih, (j + 1)*asciiw);
                if(toChange.second != 0) {
                    const Segment updatedLine = std::make_pair(paths[toChange.first][toChange.second - 1], paths[toChange.first][toChange.second]);
                    const Segment clipped = eg::tool::clip(updatedLine, lu, rd);
                    if(clipped.first.first != -1) {
                        linesPerGrid[i][j].push_back(std::make_pair(clipped, toChange));
                        addedLinesCnt[i][j]++;
                        needCalc.insert(std::make_pair(i, j));
#ifdef DEBUG_WITH_TMP_IMAGE
                        drawSegmentDirect(tmpMat2d, clipped, 255);
#endif
                    }
                }

                if(toChange.second != paths[toChange.first].size() - 1) {
                    const Segment updatedLine = std::make_pair(paths[toChange.first][toChange.second], paths[toChange.first][toChange.second + 1]);
                    const Segment clipped = eg::tool::clip(updatedLine, lu, rd);
                    if(clipped.first.first != -1) {
                        linesPerGrid[i][j].push_back(std::make_pair(clipped, std::make_pair(toChange.first, toChange.second + 1)));
                        addedLinesCnt[i][j]++;
                        needCalc.insert(std::make_pair(i, j)); // needCalc is a set. Don't worry about inserting twice~!!
#ifdef DEBUG_WITH_TMP_IMAGE
                        drawSegmentDirect(tmpMat2d, clipped, 255);
#endif
                    }
                }
            }
        }

#ifdef DEBUG_WITH_TMP_IMAGE
        for(const auto & p: needCalc) {
            int i = p.first, j = p.second;
            for(int k = 0; k < asciih; k++) {
                if(i*asciih + k < h) tmpImage(i*asciih + k, j*asciiw, 1) = 255;
                if((j + 1)*asciiw < w && k + i*asciih < h)tmpImage(k + i*asciih, (j + 1)*asciiw, 1) = 255;
            }
            for(int k = 0; k < asciiw; k++) {
                if(k + j*asciiw < w) tmpImage(i*asciih, k + j*asciiw, 1) = 255;
                if((i + 1)*asciih < h && k + j*asciiw < w)tmpImage((i + 1)*asciih, k + j*asciiw, 1) = 255;
            }
        }
        for(int i = 0; i < h; i++) {
            for(int j = 0; j < w; j++) {
                if(tmpMat2d(i, j) > 0) tmpImage(i, j, 1) = 255; // green: line drawed
            }
        }
        if(toChange.second == 0) {
            tmpMat2d.setConstant(0);
            auto t = std::make_pair(paths[toChange.first][toChange.second], paths[toChange.first][toChange.second + 1]);
            drawSegmentDirect(tmpMat2d, t, 255);
            for(int i = 0; i < h; i++) {
                for(int j = 0; j < w; j++) {
                    if(tmpMat2d(i, j) > 0) {
                        tmpImage(i, j, 0) = 135;
                        tmpImage(i, j, 1) = 206;
                        tmpImage(i, j, 2) = 235;
                    }
                }
            }
        }
        else if(toChange.second == paths[toChange.first].size() - 1) {
            tmpMat2d.setConstant(0);
            auto t = std::make_pair(paths[toChange.first][toChange.second - 1], paths[toChange.first][toChange.second]);
            drawSegmentDirect(tmpMat2d, t, 255);
            for(int i = 0; i < h; i++) {
                for(int j = 0; j < w; j++) {
                    if(tmpMat2d(i, j) > 0) {
                        tmpImage(i, j, 0) = 135;
                        tmpImage(i, j, 1) = 206;
                        tmpImage(i, j, 2) = 235;
                    }
                }
            }
        }
        else {
            tmpMat2d.setConstant(0);
            auto t = std::make_pair(paths[toChange.first][toChange.second - 1], paths[toChange.first][toChange.second]);
            drawSegmentDirect(tmpMat2d, t, 255);
            for(int i = 0; i < h; i++) {
                for(int j = 0; j < w; j++) {
                    if(tmpMat2d(i, j) < 1) continue;
                    tmpImage(i, j, 0) = 235*(tmpMat2d(i, j) > 0);
                    tmpImage(i, j, 1) = 162*(tmpMat2d(i, j) > 0);
                    tmpImage(i, j, 2) = 134*(tmpMat2d(i, j) > 0);
                }
            }
            tmpMat2d.setConstant(0);
            t = std::make_pair(paths[toChange.first][toChange.second], paths[toChange.first][toChange.second + 1]);
            drawSegmentDirect(tmpMat2d, t, 255);
            for(int i = 0; i < h; i++) {
                for(int j = 0; j < w; j++) {
                    if(tmpMat2d(i, j) < 1) continue;
                    tmpImage(i, j, 0) = 133*(tmpMat2d(i, j) > 0);
                    tmpImage(i, j, 1) = 234*(tmpMat2d(i, j) > 0);
                    tmpImage(i, j, 2) = 162*(tmpMat2d(i, j) > 0);
                }
            }
        }
        for(int i = 0; i < paths.size(); i++) {
            for(int j = 0; j < paths[i].size(); j++) {
                int x = paths[i][j].first;
                int y = paths[i][j].second;
                for(int k = 0; k < 3; k++)
                    tmpImage(x, y, k) = 255;
            }
        }
        tmpImage(paths[toChange.first][toChange.second].first - dx, paths[toChange.first][toChange.second].second - dy, 0) = 255;
        tmpImage(paths[toChange.first][toChange.second].first - dx, paths[toChange.first][toChange.second].second - dy, 1) = 110;
        tmpImage(paths[toChange.first][toChange.second].first - dx, paths[toChange.first][toChange.second].second - dy, 2) = 40;
        tmpImage(paths[toChange.first][toChange.second].first, paths[toChange.first][toChange.second].second, 0) = 110;
        tmpImage(paths[toChange.first][toChange.second].first, paths[toChange.first][toChange.second].second, 1) = 110;
        tmpImage(paths[toChange.first][toChange.second].first, paths[toChange.first][toChange.second].second, 2) = 183;

        /*
        red: line exists before and after removing lines.
        green: line drawed recently.
        blue: line removed and not redrawed.
        yellow(r+g): line exists before & after, and drawed.
        magenta(r+b): removed line and redraw.
        cyan(g+b): ??? WHy this exists?
        white: joint
        orange: pos of original vertex
        purple: pos of moved vertex
        skyblue: after line when 0 or size() - 1
        dirt, green: after line
        */

        tmpPNG.setImage(tmpImage);
        std::string oout = "main-sa-tmp-";
        std::string ext = ".png";
        tmpPNG.saveImage(oout + std::to_string(c) + ext);
        std::cout << "TMP IMAGE generated" << std::endl;
#endif

        bool forceReset = false;
        std::vector<std::pair<Index2d, double>> prevErr;
        std::vector<std::pair<Index2d, double>> prevBuff;
        // for(int i = 0; i < grCnt; i++) for(int j = 0; j < gcCnt; j++) needCalc.insert(std::make_pair(i, j));
        for(const auto & p : needCalc) { // re-calculate error per cell.
            const int i = p.first, j = p.second;
            // std::cout << i << " " << j << " needs calculate" << std::endl;
            const std::vector<SegWithOri> & linesInGrid = linesPerGrid[i][j];
            const std::vector<SegWithOri> & oriLinesInGrid = originalLinesPerGrid[i][j];
            if(linesInGrid.empty()) continue;

            double lengthSum = 0;
            for(int k = 0; k < oriLinesInGrid.size(); k++)
                lengthSum += eg::geo::euclideDist(oriLinesInGrid[k].first.first, oriLinesInGrid[k].first.second);

            if(lengthSum == 0) continue;

            double deform = 0;
            for(int k = 0; k < oriLinesInGrid.size(); k++) {
                const double lineLength = eg::geo::euclideDist(oriLinesInGrid[k].first.first, oriLinesInGrid[k].first.second);
                const double l = lineLength/lengthSum;
                const Segment before = oriLinesInGrid[k].first;
                const Dot afterFirst = getChanged(before.first, oriLinesInGrid[k].second, originalPaths, paths);
                const Dot afterSecond = getChanged(before.second, oriLinesInGrid[k].second, originalPaths, paths);
                const Segment after = std::make_pair(afterFirst, afterSecond);

                double calcedDeform;
                try{
                    calcedDeform = eg::math::calcDeform(before, after, segs, originalSegs, candidates[i][j][k]);
                    deform += l*calcedDeform;
                }
                catch (const eg::exceptions::InvalidParameter &) {
                    forceReset = true;
                    break;
                }
            }
            if(forceReset) break;

            prevErr.push_back(std::make_pair(std::make_pair(i, j), errCell(i, j)));
            prevBuff.push_back(std::make_pair(std::make_pair(i, j), buffer[i][j]));

            const std::pair<double, int> tmp = calcAISS(i, j, linesInGrid, asciih, asciiw, fCnt, distances, asciiPNGs);
            const double minVal = tmp.first;
            const int minIndex = tmp.second;

            if(minIndex == -1) {
                buffer[i][j] = ' ';
                errCell(i, j) = 0;
            }
            else {
                // std::cout << "DEFORM " << deform << " Minval " << minVal << std::endl;
                K++;
                --KMinused;
                buffer[i][j] = getAsciiFromPath(names[minIndex])[0];
                errCell(i, j) = 0.1*deform*minVal;
            }
        }

        double E = 0;
        if(!forceReset) {
            Eigen::Tensor<double, 0> tmpErr = errCell.sum();
            E = tmpErr(0)/K;
            std::cout << "E: " << E << std::endl;

            for(int i = 0; i < grCnt; i++) {
                for(int j = 0; j < gcCnt; j++) {
                    std::cout << buffer[i][j];
                }
                std::cout << "\n";
            }
            std::cout << std::endl;
            std::cout << "^^^====E: " << E << " c: " << c << " co: " << co << std::endl;
        }

        if(!forceReset && lastMinE > E) {
            std::cout << "best" << std::endl;
            prevE = E;
            lastMinE = E;
            co = 0;
            for(int i = 0; i < grCnt; i++)
                for(int j = 0; j < gcCnt; j++)
                    best[i][j] = buffer[i][j];
            c++;
        }
        else {
            double Pr, ran;
            if(!forceReset) {
                if(++co >= 5000) break;
                const double delta = prevE - E;
                const double t = 0.2*ta*std::pow(c, 0.997);
                Pr = std::exp(-std::abs(delta)/t);
                if(delta < -20) Pr = 15; // force reset
                ran = (double)mt()/INF;
                prevE = E;
                c++;
                // std::cout << "Pr " << Pr << " RAN " << ran << std::endl;
            }
            if(forceReset || (Pr >= ran)) { // goto previous state.
                std::cout << "goto prev" << std::endl;
                K += KMinused;

                updateLine(toChange, -dx, -dy, segs, paths);

                for(int i = 0; i < grCnt; i++) {
                    for(int j = 0; j < gcCnt; j++) {
                        for(int k = 0; k < addedLinesCnt[i][j]; k++) {
                            linesPerGrid[i][j].pop_back();
                        }
                    }
                }
                for(int i = 0; i < grCnt; i++) {
                    for(int j = 0; j < gcCnt; j++) {
                        for(int k = 0; k < removedLines[i][j].size(); k++) {
                            linesPerGrid[i][j].push_back(removedLines[i][j][k]);
                        }
                    }
                }
                for(const auto & p : prevErr) {
                    int i = p.first.first;
                    int j = p.first.second;
                    errCell(i, j) = p.second;
                }
                for(const auto & p : prevBuff) {
                    int i = p.first.first;
                    int j = p.first.second;
                    buffer[i][j] = p.second;
                }
            }
        }
    }

    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            std::cout << best[i][j];
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
    std::cout << "^^^====E: " << lastMinE << std::endl;

    for(int i = 0; i < grCnt; i++) {
        delete[] buffer[i];
        delete[] best[i];
    }
    delete[] best;
    delete[] buffer;
    delete[] distances;
    delete[] asciiPNGs;
    delete[] names;
}
