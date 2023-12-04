/**
 * @file main.cpp
 * @author Daniel cho
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

#define INF 987654321
#define THREAD_CNT 8

using namespace eg::imgproc;
namespace fs = std::filesystem;

std::string * names;

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
        ans[i] = 255*ans[i]; // This is necessary because it will increase error when misaligned.
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

int computeOutCode(const Dot & p, const Dot & lu, const Dot & rd) {
    const int INSIDE = 0;
    const int LEFT = 1;
    const int RIGHT = 2;
    const int BOTTOM = 4;
    const int TOP = 8;

    int code = INSIDE;

    if(p.first < lu.first)
        code |= LEFT;
    else if(p.first > rd.first)
        code |= RIGHT;
    if(p.second < lu.second)
        code |= BOTTOM;
    else if(p.second > rd.second)
        code |= TOP;

    return code;
}

/**
 * @breif do line clipping using Cohen-Sutherland algorithm
 * @see https://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland_algorithm
 */
Segment clip(const Segment & s, const Dot & lu, const Dot & rd) {
    const int INSIDE = 0;
    const int LEFT = 1;
    const int RIGHT = 2;
    const int BOTTOM = 4;
    const int TOP = 8;

    int outcode0 = computeOutCode(s.first, lu, rd);
    int outcode1 = computeOutCode(s.second, lu, rd);
    bool accept = false;

    Segment ans = s;

    while(true) {
        if(!(outcode0 | outcode1)) {
            accept = true;
            break;
        }
        else if(outcode0 & outcode1)
            break;
        else {
            int outcodeOut = (outcode1 > outcode0)?outcode1:outcode0;
            int x, y;
            int x0 = s.first.first;
            int y0 = s.first.second;
            int x1 = s.second.first;
            int y1 = s.second.second;
            int xmax = rd.first;
            int xmin = lu.first;
            int ymax = rd.second;
            int ymin = lu.second;
            if(outcodeOut & TOP) {
                x = round(x0 + (x1 - x0)*(ymax - y0)/(y1 - y0));
                y = ymax;
            }
            else if(outcodeOut & BOTTOM) {
                x = round(x0 + (x1 - x0)*(ymin - y0)/(y1 - y0));
                y = ymin;
            }
            else if(outcodeOut & RIGHT) {
                y = round(y0 + (y1 - y0)*(xmax - x0)/(x1 - x0));
                x = xmax;
            }
            else if(outcodeOut & LEFT) {
                y = round(y0 + (y1 - y0)*(xmin - x0)/(x1 - x0));
                x = xmin;
            }
            if(outcodeOut == outcode0) {
                ans.first.first = x;
                ans.first.second = y;
                outcode0 = computeOutCode(ans.first, lu, rd);
            }
            else {
                ans.second.first = x;
                ans.second.second = y;
                outcode1 = computeOutCode(ans.second, lu, rd);
            }
        }
    }
    if(accept)
        return ans;
    return std::make_pair(std::make_pair(-1, -1), std::make_pair(-1, -1));
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
    std::cout << "Converting input image gray" << std::endl;
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    std::cout << "Getting Edge of input image" << std::endl;
    t = getEdge(t, eg::edgeDetectMethod::gradient);
    t = markOutlier(t, 10);
    t = dilate(t, 1, 3);
    t = erode(t, 1, 3);
    t = dilate(t, 3, 1);
    t = erode(t, 3, 1);
    i = mat2dToImage(t);
    std::cout << "getContours" << std::endl;
    auto contours = getContours(t, eg::contourMethod::suzuki);
    Segments segs;
    std::cout << "Merging contours" << std::endl;
    Paths paths(contours.first.size());
    for(int i = 2; i < contours.first.size(); i++) {
        if(contours.first[i].size() == 0) continue;
        Path tmp = eg::trace::approxPath(contours.first[i]);
        paths[i - 2] = tmp;
        for(int j = 1; j < tmp.size(); j++)
            segs.push_back(std::make_pair(tmp[j - 1], tmp[j]));
    }
    Segments originalSegs = segs;
    Paths original = paths;

    inputImage.divideImageByLength(asciih, asciiw);
    int gcCnt = inputImage.getGridColCnt();
    int grCnt = inputImage.getGridRowCnt();

    std::cout << grCnt << "x" << gcCnt << std::endl;

    char devnull[10];


    std::cout << "not DEFORM " << eg::math::calcDeform(segs[0], segs[0], segs, originalSegs) << std::endl;

    std::cout << "Decomposing lines into grids.." << std::endl;
    std::vector<std::vector<std::vector<std::pair<Segment, std::pair<int, int>>>>> lines(grCnt);
    for(int i = 0; i < grCnt; i++)
        lines[i].resize(gcCnt);

    for(int k = 0; k < paths.size(); k++) {
        for(int l = 1; l < paths[k].size(); l++) {
            for(int i = 0; i < grCnt; i++) {
                for(int j = 0; j < gcCnt; j++) {
                    Dot lu = std::make_pair(i*asciih, j*asciiw);
                    Dot rd = std::make_pair((i + 1)*asciih, (j + 1)*asciiw);
                    Segment original = std::make_pair(paths[k][l - 1], paths[k][l]);
                    Segment tmp = clip(original, lu, rd);
                    if(tmp.first == tmp.second) continue;
                    if(tmp.first.first == -1) continue;
                    // std::cout << lu.first << "/" << lu.second << " " << rd.first << "/" << rd.second
                            //   << " " << "clip"  << tmp.first.first << " " << tmp.first.second << "/" << tmp.second.first << " " << tmp.second.second << std::endl;

                    lines[i][j].push_back(std::make_pair(tmp, std::make_pair(k, l)));
                }
            }
        }
    }
    /*
    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            for(int k = 0; k < lines[i][j].size(); k++) {
                if(lines[i][j][k].first.first == lines[i][j][k].first.second) {
                    std::cout << "SAME" << std::endl;
                }
            }
        }
    }
    */
    std::cout << "Decompose finished" << std::endl;
    std::vector<std::vector<std::vector<std::pair<Segment, std::pair<int, int>>>>> originalLines = lines;
    Mat2d decomTestMat(inputImage.info.height, inputImage.info.width);
    decomTestMat.setConstant(0);
    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            for(int k = 0; k < lines[i][j].size(); k++) {
                drawSegmentDirect(decomTestMat, lines[i][j][k].first, 255);
            }
        }
    }
    Image decomImage = mat2dToImage(decomTestMat);
    PNG decom;
    decom.openImage(inputImage.getInputPath());
    decom.setImage(decomImage);
    decom.saveImage("decomed-approx.png");

    std::cout << "Calculating initial E" << std::endl;

    int co = 0;
    int c = 0; // iteration index
    Mat2d zeros(inputImage.info.height, inputImage.info.width);
    zeros.setConstant(0);

    Mat2d sample(asciih, asciiw);

    double * distances = new double[fCnt];
    char ** buffer = new char * [grCnt];
    char ** best = new char * [grCnt];
    for(int i = 0; i < grCnt; i++) {
        buffer[i] = new char[gcCnt];
        best[i] = new char[grCnt];
    }
    double prevE = 0;
    /*std::thread ts[THREAD_CNT];
    for(int i = 0; i < THREAD_CNT; i++)
        ts[i] = std::thread(fillDist, std::ref(distances), fCnt*(i/THREAD_CNT), fCnt*((i+1)/THREAD_CNT),
                            std::ref(sample), std::ref(asciiPNGs));*/

    Mat2d errCell(grCnt, gcCnt);
    int K = grCnt*gcCnt;
    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            std::cout << "\rCalculating " << i << " " << j << std::flush;
            sample.setConstant(0);
            for(int k = 0; k < lines[i][j].size(); k++) {
                Segment tmp = lines[i][j][k].first;
                tmp.first.first -= i*asciih;
                tmp.first.second -= j*asciiw;
                tmp.second.first -= i*asciih;
                tmp.second.second -= j*asciiw;
                drawSegmentDirect(sample, tmp, 255); // setting val as 255 will increase misaligned err.
            }

            Eigen::Tensor<double, 0> tmp = sample.sum();
            if(tmp(0) < 2) {
                buffer[i][j] = ' ';
                errCell(i, j) = 0;
                --K;
                continue;
            }

            /*std::thread ts[THREAD_CNT];
            for(int k = 0; k < THREAD_CNT; k++) {
                ts[k] = std::thread(fillDist, std::ref(distances), fCnt*((double)k/THREAD_CNT), fCnt*((double)(k+1)/THREAD_CNT),
                                    std::ref(sample), std::ref(asciiPNGs));
                ts[k].join();
            }*/
            fillDist(distances, 0, fCnt, sample, asciiPNGs);
            /*
            for(int i = 0; i < THREAD_CNT; i++) {
                ts[i].join();
            }*/
            double minVal = INF;
            int minIndex = -1;
            for(int k = 0; k < fCnt; k++) {
                if(distances[k] < minVal) {
                    minVal = distances[k];
                    minIndex = k;
                }
            }
            errCell(i, j) = minVal;
            prevE += minVal;
            buffer[i][j] = getAsciiFromPath(names[minIndex])[0];
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
    double ta = prevE/(grCnt*gcCnt); // t_a
    prevE /= K;
    double lastMinE = prevE;
    std::cout << "Finished calc. E: " << prevE << std::endl;
    time_t seed = time(nullptr);
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> distribution(0, INF);
    auto mt = std::bind(distribution, engine);

    PNG tmpPNG;
    tmpPNG.openImage(inputImage.getInputPath());
    Mat2d tmpMat2d(tmpPNG.info.height, tmpPNG.info.width);

    while(true) {
        std::pair<int, int> toChange;
        int randI;
        do{
            randI = mt() % paths.size();
        }while(paths.at(randI).size() < 2);
        // std::cout << randI << "SS" << paths.at(randI).size() << std::endl;
        int randJ;
        if(paths.at(randI).size() == 1)
            randJ = 1;
        else
            randJ = (mt() % (paths[randI].size() - 1)) + 1;
        // std::cout << randI << " " << randJ << " selected." << std::endl;
        toChange = std::make_pair(randI, randJ);
        int dx, dy;
        bool lengthNotZero = false;
        do{
            lengthNotZero = true;
            dx = (mt() % asciih) - asciih/2;
            dy = (mt() % asciiw) - asciiw/2;
            if(toChange.second == 0) {
                if(paths[toChange.first][toChange.second].first + dx == paths[toChange.first][toChange.second + 1].first) {
                    if(paths[toChange.first][toChange.second].second + dx == paths[toChange.first][toChange.second + 1].second) {
                        lengthNotZero = false;
                    }
                }
            }
            else if(toChange.second == paths[toChange.first].size() - 1) {
                if(paths[toChange.first][toChange.second - 1].first + dx == paths[toChange.first][toChange.second].first) {
                    if(paths[toChange.first][toChange.second - 1].second + dx == paths[toChange.first][toChange.second].second) {
                        lengthNotZero = false;
                    }
                }
            }
            else {
                if(paths[toChange.first][toChange.second - 1].first + dx == paths[toChange.first][toChange.second].first) {
                    if(paths[toChange.first][toChange.second - 1].second + dx == paths[toChange.first][toChange.second].second) {
                        lengthNotZero = false;
                    }
                }
                if(paths[toChange.first][toChange.second].first + dx == paths[toChange.first][toChange.second + 1].first) {
                    if(paths[toChange.first][toChange.second].second + dx == paths[toChange.first][toChange.second + 1].second) {
                        lengthNotZero = false;
                    }
                }
            }
        }while(dx*dx + dy*dy > std::max(asciih, asciiw) || !lengthNotZero);
        paths.at(toChange.first).at(toChange.second).first += dx;
        paths.at(toChange.first).at(toChange.second).second += dy;

        // std::cout << "HI" << std::endl;
        int ssss = 0;
        for(int i = 0; i < toChange.first - 1; i++)
            ssss += paths.at(i).size() - 1;
        if(toChange.second == 0) {
            segs.at(ssss + toChange.second).first.first += dx;
            segs.at(ssss + toChange.second).first.second += dy;
        }
        else if(toChange.second == paths[toChange.first].size() - 1) {
            segs.at(ssss + toChange.second).second.first += dx;
            segs.at(ssss + toChange.second).second.second += dy;
        }
        else {
            segs.at(ssss + toChange.second - 1).second.first += dx;
            segs.at(ssss + toChange.second - 1).second.second += dy;
            segs.at(ssss + toChange.second).first.first += dx;
            segs.at(ssss + toChange.second).first.second += dy;
        }

        tmpMat2d.setConstant(0);
        for(int i = 0; i < grCnt; i++) {
            for(int j = 0; j < gcCnt; j++) {
                for(int k = 0; k < lines[i][j].size(); k++) {
                    drawSegmentDirect(tmpMat2d, lines[i][j][k].first, 255);
                }
            }
        }
        // tmpMat2d = drawSegments(tmpMat2d, segs, 255);
        Image tmpImage = mat2dToImage(tmpMat2d);
        tmpPNG.setImage(tmpImage);
        tmpPNG.saveImage("main-sa-tmp.png");
        std::cout << "TMP IMAGE generated" << std::endl;
        // std::cout << "HI" << std::endl;
        std::set<std::pair<int, int>> needCalc;
        for(int i = 0; i < grCnt; i++) {
            for(int j = 0; j < gcCnt; j++) {
                for(int k = 0; k < lines[i][j].size(); k++) {
                    if(lines[i][j][k].second == toChange) {
                        needCalc.insert(std::make_pair(i, j));
                        lines[i][j].erase(lines[i][j].begin() + k);
                        K--;
                        break;
                    }
                }
            }
        }

        // std::cout << "HI" << std::endl;
        for(int i = 0; i < grCnt; i++) {
            for(int j = 0; j < gcCnt; j++) {
                Dot lu = std::make_pair(i*asciih, j*asciiw);
                Dot rd = std::make_pair((i + 1)*asciih, (j + 1)*asciiw);
                Segment ttmp = std::make_pair(paths.at(toChange.first).at(toChange.second - 1), paths.at(toChange.first).at(toChange.second));
                Segment tmp = clip(ttmp, lu, rd);
                if(tmp.first.first == -1) continue;

                lines.at(i).at(j).push_back(std::make_pair(tmp, toChange));
                needCalc.insert(std::make_pair(i, j));
            }
        }

        bool forceReset = false;
        std::vector<std::pair<std::pair<int, int>, double>> prevErr;
        // std::cout << "HI" << std::endl;
        for(const auto & p : needCalc) {
            sample.setConstant(0);

            int i = p.first, j = p.second;
            if(lines[i][j].empty()) continue;
            // std::cout << "needCalc " << i << " " << j << std::endl;
            double lengthSum = 0;
            for(int k = 0; k < lines[i][j].size(); k++)
                lengthSum += eg::geo::euclideDist(lines[i][j][k].first.first, lines[i][j][k].first.second);
            // std::cout << "lengthSum: " << lengthSum << std::endl;
            if(lengthSum == 0) continue;
            double deform = 0;
            for(int k = 0; k < originalLines[i][j].size(); k++) {
                // std::cout << "calc deform of " << k << std::endl;
                // std::cout << i << " " << j << " " << k << std::endl;
                // std::cout << originalLines.at(i).at(j).at(k).first.first.first << std::endl;
                // std::cout << originalLines[i][j][k].first.first.second << std::endl;
                // std::cout << originalLines[i][j][k].first.second.first << std::endl;
                // std::cout << originalLines[i][j][k].first.second.second << std::endl;
                // std::cout << originalLines[i][j][k].second.first << std::endl;
                // std::cout << originalLines[i][j][k].second.second << std::endl;
                double l = eg::geo::euclideDist(originalLines[i][j][k].first.first, originalLines[i][j][k].first.second)/lengthSum;
                // std::cout << "HI" << std::endl;
                // std::cout << "1" << std::endl;
                Segment before = originalLines[i][j][k].first;
                Dot afterFirst = getChanged(before.first, originalLines[i][j][k].second, original, paths);
                Dot afterSecond = getChanged(before.second, originalLines[i][j][k].second, original, paths);
                Segment after = std::make_pair(afterFirst, afterSecond);
                // std::cout << "2" << std::endl;
                double tttt = eg::math::calcDeform(before, after, segs, originalSegs);
                if(tttt == 987) {
                    forceReset = true;
                    break;
                }
                deform += l*tttt;
                // std::cout << "3" << std::endl;
                Segment tmp = lines[i][j][k].first;
                tmp.first.first -= i*asciih;
                tmp.first.second -= j*asciiw;
                tmp.second.first -= i*asciih;
                tmp.second.second -= j*asciiw;
                drawSegmentDirect(sample, tmp, 255); // setting val as 255 will increase misaligned err.
            }
            if(forceReset) break;
            prevErr.push_back(std::make_pair(std::make_pair(i, j), errCell(i, j)));
            // std::cout << "sum" << std::endl;
            Eigen::Tensor<double, 0> oneCnt = sample.sum();
            if(oneCnt(0) < 2) {
                buffer[i][j] = ' ';
                errCell(i, j) = 0;
                continue;
            }
            // std::cout << "calc distance" << std::endl;
            K++;
            /*
            std::thread ts[THREAD_CNT];
            for(int k = 0; k < THREAD_CNT; k++) {
                ts[k] = std::thread(fillDist, std::ref(distances), fCnt*((double)k/THREAD_CNT), fCnt*((double)(k+1)/THREAD_CNT),
                                    std::ref(sample), std::ref(asciiPNGs));
                ts[k].join();
            }*/
            fillDist(distances, 0, fCnt, sample, asciiPNGs);
            double minVal = INF;
            int minIndex = -1;
            for(int k = 0; k < fCnt; k++) {
                if(distances[k] < minVal) {
                    minVal = distances[k];
                    minIndex = k;
                }
            }
            // std::cout << "deform " << deform << "minval " << minVal << std::endl;
            errCell(i, j) = deform*minVal;
            // if(minIndex == -1) std::cout << "MINUS" << std::endl;
            buffer[i][j] = getAsciiFromPath(names[minIndex])[0];
        }

        double E;
        if(!forceReset) {
            E = 0;
            for(int i = 0; i < grCnt; i++)
                for(int j = 0; j < gcCnt; j++)
                    E += errCell(i, j);
            // std::cout << "RAW E: " << E << std::endl;
            E /= K;
            std::cout << "E: " << E << std::endl;
        }

        if(lastMinE > E) {
            E = lastMinE;
            co = 0;
            for(int i = 0; i < grCnt; i++) {
                for(int j = 0; j < gcCnt; j++) {
                    best[i][j] = buffer[i][j];
                    std::cout << buffer[i][j];
                }
                std::cout << "\n";
            }
            std::cout << std::endl;
            std::cout << "^^^====E: " << E << " c: " << c << " co: " << co << std::endl;
        }
        else {
            double Pr, ran;
            if(!forceReset) {
                if(++co >= 5000) break;
                double delta = abs(E - prevE);
                double t = 0.2*ta*std::pow(c, 0.997);
                Pr = std::exp(-delta/t);
                ran = (double)mt()/INF;
                for(int i = 0; i < grCnt; i++) {
                    for(int j = 0; j < gcCnt; j++) {
                        std::cout << buffer[i][j];
                    }
                    std::cout << "\n";
                }
                std::cout << std::endl;
                std::cout << "^^^====E: " << E << " c: " << c << " co: " << co << std::endl;
                // std::cout << "Pr " << Pr << " RAN " << ran << std::endl;
            }
            if(forceReset || (Pr >= ran)) { // goto previous state.
                for(const auto & p : prevErr) {
                    int i = p.first.first;
                    int j = p.first.second;
                    errCell(i, j) = p.second;
                }
                paths[toChange.first][toChange.second].first -= dx;
                paths[toChange.first][toChange.second].second -= dy;
                if(toChange.second == 0) {
                    segs[ssss + toChange.second].first.first -= dx;
                    segs[ssss + toChange.second].first.second -= dy;
                }
                else if(toChange.second == paths[toChange.first].size() - 1) {
                    segs[ssss + toChange.second].second.first -= dx;
                    segs[ssss + toChange.second].second.second -= dy;
                }
                else {
                    segs[ssss + toChange.second - 1].second.first -= dx;
                    segs[ssss + toChange.second - 1].second.second -= dy;
                    segs[ssss + toChange.second].first.first -= dx;
                    segs[ssss + toChange.second].first.second -= dy;
                }
            }
        }
        if(!forceReset) {
            prevE = E;
            c++;
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
