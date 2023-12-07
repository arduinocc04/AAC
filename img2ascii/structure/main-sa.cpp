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
    int k = 0;
    for(int i = 2; i < contours.first.size(); i++) {
        if(contours.first[i].size() == 0) continue;
        Path tmp = eg::trace::approxPath(contours.first[i]);
        if(tmp.size() > 4)
            paths[k++] = tmp;
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

    std::cout << "segs size: " << segs.size() << std::endl;

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
                    Segment tmp = eg::tool::clip(original, lu, rd);
                    if(tmp.first == tmp.second) continue;
                    if(tmp.first.first == -1) continue;

                    lines[i][j].push_back(std::make_pair(tmp, std::make_pair(k, l)));
                }
            }
        }
    }

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
    int c = 1; // iteration index

    Mat2d sample(asciih, asciiw);

    double * distances = new double[fCnt];
    char ** buffer = new char * [grCnt];
    char ** best = new char * [grCnt];
    for(int i = 0; i < grCnt; i++) {
        buffer[i] = new char[gcCnt];
        best[i] = new char[gcCnt];
    }
    for(int i = 0; i < grCnt; i++)
        for(int j = 0; j < gcCnt; j++)
            buffer[i][j] = ' ';
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
                drawSegmentDirect(sample, tmp, 1);
            }

            Eigen::Tensor<double, 0> tmp = sample.sum();
            if(tmp(0) < 4) {
                buffer[i][j] = ' ';
                errCell(i, j) = 0;
                --K;
                continue;
            }

            std::thread ts[THREAD_CNT];
            for(int k = 0; k < THREAD_CNT; k++) {
                ts[k] = std::thread(fillDist, std::ref(distances), fCnt*((double)k/THREAD_CNT), fCnt*((double)(k+1)/THREAD_CNT),
                                    std::ref(sample), std::ref(asciiPNGs));
                ts[k].join();
            }
            // fillDist(distances, 0, fCnt, sample, asciiPNGs);
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
    double ta = prevE/(grCnt*gcCnt);
    if(K > 0) prevE /= K;
    double lastMinE = prevE;
    std::cout << "Finished calc. E: " << prevE << std::endl;
    time_t seed = time(nullptr);
    std::cout << "seed: " << seed << std::endl;
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> distribution(0, INF);
    auto mt = std::bind(distribution, engine);

    PNG tmpPNG;
    tmpPNG.openImage(inputImagePath);
    Mat2d tmpMat2d(asciih*grCnt, asciiw*gcCnt);

    while(true) {
        std::pair<int, int> toChange;
        int randI;
        std::cout << paths.size() << std::endl;
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
            if(toChange.second == 0) {
                if(paths[toChange.first][toChange.second].first + dx == paths[toChange.first][toChange.second + 1].first) {
                    if(paths[toChange.first][toChange.second].second + dy == paths[toChange.first][toChange.second + 1].second) {
                        lengthNotZero = false;
                    }
                }
            }
            else if(toChange.second == paths[toChange.first].size() - 1) {
                if(paths[toChange.first][toChange.second - 1].first == paths[toChange.first][toChange.second].first + dx) {
                    if(paths[toChange.first][toChange.second - 1].second == paths[toChange.first][toChange.second].second + dy) {
                        lengthNotZero = false;
                    }
                }
            }
            else {
                if(paths[toChange.first][toChange.second - 1].first == paths[toChange.first][toChange.second].first + dx) {
                    if(paths[toChange.first][toChange.second - 1].second == paths[toChange.first][toChange.second].second + dy) {
                        lengthNotZero = false;
                    }
                }
                if(paths[toChange.first][toChange.second].first + dx == paths[toChange.first][toChange.second + 1].first) {
                    if(paths[toChange.first][toChange.second].second + dy == paths[toChange.first][toChange.second + 1].second) {
                        lengthNotZero = false;
                    }
                }
            }
        }while(dx*dx + dy*dy > std::max(asciih, asciiw)*std::max(asciih, asciiw) || !lengthNotZero);
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
            segs[ssss + toChange.second - 1].second.first += dx;
            segs[ssss + toChange.second - 1].second.second += dy;
        }
        else {
            segs.at(ssss + toChange.second - 1).second.first += dx;
            segs.at(ssss + toChange.second - 1).second.second += dy;
            segs.at(ssss + toChange.second).first.first += dx;
            segs.at(ssss + toChange.second).first.second += dy;
        }

        int h = inputImage.info.height;
        int w = inputImage.info.width;
        Image tmpImage(h, w, 4);
        tmpMat2d.setConstant(0);

        std::vector<std::vector<std::vector<std::pair<Segment, std::pair<int, int>>>>> removedLines(grCnt);
        for(int i = 0; i < grCnt; i++) removedLines[i].resize(gcCnt);

        std::set<std::pair<int, int>> needCalc;
        for(int i = 0; i < grCnt; i++) {
            for(int j = 0; j < gcCnt; j++) {
                std::vector<int> toRemove;
                for(int k = 0; k < lines[i][j].size(); k++) {
                    if(lines[i][j][k].second == toChange || (lines[i][j][k].second.first == toChange.first && toChange.second == 0 && lines[i][j][k].second.second == 1)) {
                        needCalc.insert(std::make_pair(i, j));
                        drawSegmentDirect(tmpMat2d, lines[i][j][k].first, 255);
                        removedLines[i][j].push_back(lines[i][j][k]);
                        toRemove.push_back(k);
                        K--;
                        break;
                    }
                    else if(lines[i][j][k].second.first == toChange.first && (lines[i][j][k].second.second == toChange.second + 1 || (toChange.second == lines[i][j][k].second.second))) {
                        needCalc.insert(std::make_pair(i, j));
                        drawSegmentDirect(tmpMat2d, lines[i][j][k].first, 255);
                        removedLines[i][j].push_back(lines[i][j][k]);
                        toRemove.push_back(k);
                        K--;
                        break;
                    }
                }
                for(int k = 0; k < toRemove.size(); k++) {
                    lines[i][j].erase(lines[i][j].begin() + toRemove[k] - k);
                }
            }
        }

        for(int i = 0; i < h; i++) {
            for(int j = 0; j < w; j++) {
                tmpImage(i, j, 3) = 255;
                tmpImage(i, j, 2) = 255*(tmpMat2d(i, j) > 0); // blue: line removed
            }
        }
        tmpMat2d.setConstant(0);
        for(int i = 0; i < grCnt; i++) {
            for(int j = 0; j < gcCnt; j++) {
                for(int k = 0; k < lines[i][j].size(); k++) {
                    drawSegmentDirect(tmpMat2d, lines[i][j][k].first, 127);
                }
            }
        }
        for(int i = 0; i < h; i++) {
            for(int j = 0; j < w; j++) {
                tmpImage(i, j, 3) = 255;
                tmpImage(i, j, 0) = 255*(tmpMat2d(i, j) > 0); //red: line exist
            }
        }
        tmpMat2d.setConstant(0);
        // std::cout << "HI" << std::endl;
        std::vector<std::vector<int>> addedLinesCnt(grCnt);
        for(int i = 0; i < grCnt; i++) addedLinesCnt[i].resize(gcCnt);
        for(int i = 0; i < grCnt; i++)
            for(int j = 0; j < gcCnt; j++)
                addedLinesCnt[i][j] = 0;

        for(int i = 0; i < grCnt; i++) {
            for(int j = 0; j < gcCnt; j++) {
                Dot lu = std::make_pair(i*asciih, j*asciiw);
                Dot rd = std::make_pair((i + 1)*asciih, (j + 1)*asciiw);
                if(toChange.second != 0) {
                    Segment ttmp = std::make_pair(paths[toChange.first][toChange.second - 1], paths[toChange.first][toChange.second]);
                    Segment tmp = eg::tool::clip(ttmp, lu, rd);
                    if(tmp.first.first != -1) {
                        lines[i][j].push_back(std::make_pair(tmp, toChange));
                        addedLinesCnt[i][j]++;
                        needCalc.insert(std::make_pair(i, j));
                        drawSegmentDirect(tmpMat2d, tmp, 255);
                    }
                }

                if(toChange.second != paths[toChange.first].size() - 1) {
                    Segment ttmp = std::make_pair(paths[toChange.first][toChange.second], paths[toChange.first][toChange.second + 1]);
                    Segment tmp = eg::tool::clip(ttmp, lu, rd);
                    if(tmp.first.first != -1) {
                        lines[i][j].push_back(std::make_pair(tmp, std::make_pair(toChange.first, toChange.second + 1)));
                        addedLinesCnt[i][j]++;
                        needCalc.insert(std::make_pair(i, j)); // needCalc is set.
                        drawSegmentDirect(tmpMat2d, tmp, 255);
                    }
                }
            }
        }
        for(int i = 0; i < h; i++) {
            for(int j = 0; j < w; j++) {
                tmpImage(i, j, 1) = 255*(tmpMat2d(i, j) > 0); // green: line drawed
            }
        }
        if(toChange.second == 0) {
            tmpMat2d.setConstant(0);
            auto t = std::make_pair(paths[toChange.first][toChange.second], paths[toChange.first][toChange.second + 1]);
            drawSegmentDirect(tmpMat2d, t, 255);
            for(int i = 0; i < h; i++) {
                for(int j = 0; j < w; j++) {
                    if(tmpMat2d(i, j) < 1) continue;
                    tmpImage(i, j, 0) = 135*(tmpMat2d(i, j) > 0);
                    tmpImage(i, j, 1) = 206*(tmpMat2d(i, j) > 0);
                    tmpImage(i, j, 2) = 235*(tmpMat2d(i, j) > 0);
                }
            }
        }
        else if(toChange.second == paths[toChange.first].size() - 1) {
            tmpMat2d.setConstant(0);
            auto t = std::make_pair(paths[toChange.first][toChange.second - 1], paths[toChange.first][toChange.second]);
            drawSegmentDirect(tmpMat2d, t, 255);
            for(int i = 0; i < h; i++) {
                for(int j = 0; j < w; j++) {
                    if(tmpMat2d(i, j) < 1) continue;
                    tmpImage(i, j, 0) = 135*(tmpMat2d(i, j) > 0);
                    tmpImage(i, j, 1) = 206*(tmpMat2d(i, j) > 0);
                    tmpImage(i, j, 2) = 235*(tmpMat2d(i, j) > 0);
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
                // std::cout << "2" << std::endl
                double tttt;
                try{
                    tttt = eg::math::calcDeform(before, after, segs, originalSegs);
                }
                catch (const eg::exceptions::InvalidParameter &) {
                    forceReset = true;
                    break;
                }
                deform += l*tttt;
                // std::cout << "3" << std::endl;
            }
            if(forceReset) break;

            for(int k = 0; k < lines[i][j].size(); k++)
                drawSegmentDirect(sample, lines[i][j][k].first, 1);

            prevErr.push_back(std::make_pair(std::make_pair(i, j), errCell(i, j)));
            // std::cout << "sum" << std::endl;
            Eigen::Tensor<double, 0> oneCnt = sample.sum();
            if(oneCnt(0) < 4) {
                buffer[i][j] = ' ';
                errCell(i, j) = 0;
                continue;
            }
            // std::cout << "calc distance" << std::endl;
            K++;
            std::thread ts[THREAD_CNT];
            for(int k = 0; k < THREAD_CNT; k++) {
                ts[k] = std::thread(fillDist, std::ref(distances), fCnt*((double)k/THREAD_CNT), fCnt*((double)(k+1)/THREAD_CNT),
                                    std::ref(sample), std::ref(asciiPNGs));
                ts[k].join();
            }
            // fillDist(distances, 0, fCnt, sample, asciiPNGs);
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
            // E = lastMinE;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            lastMinE = E;
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
            std::cout << "best" << std::endl;
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
                //todo: update lines.... I didn't do it!!!!!!!!!!!!!!!!
                std::cout << "goto prev" << std::endl;
                for(int i = 0; i < grCnt; i++) {
                    for(int j = 0; j < gcCnt; j++) {
                        for(int k = 0; k < addedLinesCnt[i][j]; k++) {
                            lines[i][j].pop_back();
                        }
                    }
                }
                for(int i = 0; i < grCnt; i++) {
                    for(int j = 0; j < gcCnt; j++) {
                        for(int k = 0; k < removedLines[i][j].size(); k++) {
                            lines[i][j].push_back(removedLines[i][j][k]);
                        }
                    }
                }
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
                    segs[ssss + toChange.second - 1].second.first -= dx;
                    segs[ssss + toChange.second - 1].second.second -= dy;
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
