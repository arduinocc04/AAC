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

#include "egLoader.hpp"
#include "egProcessing.hpp"
#include "egMath.hpp"
#include "egTrace.hpp"

#define INF 987654321
#define THREAD_CNT 15

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
    for(int i = 0; i < contours.first.size(); i++) {
        if(contours.first[i].size() == 0) continue;
        Segments tmp = eg::trace::decomposePathToSegments(contours.first[i], eg::pathDecomMethod::greedy);
        for(int j = 0; j < tmp.size(); j++)
            segs.push_back(tmp[j]);
    }
    Segments original = segs;

    inputImage.divideImageByLength(asciih, asciiw);
    int gcCnt = inputImage.getGridColCnt();
    int grCnt = inputImage.getGridRowCnt();

    std::cout << gcCnt << "x" << grCnt << std::endl;

    char devnull[10];

    std::mt19937 mt(time(nullptr));

    double lastMinE = INF;
    double prevE = INF;

    int co = 0;
    double ta = -1; // t_a
    int c = 0; // iteration index
    Mat2d zeros(inputImage.info.height, inputImage.info.width);
    zeros.setConstant(0);
    double * distances = new double[fCnt];
    while(true) {
        int K = grCnt*gcCnt;
        int toChange = mt() % segs.size();
        Segment before = segs[toChange];
        double E;
        if(c) {
            segs[toChange].first.first += (mt() % (2*asciih)) - asciih;
            segs[toChange].first.second += (mt() % (2*asciiw)) - asciiw;
            segs[toChange].second.first += (mt() % (2*asciih)) - asciih;
            segs[toChange].second.second += (mt() % (2*asciiw)) - asciiw;
            E = eg::math::calcDeform(before, segs[toChange], segs);
        }
        else
            E = 1;

        double sum = 0;
        Mat2d drawn = drawSegments(zeros, segs, 1);
        Mat2d tttttmp = 255*drawn;
        Image tmp = mat2dToImage(tttttmp);
        std::cout << "TMP IMAGE GEN" << std::endl;
        inputImage.setImage(tmp);
        inputImage.saveImage("asdf.png");
        tmp = mat2dToImage(drawn);
        inputImage.setImage(tmp);
        for(int i = 0; i < grCnt; i++) {
            for(int j = 0; j < gcCnt; j++) {
                Image raw = inputImage.getImageAtGrid(i, j);
                Mat2d sample = cvtGray(raw, eg::grayCvtMethod::mean);

                Eigen::Tensor<double, 0> tmp = sample.sum();
                if(tmp(0) < 2) {
                    std::cout << " ";
                    --K;
                    continue;
                }

                sample = inflate(sample, asciih, asciiw);
                std::thread ts[THREAD_CNT];
                for(int i = 0; i < THREAD_CNT; i++) {
                    ts[i] = std::thread(fillDist, std::ref(distances), fCnt*(i/THREAD_CNT), fCnt*((i+1)/THREAD_CNT), std::ref(sample), std::ref(asciiPNGs));
                    ts[i].join();
                }
                double minVal = INF;
                int minIndex = -1;
                for(int k = 0; k < fCnt; k++) {
                    if(distances[k] < minVal) {
                        minVal = distances[k];
                        minIndex = k;
                    }
                }
                sum += minVal;
                E += minVal;
                std::string ascii = getAsciiFromPath(names[minIndex]);
                std::cout << ascii;
            }
            std::cout << std::endl;
        }
        if(c == 0)
            ta = sum/(grCnt*gcCnt);

        std::cout << "E: " << E << std::endl;
        if(lastMinE > E) {
            E = lastMinE;
            co = 0;
        }
        else {
            if(++co >= 50) break;
            double delta = abs(E - prevE);
            double t = 0.2*ta*std::pow(c, 0.997);
            double Pr = std::exp(-delta/t);
            double ran = (double)mt()/mt.max();
            if(Pr >= ran)
                segs[toChange] = before;
        }
        prevE = E;
        c++;
    }
    delete[] distances;
    delete[] asciiPNGs;
    delete[] names;
}
