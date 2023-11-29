/**
 * @file main.cpp
 * @author Daniel cho
 * @date 2023.11.17
 * @version 0.0.1
 */
#include <iostream>
#include <filesystem>

#include "egLoader.hpp"
#include "egProcessing.hpp"
#include "egMath.hpp"

#define DIST_METHOD 2 // 0: RMSE
#define PRINT_INPUT_IMAGE 0

#define INF 987654321

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
        ans[i] = binary(ans[i], 70);
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
    int asciih = asciiPNGs[0].dimensions()[0];
    int asciiw = asciiPNGs[0].dimensions()[1];

    eg::PNG inputImage;
    std::cout << "Opening input image" << std::endl;
    inputImage.openImage(inputImagePath);
    Image i = inputImage.copy();
    std::cout << "Converting input image gray" << std::endl;
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    std::cout << "Getting Edge of input image" << std::endl;
    t = getEdge(t, eg::edgeDetectMethod::gradient);
    //t = markOutlier(t, 10);
    for(int i = 0; i < 10; i++) {
        std::cout << i + 1 << "/10 blurring input image" << std::endl;
        t = blur(t, eg::blurMethod::gaussian);
    }
    std::cout << "Get Binary of input image" << std::endl;
    t = binary(t, 10);
    //t = getContours(t, -1);
    //t = markOutlier(t, 1);
    i = mat2dToImage(t);
    inputImage.setImage(i);
    inputImage.divideImageByLength(asciih, asciiw);

    int gcCnt = inputImage.getGridColCnt();
    int grCnt = inputImage.getGridRowCnt();

    std::cout << gcCnt << "x" << grCnt << std::endl;

    char null[10];

    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            Image raw = inputImage.getImageAtGrid(i, j);
            Mat2d sample = cvtGray(raw, eg::grayCvtMethod::mean);
            //sample = getEdge(sample, eg::edgeDetectMethod::gradient);
            //sample = getContours(sample, -1);
            //sample = markOutlier(sample, 1);

            if(PRINT_INPUT_IMAGE) {
                std::cout << "===Printing input image at grid " << i << " " << j << std::endl;
                print(sample);
                std::cin >> null;
            }

            Eigen::Tensor<double, 0> tmp = sample.sum();
            if(tmp(0) < asciih*asciiw/10) {
                std::cout << " ";
                continue;
            }
            double minVal = INF;
            int minIndex = -1;
            for(int k = 0; k < fCnt; k++) {
                Mat2d tmp = inflate(sample, asciih, asciiw);
                double dist;
                switch(DIST_METHOD) {
                    case 0:
                        dist = eg::math::compareMat2d(tmp, asciiPNGs[k], eg::matCmpMethod::rmse);
                        break;
                    case 1:
                        dist = eg::math::compareMat2d(tmp, asciiPNGs[k], eg::matCmpMethod::shape);
                        break;
                    case 2:
                        dist = eg::math::compareMat2d(tmp, asciiPNGs[k], eg::matCmpMethod::logpolar);
                        break;
                    default:
                        throw eg::exceptions::InvalidParameter();
                }
                if(dist < minVal) {
                    minVal = dist;
                    minIndex = k;
                }
            }
            std::string ascii = getAsciiFromPath(names[minIndex]);
            std::cout << ascii;
        }
        std::cout << std::endl;
    }
    delete[] asciiPNGs;
    delete[] names;
}
