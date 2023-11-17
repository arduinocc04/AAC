/**
 * @file main.cpp
 * @author Daniel cho
 * @date 2023.11.17
 * @version 0.0.1
 */
#include <iostream>
#include <filesystem>

#include "PNGLoader.hpp"

#define DIST_METHOD 0 // 0: RMSE
#define PRINT_INPUT_IMAGE 0
namespace fs = std::filesystem;

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

eg::PNG * getAllImages(std::string path, int fCnt) {
    eg::PNG * ans = new eg::PNG[fCnt];

    int i = 0;
    for(const auto & entry : fs::directory_iterator(path)) {
        ans[i].openImage(entry.path());
        ans[i].cvtGray(eg::grayCvtMethod::mean);
        ans[i].binary(70);
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
    eg::PNG * asciiPNGs = getAllImages(asciiImagesPath, fCnt);
    std::cout << "Opened ascii images" << std::endl;
    int asciih = asciiPNGs[0].getPlayground()->dimensions()[0];
    int asciiw = asciiPNGs[0].getPlayground()->dimensions()[1];

    eg::PNG inputImage;
    std::cout << "Opening input image" << std::endl;
    inputImage.openImage(inputImagePath);
    std::cout << "Converting input image gray" << std::endl;
    inputImage.cvtGray(eg::grayCvtMethod::mean);
    std::cout << "Getting Edge of input image" << std::endl;
    inputImage.getEdge(eg::edgeDetectMethod::gradient);
    for(int i = 0; i < 10; i++) {
        std::cout << i + 1 << "/10 blurring input image" << std::endl;
        inputImage.blur(eg::blurMethod::gaussian);
    }

    std::cout << "Get Binary of input image" << std::endl;
    inputImage.binary(10);
    inputImage.dividePlaygroundByLength(asciih, asciiw);

	if(PRINT_INPUT_IMAGE)
		print(*inputImage.getPlayground());

    int gcCnt = inputImage.getGridColCnt();
    int grCnt = inputImage.getGridRowCnt();

    std::cout << gcCnt << "x" << grCnt << std::endl;

    for(int i = 0; i < grCnt; i++) {
        for(int j = 0; j < gcCnt; j++) {
            Eigen::Tensor<double, 2> sample = inputImage.getPlaygroundAtGrid(i, j);

			Eigen::Tensor<double, 0> tmp = sample.sum();
			if(tmp(0) < 255*5) {
				std::cout << " ";
				continue;
			}
            double minVal = 987654321;
            int minIndex = -1;
            for(int k = 0; k < fCnt; k++) {
                Eigen::Tensor<double, 2> tmp = eg::math::inflate(sample, asciih, asciiw);
                double dist;
                switch(DIST_METHOD) {
                    case 0:
                        dist = eg::math::rmse(tmp, *(asciiPNGs[k].getPlayground()));
                        break;
                    default:
                        throw eg::exceptions::InvalidParameter();
				}
                if(dist < minVal) {
                    minVal = dist;
                    minIndex = k;
                }
            }
            std::string ascii = getAsciiFromPath(asciiPNGs[minIndex].getInputPath());
            std::cout << ascii;
        }
        std::cout << std::endl;
    }
    delete[] asciiPNGs;
}
