#include <iostream>
#include <string>
#include <filesystem>
#include "egLoader.hpp"
#include "egMath.hpp"
#include "egProcessing.hpp"
#include <cmath>
#include <thread>

using namespace eg::imgproc;

#define THREAD_CNT 8
#define MAX_ASCIIS_CNT 15'000

int p[MAX_ASCIIS_CNT];
std::string names[MAX_ASCIIS_CNT];

int find(int n) {
    if(n == p[n]) return n;
    return p[n] = find(p[n]);
}

void merge(int n1, int n2) {
    n1 = find(n1);
    n2 = find(n2);
    p[n1] = n2;
}

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

Mat2d * getAllImages(std::string path, int fCnt) {
    Mat2d * ans = new Mat2d[fCnt];

    int i = 0;
    eg::PNG t;
    for(const auto & entry : fs::directory_iterator(path)) {
        t.openImage(entry.path());
        names[i] = entry.path();
        Image ti = t.copy();
        ans[i] = cvtGray(ti, eg::grayCvtMethod::mean);
        ans[i] = binary(ans[i], 70);
        i += 1;
    }
    std::cout << std::endl;

    return ans;
}

void fillRmse(Eigen::Tensor<double, 2> & rmses, int iStart, int iEnd, int jEnd, Mat2d * pngs) {
    for(int i = 0; i < iEnd; i++)
        for(int j = i + 1; j < jEnd; j++)
            rmses(i, j) = eg::math::rmse(pngs[i], pngs[j]);
}

int main(int argc, char * argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);

    if(argc != 3) {
        std::cout << "Use program properly!" << std::endl;
        return -1;
    }

    std::string path = argv[1];
    int fCnt = getFileCount(path);

    std::cout << "Handling " << fCnt << " images." << std::endl;

    Mat2d * pngs = getAllImages(path, fCnt);

    Eigen::Tensor<double, 2> rmses(fCnt, fCnt);

    const double thres = std::stod(argv[2]);

    std::thread ts[THREAD_CNT];
    for(int i = 0; i < THREAD_CNT; i++) {
        ts[i] = std::thread(fillRmse, std::ref(rmses), fCnt*(i/THREAD_CNT), fCnt*((i+1)/THREAD_CNT), fCnt, std::ref(pngs));
        ts[i].join();
    }
    std::cout << std::endl;

    std::cout << "MAX: " << rmses.maximum() << " MIN: " << rmses.minimum() << std::endl;

    for(int thres = 200; thres < 4000; thres += 200) {
        for(int i = 0; i < MAX_ASCIIS_CNT; i++) p[i] = i;
        for(int i = 0; i < fCnt; i++)
            for(int j = i + 1; j < fCnt; j++)
                if(rmses(i, j) < thres) merge(i, j);

        int c = 0;
        std::cout << "========" << thres << "========" << std::endl;
        for(int i = 0; i < fCnt; i++) {
            if(p[i] == i) {
                std::cout << names[i] << std::endl;
                c++;
            }
        }

        std::cout << c << std::endl;
    }
    delete[] pngs;
}
