/**
 * @file extractCenterline.cpp
 * @author Daniel Cho
 * @version 0.0.1
 */
#include <iostream>
#include <vector>

#include "egLoader.hpp"
#include "egProcessing.hpp"

using namespace eg::imgproc;

std::vector<int> p;

int dx[] = {0, 1, 1, 1, 0, -1, -1, -1}, dy[] = {1, 1, 0, -1, -1, -1, 0, 1};

int find(int n) {
    if(n == p[n]) return n;
    return p[n] = find(p[n]);
}

void merge(int n1, int n2) {
    n1 = find(n1);
    n2 = find(n2);
    p[n1] = n2;
}

int main(int argc, char * argv[]) {
    if(argc != 2) {
        std::cout << "Use program properly!" << std::endl;
        return -1;
    }

    std::string inputPath = argv[1];

    eg::PNG png;
    png.openImage(inputPath);
    Image i = png.copy();
    Mat2d t = cvtGray(i, eg::grayCvtMethod::mean);
    i = mat2dToImage(t);
    png.setImage(i);
    png.saveImage("grrray-" + inputPath);
    t = getEdge(t, eg::edgeDetectMethod::gradient);
    i = mat2dToImage(t);
    png.setImage(i);
    png.saveImage("eddge-" + inputPath);
    for(int i = 0; i < 10; i++)
        t = blur(t, eg::blurMethod::gaussian);
    i = mat2dToImage(t);
    png.setImage(i);
    png.saveImage("tmp-" + inputPath);
    t = extractCenterline(t, eg::centerlineMethod::grassfire);

    p.resize(png.info.width*png.info.height);
    for(int i = 0; i < p.size(); i++)
        p[i] = i;
    Mat2d m = t;
    std::vector<int> crit(png.info.width*png.info.height);
    for(int i = 0; i < crit.size(); i++) crit[i] = 0;

    for(int i = 0; i < png.info.height; i++) {
        for(int j = 0; j < png.info.width; j++) {
            if(m(i, j) == 0) continue;
            for(int di = 0; di < 8; di++) {
                for(int dj = 0; dj < 8; dj++) {
                    int tmpI = i + dx[di], tmpJ = j + dy[dj];
                    if(tmpI < 0 || tmpI >= png.info.height || tmpJ < 0 || tmpJ >= png.info.width) continue;
                    if(m(tmpI, tmpJ) != 0) merge(tmpI*png.info.width + tmpJ, i*png.info.width + j);
                }
            }
        }
    }
    for(int i = 0; i < png.info.height; i++) {
        for(int j = 0; j < png.info.width; j++) {
            int t = find(i*png.info.width + j);
            std::cout << i << " " << j << " " << t << " " << crit[t] << " " << m(i, j) << std::endl;
            if(crit[t] < m(i, j))
                crit[t] = m(i, j);
        }
    }

    for(int i = 0; i < png.info.height; i++) {
        for(int j = 0; j < png.info.width; j++) {
            if(!m(i, j) || m(i, j) + 2 < crit[find(i*png.info.width + j)])
                m(i, j) = 0;
            else
                m(i, j) = 255;
        }
    }
    i = mat2dToImage(m);
    png.setImage(i);
    png.saveImage("extCen-" + inputPath);
}
