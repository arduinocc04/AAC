#include <iostream>
#include <vector>
#include <algorithm>

#include "egGeometry.hpp"

std::vector<std::pair<int, int>> dots;

int main() {
    dots.push_back(std::make_pair(0, 0));
    dots.push_back(std::make_pair(0, 4));
    dots.push_back(std::make_pair(2, 2));
    dots.push_back(std::make_pair(4, 0));
    dots.push_back(std::make_pair(4, 4));
    dots.push_back(std::make_pair(1, 1));
    dots.push_back(std::make_pair(2, 5));

    std::sort(dots.begin(), dots.end());

    std::vector<std::pair<int, int>> hull = eg::geo::getConvexHull(dots);

    for(int i = 0; i < hull.size(); i++) {
        std::cout << hull[i].first << " " << hull[i].second << std::endl;
    }
}
