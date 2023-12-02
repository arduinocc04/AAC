/**
 * @file egTools.cpp
 * @author Daniel Cho
 * @date 2023.12.1
 * @version 0.0.1
 */
#include "egTool.hpp"

Dots eg::tool::merge(const Paths & a) {
    Dots p;
    for(int i = 0; i < a.size(); i++)
        for(int j = 0; j < a[i].size(); j++)
            p.push_back(a[i][j]);
    return p;
}


