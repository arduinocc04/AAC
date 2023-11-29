/**
 * @file trace.cpp
 * @author Daniel Cho
 * @date 2023.11.29
 * @version 0.0.1
 */
#include <vector>
#include <queue>
#include <iostream>

#include "egGeometry.hpp"
#include "egMath.hpp"
#include "egProcessing.hpp"
#include "egTrace.hpp"

#define x first
#define y second

Segments decomposePathGreedy(Path & a) {
    Segments ans;
    int startIdx = 0;
    for(int i = 1; i < a.size(); i++) {
        Segment tmp;
        tmp.first = a[startIdx];
        tmp.second = a[i];
        bool flag = true;
        for(int j = startIdx + 1; j < i; j++) {
            if(eg::geo::distSegDot(tmp, a[j]) > 1.42) {
                flag = false;
                break;
            }
        }
        if(!flag) {
            tmp.second = a[i - 1];
            ans.push_back(tmp);
            startIdx = i;
        }
    }
    ans.push_back(std::make_pair(a[startIdx], a[a.size() - 1]));
    return ans;
}

Segments decomposePathPotrace(Path & a) {
}

Segments eg::trace::decomposePathToSegments(Path & a, int method) {
    if(!a.size())
        throw eg::exceptions::InvalidParameter();
    switch(method) {
        case eg::pathDecomMethod::greedy:
            return decomposePathGreedy(a);
        case eg::pathDecomMethod::potrace:
            return decomposePathPotrace(a);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

