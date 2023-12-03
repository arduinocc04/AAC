/**
 * @file trace.cpp
 * @author Daniel Cho
 * @date 2023.11.29
 * @version 0.0.1
 */
#include <vector>
#include <queue>

#include "egGeometry.hpp"
#include "egMath.hpp"
#include "egProcessing.hpp"
#include "egTrace.hpp"

#define x first
#define y second

Segments decomposePathGreedy(const Path & a) {
    Segments ans;
    int startIdx = 0;
    const double THRESHOLD = 5;
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
            if(eg::geo::euclideDist(tmp.first, tmp.second) >= THRESHOLD)
            ans.push_back(tmp);
            startIdx = i;
        }
    }
    if(eg::geo::euclideDist(a[startIdx], a[a.size() - 1]) >= THRESHOLD)
        ans.push_back(std::make_pair(a[startIdx], a[a.size() - 1]));
    return ans;
}

Segments eg::trace::decomposePathToSegments(const Path & a, int method) {
    if(!a.size())
        throw eg::exceptions::InvalidParameter();
    switch(method) {
        case eg::pathDecomMethod::greedy:
            return decomposePathGreedy(a);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

