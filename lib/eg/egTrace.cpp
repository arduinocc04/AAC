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

/**
 * @todo Implement correctly. Segment a[a.size() - 1] -- a[0] may needed.
 */
Segments approxPathGreedy(const Path & a) {
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

Path eg::trace::approxPath(const Path & a) {
    Path ans;
    int startIdx = 0;
    const double THRESHOLD = 4;
    if(a.size() <= 1) return ans;
    for(int i = 1; i < a.size(); i++) {
        bool flag = true;
        Segment tmp = std::make_pair(a[startIdx], a[i]);
        for(int j = startIdx + 1; j < i; j++) {
            if(eg::geo::distSegDot(tmp, a[j]) > 1.42) {
                flag = false;
                break;
            }
        }
        if(!flag) {
            tmp.second = a[i - 1];
            if(eg::geo::euclideDist(tmp.first, tmp.second) >= THRESHOLD)
                ans.push_back(a[startIdx]);
            startIdx = i;
        }
    }
    if(eg::geo::euclideDist(a[startIdx], a[a.size() - 1]) >= THRESHOLD) {
        ans.push_back(a[startIdx]);
        ans.push_back(a[a.size() - 1]);
    }
    if(eg::geo::euclideDist(a[a.size() - 1], a[0]) < 1.5) {
        ans.pop_back();
        ans.push_back(a[0]);
    }
    return ans;
}

Segments eg::trace::approxPathToSegments(const Path & a, int method) {
    if(!a.size())
        throw eg::exceptions::InvalidParameter();
    switch(method) {
        case eg::pathDecomMethod::greedy:
            return approxPathGreedy(a);
        default:
            throw eg::exceptions::InvalidParameter();
    }
}

