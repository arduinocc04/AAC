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

int computeOutCode(const Dot & p, const Dot & lu, const Dot & rd) {
    const int INSIDE = 0;
    const int LEFT = 1;
    const int RIGHT = 2;
    const int BOTTOM = 4;
    const int TOP = 8;

    int code = INSIDE;

    if(p.first < lu.first)
        code |= LEFT;
    else if(p.first > rd.first)
        code |= RIGHT;
    if(p.second < lu.second)
        code |= BOTTOM;
    else if(p.second > rd.second)
        code |= TOP;

    return code;
}

/**
 * @breif do line clipping using Cohen-Sutherland algorithm
 * @see https://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland_algorithm
 */
Segment eg::tool::clip(const Segment & s, const Dot & lu, const Dot & rd) {
    const int INSIDE = 0;
    const int LEFT = 1;
    const int RIGHT = 2;
    const int BOTTOM = 4;
    const int TOP = 8;

    int outcode0 = computeOutCode(s.first, lu, rd);
    int outcode1 = computeOutCode(s.second, lu, rd);
    bool accept = false;

    Segment ans = s;

    while(true) {
        if(!(outcode0 | outcode1)) {
            accept = true;
            break;
        }
        else if(outcode0 & outcode1)
            break;
        else {
            int outcodeOut = (outcode1 > outcode0)?outcode1:outcode0;
            int x, y;
            int x0 = s.first.first;
            int y0 = s.first.second;
            int x1 = s.second.first;
            int y1 = s.second.second;
            int xmax = rd.first;
            int xmin = lu.first;
            int ymax = rd.second;
            int ymin = lu.second;
            if(outcodeOut & TOP) {
                x = round(x0 + (x1 - x0)*(ymax - y0)/(y1 - y0));
                y = ymax;
            }
            else if(outcodeOut & BOTTOM) {
                x = round(x0 + (x1 - x0)*(ymin - y0)/(y1 - y0));
                y = ymin;
            }
            else if(outcodeOut & RIGHT) {
                y = round(y0 + (y1 - y0)*(xmax - x0)/(x1 - x0));
                x = xmax;
            }
            else if(outcodeOut & LEFT) {
                y = round(y0 + (y1 - y0)*(xmin - x0)/(x1 - x0));
                x = xmin;
            }
            if(outcodeOut == outcode0) {
                ans.first.first = x;
                ans.first.second = y;
                outcode0 = computeOutCode(ans.first, lu, rd);
            }
            else {
                ans.second.first = x;
                ans.second.second = y;
                outcode1 = computeOutCode(ans.second, lu, rd);
            }
        }
    }
    if(accept)
        return ans;
    return std::make_pair(std::make_pair(-1, -1), std::make_pair(-1, -1));
}
