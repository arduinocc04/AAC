/**
 * @file egTools.cpp
 * @author Daniel Cho
 * @date 2023.12.1
 * @version 0.0.1
 */
#include <queue>

#include "egTool.hpp"

Dots eg::tool::merge(const Paths & a) {
    Dots p;
    for(int i = 0; i < a.size(); i++)
        for(int j = 0; j < a[i].size(); j++)
            p.push_back(a[i][j]);
    return p;
}

Paths eg::tool::Mat2dToBorders(const Eigen::Tensor<int, 2> & ans, int nbd, int h, int w) {
    const int dx8ccw[] = {-1, -1, 0, 1, 1, 1, 0, -1};
    const int dy8ccw[] = {0, -1, -1, -1, 0, 1, 1, 1};
    Paths borders;
    std::vector<bool> calced(nbd);
    std::vector<std::vector<bool>> used(h);
    for(int i = 0; i < h; i++)
        used[i].resize(w);
    calced[0] = calced[1] = true;

    for(int i = 0; i < h; i++) {
        for(int j = 0; j < w; j++) {
            int cnbd = abs(ans(i, j));
            if(!calced[cnbd]) {
                for(int k = 0; k < h; k++)
                    for(int l = 0; l < w; l++)
                        used[k][l] = false;
                std::queue<Dot> q;

                used[i][j] = true;
                calced[cnbd] = true;
                int last;
                Path res;
                res.push_back(std::make_pair(i, j));
                for(int k = 3; k < 7; k++) { // go below or right
                    int tx = i + dx8ccw[k];
                    int ty = j + dy8ccw[k];
                    if(tx < 0 || tx >= h || ty < 0 || ty >= w)
                        continue;
                    if(abs(ans(tx, ty)) == cnbd) {
                        q.push(std::make_pair(tx, ty));
                        used[tx][ty] = true;
                        last = k;
                        break;
                    }
                }
                while(q.size()) {
                    Dot tmp = q.front();
                    res.push_back(tmp);
                    q.pop();
                    for(int k = 0; k < 8; k++) {
                        int tx = tmp.first + dx8ccw[(k + 3) % 8];
                        int ty = tmp.second + dy8ccw[(k + 3) % 8];
                        if(tx < 0 || tx >= h || ty < 0 || ty >= w || used[tx][ty])
                            continue;
                        if(abs(ans(tx, ty)) == cnbd) {
                            q.push(std::make_pair(tx, ty));
                            used[tx][ty] = true;
                            break;
                        }
                    }
                }

                std::reverse(res.begin(), res.end());// this will change order to clock wise.

                while(q.size()) q.pop();
                bool flag = false;
                for(int k = last + 1; k < 7; k++) {
                    int tx = i + dx8ccw[k];
                    int ty = j + dy8ccw[k];
                    if(tx < 0 || tx >= h || ty < 0 || ty >= w)
                        continue;
                    if(abs(ans(tx, ty)) == cnbd && !used[tx][ty]) {
                        flag = true;
                        q.push(std::make_pair(tx, ty));
                        used[tx][ty] = true;
                        break;
                    }
                }
                if(flag) { 
                    while(q.size()) {
                        Dot tmp = q.front();
                        res.push_back(tmp);
                        q.pop();
                        for(int k = 0; k < 8; k++) {
                            int tx = tmp.first + dx8ccw[(k + 3) % 8];
                            int ty = tmp.second + dy8ccw[(k + 3) % 8];
                            if(tx < 0 || tx >= h || ty < 0 || ty >= w || used[tx][ty])
                                continue;
                            if(abs(ans(tx, ty)) == cnbd) {
                                q.push(std::make_pair(tx, ty));
                                used[tx][ty] = true;
                                break;
                            }
                        }
                    }
                }
                // std::reverse(res.begin(), res.end());// this will change order to counter clock wise.
                borders.push_back(res);
            }
        }
    }
    return borders;
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


Segment eg::tool::clip(const Segment & s, const std::pair<int, int> & lu, const std::pair<int, int> & rd) {
    const std::pair<int, int> real_rd = std::make_pair(rd.first - 1, rd.second - 1); // we need this because [lu, rd)
    const int INSIDE = 0;
    const int LEFT = 1;
    const int RIGHT = 2;
    const int BOTTOM = 4;
    const int TOP = 8;

    int outcode0 = computeOutCode(s.first, lu, real_rd);
    int outcode1 = computeOutCode(s.second, lu, real_rd);
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
            int xmax = real_rd.first;
            int xmin = lu.first;
            int ymax = real_rd.second;
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
                outcode0 = computeOutCode(ans.first, lu, real_rd);
            }
            else {
                ans.second.first = x;
                ans.second.second = y;
                outcode1 = computeOutCode(ans.second, lu, real_rd);
            }
        }
    }
    if(accept)
        return ans;
    return std::make_pair(std::make_pair(-1, -1), std::make_pair(-1, -1));
}
