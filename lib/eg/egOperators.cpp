/**
 * @file egOperators.cpp
 * @author Daniel Cho
 * @date 2023.11.29
 * @version 0.0.1
 */
#ifndef __EGTYPES_H
#include "egTypes.hpp"
#endif
#define __EGTYPES_H
namespace eg {
Dot operator-(const Dot & a, const Dot & b) {
    Dot res;
    res.first = a.first - b.first;
    res.second = a.second - b.second;
    return res;
}

Dot operator+(const Dot & a, const Dot & b) {
    Dot res;
    res.first = a.first + b.first;
    res.second = a.second + b.second; return res;
}

Dot operator+(const Dot & a) {
    return a;
}

Dot operator-(const Dot & a) {
    Dot res;
    res.first = -a.first;
    res.second = -a.second;
    return res;
}
}
