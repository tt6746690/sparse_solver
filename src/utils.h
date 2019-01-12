#ifndef __UTILS_H__
#define __UTILS_H__

#include <chrono>
#include <string>
#include <fstream>
#include <vector>

using time_point_t = std::chrono::high_resolution_clock::time_point;
using sec_t = std::chrono::duration<double>;

inline time_point_t now() {
    return std::chrono::high_resolution_clock::now();
}

inline decltype(auto) diff(time_point_t t1, time_point_t t2) {
    return std::chrono::duration_cast<sec_t>(t2-t1).count();
}

inline void filetovec(const std::string& file, std::vector<double>& v) {
    v.clear();
    std::ifstream in(file);
    std::string l;
    while(std::getline(in, l))
        v.push_back(std::stod(l));
};

inline void vectofile(const std::string& file, const std::vector<double>& v) {
    std::ofstream out(file);
    for (auto x : v) 
        out << x << '\n';
};


#endif