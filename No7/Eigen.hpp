#pragma once
#include <random>
#include "CG_algorithm.hpp"

namespace Eigen{

const double EPS = 1e-5;
const int MAX = 100;

std::pair<double, std::vector<double>> PowerIteration(SparseMatrix& csr, bool inv = false){
    std::vector<double> x(csr.m); // 注意 : 初期値ベクトルが固有ベクトルになってる場合は挙動が違う
    std::uniform_real_distribution<double> mt(0.0, 10000.0);
    std::default_random_engine eng;
    for(int i = 0; i < csr.m; i++) x[i] = mt(eng);
    double lambda;
    for(int i = 0; i < MAX; i++){
        Calculation::regularlize(x);
        std::vector<double> y;
        if(inv) y = CG_Algorithm::solve(csr, x); // 逆行列が与えられた場合
        else y = Calculation::multiply(csr, x);
        // because of regularlization, dot(x, x) = 1 
        lambda = Calculation::dot(y, x); 
        double tmp = Calculation::dot(y, y);
        if(std::abs(tmp - lambda * lambda) < EPS * EPS) break;
        std::swap(x, y);
    }
    return std::make_pair(lambda, x);
}

// mu = 0 として絶対値最小を求める
// 1/lambda が返ることに注意
std::pair<double, std::vector<double>> InverseIteration(SparseMatrix& csr){
    return PowerIteration(csr, true);
}

}
