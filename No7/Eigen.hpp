#pragma once
#include <random>
#include "CG_algorithm.hpp"

namespace Eigen{

const double EPS = 1e-6;
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

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> QR_decomposition(std::vector<std::vector<double>>& A){
    int n = A.size();
    std::vector<std::vector<double>> R = A;
    std::vector<std::vector<double>> Q = Calculation::identity(n);
    for(int k = 0; k < n - 1; k++){
        double tmp = 0.0;
        for(int i = k; i < n; i++) tmp += R[i][k] * R[i][k];
        tmp = std::sqrt(tmp);
        if(R[k][k] > 0) tmp = -tmp;

        std::vector<double> x(n);
        for(int i = k; i < n; i++) x[i] = R[i][k] - (i == k ? tmp : 0.0);
        Calculation::regularlize(x);

        // Householder Matrix
        auto H = Calculation::identity(n);
        for(int i = k; i < n; i++) for(int j = k; j < n; j++) H[i][j] -= 2*x[i]*x[j];
        R = Calculation::matdot(H, R);
        Q = Calculation::matdot(Q, Calculation::transpose(H));
    }
    return std::make_pair(Q, R);
}

std::vector<double> QR(std::vector<std::vector<double>> A){
    int n = A.size();
    std::vector<double> lambda(n);
    for(int i = 0; i < MAX; i++){
        auto[Q, R] = QR_decomposition(A);
        A = Calculation::matdot(R, Q);
        double e = 0;
        for(int j = 0; j < n; j++) e += std::abs(lambda[j] - A[j][j]), lambda[j] = A[j][j];
        if(e / n < EPS) break;
    }
    return lambda;
}

}
