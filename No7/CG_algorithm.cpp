#include "distr/load.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <functional>

namespace CG_Algorithm{

const double EPS = 1e-5;

template<typename T>
std::vector<T> multiply(const SparseMatrix& A, const std::vector<T>& x){
    assert(A.n == (int)x.size());
    std::vector<T> y((int)x.size());
    for(int i = 0; i < A.m; i++){
        for(int j = A.rowptr[i]; j < A.rowptr[i + 1]; j++){
            y[i] += A.val[j] * x[A.colind[j]];
        }
    }
    return y;
}

template<typename T>
T dot(const std::vector<T>& x, const std::vector<T>& y){
    assert(x.size() == y.size());
    T ans = 0;
    for(int i = 0; i < (int)x.size(); i++) ans += x[i] * y[i];
    return ans;
}

template<typename T>
double abs(const std::vector<T>& x){
    double ans = 0;
    for(int i = 0; i < (int)x.size(); i++) ans += x[i] * x[i];
    return std::sqrt(ans);
}

// return x s.t. Ax = b
template<typename T>
std::vector<T> solve(SparseMatrix& A, std::vector<T>& b){
    std::vector<T> x((int)b.size(), 0);
    std::vector<T> r = b;
    std::vector<T> p = b;
    std::vector<T> nr((int)b.size());

    while(1){
        std::vector<T> y = multiply(A, p);
        double alpha = dot(r, r) / dot(p, y);
        for(int i = 0; i < (int)x.size(); i++) x[i] += alpha*p[i];
        for(int i = 0; i < (int)r.size(); i++) nr[i] = r[i] - alpha*y[i];
        if(abs(r) < EPS) break;
        double beta = dot(nr, nr) / dot(r, r);
        std::swap(nr, r);
        for(int i = 0; i < (int)p.size(); i++) p[i] = r[i] + beta*p[i];
    }
    return x;
}
}

const double pi = std::acos(-1);
double f(double x, double y){return 0;}
double up(double y){return 0;}
double left(double x){return std::sin(pi*x);}
double down(double y){return 0;}
double right(double x){return std::min(x, 1 - x);}

// make vector b s.t. Ax = b for poisson's equation
std::vector<double> make_vector(int grid_x, int grid_y, double delta_x, double delta_y){
    std::vector<double> b(grid_x * grid_y);
    double cur_x = 0;
    double cur_y = 0;
    for(int i = 0; i < grid_x; i++){
        cur_x += delta_x;
        cur_y = 0;
        for(int j = 0; j < grid_y; j++){
            cur_y += delta_y;

            double val = f(cur_x, cur_y);
            if(i == 0) val += up(cur_y);
            if(j == 0) val += left(cur_x);
            if(i == grid_x - 1) val += down(cur_y);
            if(j == grid_y - 1) val += right(cur_x);

            b[i*grid_y + j] = -val;
        }
    }
    return b;
}

using namespace std;

int main(){
    const double GRID_SIDE = 100; // 1辺
    const double delta = 1.0 / GRID_SIDE;

    SparseMatrix csr;
    csr.stencil5(GRID_SIDE);

    auto b = make_vector(GRID_SIDE, GRID_SIDE, delta, delta);
    auto x = CG_Algorithm::solve(csr, b);

    cout << GRID_SIDE << " " << GRID_SIDE << '\n';

    for(int i = 0; i < GRID_SIDE; i++){
        for(int j = 0; j < GRID_SIDE; j++){
            cout << x[i*GRID_SIDE + j] << " ";
        }
        cout << '\n';        
    }
}

