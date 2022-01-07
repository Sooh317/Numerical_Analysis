#include "CG_algorithm.hpp"
#include "Eigen.hpp"
using namespace std;

void kadai_A(){
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

void kadai_B(){
    SparseMatrix csr;
    csr.stencil5(2);

    auto [l, x] = Eigen::PowerIteration(csr);
    auto [s, y] = Eigen::InverseIteration(csr);
    auto z = Calculation::multiply(csr, x);
    cout << l << endl;
    for(int i = 0; i < (int)x.size(); i++){
        cout << x[i] << ' ';
    }
    cout << endl;
    cout << 1/s << endl;
    for(int i = 0; i < (int)y.size(); i++){
        cout << y[i] << ' ';
    }
    cout << endl;
}

int main(){
    // kadai_A();
    kadai_B();
}
