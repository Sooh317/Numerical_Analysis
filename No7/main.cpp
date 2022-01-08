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
    SparseMatrix csr("Kuu.mtx");

    auto [l, x] = Eigen::PowerIteration(csr);
    auto [s, y] = Eigen::InverseIteration(csr);
    auto z = Calculation::multiply(csr, x);
    cout << "absmax: " << l << '\n';
    cout << "absmin: " << 1/s << '\n';
    cout << "cond:   " << abs(l*s) << '\n';
    // result 
    // absmax: 53.8028
    // absmin: 0.00343204
    // cond:   15676.6
}

void kadai_C(){
    vector<vector<double>> hilbert(6, vector<double>(6));
    for(int i = 0; i < 6; i++) for(int j = 0; j < 6; j++) hilbert[i][j] = 1.0 / (i + j + 1);
    auto res = Eigen::QR(hilbert);
    for(int i = 0; i < res.size(); i++) cout << res[i] << endl;
    // result
    // [1.6189, 0.242361, 0.0163215, 0.000615748, 1.25708e-05, 1.0828e-07]
}

int main(){
    // kadai_A();
    kadai_B();
    // kadai_C();
}
