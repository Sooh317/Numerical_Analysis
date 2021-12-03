#include <random>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <numeric>
#include <algorithm>

template<typename T>
struct Matrix{
    std::vector<std::vector<T>> A;

    Matrix(int _n, int _m):A(_n, std::vector<T>(_m, 0)){}
    Matrix(std::vector<std::vector<T>> &_A){A = _A;}

    int height()const{return A.size();}
    int width()const{return A[0].size();}

    std::vector<T> multiply_vector(const std::vector<T>& b){
        int n = height(), m = width();
        std::vector<T> res(n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                res[i] += A[i][j] * b[j];
            }
        }
        return res;
    }

    std::vector<T> solve_linear_system(const std::vector<T>& b1){
        int n = b1.size();
        auto [LU, order] = LU_decomposition();
        std::vector<T> b(n);
        for(int i = 0; i < n; i++) b[i] = b1[order[i]];
        std::vector<T> y(n), x(n);
        for(int i = 0; i < n; i++){
            y[i] = b[i];
            for(int j = 0; j < i; j++){
                y[i] -= LU[i][j] * y[j];
            }
        }
        for(int i = n - 1; i >= 0; i--){
            x[i] = y[i];
            for(int j = n - 1; j > i; j--){
                x[i] -= LU[i][j] * x[j];
            }
            x[i] /= LU[i][i];
        }
        return x;
    }

    // [LU, 置換]を返す
    std::pair<Matrix, std::vector<int>> LU_decomposition(){ // Lの対角部が1、LUの狭義下三角がL
        assert(height() == width());
        Matrix LU = Matrix(*this);

        int n = height();
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);

        for(int k = 0; k < n; k++){
            for(int i = k + 1; i < n; i++){
                if(std::abs(LU[i][k]) > std::abs(LU[k][k])){
                    std::swap(LU[i], LU[k]);
                    std::swap(order[i], order[k]);
                }
            }
            for(int i = k + 1; i < n; i++){
                LU[i][k] /= LU[k][k];
            }
            for(int i = k + 1; i < n; i++){
                for(int j = k + 1; j < n; j++){
                    LU[j][i] -= LU[j][k] * LU[k][i];
                }
            }
        }
        return std::make_pair(LU, order);
    }

    Matrix &operator*=(const Matrix &B){
        int n = height(), m = width(), l = B.width();
        assert(m == B.height());
        std::vector<std::vector<T>> C(n, std::vector<T>(m, 0));
        for(int i = 0; i < n; i++) for(int j = 0; j < l; j++){
            for(int k = 0; k < m; k++){
                C[i][j] += ((*this)[i][k] * B[k][j]);
            }
        }
        A.swap(C);
        return (*this);
    }

    Matrix operator*(const Matrix &B)const{return (Matrix(*this) *= B);}
    inline std::vector<T> &operator[](int k){ return A.at(k);}
    inline const std::vector<T> &operator[](int k)const{ return A.at(k);}
};

using namespace std;

void check(Matrix<double> LU, vector<int> ord){
    int n = LU.height();

    Matrix<double> L(n, n), U(n, n);
    for(int i = 0; i < n; i++) L[i][i] = 1;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i > j) L[i][j] = LU[i][j];
            else U[i][j] = LU[i][j];
        }
    }

    cout << "L*U = \n";
    auto B = L*U;

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            cout << B[i][j] << " ";
        }
        cout << endl;
    }
}

void test(){
    std::random_device rnd;
    std::default_random_engine eng(rnd());
    std::uniform_real_distribution<double> dist(0.0, 1000.0);

    for(int n = 1; n <= 10; n++){
        vector<vector<double>> hilbert(n, vector<double>(n));
        for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) hilbert[i][j] = 1.0 / (i + j + 1);
        Matrix<double> A(hilbert);

        vector<double> x(n);
        for(int i = 0; i < n; i++) x[i] = dist(eng);
        vector<double> b = A.multiply_vector(x); // Ax = b
        vector<double> x1 = A.solve_linear_system(b);

        vector<double> error(n);

        cout << "error when n = " << n << endl;
        for(int i = 0; i < n; i++){
            error[i] = x1[i] - x[i];
            cout << error[i] << endl;
        }
        cout << endl;
    }
}

int main(){
    cout << fixed << setprecision(20);
    test();
}
