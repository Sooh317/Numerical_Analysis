#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <algorithm>

const int S = 1000; 
const int N = S + 2; // 境界も含めた大きさ
const int ITER = 100;
float dp[N*N], ndp[N*N];
float r = 0.1;
float rr = 1 - 4*r;


int main(){
    double elapsed;
    struct timespec t1, t2;

    for(int i = 0; i < N*N; i++){
        if(i < N || i >= N*(N-1) || i % N == 0 || i % N == N - 1) continue;
        dp[i] = 1.0;
    }

    clock_gettime(CLOCK_REALTIME, &t1);
    for(int _ = 0; _ < ITER; _++){
        for(int i = 1; i < N - 1; i++){
            for(int j = 1; j < N - 1; j++){
                ndp[i*N + j] = rr*dp[i*N + j] + r*(dp[(i-1)*N + j] + dp[(i+1)*N + j] + dp[i*N + j - 1] + dp[i*N + j + 1]);
            }
        }
        std::swap(dp, ndp);
    }
    clock_gettime(CLOCK_REALTIME, &t2);

    elapsed = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)*(1.0e-9);
    printf("elapsed without OpenMP : %lf\n", elapsed);


    #pragma omp parallel
    for(int i = 0; i < N*N; i++){
        if(i < N || i >= N*(N-1) || i % N == 0 || i % N == N - 1) continue;
        dp[i] = 1.0;
    }

    clock_gettime(CLOCK_REALTIME, &t1);
    for(int _ = 0; _ < ITER; _++){
        #pragma omp parallel for
        for(int i = 1; i < N - 1; i++){
            for(int j = 1; j < N - 1; j++){
                ndp[i*N + j] = rr*dp[i*N + j] + r*(dp[(i-1)*N + j] + dp[(i+1)*N + j] + dp[i*N + j - 1] + dp[i*N + j + 1]);
            }
        }
        std::swap(dp, ndp);
    }
    clock_gettime(CLOCK_REALTIME, &t2);

    elapsed = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)*(1.0e-9);
    printf("elapsed with OpenMP : %lf\n", elapsed);
}
