#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <algorithm>

constexpr int S = 1 << 12; 
constexpr int N = S + 2; // 境界も含めた大きさ
constexpr int ITER = 100;
float dp[2][N*N];
float dp1[2][N][N];
float r = 0.1;
float rr = 1 - 4*r;


int main(){
    double elapsed;
    struct timespec t1, t2;
    omp_set_num_threads(19);

    for(int i = 0; i < N*N; i++){
        if(i < N || i >= N*(N-1) || i % N == 0 || i % N == N - 1) continue;
        dp[0][i] = 1.0;
    }

    clock_gettime(CLOCK_REALTIME, &t1);
    for(int k = 0; k < ITER; k++){
        for(int i = 1; i < N - 1; i++){
            for(int j = 1; j < N - 1; j++){
                int idx = i*N + j;
                int bit = k & 1;
                dp[bit ^ 1][idx] = rr*dp[bit][idx] + r*(dp[bit][idx - N] + dp[bit][idx + N] + dp[bit][idx - 1] + dp[bit][idx + 1]);
            }
        }
    }
    clock_gettime(CLOCK_REALTIME, &t2);

    elapsed = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)*(1.0e-9);
    printf("elapsed without OpenMP : %lf\n", elapsed);

    #pragma omp parallel
    for(int i = 0; i < N*N; i++){
        if(i < N || i >= N*(N-1) || i % N == 0 || i % N == N - 1) continue;
        dp[0][i] = 1.0;
    }

    clock_gettime(CLOCK_REALTIME, &t1);

    for(int k = 0; k < ITER; k++){
        #pragma omp parallel for
        for(int i = 1; i < N - 1; i++){
            for(int j = 1; j < N - 1; j++){
                int idx = i*N + j;
                int bit = k & 1;
                dp[bit ^ 1][idx] = rr*dp[bit][idx] + r*(dp[bit][idx - N] + dp[bit][idx + N] + dp[bit][idx - 1] + dp[bit][idx + 1]);
            }
        }
    }
    clock_gettime(CLOCK_REALTIME, &t2);

    elapsed = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)*(1.0e-9);
    printf("elapsed with OpenMP (1-d array) : %lf\n", elapsed);


    for(int i = 1; i < N - 1; i++){
        for(int j = 1; j < N - 1; j++){
            dp1[0][i][j] = 1.0;
        }
    }

    clock_gettime(CLOCK_REALTIME, &t1);
    for(int k = 0; k < ITER; k++){
        #pragma omp parallel for
        for(int i = 1; i < N - 1; i++){
            for(int j = 1; j < N - 1; j++){
                int bit = k & 1;
                dp1[bit ^ 1][i][j] = rr*dp1[bit][i][j] + r*(dp1[bit][i - 1][j] + dp1[bit][i + 1][j] + dp1[bit][i][j - 1] + dp1[bit][i][j + 1]);
            }
        }
    }
    clock_gettime(CLOCK_REALTIME, &t2);

    elapsed = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)*(1.0e-9);
    printf("elapsed with OpenMP (2-d array) : %lf\n", elapsed);

    // 確認用
/*
    printf("1-d\n");
    for(int i = 1; i <= 10; i++){
        for(int j = 1; j <= 10; j++){
            printf("%f ", dp[ITER & 1][i*N + j]);
        }
        printf("\n");
    }


    printf("2-d\n");
    for(int i = 1; i <= 10; i++){
        for(int j = 1; j <= 10; j++){
            printf("%f ", dp1[ITER & 1][i][j]);
        }
        printf("\n");
    }
*/
}
