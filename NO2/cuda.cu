#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

const int ITER = 100;
const int S = 10;

double elapsed(){
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	return ts.tv_sec + ts.tv_nsec*1.0e-9;
}

__global__
void kernel(int n, float r, float* a){
    int NN = n*n;
    float rr = 1 - 4*r;
    for(int i = 0; i < ITER; i++){
        float *dp = a + (i & 1)*NN, *ndp = a + (1 - (i & 1))*NN;
	    int idx = blockIdx.x * blockDim.x + threadIdx.x;
	    while(idx < NN){
            bool dame = (idx - n < 0 || idx % n == 0 || idx % n == n - 1 || idx + n >= NN);
            if(!dame) ndp[idx] = rr*dp[idx] + r*(dp[idx - n] + dp[idx + n] + dp[idx - 1] + dp[idx + 1]);
		    idx += blockDim.x * gridDim.x;
	    }
        __syncthreads();
    }
}

/*
__global__
void kernel(int N, float r, float rr, float *u){
    int i = (threadIdx.x + 1) * 12;
    for(int t=1; t<=100; t++){
        int k = (t % 2) * N * N;
        int kk = (1 - (t % 2)) * N * N;
        for(int j=1; j<N-1; j++) u[k + i + j] = rr * u[kk + i + j] + r * (u[kk + i - N + j] + u[kk + i + j - 1] + u[kk + i + N + j] + u[kk + i + j + 1]);
        __syncthreads();
    } 
}
*/


int main(){
    const int N = S + 2;
	int size = (2*N*N) * sizeof(float);
	float* dp_host = (float*)malloc(size);
    const float r = 0.1f;

	for(int i = 0; i < N*N; i++){
        if(i < N || i >= N*(N-1) || i % N == 0 || i % N == N - 1) continue;
        dp_host[i] = 1.0;
    }
		
	float *dp_dev;
	cudaMalloc((void**)&dp_dev,size);

	double t0 = elapsed();
	cudaMemcpy(dp_dev,dp_host,size,cudaMemcpyHostToDevice);

	double t1 = elapsed();
	
	kernel<<<1,16>>>(N, r, dp_dev);
	cudaDeviceSynchronize();

	double t2 = elapsed();

	cudaMemcpy(dp_host,dp_dev,size,cudaMemcpyDeviceToHost);

	double t3 = elapsed();
	for(int i = 0; i < N; i++){
	    for(int j = 0; j < N; j++){
            printf("%f ", dp_host[(ITER & 1)*(N*N) + i*N + j]);
        }
        printf("\n");
    }

	printf("H2D : %f\nCOMP: %f\nD2H : %f\n",t1-t0,t2-t1,t3-t2);

	cudaFree(dp_dev);
	free(dp_host);
	return 0;
}
