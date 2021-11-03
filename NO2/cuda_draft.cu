#include <time.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

const int ITER = 100;
const long S = 1 << 12;
const int BLOCK_NUM = 256;
const int THREAD_NUM = 1024;

double elapsed(){
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	return ts.tv_sec + ts.tv_nsec*1.0e-9;
}

__global__
void kernel(long n, long NN, int width, int height, float r, float rr, float* a){
    float *dp = a + (i & 1)*NN, *ndp = a + (1 - (i & 1))*NN;
	long base = n + 1;
	long idx = base + (blockIdx.x / width) * n + (blockIdx.x % width) * blockDim.x + threadIdx.x;
	while(idx < NN){
		ndp[idx] = rr*dp[idx] + r*(dp[idx - n] + dp[idx + n] + dp[idx - 1] + dp[idx + 1]);
	    idx += n * height;
	}
}


int main(){
    const long N = S + 2;
	const long NN = N*N;
	const long width = S / THREAD_NUM;
	const int height = BLOCK_NUM / (S / THREAD_NUM);
	int size = (2*N*N) * sizeof(float);
	float* dp_host = (float*)malloc(size);
    const float r = 0.1f;
	const float rr = 1.0 - 4*r;

	for(int i = 0; i < N*N; i++){
        if(i < N || i >= N*(N-1) || i % N == 0 || i % N == N - 1) continue;
        dp_host[i] = 1.0;
    }
		
	float *dp_dev;
	cudaMalloc((void**)&dp_dev,size);

	double t0 = elapsed();
	cudaMemcpy(dp_dev,dp_host,size,cudaMemcpyHostToDevice);

	double t1 = elapsed();
	
	for(int i = 0; i < ITER; i++){
		kernel<<<BLOCK_NUM,THREAD_NUM>>>(N, NN, width, height, r, rr, dp_dev);
	}
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