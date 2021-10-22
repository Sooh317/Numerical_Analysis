all: omp cuda

omp: omptest.c
	gcc -O3 -fopenmp -o $@ $^

cuda: cudatest.cu
	nvcc -o $@ $^

clean:
	rm omp cuda