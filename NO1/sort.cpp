// compile like this :  mpic++ -O3 main.cpp
// node number -> 4
// task number per node -> 8
#pragma GCC optimize("unroll-loops")
#include <mpi.h>
#include <stdio.h>
#include <algorithm>
#include <functional>
#include <stdlib.h>

constexpr int PROC = 1<<5;
constexpr int procp = 5;
constexpr int NUM = 1 << 23;
constexpr int nump = 23;
constexpr int part = NUM / PROC;
constexpr int partp = nump - procp;
int myid, numproc, myp;
int buf[NUM];

void internal(int i, int j, int sz, int myid){
    int base = myid << partp;
    if(base & (sz - 1)) return;
    int d = 1 << (i - j);
    for(int k = base; k < base + sz; k++){
        int up = (((k >> i) & 2) == 0 ? 1 : 0);
        int idx = k - base;
        if((k & d) == 0 && (buf[idx] > buf[idx | d]) == up){
            int t = buf[idx];
            buf[idx] = buf[idx | d];
            buf[idx | d] = t;
        }
    }
}

int gather(int d, int myid){
    int width = 1 << d;
    int num = width >> partp;
    int np = d - partp;
    if((myid & (num - 1)) == 0){
        int tmp = part;
        for(int i = 1; i < num; i++){
            MPI_Recv(buf + tmp, part, MPI_INT, myid + i,  myid + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            tmp += part;
        }
    }
    else{
        MPI_Send(buf, part, MPI_INT, (myid >> np) << np, myid, MPI_COMM_WORLD);       
    }
    return width;
}

void distribute(int width){
    int d = width >> 1; // next width
    int num = d >> partp; // num parts compose one set
    if((myp & (width - 1)) == 0){
        MPI_Send(buf + d, d, MPI_INT, myid + num, myid, MPI_COMM_WORLD);
    }
    else if((((myid + num) << partp) & (width - 1)) == 0){
        MPI_Recv(buf, d, MPI_INT, myid - num, myid - num, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char** argv){
    double t1, t2, elapsed;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    myp = myid * part;

    if(myid == 0){
        for(int i = 0; i < NUM; i++) buf[i] = rand();
    }
    MPI_Scatter(buf, part, MPI_INTEGER, buf, part, MPI_INTEGER, 0, MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();

    if(myid & 1) std::sort(buf, buf + part, std::greater<int>());
    else std::sort(buf, buf + part);

    for(int i = partp; i < nump; i++){
        int width = gather(i + 1, myid);
        for(int j = 0; j <= i; j++){
            internal(i, j, width, myid);
            if(width > part){
                distribute(width);
                width >>= 1;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();

    elapsed = t2 - t1;
    MPI_Gather(buf, part, MPI_INT, buf, part, MPI_INT, 0, MPI_COMM_WORLD);

    if(myid == 0){
        printf("%lf\n", elapsed);
        for(int i = 0; i < NUM - 1; i++){
            if(buf[i] > buf[i + 1]){
                printf("fail\n");
                break;
            }
        }
    }

    MPI_Finalize();
    return 0;
}
