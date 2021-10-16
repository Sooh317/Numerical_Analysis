#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define PROC 128
#define NUM (1<<20)

const int part = NUM / PROC;

void internal(int* a, int i, int j, int sz, int myid){
    int base = myid * part;
    if(base % sz) return;
    int d = 1 << (i - j);
    // if(myid == 0){
    //     printf("a[0] : %d, a[1] : %d", a[0], a[1]);
    // }
    for(int k = base; k < base + sz; k++){
        int up = (((k >> i) & 2) == 0 ? 1 : 0);
        int idx = k - base;
        if((k & d) == 0 && (a[idx] > a[idx | d]) == up){
            int t = a[idx];
            a[idx] = a[idx | d];
            a[idx | d] = t;
        }
    }
    // if(myid == 0){
    //     printf("a[0] : %d, a[1] : %d", a[0], a[1]);
    // }
}

// gatherがちゃんと動いているか確認
int gather(int d, int myid, int* buf){
    int width = 1 << d;
    if(width <= part) return part;
    int num = width / part;
    if(myid % num == 0){ 
        //printf("myid : %d is gathering from %d\n", myid, myid + 1);
        for(int i = 1; i < num; i++){
            MPI_Recv(buf + i*part, part, MPI_INT, myid + i,  myid + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else{
        //printf("myid : %d is sending to %d\n", myid, myid / num * num);
        MPI_Send(buf, part, MPI_INT, myid / num * num, myid, MPI_COMM_WORLD);
    }
    return width;
}

void distribute(int myid, int width, int ij, int* buf){
    //printf("IN DISTRIBUTE ||||| width : %d", width);
    int d = width / 2; // next width
    int num = d / part; // num parts compose one set
    if(myid * part % width == 0){
        MPI_Send(buf + d, d, MPI_INT, myid + num, myid, MPI_COMM_WORLD);
    }
    else if((myid + num) * part % width == 0){
        MPI_Recv(buf, d, MPI_INT, myid - num, myid - num, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char** argv){
	int myid, numproc;
    int t1, t2, elapsed, lg = 0, tmp = NUM;
    while(tmp != 1) lg++, tmp >>= 1;

    //printf("log : %d\n", lg);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int buf[NUM];

	if(myid == 0){
        for(int i = 0; i < NUM; i++){
            buf[i] = rand();
            //printf("%d ", buf[i]);
        }
        //printf("\n");
    }
    MPI_Scatter(buf, NUM / PROC, MPI_INTEGER, buf, NUM / PROC, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    
    
    for(int i = 0; i < lg; i++){
        int width = gather(i + 1, myid, buf);
        for(int j = 0; j <= i; j++){
            internal(buf, i, j, width, myid);
            if(width > part){
                distribute(myid, width, i - j, buf);
                width >>= 1;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();

    elapsed = t2 - t1;
    MPI_Gather(buf, part, MPI_INT, buf, part, MPI_INT, 0, MPI_COMM_WORLD);

    if(myid == 0){
        printf("%d\n", elapsed);
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
