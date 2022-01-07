#include "load.h"

#include <stdlib.h>

// grid size = n*n
void SparseMatrix::stencil5(int n){
	int nrow = n*n;
	int nnz = nrow + 4 * (n-1) * n;

	if(val) free(val);
	if(rowptr) free(rowptr);
	if(colind) free(colind);
	this->m = nrow;
	this->n = nrow;
	val = (double*)malloc(nnz*sizeof(double));
	colind = (int*)malloc(nnz*sizeof(int));
	rowptr = (int*)malloc((nrow+1)*sizeof(int));
	rowptr[nrow] = nnz;

	int cnz = 0;
	for(int i = 0; i < nrow; i++){
		rowptr[i] = cnz;
		if(i >= n){ // -n
			val[cnz] = 1.0;
			colind[cnz] = i - n;
			cnz++;
		}
		if((i % n) != 0){ // -1
			val[cnz] = 1.0;
			colind[cnz] = i-1;
			cnz++;
		}
		// 0
		val[cnz] = -4.0;
		colind[cnz] = i;
		cnz++;
		if((i % n) != n-1){ // +1
			val[cnz] = 1.0;
			colind[cnz] = i+1;
			cnz++;
		}
		if(i < nrow - n){ // +n
			val[cnz] = 1.0;
			colind[cnz] = i + n;
			cnz++;
		}
	}
}

