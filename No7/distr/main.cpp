#include <stdio.h>
#include "load.h"

int main(){
//	SparseMatrix csr("../msc00726/msc00726.mtx");
	SparseMatrix csr;
	csr.stencil5(100);

	fprintf(stderr,"#row = %d\n",csr.m);
	csr.GenerateBitmap("csr.bmp");
	return 0;
}
