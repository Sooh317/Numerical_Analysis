#include "load.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <map>

#include "bitmap.h"

SparseMatrix::SparseMatrix(){
	m = n = 0;
	rowptr = NULL;
	colind = NULL;
	val = NULL;
}

// read MM format
SparseMatrix::SparseMatrix(const char* filename){
	SparseMatrix();
	tCSRLoadFromMM(filename);
	if(m == 0) return;
	Expand();
}
SparseMatrix::~SparseMatrix(){
	if(rowptr) free(rowptr);
	if(colind) free(colind);
	if(val) free(val);
}

void SparseMatrix::Dump(){
	printf("m = %d, n = %d, nnz = %d\n",m,n,rowptr[m]);
}

void SparseMatrix::GenerateBitmap(const char* filename){
	if(m < 2000 || n < 2000){
		fprintf(stderr,"Generate Bmp : m = %d, nnz = %d, h = 1, w = 1\n",m,rowptr[m]);
		Bitmap bitmap(m,n);
		for(int i = 0; i < m; i++){
			for(int j = rowptr[i]; j < rowptr[i+1]; j++){
				unsigned x = i;
				unsigned y = colind[j];
				bitmap.set1(x,y);
			}
		}
		bitmap.write(filename);
	}
	else{
		int h = m / 1000;
		int w = n / 1000;
		fprintf(stderr,"Generate Bmp : m = %d, nnz = %d, h = %d, w = %d\n",m,rowptr[m],h,w);
		Bitmap bitmap(1000,1000);
		for(int i = 0; i < m; i++){
			for(int j = rowptr[i]; j < rowptr[i+1]; j++){
				unsigned x = i / h;
				unsigned y = colind[j] / w;
				bitmap.set1(x,y);
			}
		}
		bitmap.write(filename);
	}
}

int SparseMatrix::tCSRLoadFromMM(const char* filename){
	static char buf[4096];
	FILE* fp;
	fp = fopen(filename,"r");
	if(fp == NULL){
		fprintf(stderr,"fail at opening file\n");
		return -1;
	}
	// skip comment line
	while(1){
		if(fgets(buf,4096,fp) == NULL){
			fclose(fp);
			return -1;
		}
		if(buf[0] != '%'){
			break;
		}
	}
	int nz;
	sscanf(buf,"%d %d %d",&n,&m,&nz);
	val = (double*)malloc(nz*sizeof(double));
	colind = (int*)malloc(nz*sizeof(int));
	rowptr = (int*)malloc((m+1)*sizeof(int));
	rowptr[m] = nz;
	int tmpR = -1;
	for(int i = 0; i < nz; i++){
		if(i % 1000000 == 0) fprintf(stderr,"%d / %d\n",i,nz);
		int x,y;
		double f;
		if(fscanf(fp,"%d %d %lf",&y,&x,&f) != 3){
			fclose(fp);
			m = n = 0;
			return -1;
		}
		val[i] = f;
		colind[i] = y-1;
		if(x-1 != tmpR){
			while(tmpR != x-1){
				tmpR++;
				rowptr[tmpR] = 0;
			}
			rowptr[x-1] = i;
			tmpR = x-1;
		}
	}
	fclose(fp);
	return 0;
}

int SparseMatrix::Expand(){
	std::vector<std::vector<std::pair<int,double> > > vlist(n);
	for(int i = 0; i < m; i++){
		for(int j = rowptr[i]; j < rowptr[i+1]; j++){
			vlist[colind[j]].push_back(std::make_pair(i,val[j]));
			if(i != colind[j]){
				vlist[i].push_back(std::make_pair(colind[j],val[j]));
			}
		}
	}
	int nnz = 0;
	for(int i = 0; i < m; i++){
		std::sort(vlist[i].begin(),vlist[i].end());
		nnz += vlist[i].size();
	}
	bool both = false;
	for(int i = 0; i < m; i++){
		for(int j = 1; j < vlist[i].size(); j++){
			if(vlist[i][j-1] == vlist[i][j]){
				both = true;
				break;
			}
		}
	}
	if(both){ // upper and lower part exist
		return 0;
	}
	// else upper or lower only
	int newN = m;
	int newM = n;
	free(rowptr);
	free(colind);
	free(val);
	rowptr = (int*)malloc((newM+1)*sizeof(int));
	colind = (int*)malloc(nnz*sizeof(int));
	val = (double*)malloc(nnz*sizeof(double));
	if(rowptr == NULL) {
		fprintf(stderr,"malloc error\n");
		return -1;
	}
	int cntr = 0;
	rowptr[newM] = nnz;
	for(int i = 0; i < newM; i++){
		rowptr[i] = cntr;
		for(int j = 0; j < vlist[i].size(); j++){
			colind[cntr] = vlist[i][j].first;
			val[cntr] = vlist[i][j].second;
			cntr++;
		}
	}
	if(cntr != nnz) {
		fprintf(stderr,"not match nnz\n");
		return -1;
	}
	m = newM;
	n = newN;
	return 0;
}


