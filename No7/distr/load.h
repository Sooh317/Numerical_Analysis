#ifndef LOAD_H
#define LOAD_H

class SparseMatrix{
public:
	int m,n;
	int* rowptr; // size = m+1
	int* colind; // size = rowptr[m]
	double* val; // size = rowptr[m]

	SparseMatrix();
	SparseMatrix(const char* filename); // read MM format (symmetric)
	~SparseMatrix();

	void stencil5(int size); // (grid size) = size * size

	void Dump();
	void GenerateBitmap(const char* filename);

	// subroutine
	int tCSRLoadFromMM(const char* filename);
	int Expand();
};


#endif
