all: DLCA FracDimDLCA Percolates

DLCA: DLCA.c
	gcc DLCA.c -o DLCA -lm

FracDimDLCA: FracDimDLCA.c
	gcc FracDimDLCA.c -o FracDimDLCA -lm

Percolates: PercolatesDLCA.c
	gcc PercolatesDLCA.c -o Percolates 
