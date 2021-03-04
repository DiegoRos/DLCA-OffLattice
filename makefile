all: DLCA FracDimDLCA Percolates

DLCA: DLCA.c
	gcc DLCA.c -D DEBUG -o DLCA

FracDimDLCA: FracDimDLCA.c
	gcc FracDimDLCA.c -o FracDimDLCA

Percolates: PercolatesDLCA.c
	gcc PercolatesDLCA.c -o Percolates 
