all: DLCA FracDimDLCA Percolates

DLCA: DLCA.c
	gcc DLCA.c -o DLCA.exe -lm

FracDimDLCA: FracDimDLCA.c
	gcc FracDimDLCA.c -o FracDimDLCA.exe -lm

Percolates: PercolatesDLCA.c
	gcc PercolatesDLCA.c -o Percolates.exe

Debug: DLCA.c
	gcc DLCA.c -D DEBUG -o DLCA -lm
	gcc -g -D DEBUG DLCA.c -lm

Rg: DLCA.c
	gcc DLCA.c -D RGINFO -o DLCA.exe

