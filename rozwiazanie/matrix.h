#ifndef _MATRIX_
#define _MATRIX_

#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

int CreateMatrix( double*** pTab, int nSize ); //!!!!!! 0 error 1 ok
void DeleteMatrix( double*** pTab, int nSize );
void InverseMatrix( double** pInv, double **pTab, int nSize, double det ); //odwracanie macierzy input:pTab,nSize,det || output: pInv
double Det( double** pTab, int nSize ); //wyznacznik
void LayoutEqu( double** pInv, double* pB, double* pRes, int nSize );
void PrintMatrix( double** pTab, int nSize );

#endif // !_MATRIX_