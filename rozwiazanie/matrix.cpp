#include "matrix.h"
void Complement( double** pTabO, double** pTabI, int nRow, int nCol, int nDim ); // 2x continue
void ComplMatrix( double** pTabD, double** pTab, int nDim );
void TransMatrix( double** pTab, int nDim );


//====================================
int CreateMatrix( double*** pTab, int nSize )
{
	try {
		*pTab = ( double** )malloc( nSize * sizeof( double* ) );
		if( !pTab )
			return 0;
		for ( int i = 0; i < nSize; i++ )
		{
			( *pTab )[ i ] = ( double* )malloc( nSize * sizeof( double ) );
			if ( !(pTab[ i ]) )
				return 0;
			memset( ( *pTab )[ i ], 0, nSize * sizeof( double ) );
		}
	}
	catch ( ... )
	{
		return 0;
	}
	return 1;
}




//====================================
void DeleteMatrix( double*** pTab, int nSize )
{
	double **pPom = *pTab;
	for ( int i = 0; i < nSize; i++ )
		free( *pPom++ );
	free( *pTab );
}




//====================================
void InverseMatrix( double** pInv, double **pTab, int nSize, double det )
{
	ComplMatrix( pInv, pTab, nSize );
	//PrintMatrix( pInv, nSize );
	TransMatrix( pInv, nSize );
	//PrintMatrix( pInv, nSize );
	for ( int i = 0; i < nSize; i++ )
	{
		double *pV = *pInv++;
		for ( int j = 0; j < nSize; j++ )
			*pV++ /= det;
	}
}




//====================================
double Det( double** pTab, int nSize )
{
	if ( nSize == 1 )
		return **pTab;
	if ( nSize == 2 )
	{
		return pTab[ 0 ][ 0 ] * pTab[ 1 ][ 1 ] - pTab[ 0 ][ 1 ] * pTab[ 1 ][ 0 ];
	}
	else
	{
		double res = 0;
		double **pPomTab = NULL;
		CreateMatrix( &pPomTab, nSize - 1 );
		double *pV = *pTab;
		for ( int i = 0; i < nSize; i++ )
		{
			Complement( pPomTab, pTab, 0, i, nSize );
			res += ( i % 2 ? -1 : 1 ) * *pV++ * Det( pPomTab, nSize - 1 ); // rozwinêcie Laplace'a
		}
		DeleteMatrix( &pPomTab, nSize - 1 );
		return res;
	}
}




//====================================
void LayoutEqu( double** pInv, double* pB, double* pRes, int nSize )
{
	for ( int i = 0; i < nSize; i++ )
	{
		double *pInvIt = *pInv++;
		double *pBIt = pB;
		for ( int j = 0; j < nSize; j++ )
			*pRes += *pInvIt++ * *pBIt++;
		pRes++;
	}
}




//====================================
void PrintMatrix( double** pTab, int nSize )
{
	for ( int i = 0; i < nSize; i++ )
	{
		double* pPom = *pTab++;
		for ( int j = 0; j < nSize; j++ )
			printf( "%lf\t", *pPom++ );
		printf( "\n" );
	}
	printf( "\n" );
}




//====================================
void Complement( double** pTabO, double** pTabI, int nRow, int nCol, int nDim )
{
	for ( int i = 0; i < nDim; i++ )
	{
		if ( i == nRow )
		{
			pTabI++;
			continue;
		}

		double *pI = *pTabI++;
		double *pO = *pTabO++;
		for ( int j = 0; j < nDim; j++ )
		{
			if ( j == nCol )
			{
				pI++;
				continue;
			}
			*pO++ = *pI++;
		}

	}
}




//====================================
void ComplMatrix( double** pTabD, double** pTab, int nDim )
{
	double **pPomTab = NULL;
	CreateMatrix( &pPomTab, nDim - 1 );
	double **pTabCopy = pTab;
	for ( int i = 0; i < nDim; i++ )
	{
		double *pI = *pTabCopy++;
		double *pO = *pTabD++;
		for ( int j = 0; j < nDim; j++ )
		{
			Complement( pPomTab, pTab, i, j, nDim );
			*pO++ = ( ( i + j ) % 2 ? -1 : 1 ) * Det( pPomTab, nDim - 1 );
		}
	}
	DeleteMatrix( &pPomTab, nDim - 1 );
}




//====================================
void TransMatrix( double** pTab, int nDim )
{
	for ( int i = 0; i < nDim; i++ )
	{
		for ( int j = i + 1; j < nDim; j++ )
		{
			double pom = pTab[ i ][ j ];
			pTab[ i ][ j ] = pTab[ j ][ i ];
			pTab[ j ][ i ] = pom;
		}
	}
}