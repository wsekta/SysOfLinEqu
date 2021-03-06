//rozwiązać układ równań 6x6
#include "matrix.h"
#include <math.h>

#define _DEBUG_

void ReadData( FILE* fin, double**pMatrix, double* b, int nDim );

int main( int argc, char** argv )
{
	if ( argc != 2 )
	{
		//mesg Usage: rozwiazanie.exe <input_file>
		printf( "Aby uzyc programu wlasciwie, wpisz komede: rozwiazanie.exe <input_file>" );
		return 1;
	}


	//otworzyć plik argv[1]
	FILE *pFile = fopen( argv[ 1 ], "r" );
	if ( pFile == NULL )
	{
		printf( "Nie udalo sie otworzyc pliku %s", argv[ 1 ] );
		return 1;
	}

	int nDim = 0;
	fscanf( pFile, "%d", &nDim );
	if ( nDim < 1 )
	{
		perror( "main: Wymiar problemu niedodatni" );
		return 1;
	}

	double** pMx = NULL;
	double* pB = ( double* )malloc( nDim * sizeof( double ) );
	if ( !pB )
	{
		perror( "main: Blad alokowania pamieci" );
		return 1;
	}
	memset( pB, 0, nDim * sizeof( double ) );

	
	//wykreować macierz równań
	if ( !CreateMatrix( &pMx, nDim ) )
	{
		perror( "main: Blad alokowania pamieci" );
		return 1;
	}

	
	//wczytać dane ( ReadData() )
	ReadData( pFile, pMx, pB, nDim );

	
	//oblicz wyznacznik
	double fDet = Det( pMx, nDim );
	if ( fabs(fDet) < 1e-10 )
	{
		printf( "main: wyznacznik równy zero, nie mozna odwrocic macierzy!" );
		return 1;
	}
	
	//obróć macierz
	double **pInvMx = NULL;
	if ( !CreateMatrix( &pInvMx, nDim ) )
	{
		perror( "main: Blad alokowania pamieci" );
		return 1;
	}
	InverseMatrix( pInvMx, pMx, nDim, fDet );
	
	//rozwiązać układ równań
	double *pRes = ( double* )malloc( nDim * sizeof( double ) );
	if ( !pRes )
	{
		perror( "main: Blad alokowania pamieci" );
		return 1;
	}
	memset( pRes, 0, nDim * sizeof( double ) );
	LayoutEqu( pInvMx, pB, pRes, nDim );


#ifdef _DEBUG_
	//wydruk kontrolny
	printf( "wymiar problemu:\t%d\n", nDim );
	printf( "\nmacierz:\n" );
	PrintMatrix( pMx, nDim );
	printf( "\nwyrazy wolne:\n" );
	double *pPom = pB;
	for ( int i = 0; i < nDim; i++ )
		printf( "%lf\t", *pPom++ );
	printf( "\n\nwyznacznik:\t%lf\n", fDet );
	printf( "\nmacierz odwrotna:\n" );
	PrintMatrix( pInvMx, nDim );
	printf( "\nrozwiazania:\n" );
	pPom = pRes;
	for ( int i = 0; i < nDim; i++ )
		printf( "x%d = %lf\t", i, *pPom++ );
	printf( "\n\n\n" );
#endif // _DEBUG_

	
	//zwolnić pamięć
	DeleteMatrix( &pMx, nDim );
	DeleteMatrix( &pInvMx, nDim );
	free( pB );
	free( pRes );
	
	return 0;
}

void ReadData( FILE* fin, double**pMatrix, double* b, int nDim )
{
	for ( int i = 0; i < nDim; i++ )
	{
		double *pPom = *pMatrix++;
		for ( int j = 0; j < nDim; j++ )
		{
			fscanf( fin, "%lf", pPom++ );
		}
		fscanf( fin, "%lf", b++ );
	}
}

/*
UWAGI:
	-wydruk na warunkowej kompilacji
	-parametrem programu jest nazwa pliku z danymi
		-pierwsza linia rozmiar
		-następne dane

		3
		1 2 3 1
		2 0 0 -1
		-1 -1 -1 1
		//det =
	-wykreować dynamicznie macierz i wektor niewiadomych, na koniec zwolnić pamięć
	-tylko standardowe i/o
	-warning dla fopen (_CRT_SECURE_NO_WARNINGS)
*/