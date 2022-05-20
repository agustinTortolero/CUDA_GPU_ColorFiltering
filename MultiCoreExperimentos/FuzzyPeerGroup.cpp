#include "stdafx.h"
#include <omp.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>    // std::sort

#define K	1024
#define d	.95
#define q	1

#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b) 


using namespace std;

void PeerGroup_Multi(unsigned char* d_Pout, const unsigned char* d_Pin, unsigned char* Noise, int n, int m, int channels, int nThreads)
{
	int Col = 0, Row = 0, c = 0, i = 0, j = 0;
	unsigned char array1[9];


	int x = 0, posicion[9], hold2 = 0, F = 0;
	unsigned char vectR[9], vectG[9], vectB[9];
	float disteucl = 0.0, disteucl1[9], hold;
	float val1 = 0.0, val2 = 0.0, val3 = 0.0, arriva = 0.0, abajo = 0.0, dist_M = 0.0;

	int  P = 0;


	for (Row = 1; Row < n - 1; Row++) {

#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F ,x,vectR,vectG,vectB,val1,val2,val3,arriva,abajo,dist_M,P) shared(d_Pout, d_Pin,Noise,n, m,channels,Row ) schedule(static)
		for (Col = 1; Col < m - 1; Col++) {
			int tid = omp_get_thread_num();
			//hacer el arreglo
			F = 0;
			P = 0;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					vectR[F] = d_Pin[((Row + i) * n + (Col + j)) * 3 + 0];
					vectG[F] = d_Pin[((Row + i) * n + (Col + j)) * 3 + 1];
					vectB[F] = d_Pin[((Row + i) * n + (Col + j)) * 3 + 2];

					posicion[F] = F;
					F++;
				}
			}

			for (F = 0; F <= 8; F++) {
				arriva = min(vectR[F], vectR[4]) + K;
				abajo = max(vectR[F], vectR[4]) + K;
				val1 = arriva / abajo;

				arriva = min(vectG[F], vectG[4]) + K;
				abajo = max(vectG[F], vectG[4]) + K;
				val2 = arriva / abajo;

				arriva = min(vectB[F], vectB[4]) + K;
				abajo = max(vectB[F], vectB[4]) + K;
				val3 = arriva / abajo;

				dist_M = min(min(val1, val2), val3);
				if (dist_M>d)	P++;
			}

			if (P <= (q + 1)) {
				Noise[(Row * n + Col)] = 255;
			}
			else {
				Noise[(Row * n + Col)] = 0;
			}

		}
	}



	////Filtrado
	float sumR = 0.0, sumG = 0.0, sumB = 0.0;
	int Div = 0;
	for (Row = 1; Row < n - 1; Row++) {

#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F ,x,vectR,vectG,vectB,sumR,sumG,sumB,Div) shared(d_Pout, d_Pin,Noise,n, m,channels,Row ) schedule(static)
		for (Col = 1; Col < m - 1; Col++) {
			//printf("%d     ",nThreads);

			sumR = 0.0, sumG = 0.0, sumB = 0.0;
			if (Noise[(Row * n + Col)] == 255) {
				Div = 0;
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {

						if (Noise[((Row + i) * n + (Col + j))] == 0) {//solo los que no son Noise
							Div++;
							sumR += d_Pin[((Row + i) * n + (Col + j)) * 3 + 0];
							sumG += d_Pin[((Row + i) * n + (Col + j)) * 3 + 1];
							sumB += d_Pin[((Row + i) * n + (Col + j)) * 3 + 2];
						}
					}
				}


				d_Pout[((Row*n) + Col) * 3 + 0] = sumR / Div;
				d_Pout[((Row*n) + Col) * 3 + 1] = sumG / Div;
				d_Pout[((Row*n) + Col) * 3 + 2] = sumB / Div;

				//printf("ruido");
				//d_Pout[((Row*n) + Col)*3 + 0]=0;
				//d_Pout[((Row*n) + Col)*3 + 1]=0;
				//d_Pout[((Row*n) + Col)*3 + 2]=0;

			}//fin de if
			else {
				d_Pout[((Row*n) + Col) * 3 + 0] = d_Pin[((Row*n) + Col) * 3 + 0];
				d_Pout[((Row*n) + Col) * 3 + 1] = d_Pin[((Row*n) + Col) * 3 + 1];
				d_Pout[((Row*n) + Col) * 3 + 2] = d_Pin[((Row*n) + Col) * 3 + 2];

			}



		}

	}//for rows


}//final


void Detection_FuzzyMetric(unsigned char* Noise, const unsigned char* d_Pin, int n, int m, int channels, int nThreads)
{
	int Col = 0, Row = 0, c = 0, i = 0, j = 0;
	unsigned char array1[9];



	unsigned char vectR[9], vectG[9], vectB[9];
	float val1 = 0.0, val2 = 0.0, val3 = 0.0, arriva = 0.0, abajo = 0.0, dist_M = 0.0;

	int  P = 0, F = 0;


	for (Row = 1; Row < n - 1; Row++) {

#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F,vectR,vectG,vectB,val1,val2,val3,arriva,abajo,dist_M,P) shared( d_Pin,Noise,n, m,channels,Row ) schedule(static)
		for (Col = 1; Col < m - 1; Col++) {
			//hacer el arreglo
			P = 0; F = 0;
			for (i = -1; i <= 1; i++) {
				for (j = -1; j <= 1; j++) {
					vectR[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
					vectG[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
					vectB[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
					F++;

				}
			}

			for (F = 0; F <= 8; F++) {
				arriva = min(vectR[F], vectR[4]) + K;
				abajo = max(vectR[F], vectR[4]) + K;
				val1 = arriva / abajo;

				arriva = min(vectG[F], vectG[4]) + K;
				abajo = max(vectG[F], vectG[4]) + K;
				val2 = arriva / abajo;

				arriva = min(vectB[F], vectB[4]) + K;
				abajo = max(vectB[F], vectB[4]) + K;
				val3 = arriva / abajo;

				dist_M = min(min(val1, val2), val3);
				if (dist_M>d)	P++;
			}

			if (P <= (q + 1)) {
				Noise[(Row * m + Col)] = 255;
			}
			else {
				Noise[(Row * m + Col)] = 0;
			}

		}
	}

}

void AMF_Filtering(unsigned char* d_Pout, const unsigned char* d_Pin, unsigned char* Noise, int n, int m, int channels, int nThreads)
{

	int Col = 0, Row = 0, c = 0, i = 0, j = 0;
	unsigned char array1[9];
	int F = 0;
	unsigned char vectR[9], vectG[9], vectB[9];

	int  P = 0;

	float sumR = 0.0, sumG = 0.0, sumB = 0.0;
	float Div = 0;
	for (Row = 1; Row < n - 1; Row++) {

		#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F,vectR,vectG,vectB,sumR,sumG,sumB,Div) shared(d_Pout, d_Pin,Noise,n, m,channels,Row ) schedule(static)
		for (Col = 1; Col < m - 1; Col++) {
			//printf("%d     ",nThreads);

			sumR = 0.0, sumG = 0.0, sumB = 0.0;
			if (Noise[(Row * m + Col)] == 255) {
				Div = 0;
				for (i = -1; i <= 1; i++) {
					for (j = -1; j <= 1; j++) {

						if (Noise[((Row + i) * m + (Col + j))] == 0) {//solo los que no son Noise
							Div++;
							sumR += d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
							sumG += d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
							sumB += d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
						}
					}
				}


				d_Pout[((Row*m) + Col) * 3 + 0] = sumR / Div;
				d_Pout[((Row*m) + Col) * 3 + 1] = sumG / Div;
				d_Pout[((Row*m) + Col) * 3 + 2] = sumB / Div;

				//printf("ruido");
				//d_Pout[((Row*n) + Col)*3 + 0]=0;
				//d_Pout[((Row*n) + Col)*3 + 1]=0;
				//d_Pout[((Row*n) + Col)*3 + 2]=0;

			}//fin de if
			else {
				d_Pout[((Row*m) + Col) * 3 + 0] = d_Pin[((Row*m) + Col) * 3 + 0];
				d_Pout[((Row*m) + Col) * 3 + 1] = d_Pin[((Row*m) + Col) * 3 + 1];
				d_Pout[((Row*m) + Col) * 3 + 2] = d_Pin[((Row*m) + Col) * 3 + 2];

			}



		}


	}

}
void AMF_Filtering_ExpWindow(unsigned char* d_Pout, const unsigned char* d_Pin, unsigned char* Noise, int n, int m, int channels, int nThreads)
{

	int Col = 0, Row = 0, c = 0, i = 0, j = 0, i_exp = 1, j_exp = 1;
	unsigned char array1[9];
	int F = 0;
	unsigned char vectR[9], vectG[9], vectB[9];
	float sumR = 0.0, sumG = 0.0, sumB = 0.0;
	int Div = 0;
	for (Row = 9; Row < n - 9; Row++) {

		#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F,vectR,vectG,vectB,sumR,sumG,sumB,Div) shared(d_Pout, d_Pin,Noise,n, m,channels,Row ) schedule(static)
		for (Col = 9; Col < m - 9; Col++) {
			//printf("%d     ",nThreads);
			if (Noise[(Row * m + Col)] == 255) {
				sumR = 0.0, sumG = 0.0, sumB = 0.0;

				Div = 0;
				for (i = -1; i <= 1; i++) {
					for (j = -1; j <= 1; j++) {
						if (Noise[((Row + i) * m + (Col + j))] == 0) {//solo los que no son Noise
							Div++;
							sumR += d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
							sumG += d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
							sumB += d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
						}
					}
				}

				//si Div==0 esque todos son ruido, se expande la ventana en 2 para ser de 5x5
				if (Div == 0) {
					sumR = 0.0, sumG = 0.0, sumB = 0.0;
					for (i = -2; i <= 2; i++) {
						for (j = -2; j <= 2; j++) {
							if (Noise[((Row + i) * m + (Col + j))] == 0) {//solo los que no son Noise
								Div++;
								sumR += d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
								sumG += d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
								sumB += d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
							}
						}
					}
				}

				if (Div == 0) {
					sumR = 0.0, sumG = 0.0, sumB = 0.0;
					for (i = -3; i <= 3; i++) {
						for (j = -3; j <= 3; j++) {
							if (Noise[((Row + i) * m + (Col + j))] == 0) {//solo los que no son Noise
								Div++;
								sumR += d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
								sumG += d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
								sumB += d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
							}
						}
					}

				}
				if (Div == 0) {
					Div = 7 * 7;
					printf("No se cumplio ninguna \n");
				}

				d_Pout[((Row*m) + Col) * 3 + 0] = sumR / Div;
				d_Pout[((Row*m) + Col) * 3 + 1] = sumG / Div;
				d_Pout[((Row*m) + Col) * 3 + 2] = sumB / Div;

			}//fin de if
			else {
				d_Pout[((Row*m) + Col) * 3 + 0] = d_Pin[((Row*m) + Col) * 3 + 0];
				d_Pout[((Row*m) + Col) * 3 + 1] = d_Pin[((Row*m) + Col) * 3 + 1];
				d_Pout[((Row*m) + Col) * 3 + 2] = d_Pin[((Row*m) + Col) * 3 + 2];

			}



		}


	}

}

void VMF_Filtering(unsigned char* d_Pout, const unsigned char* d_Pin, unsigned char* Noise, int n, int m, int channels, int nThreads)
{
	int Col = 0, Row = 0, i = 0, j = 0;
	unsigned char array1[9];
	int F = 0;
	unsigned char vectR[9], vectG[9], vectB[9];

	unsigned int  P = 0, posicion[9];
	float disteucl = 0, disteucl1[9], hold;
	float arrayFiltradoR[9], arrayFiltradoG[9], arrayFiltradoB[9];
	float mn, mx;
	int posMin = 0;


	unsigned int c = 0;

	for (Row = 1; Row < n - 1; Row++) {
#pragma omp parallel for  num_threads(nThreads) private(Col, i, j, F, disteucl1, disteucl, vectR, vectG, vectB, hold, posicion, posMin, mn, mx,c,arrayFiltradoR,arrayFiltradoG,arrayFiltradoB) shared(d_Pout, d_Pin,Noise,n, m,channels,Row ) schedule(static)
		for (Col = 1; Col < m - 1; Col++) {             
			//printf("%d     ",nThreads);

			if (Noise[(Row * m + Col)] == 255) {
				c = 0;
				F = 0;
				for (i = -1; i <= 1; i++) {
					for (j = -1; j <= 1; j++) {
						posicion[F] = 0;

						if (Noise[((Row + i) * m + (Col + j))] == 0) {//solo los que no son Noise

							arrayFiltradoR[c] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
							arrayFiltradoG[c] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
							arrayFiltradoB[c] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
							posicion[c] = c;
							c++;
						}
						F++;
					}
				}

				if (c == 0)	c = 1; //Todos son Noise, se hace c=1 como seguro
								   //if (c>= 1)	c = 9; //si todos son ruido se aplica VMF a la ventana

				for (i = 0; i <= c - 1; i++) {
					disteucl = 0;
					for (j = 0; j <= c - 1; j++) {
						float distR = abs(arrayFiltradoR[i] - arrayFiltradoR[j]);
						float distG = abs(arrayFiltradoG[i] - arrayFiltradoG[j]);
						float distB = abs(arrayFiltradoB[i] - arrayFiltradoB[j]);
						disteucl += distR + distB + distG;

					}
					disteucl1[i] = disteucl;
				}

				mn = disteucl1[1];
				mx = disteucl1[1];

				posMin = 0;

				for (i = 0; i <= c - 1; i++)
				{
					if (mn>disteucl1[i])
					{
						mn = disteucl1[i];
						posMin = posicion[i];
					}
					else if (mx<disteucl1[i])
					{

					}
				}


				d_Pout[(Row * m + Col) * 3 + 0] = arrayFiltradoR[posMin];
				d_Pout[(Row * m + Col) * 3 + 1] = arrayFiltradoG[posMin];
				d_Pout[(Row * m + Col) * 3 + 2] = arrayFiltradoB[posMin];

			}//fin de if
			else {
				d_Pout[((Row*m) + Col) * 3 + 0] = d_Pin[((Row*m) + Col) * 3 + 0];
				d_Pout[((Row*m) + Col) * 3 + 1] = d_Pin[((Row*m) + Col) * 3 + 1];
				d_Pout[((Row*m) + Col) * 3 + 2] = d_Pin[((Row*m) + Col) * 3 + 2];

			}
		}
	}
}

void VMF_Filtering_Single(unsigned char* d_Pout, const unsigned char* d_Pin, unsigned char* Noise, int n, int m, int channels)
{

	int Col = 0, Row = 0, i = 0, j = 0;
	unsigned char array1[9];
	int F = 0;
	unsigned char vectR[9], vectG[9], vectB[9];

	unsigned int  P = 0, posicion[9];
	float disteucl = 0, disteucl1[9], hold;
	float arrayFiltradoR[9], arrayFiltradoG[9], arrayFiltradoB[9];
	float mn, mx;
	int posMin = 0;


	unsigned int c = 0;

	for (Row = 1; Row < n - 1; Row++) {
		for (Col = 1; Col < m - 1; Col++) {
			//printf("%d     ",nThreads);

			if (Noise[(Row * m + Col)] == 255) {
				c = 0;
				F = 0;
				for (i = -1; i <= 1; i++) {
					for (j = -1; j <= 1; j++) {
						posicion[F] = 0;

						if (Noise[((Row + i) * m + (Col + j))] == 0) {//solo los que no son Noise

							arrayFiltradoR[c] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
							arrayFiltradoG[c] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
							arrayFiltradoB[c] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
							posicion[c] = c;
							c++;
						}
						F++;
					}
				}

				if (c == 0)	c = 1; //Todos son Noise, se hace c=1 como seguro
								   //if (c>= 1)	c = 9; //si todos son ruido se aplica VMF a la ventana

				for (i = 0; i <= c - 1; i++) {
					disteucl = 0;
					for (j = 0; j <= c - 1; j++) {
						float distR = abs(arrayFiltradoR[i] - arrayFiltradoR[j]);
						float distG = abs(arrayFiltradoG[i] - arrayFiltradoG[j]);
						float distB = abs(arrayFiltradoB[i] - arrayFiltradoB[j]);
						disteucl += distR + distB + distG;

					}
					disteucl1[i] = disteucl;
				}

				mn = disteucl1[1];
				mx = disteucl1[1];

				posMin = 0;

				for (i = 0; i <= c - 1; i++)
				{
					if (mn>disteucl1[i])
					{
						mn = disteucl1[i];
						posMin = posicion[i];
					}
					else if (mx<disteucl1[i])
					{

					}
				}


				d_Pout[(Row * m + Col) * 3 + 0] = arrayFiltradoR[posMin];
				d_Pout[(Row * m + Col) * 3 + 1] = arrayFiltradoG[posMin];
				d_Pout[(Row * m + Col) * 3 + 2] = arrayFiltradoB[posMin];

			}//fin de if
			else {
				d_Pout[((Row*m) + Col) * 3 + 0] = d_Pin[((Row*m) + Col) * 3 + 0];
				d_Pout[((Row*m) + Col) * 3 + 1] = d_Pin[((Row*m) + Col) * 3 + 1];
				d_Pout[((Row*m) + Col) * 3 + 2] = d_Pin[((Row*m) + Col) * 3 + 2];

			}
		}
	}
}

void FiltradoPropuesta(unsigned char* d_Pout, const unsigned char* d_Pin, unsigned char* Noise, int n, int m, int channels, int nThreads)
{
	int Col = 0, Row = 0, i = 0, j = 0;
	unsigned char array1[9];
	int F = 0;
	unsigned char vectR[9], vectG[9], vectB[9];

	unsigned int  P = 0, posicion[9], contador_nNoise = 0;
	float disteucl = 0, disteucl1[9], hold;
	float arrayFiltradoR[9], arrayFiltradoG[9], arrayFiltradoB[9];
	float mn, mx;
	int posMin = 0;
	float sumR = 0.0, sumG = 0.0, sumB = 0.0;
	int Div = 0;


	unsigned int c = 0;

	for (Row = 1; Row < n - 1; Row++) {
#pragma omp parallel for  num_threads(nThreads) private(Col, i, j, F, disteucl1, disteucl, vectR, vectG, vectB, hold, posicion, posMin, mn, mx,c) shared(d_Pout, d_Pin,Noise,n, m,channels,Row ) schedule(static)

		for (Col = 1; Col < m - 1; Col++) {
			if (Noise[(Row * m + Col)] == 255) {
				c = 0;
				for (i = -1; i <= 1; i++) {
					for (j = -1; j <= 1; j++) {
						vectR[c] = d_Pin[((Row + i) * n + (Col + j)) * 3 + 0];
						vectG[c] = d_Pin[((Row + i) * n + (Col + j)) * 3 + 1];
						vectB[c] = d_Pin[((Row + i) * n + (Col + j)) * 3 + 2];

						posicion[c] = c;
						c++;
					}
				}
				for (i = 0; i <= 8; i++) {
					disteucl = 0;
					disteucl1[i] = 0;
					for (j = 0; j <= 8; j++) {
						float distR = abs(vectR[i] - vectR[j]);
						float distG = abs(vectG[i] - vectG[j]);
						float distB = abs(vectB[i] - vectB[j]);
						disteucl += distR + distB + distG;
					}
					disteucl1[i] = disteucl;
				}
				mn = disteucl1[0];
				mx = disteucl1[0];
				posMin = 0;
				for (i = 0; i <= 8; i++) {
					if (mn > disteucl1[i]) {
						mn = disteucl1[i];
						posMin = posicion[i];
					}
				}
				d_Pout[(Row * m + Col) * 3 + 0] = vectR[posMin];
				d_Pout[(Row * m + Col) * 3 + 1] = vectG[posMin];
				d_Pout[(Row * m + Col) * 3 + 2] = vectB[posMin];
			}//fin de if (Noise[(Row * m + Col)] == 255)
			else {
				d_Pout[((Row*m) + Col) * 3 + 0] = d_Pin[((Row*m) + Col) * 3 + 0];
				d_Pout[((Row*m) + Col) * 3 + 1] = d_Pin[((Row*m) + Col) * 3 + 1];
				d_Pout[((Row*m) + Col) * 3 + 2] = d_Pin[((Row*m) + Col) * 3 + 2];

			}
		}
	}
}


#define maxCUDA( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define minCUDA( a, b ) ( ((a) < (b)) ? (a) : (b) )

inline void s(unsigned char* a, unsigned char*b)
{
	int tmp;
	if (*a>*b) {//si a es mayor a b, se intercambian a y b.
		tmp = *b;
		*b = *a;
		*a = tmp;
	}
}

#define min3(a,b,c) s(a, b); s(a,c);
#define max3(a,b,c) s(b, c); s(a,c);

#define minmax3(a,b,c)			max3(a, b, c); s(a,b);
#define minmax4(a,b,c,d)		s(a, b); s(c,d);s(a, c); s(b,d);
#define minmax5(a,b,c,d,e)		s(a, b); s(c,d);min3(a,c,e);max3(b,d,e);

#define minmax6(a,b,c,d,e,f)	s(a,d);s(b,e);s(c,f);min3(a,b,c);max3(d,e,f);


void MPGFMF(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads)
{
	int Col = 0, Row = 0;

	int x = 0, F = 0, i, j;
	unsigned char vectR[9], vectG[9], vectB[9],Noise=0;
	
	float arriva = 0.0, abajo = 0.0, val1, val2, val3, dist_M = 0;
	unsigned int P = 0, hold=0;


	for (Row = 1; Row < n - 1; Row++) {
#pragma omp parallel for  num_threads(nThreads) private(Col, i, j, F, vectR, vectG, vectB, hold, Noise, P, x, arriva, abajo, val1, val2, val3,dist_M) shared(d_Pout, d_Pin,n, m,channels,Row ) schedule(static)

		for (Col = 1; Col < m - 1; Col++) {
			F = 0; Noise = 0; P = 0;
			for ( i = -1; i <= 1; i++) {
				for ( j = -1; j <= 1; j++) {
					
					vectR[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
					vectG[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
					vectB[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
					
					F++;
				}
			}
			
			for (F = 0; F <= 8; F++) {
				
				arriva = minCUDA(vectR[F], vectR[4]) + K;
				abajo = maxCUDA(vectR[F], vectR[4]) + K;
				val1 = arriva / abajo;

				arriva = minCUDA(vectG[F], vectG[4]) + K;
				abajo = maxCUDA(vectG[F], vectG[4]) + K;
				val2 = arriva / abajo;

				arriva = minCUDA(vectB[F], vectB[4]) + K;
				abajo = maxCUDA(vectB[F], vectB[4]) + K;
				val3 = arriva / abajo;

				dist_M = minCUDA(minCUDA(val1, val2), val3);
				
				if (dist_M>d)	P++;
				
			}

			if (P <= (q + 1)) {
				Noise = 255;
			}
			else {
				Noise = 0;
			}
			if (Noise == 255) {
				
		
			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 7; x++) {
					if (vectR[x] > vectR[x + 1]) {
						hold = vectR[x];
						vectR[x] = vectR[x + 1];
						vectR[x + 1] = hold;
					}
				}
			}

			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 7; x++) {
					if (vectG[x] > vectG[x + 1]) {
						hold = vectG[x];
						vectG[x] = vectG[x + 1];
						vectG[x + 1] = hold;
					}
				}
			}
			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 7; x++) {
					if (vectB[x] > vectB[x + 1]) {
						hold = vectB[x];
						vectB[x] = vectB[x + 1];
						vectB[x + 1] = hold;
					}
				}
			}

			
			d_Pout[(Row * n + Col) * 3 + 0] = vectR[4];
			d_Pout[(Row * n + Col) * 3 + 1] = vectG[4];
			d_Pout[(Row * n + Col) * 3 + 2] = vectB[4];
			

			}
			else {
				
				d_Pout[((Row*m) + Col) * 3 + 0] = vectR[4];
				d_Pout[((Row*m) + Col) * 3 + 1] = vectG[4];
				d_Pout[((Row*m) + Col) * 3 + 2] = vectB[4];
				
			}
			
		}
	}
}





