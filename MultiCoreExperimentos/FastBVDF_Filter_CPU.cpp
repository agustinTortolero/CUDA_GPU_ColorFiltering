#include "stdafx.h"
#include <omp.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>    // std::sort

struct pixelUnit
{
	float disteucl1;
	unsigned char R;
	unsigned char G;
	unsigned char B;
	float RUni;
	float GUni;
	float BUni;
	int index;

};

int compare(const struct pixelUnit *a, const struct pixelUnit *b)
{
	if (a->disteucl1 < b->disteucl1) return -1;
	if (a->disteucl1 == b->disteucl1) return 0;
	if (a->disteucl1 > b->disteucl1) return 1;
}




void FastBVDF_CPU_Multi(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads)
{
	int x = 0, posicion[9], F = 0;
	float vectR[9], vectG[9], vectB[9];
	float valMag, disteucl;//, disteucl1[9], hold, hold2;
	int sizeFloat = n*m* sizeof(float)* channels;
	float *VectorUnitImg;	VectorUnitImg = (float *)malloc(sizeFloat);
	int Row, Col, i, j;

	struct pixelUnit pixeles[9];
	struct pixelUnit pixel;

	for (Row = 2; Row <= m - 2; Row++) {
#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F ,x, pixeles, disteucl,valMag) shared(d_Pout, d_Pin,VectorUnitImg,n, m,channels,Row ) schedule(static)
		//hacer el arreglo
		for (Col = 2; Col < m - 2; Col++) {
			F = 0;
			valMag = 0;
			for ( i = -1; i <= 1; i++) {
				for ( j = -1; j <= 1; j++) {

					pixeles[F].R = d_Pin[((Row + i) * m + (Col + j)) * channels + 0];
					pixeles[F].G = d_Pin[((Row + i) * m + (Col + j)) * channels + 1];
					pixeles[F].B = d_Pin[((Row + i) * m + (Col + j)) * channels + 2];

					/* tenia esto, esta mal el parentesis*/
					//valMag = sqrt(pixeles[F].R*pixeles[F].R) + (pixeles[F].G*pixeles[F].G) + (pixeles[F].B*pixeles[F].B);

					
					if (pixeles[F].R == 0 && pixeles[F].G == 0 && pixeles[F].B == 0){
						pixeles[F].RUni = 10;
						pixeles[F].GUni = 10;
						pixeles[F].BUni = 10;
					
					}
					else{
						valMag = sqrt( (pixeles[F].R*pixeles[F].R) + (pixeles[F].G*pixeles[F].G) + (pixeles[F].B*pixeles[F].B) );
						pixeles[F].RUni = pixeles[F].R / valMag;
						pixeles[F].GUni = pixeles[F].G / valMag;
						pixeles[F].BUni = pixeles[F].B / valMag;
					}




					//unsigned char Prueba = (unsigned char)(255*(pixel.B/valMag));
					F++;
				}
			}

			disteucl = 0;
			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 8; x++) {
					//disteucl += abs(vectB[F]-vectB[x])+abs(vectG[F]-vectG[x])+abs(vectR[F]-vectR[x]);
					////Estaba este
					//disteucl += sqrt(pow(pixeles[F].RUni - pixeles[x].RUni, 2)
						//+ pow(pixeles[F].GUni - pixeles[x].GUni, 2)
						//+ pow(pixeles[F].BUni - pixeles[x].BUni, 2));

					/////
					//disteucl += pow(pixeles[F].RUni - pixeles[x].RUni, 2)
					//+ pow(pixeles[F].GUni - pixeles[x].GUni, 2)
					//+ pow(pixeles[F].BUni - pixeles[x].BUni, 2);
					disteucl += sqrt( (pixeles[F].RUni - pixeles[x].RUni) * (pixeles[F].RUni - pixeles[x].RUni)
						+ (pixeles[F].GUni - pixeles[x].GUni) * (pixeles[F].GUni - pixeles[x].GUni)
						+ (pixeles[F].BUni - pixeles[x].BUni)*(pixeles[F].BUni - pixeles[x].BUni)    );

				}
				pixeles[F].disteucl1 = disteucl;
				disteucl = 0;
			}

			qsort(pixeles, 9, sizeof(struct pixelUnit), (int(*)(const void*, const void*)) compare);

			d_Pout[(Row * m + Col) * channels + 0] = pixeles[0].R;
			d_Pout[(Row * m + Col) * channels + 1] = pixeles[0].G;
			d_Pout[(Row * m + Col) * channels + 2] = pixeles[0].B;


		}
	}
	free(VectorUnitImg);

}
