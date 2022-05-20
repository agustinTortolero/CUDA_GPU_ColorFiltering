#include "stdafx.h"
#include <omp.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>    // std::sort

struct pixelFastDDF
{
	float disteucl1;
	unsigned char R;
	unsigned char G;
	unsigned char B;
	float RUni;
	float GUni;
	float BUni;

	float PixelMag;


	int index;

};

int compare(const struct pixelFastDDF *a, const struct pixelFastDDF *b)
{
	if (a->disteucl1 < b->disteucl1) return -1;
	if (a->disteucl1 == b->disteucl1) return 0;
	if (a->disteucl1 > b->disteucl1) return 1;
}


////Con qsort
void FastDDF_CPU_Single(unsigned char* d_Pout, unsigned char* d_Pin, int n, int m, int channels)
{
	int x = 0, posicion[9], F = 0;
	float vectR[9], vectG[9], vectB[9];
	float valMag, disteucl, distMag;//, disteucl1[9], hold, hold2;
	int sizeFloat = n*m* sizeof(float)* channels;

	float *VectorUnitImg;	VectorUnitImg = (float *)malloc(sizeFloat);
	float *MagImg;	MagImg = (float *)malloc(sizeFloat);

	struct pixelFastDDF pixeles[9];
	struct pixelFastDDF pixel;

	for (int Col = 2; Col <= n - 2; Col++) {
		for (int Row = 2; Row <= m - 2; Row++) {
			pixel.R = d_Pin[((Row)* n + (Col)) * channels + 0];
			pixel.G = d_Pin[((Row)* n + (Col)) * channels + 1];
			pixel.B = d_Pin[((Row)* n + (Col)) * channels + 2];

			//la siguiente linea se puede hacer mas rapido
			MagImg[((Row)* n + (Col)) * channels + 0] = sqrt(pixel.R*pixel.R + pixel.G*pixel.G + pixel.B*pixel.B);


			VectorUnitImg[((Row)* n + (Col)) * channels + 0] = (pixel.R / MagImg[((Row)* n + (Col)) * channels + 0]);
			VectorUnitImg[((Row)* n + (Col)) * channels + 1] = (pixel.G / MagImg[((Row)* n + (Col)) * channels + 0]);// para visualizarla multiplciar por 255
			VectorUnitImg[((Row)* n + (Col)) * channels + 2] = (pixel.B / MagImg[((Row)* n + (Col)) * channels + 0]);

		}
	}

	for (int Col = 2; Col <= n - 2; Col++) {
		for (int Row = 2; Row <= m - 2; Row++) {

			//hacer el arreglo
			F = 0;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {

					pixeles[F].R = d_Pin[((Row + i) * n + (Col + j)) * channels + 0];
					pixeles[F].G = d_Pin[((Row + i) * n + (Col + j)) * channels + 1];
					pixeles[F].B = d_Pin[((Row + i) * n + (Col + j)) * channels + 2];


					pixeles[F].RUni = (VectorUnitImg[((Row + i) * n + (Col + j)) * channels + 0]);
					pixeles[F].GUni = (VectorUnitImg[((Row + i) * n + (Col + j)) * channels + 1]);
					pixeles[F].BUni = (VectorUnitImg[((Row + i) * n + (Col + j)) * channels + 2]);

					pixeles[F].PixelMag = (MagImg[((Row + i)* n + (Col + j)) * channels + 0]);


					//unsigned char Prueba = (unsigned char)(255*(pixel.B/valMag));
					pixeles[F].index = F;
					F++;
				}
			}

			disteucl = 0; distMag = 0;
			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 8; x++) {
					//disteucl += abs(vectB[F]-vectB[x])+abs(vectG[F]-vectG[x])+abs(vectR[F]-vectR[x]);
					//disteucl += sqrt( (pixeles[F].RUni-pixeles[x].RUni)*(pixeles[F].RUni-pixeles[x].RUni)
					//				 +(pixeles[F].GUni-pixeles[x].GUni)*(pixeles[F].GUni-pixeles[x].GUni)
					//			 +(pixeles[F].BUni-pixeles[x].BUni)*(pixeles[F].GUni-pixeles[x].GUni) );
					disteucl += sqrt(pow(pixeles[F].RUni - pixeles[x].RUni, 2)
						+ pow(pixeles[F].GUni - pixeles[x].GUni, 2)
						+ pow(pixeles[F].BUni - pixeles[x].BUni, 2));

					/// cambiar la siguiente linea por el VMF
					//distMag += sqrt(pow(pixeles[F].PixelMag - pixeles[x].PixelMag, 2));
					distMag += sqrt(pow(pixeles[F].R - pixeles[x].R, 2)
						+ pow(pixeles[F].G - pixeles[x].G, 2)
						+ pow(pixeles[F].B - pixeles[x].B, 2));
				}

				pixeles[F].disteucl1 = disteucl*distMag;
				disteucl = 0; distMag = 0;
			}

			qsort(pixeles, 9, sizeof(struct pixelFastDDF), (int(*)(const void*, const void*)) compare);

			d_Pout[(Row * n + Col) * channels + 0] = pixeles[0].R;
			d_Pout[(Row * n + Col) * channels + 1] = pixeles[0].G;
			d_Pout[(Row * n + Col) * channels + 2] = pixeles[0].B;


		}
	}
	free(VectorUnitImg);

}


void FastDDF_CPU_Multi(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads)
{
	int x = 0, posicion[9], F = 0;
	float vectR[9], vectG[9], vectB[9];
	float valMag, disteucl, distMag;//, disteucl1[9], hold, hold2;
	int sizeFloat = n*m* sizeof(float)* channels;
	float *VectorUnitImg;	VectorUnitImg = (float *)malloc(sizeFloat);
	int Row, Col, i, j;

	struct pixelFastDDF pixeles[9];
	struct pixelFastDDF pixel;





	for (Row = 2; Row <= m - 2; Row++) {

#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F ,x, pixeles, disteucl,distMag) shared(d_Pout, d_Pin,n, m,channels,Row ) schedule(static)
		//hacer el arreglo
		for (Col = 1; Col < m - 1; Col++) {
			F = 0;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {

					pixeles[F].R = d_Pin[((Row + i) * n + (Col + j)) * channels + 0];
					pixeles[F].G = d_Pin[((Row + i) * n + (Col + j)) * channels + 1];
					pixeles[F].B = d_Pin[((Row + i) * n + (Col + j)) * channels + 2];

					valMag = sqrt(pixeles[F].R*pixeles[F].R) + (pixeles[F].G*pixeles[F].G) + (pixeles[F].B*pixeles[F].B);

					pixeles[F].RUni = pixeles[F].R / valMag;
					pixeles[F].GUni = pixeles[F].G / valMag;
					pixeles[F].BUni = pixeles[F].B / valMag;

					//unsigned char Prueba = (unsigned char)(255*(pixel.B/valMag));
					pixeles[F].index = F;
					F++;
				}
			}

			disteucl = 0;
			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 8; x++) {
					//disteucl += abs(vectB[F]-vectB[x])+abs(vectG[F]-vectG[x])+abs(vectR[F]-vectR[x]);
					//disteucl += sqrt( (pixeles[F].RUni-pixeles[x].RUni)*(pixeles[F].RUni-pixeles[x].RUni)
					//				 +(pixeles[F].GUni-pixeles[x].GUni)*(pixeles[F].GUni-pixeles[x].GUni)
					//			 +(pixeles[F].BUni-pixeles[x].BUni)*(pixeles[F].GUni-pixeles[x].GUni) );
					disteucl += sqrt(pow(pixeles[F].RUni - pixeles[x].RUni, 2)
						+ pow(pixeles[F].GUni - pixeles[x].GUni, 2)
						+ pow(pixeles[F].BUni - pixeles[x].BUni, 2));

					distMag += sqrt(pow(pixeles[F].R - pixeles[x].R, 2)
						+ pow(pixeles[F].G - pixeles[x].G, 2)
						+ pow(pixeles[F].B - pixeles[x].B, 2));
				}

				pixeles[F].disteucl1 = disteucl*distMag;
				disteucl = 0; distMag = 0;
			}

			qsort(pixeles, 9, sizeof(struct pixelFastDDF), (int(*)(const void*, const void*)) compare);

			d_Pout[(Row * n + Col) * channels + 0] = pixeles[0].R;
			d_Pout[(Row * n + Col) * channels + 1] = pixeles[0].G;
			d_Pout[(Row * n + Col) * channels + 2] = pixeles[0].B;


		}
	}
	free(VectorUnitImg);

}

