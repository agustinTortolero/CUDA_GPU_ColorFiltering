#include "stdafx.h"
#include <omp.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>    // std::sort

struct pixel
{
	float alpha;
	unsigned char R;
	unsigned char G;
	unsigned char B;
	int index;

};

int compareDDF(const struct pixel *a, const struct pixel *b)
{
	if (a->alpha < b->alpha) return -1;
	if (a->alpha == b->alpha) return 0;
	if (a->alpha > b->alpha) return 1;
}


////Con qsort
void DDF_CPU_Single(unsigned char* d_Pout, unsigned char* d_Pin, int n, int m, int channels)
{
	int x = 0, posicion[9], F = 0;
	float vectR[9], vectG[9], vectB[9];
	float arriva = 0.0, abajo = 0.0, valAngulo = 0.0, valDistL2 = 0.0, r = 0.0;

	struct pixel pixeles[9];

	for (int Col = 2; Col <= n - 2; Col++) {
		for (int Row = 2; Row <= m - 2; Row++) {

			//hacer el arreglo
			F = 0;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					pixeles[F].R = d_Pin[((Row + i) * n + (Col + j)) * channels + 0];
					pixeles[F].G = d_Pin[((Row + i) * n + (Col + j)) * channels + 1];
					pixeles[F].B = d_Pin[((Row + i) * n + (Col + j)) * channels + 2];

					pixeles[F].index = F;
					F++;
				}
			}

			valAngulo = 0.0, valDistL2 = 0.0;
			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 8; x++) {

					arriva = (pixeles[F].R * pixeles[x].R) + (pixeles[F].G * pixeles[x].G) + (pixeles[F].B * pixeles[x].B);
					abajo = sqrt((pixeles[F].R*pixeles[F].R) + (pixeles[F].G*pixeles[F].G) + (pixeles[F].B*pixeles[F].B)) * sqrt((pixeles[x].R*pixeles[x].R) + (pixeles[x].G*pixeles[x].G) + (pixeles[x].B*pixeles[x].B));
					valAngulo += acos(arriva / abajo); ///le quite el abs y no se aprecia cambio en la calidad del filtrado			 

					valDistL2 += sqrt(pow(pixeles[F].R - pixeles[x].R, 2)
						+ pow(pixeles[F].G - pixeles[x].G, 2)
						+ pow(pixeles[F].B - pixeles[x].B, 2));
				}
				pixeles[F].alpha = valAngulo*valDistL2;
				valAngulo = 0.0, valDistL2 = 0.0;
			}

			qsort(pixeles, 9, sizeof(struct pixel), (int(*)(const void*, const void*)) compareDDF);

			d_Pout[(Row * n + Col) * channels + 0] = pixeles[0].R;
			d_Pout[(Row * n + Col) * channels + 1] = pixeles[0].G;
			d_Pout[(Row * n + Col) * channels + 2] = pixeles[0].B;


		}
	}


}

void DDF_CPU_Multi(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads)
{
	int Col = 0, Row = 0, c = 0, d = 0, i = 0, j = 0;
	unsigned char array1[9];


	int x = 0, posicion[9], F = 0;//variables auxiliares para bubbleSort
	float vectR[9], vectG[9], vectB[9];
	float arriva = 0.0, abajo = 0.0, valAngulo = 0.0, valDistL2 = 0.0;

	struct pixel pixeles[9];

	for (Row = 1; Row < n - 1; Row++) {
#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F ,x, pixeles, arriva, abajo,valAngulo,valDistL2) shared(d_Pout, d_Pin,n, m,channels,Row ) schedule(static)

		for (Col = 1; Col < m - 1; Col++) {
			int tid = omp_get_thread_num();
			//hacer el arreglo
			F = 0;
			for (i = -1; i <= 1; i++) {
				for (j = -1; j <= 1; j++) {
					pixeles[F].R = d_Pin[((Row + i) * n + (Col + j)) * channels + 0];
					pixeles[F].G = d_Pin[((Row + i) * n + (Col + j)) * channels + 1];
					pixeles[F].B = d_Pin[((Row + i) * n + (Col + j)) * channels + 2];

					pixeles[F].index = F;
					F++;
				}
			}

			valAngulo = 0;
			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 8; x++) {

					arriva = (pixeles[F].R * pixeles[x].R) + (pixeles[F].G * pixeles[x].G) + (pixeles[F].B * pixeles[x].B);
					abajo = sqrt((pixeles[F].R*pixeles[F].R) + (pixeles[F].G*pixeles[F].G) + (pixeles[F].B*pixeles[F].B)) * sqrt((pixeles[x].R*pixeles[x].R) + (pixeles[x].G*pixeles[x].G) + (pixeles[x].B*pixeles[x].B));
					valAngulo += fabs(acos(arriva / abajo));

					valDistL2 += sqrt(pow(pixeles[F].R - pixeles[x].R, 2)
						+ pow(pixeles[F].G - pixeles[x].G, 2)
						+ pow(pixeles[F].B - pixeles[x].B, 2));


				}
				pixeles[F].alpha = valAngulo*valDistL2;
				valAngulo = 0; valDistL2 = 0;
			}

			qsort(pixeles, 9, sizeof(struct pixel), (int(*)(const void*, const void*)) compareDDF);

			d_Pout[(Row * n + Col) * channels + 0] = pixeles[0].R;
			d_Pout[(Row * n + Col) * channels + 1] = pixeles[0].G;
			d_Pout[(Row * n + Col) * channels + 2] = pixeles[0].B;

		}
	}
}