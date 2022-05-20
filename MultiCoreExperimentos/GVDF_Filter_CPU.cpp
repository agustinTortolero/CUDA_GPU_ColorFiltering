#include "stdafx.h"
#include <omp.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>    // std::sort



float Magnitud_L1(unsigned char* VectR, unsigned char* VectG, unsigned char* VectB, unsigned int i, unsigned int j) {

	float distR = abs(VectR[i] - VectR[j]);
	float distG = abs(VectG[i] - VectG[j]);
	float distB = abs(VectB[i] - VectB[j]);
	return distR + distB + distG;

}

float Angulo(unsigned char* VectR, unsigned char* VectG, unsigned char* VectB, unsigned int i, unsigned int j)
{

	float VectR_i_Cuadrado = VectR[i] * VectR[i];
	float VectG_i_Cuadrado = VectG[i] * VectG[i];
	float VectB_i_Cuadrado = VectB[i] * VectB[i];

	float VectR_j_Cuadrado = VectR[j] * VectR[j];
	float VectG_j_Cuadrado = VectG[j] * VectG[j];
	float VectB_j_Cuadrado = VectB[j] * VectB[j];


	float arriva = (VectR[i] * VectR[j]) + (VectG[i] * VectG[j]) + (VectB[i] * VectB[j]);
	float abajo = sqrt(VectR_i_Cuadrado + VectG_i_Cuadrado + VectB_i_Cuadrado)*sqrt(VectR_j_Cuadrado + VectG_j_Cuadrado + VectB_j_Cuadrado);
	return acos(arriva / abajo);

}


void GVDF_Filter_CPU_Multi(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels)
{
	int Col = 0, Row = 0, c = 0, d = 0, i = 0, j = 0;
	unsigned char array1[9];


	int x = 0, posicion[9], hold2 = 0, F = 0;
	unsigned char vectR[9], vectG[9], vectB[9];
	float disteucl = 0.0, disteucl1[9], distMag[9], hold;
	float D[40];
	float mn, mx;
	int posMin = 0;
	unsigned int posAux_i = 0, posAux_j = 0;

	const unsigned int nElementos = 3;

	for (Row = 1; Row < n - 1; Row++) {
		//#pragma omp parallel for  num_threads(8) private(Col ,i, j, F ,x,disteucl1,disteucl,vectR,vectG,vectB,hold,hold2,posicion,posMin,mn,mnx) shared(d_Pout, d_Pin,n, m,channels,Row ) schedule(static)

		for (Col = 1; Col < m - 1; Col++) {
			int tid = omp_get_thread_num();
			//hacer el arreglo
			F = 0;

			for (i = -1; i <= 1; i++) {
				for (j = -1; j <= 1; j++) {
					vectR[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
					vectG[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
					vectB[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];

					posicion[F] = F;
					F++;
				}
			}			//D[0]=Magnitud(vectR, vectG, vectB, i, j//i==0 y j==0 no se hace
			D[0] = (Angulo(vectR, vectG, vectB, 0, 1));
			D[1] = (Angulo(vectR, vectG, vectB, 0, 2));
			D[2] = (Angulo(vectR, vectG, vectB, 0, 3));
			D[3] = (Angulo(vectR, vectG, vectB, 0, 4));
			D[4] = (Angulo(vectR, vectG, vectB, 0, 5));
			D[5] = (Angulo(vectR, vectG, vectB, 0, 6));
			D[6] = (Angulo(vectR, vectG, vectB, 0, 7));
			D[7] = (Angulo(vectR, vectG, vectB, 0, 8));
			disteucl1[0] = D[0] + D[1] + D[2] + D[3] + D[4] + D[5] + D[6] + D[7];

			//i=1,j=0 ya esta es D[0]
			//i=1,j=1 No se hace
			D[8] = (Angulo(vectR, vectG, vectB, 1, 2));
			D[9] = (Angulo(vectR, vectG, vectB, 1, 3));
			D[10] = (Angulo(vectR, vectG, vectB, 1, 4));
			D[11] = (Angulo(vectR, vectG, vectB, 1, 5));
			D[12] = (Angulo(vectR, vectG, vectB, 1, 6));
			D[13] = (Angulo(vectR, vectG, vectB, 1, 7));
			D[14] = (Angulo(vectR, vectG, vectB, 1, 8));
			disteucl1[1] = D[0] + D[8] + D[9] + D[10] + D[11] + D[12] + D[13] + D[14];

			//i=2,j=0 ya esta es D[1]
			//i=2,j=1 ya esta es D[8]
			//i=2,j=2 No se hace
			D[15] = (Angulo(vectR, vectG, vectB, 2, 3));
			D[16] = (Angulo(vectR, vectG, vectB, 2, 4));
			D[17] = (Angulo(vectR, vectG, vectB, 2, 5));
			D[18] = (Angulo(vectR, vectG, vectB, 2, 6));
			D[19] = (Angulo(vectR, vectG, vectB, 2, 7));
			D[20] = (Angulo(vectR, vectG, vectB, 2, 8));
			disteucl1[2] = D[1] + D[8] + D[15] + D[16] + D[17] + D[18] + D[19] + D[20];

			//i=3,j=0 ya esta es D[2]
			//i=3,j=1 ya esta es D[9]
			//i=3,j=2 ya esta es D[15]
			//i=3,j=3 No se hace
			D[21] = (Angulo(vectR, vectG, vectB, 3, 4));
			D[22] = (Angulo(vectR, vectG, vectB, 3, 5));
			D[23] = (Angulo(vectR, vectG, vectB, 3, 6));
			D[24] = (Angulo(vectR, vectG, vectB, 3, 7));
			D[25] = (Angulo(vectR, vectG, vectB, 3, 8));
			disteucl1[3] = D[2] + D[9] + D[15] + D[21] + D[22] + D[23] + D[24] + D[25];

			//i=4,j=0 ya esta es D[3]
			//i=4,j=1 ya esta es D[10]
			//i=4,j=2 ya esta es D[16]
			//i=4,j=3 ya esta es D[21]
			//i=4,j=4 No se hace
			D[26] = (Angulo(vectR, vectG, vectB, 4, 5));
			D[27] = (Angulo(vectR, vectG, vectB, 4, 6));
			D[28] = (Angulo(vectR, vectG, vectB, 4, 7));
			D[29] = (Angulo(vectR, vectG, vectB, 4, 8));
			disteucl1[4] = D[3] + D[10] + D[16] + D[21] + D[26] + D[27] + D[28] + D[29];

			//i=5,j=0 ya esta es D[4]
			//i=5,j=1 ya esta es D[11]
			//i=5,j=2 ya esta es D[17]
			//i=5,j=3 ya esta es D[22]
			//i=5,j=4 ya esta es D[26]
			//i=5,j=5 No se hace
			D[30] = (Angulo(vectR, vectG, vectB, 5, 6));
			D[31] = (Angulo(vectR, vectG, vectB, 5, 7));
			D[32] = (Angulo(vectR, vectG, vectB, 5, 8));
			disteucl1[5] = D[4] + D[11] + D[17] + D[22] + D[26] + D[30] + D[31] + D[32];

			//i=6,j=0 ya esta es D[5]
			//i=6,j=1 ya esta es D[12]
			//i=6,j=2 ya esta es D[18]
			//i=6,j=3 ya esta es D[23]
			//i=6,j=4 ya esta es D[27]
			//i=6,j=5 ya esta es D[30]
			//i=6,j=6 No se hace
			D[33] = (Angulo(vectR, vectG, vectB, 6, 7));
			D[34] = (Angulo(vectR, vectG, vectB, 6, 8));
			disteucl1[6] = D[5] + D[12] + D[18] + D[23] + D[27] + D[30] + D[33] + D[34];

			//i=7,j=0 ya esta es D[6]
			//i=7,j=1 ya esta es D[13]
			//i=7,j=2 ya esta es D[19]
			//i=7,j=3 ya esta es D[24]
			//i=7,j=4 ya esta es D[28]
			//i=7,j=5 ya esta es D[31]
			//i=7,j=6 ya esta es D[33]
			//i=7,j=7 No se hace
			D[35] = (Angulo(vectR, vectG, vectB, 7, 8));
			disteucl1[7] = D[6] + D[13] + D[19] + D[24] + D[28] + D[31] + D[33] + D[35];

			//i=8,j=0 ya esta es D[7]
			//i=8,j=1 ya esta es D[14]
			//i=8,j=2 ya esta es D[20]
			//i=8,j=3 ya esta es D[25]
			//i=8,j=4 ya esta es D[29]
			//i=8,j=5 ya esta es D[32]
			//i=8,j=6 ya esta es D[34]
			//i=8,j=7 ya esta es D[35]
			//i=8,j=8 No se hace
			disteucl1[8] = D[7] + D[14] + D[20] + D[25] + D[29] + D[32] + D[34] + D[35];


			for (F = 0; F <= 8; F++) {
				for (x = 0; x <= 7; x++) {
					if (disteucl1[x] > disteucl1[x + 1]) {
						hold = disteucl1[x];
						hold2 = posicion[x];
						disteucl1[x] = disteucl1[x + 1];
						posicion[x] = posicion[x + 1];
						disteucl1[x + 1] = hold;
						posicion[x + 1] = hold2;
					}
				}
			}

			for (j = 0; j <= nElementos - 1; j++) {
				distMag[j] = 0;
				//posAux_j = posicion[j];
				posAux_j = j;
				posicion[j] = j;
				for (i = 0; i <= nElementos - 1; i++) {
					posAux_i = posicion[i];

					//distMag[j] = (Magnitud_L1(vectR, vectG, vectB, posAux_i, posAux_j)) + distMag[j];
					distMag[j] = (Magnitud_L1(vectR, vectG, vectB, i, j)) + distMag[j];

				}

			}




			mn = distMag[0];
			mx = distMag[0];

			posMin = 0;

			for (i = 0; i <= nElementos - 1; i++)
			{
				if (mn>distMag[i])
				{
					mn = distMag[i];
					posMin = posicion[i];
				}
				else if (mx<distMag[i])
				{

				}
			}




			d_Pout[(Row * m + Col) * 3 + 0] = vectR[posMin];
			d_Pout[(Row * m + Col) * 3 + 1] = vectG[posMin];
			d_Pout[(Row * m + Col) * 3 + 2] = vectB[posMin];

			/*
			d_Pout[(Row * m + Col) * 3 + 0] = (vectR[posicion[0]] + vectR[posicion[1]] + vectR[posicion[2]])/3;
			d_Pout[(Row * m + Col) * 3 + 1] = (vectG[posicion[0]] + vectG[posicion[1]] + vectG[posicion[2]])/3;
			d_Pout[(Row * m + Col) * 3 + 2] = (vectB[posicion[0]] + vectB[posicion[1]] + vectB[posicion[2]])/3;
			*/
		}
	}


}