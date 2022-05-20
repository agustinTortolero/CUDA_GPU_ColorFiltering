#include "stdafx.h"
#include <omp.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>    // std::sort

#define		NormalizeAng	2*(255*255)

float MagnitudL1(unsigned char* VectR, unsigned char* VectG, unsigned char* VectB, unsigned int i, unsigned int j) {

	float distR = (VectR[i] - VectR[j]);
	float distG = (VectG[i] - VectG[j]);
	float distB = (VectB[i] - VectB[j]);

	return sqrt((distR)*(distR)+(distG)*(distG)+(distB)*(distB));
	//return distR + distB + distG;

}
//Gran
float S_shape(float Nabla, unsigned int a, unsigned int b) {

	if (Nabla <= a)		return 0;

	if (a <= Nabla && Nabla <= ((a + b) / 2)) {
		float aux = (Nabla - a) / (b - a);
		return 2 * aux*aux;
	}

	if (((a + b) / 2) <= Nabla && Nabla <= b) {
		float aux = ((Nabla - b) / (b - a));
		return 1 - (2 * aux*aux);
	}

	if (Nabla >= b)		return 1;

}
//Peque
float Z_shape(float Nabla, unsigned int a, unsigned int b) {

	if (Nabla <= a)		return 1;
	if (a <= Nabla && Nabla <= ((a + b) / 2)) {
		float aux = (Nabla - a) / (b - a);
		return 1 - (2 * aux*aux);
	}
	if (((a + b) / 2) <= Nabla && Nabla <= b) {
		float aux = (Nabla - b) / (b - a);
		return 2 * aux*aux;
	}
	if (Nabla >= b)		return 0;
}


void FTSCF_CPU_Multi(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads)
{
	int Col = 0, Row = 0, c = 0, d = 0, i = 0, j = 0;

	int x = 0, posicion[9], hold2 = 0, F = 0;
	unsigned char vectR[25], vectG[25], vectB[25];
	float disteucl = 0.0, disteucl1[9], hold;
	float D[45];
	float mn, mx;
	int posMin = 0;

	float uGran[3], uPeque[3], rs[9], r = 0;
	const unsigned char a = 10, b = 60;

	posicion[0] = 6; posicion[1] = 7; posicion[2] = 8; posicion[3] = 11; posicion[4] = 12;
	posicion[5] = 13; posicion[6] = 16; posicion[7] = 17; posicion[8] = 18;



	for (Row = 3; Row < n - 3; Row++) {
		//#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F ,x,disteucl1,disteucl,vectR,vectG,vectB,hold,hold2,posicion,posMin,mn,mnx) shared(d_Pout, d_Pin,n, m,channels,Row ) schedule(static)

		for (Col = 3; Col < m - 3; Col++) {
			//int tid = omp_get_thread_num();
			//hacer el arreglo
			F = 0;

			for (i = -2; i <= 2; i++) {
				for (j = -2; j <= 2; j++) {
					vectR[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
					vectG[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
					vectB[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];


					F++;
				}
			}

			//NW
			D[0] = (MagnitudL1(vectR, vectG, vectB, 12, 6));
			uGran[0] = S_shape(D[0], a, b);

			D[1] = (MagnitudL1(vectR, vectG, vectB, 12, 8));
			uGran[1] = S_shape(D[1], a, b);

			D[2] = (MagnitudL1(vectR, vectG, vectB, 12, 16));
			uGran[2] = S_shape(D[2], a, b);

			D[3] = (MagnitudL1(vectR, vectG, vectB, 16, 10));
			uPeque[0] = Z_shape(D[3], a, b);

			D[4] = (MagnitudL1(vectR, vectG, vectB, 2, 8));
			uPeque[1] = Z_shape(D[4], a, b);

			rs[0] = uGran[0] * uGran[1] * uGran[2] * uPeque[0] * uPeque[1];

			//N
			D[5] = (MagnitudL1(vectR, vectG, vectB, 12, 7));
			uGran[0] = S_shape(D[5], a, b);

			D[6] = (MagnitudL1(vectR, vectG, vectB, 12, 13));
			uGran[1] = S_shape(D[6], a, b);

			D[7] = (MagnitudL1(vectR, vectG, vectB, 12, 11));
			uGran[2] = S_shape(D[7], a, b);

			D[8] = (MagnitudL1(vectR, vectG, vectB, 11, 6));
			uPeque[0] = Z_shape(D[8], a, b);

			D[9] = (MagnitudL1(vectR, vectG, vectB, 8, 13));
			uPeque[1] = Z_shape(D[9], a, b);

			rs[1] = uGran[0] * uGran[1] * uGran[2] * uPeque[0] * uPeque[1];

			//NE
			//D[10] = (MagnitudL1(vectR, vectG, vectB, 12, 8));
			// es D[1]
			uGran[0] = S_shape(D[1], a, b);

			//D[11] = (MagnitudL1(vectR, vectG, vectB, 12, 6));
			// es D[0]
			uGran[1] = S_shape(D[0], a, b);

			D[10] = (MagnitudL1(vectR, vectG, vectB, 12, 18));
			uGran[2] = S_shape(D[10], a, b);

			D[11] = (MagnitudL1(vectR, vectG, vectB, 18, 14));
			uPeque[0] = Z_shape(D[11], a, b);

			D[12] = (MagnitudL1(vectR, vectG, vectB, 6, 2));
			uPeque[1] = Z_shape(D[12], a, b);

			rs[2] = uGran[0] * uGran[1] * uGran[2] * uPeque[0] * uPeque[1];

			//E			
			//D[15] = (MagnitudL1(vectR, vectG, vectB, 12, 13));
			//es D[6]
			uGran[0] = S_shape(D[6], a, b);

			//D[16] = (MagnitudL1(vectR, vectG, vectB, 12, 7));
			//es D[5]
			uGran[1] = S_shape(D[5], a, b);

			D[13] = (MagnitudL1(vectR, vectG, vectB, 12, 17));
			uGran[2] = S_shape(D[13], a, b);

			D[14] = (MagnitudL1(vectR, vectG, vectB, 7, 8));
			uPeque[0] = Z_shape(D[14], a, b);

			D[15] = (MagnitudL1(vectR, vectG, vectB, 17, 18));
			uPeque[1] = Z_shape(D[15], a, b);

			rs[3] = uGran[0] * uGran[1] * uGran[2] * uPeque[0] * uPeque[1];

			//SE
			//D[20] = (MagnitudL1(vectR, vectG, vectB, 12, 18));
			//es D[10]
			uGran[0] = S_shape(D[10], a, b);

			//es D[2]
			//D[21] = (MagnitudL1(vectR, vectG, vectB, 12, 16));
			uGran[1] = S_shape(D[2], a, b);

			//es D[1]
			//D[22] = (MagnitudL1(vectR, vectG, vectB, 12, 8));
			uGran[2] = S_shape(D[1], a, b);

			D[16] = (MagnitudL1(vectR, vectG, vectB, 16, 22));
			uPeque[0] = Z_shape(D[16], a, b);

			D[17] = (MagnitudL1(vectR, vectG, vectB, 8, 14));
			uPeque[1] = Z_shape(D[17], a, b);

			rs[4] = uGran[0] * uGran[1] * uGran[2] * uPeque[0] * uPeque[1];

			//S
			//D[18] = (MagnitudL1(vectR, vectG, vectB, 12, 17));
			//es D[13]
			uGran[0] = S_shape(D[13], a, b);
			//es D[7]
			//D[26] = (MagnitudL1(vectR, vectG, vectB, 12, 11));
			uGran[1] = S_shape(D[7], a, b);
			//es D[6]
			//D[27] = (MagnitudL1(vectR, vectG, vectB, 12, 13));
			uGran[2] = S_shape(D[6], a, b);

			D[18] = (MagnitudL1(vectR, vectG, vectB, 11, 16));
			uPeque[0] = Z_shape(D[18], a, b);

			D[19] = (MagnitudL1(vectR, vectG, vectB, 13, 18));
			uPeque[1] = Z_shape(D[19], a, b);

			rs[5] = uGran[0] * uGran[1] * uGran[2] * uPeque[0] * uPeque[1];

			//SW
			//es D[2]
			//D[30] = (MagnitudL1(vectR, vectG, vectB, 12, 16));
			uGran[0] = S_shape(D[2], a, b);
			//es D[0]
			//D[31] = (MagnitudL1(vectR, vectG, vectB, 12, 6));
			uGran[1] = S_shape(D[0], a, b);
			//es D[10]
			//D[32] = (MagnitudL1(vectR, vectG, vectB, 12, 18));
			uGran[2] = S_shape(D[10], a, b);

			D[20] = (MagnitudL1(vectR, vectG, vectB, 6, 10));
			uPeque[0] = Z_shape(D[20], a, b);

			D[21] = (MagnitudL1(vectR, vectG, vectB, 18, 22));
			uPeque[1] = Z_shape(D[21], a, b);

			rs[6] = uGran[0] * uGran[1] * uGran[2] * uPeque[0] * uPeque[1];

			//W
			//Es D[7]
			//D[35] = (MagnitudL1(vectR, vectG, vectB, 12, 11));
			uGran[0] = S_shape(D[7], a, b);
			//Es D[5]
			//D[36] = (MagnitudL1(vectR, vectG, vectB, 12, 7));
			uGran[1] = S_shape(D[5], a, b);
			//es D[13]
			//D[37] = (MagnitudL1(vectR, vectG, vectB, 12, 17));
			uGran[2] = S_shape(D[13], a, b);

			D[21] = (MagnitudL1(vectR, vectG, vectB, 6, 7));
			uPeque[0] = Z_shape(D[21], a, b);

			D[22] = (MagnitudL1(vectR, vectG, vectB, 16, 17));
			uPeque[1] = Z_shape(D[22], a, b);

			rs[7] = uGran[0] * uGran[1] * uGran[2] * uPeque[0] * uPeque[1];

			mn = rs[0];
			r = rs[0];

			for (i = 0; i <= 7; i++)
			{
				if (r<rs[i])
				{
					r = rs[i];

				}
			}

			//Filtro VMF
			const float THS = .2;
			if (r > THS) {

				D[23] = (MagnitudL1(vectR, vectG, vectB, 6, 8));
				D[24] = (MagnitudL1(vectR, vectG, vectB, 6, 13));
				D[25] = (MagnitudL1(vectR, vectG, vectB, 6, 16));
				D[26] = (MagnitudL1(vectR, vectG, vectB, 6, 17));
				D[27] = (MagnitudL1(vectR, vectG, vectB, 6, 18));

				disteucl1[0] = D[0] + D[8] + D[21] + D[23] + D[24] + D[25] + D[26] + D[27];

				D[28] = (MagnitudL1(vectR, vectG, vectB, 7, 11));
				D[29] = (MagnitudL1(vectR, vectG, vectB, 7, 13));
				D[30] = (MagnitudL1(vectR, vectG, vectB, 7, 16));
				D[31] = (MagnitudL1(vectR, vectG, vectB, 7, 17));
				D[32] = (MagnitudL1(vectR, vectG, vectB, 7, 18));
				disteucl1[1] = D[21] + D[5] + D[14] + D[28] + D[29] + D[30] + D[31] + D[32];

				//es D[26] D[33] = (MagnitudL1(vectR, vectG, vectB, 8, 6));
				D[33] = (MagnitudL1(vectR, vectG, vectB, 8, 11));
				D[34] = (MagnitudL1(vectR, vectG, vectB, 8, 16));
				D[35] = (MagnitudL1(vectR, vectG, vectB, 8, 17));
				D[36] = (MagnitudL1(vectR, vectG, vectB, 8, 18));
				disteucl1[2] = D[14] + D[1] + D[9] + D[26] + D[33] + D[34] + D[35] + D[36];

				//es D[28]  D[37] = (MagnitudL1(vectR, vectG, vectB, 11, 7));
				//es D[33]     D[38] = (MagnitudL1(vectR, vectG, vectB, 11, 8)); 
				D[37] = (MagnitudL1(vectR, vectG, vectB, 11, 13));
				D[38] = (MagnitudL1(vectR, vectG, vectB, 11, 17));
				D[39] = (MagnitudL1(vectR, vectG, vectB, 11, 18));
				disteucl1[3] = D[7] + D[8] + D[18] + D[28] + D[33] + D[37] + D[38] + D[39];

				//Central ya estan todas las d calculadas
				disteucl1[4] = D[0] + D[5] + D[1] + D[7] + D[6] + D[2] + D[13] + D[10];

				D[40] = (MagnitudL1(vectR, vectG, vectB, 13, 16));
				D[41] = (MagnitudL1(vectR, vectG, vectB, 13, 17));
				disteucl1[5] = D[6] + D[19] + D[9] + D[24] + D[29] + D[37] + D[40] + D[41];

				D[42] = (MagnitudL1(vectR, vectG, vectB, 16, 18));
				disteucl1[6] = D[18] + D[2] + D[22] + D[25] + D[30] + D[34] + D[40] + D[42];

				disteucl1[7] = D[22] + D[13] + D[15] + D[26] + D[31] + D[35] + D[38] + D[41];

				disteucl1[8] = D[19] + D[10] + D[15] + D[27] + D[32] + D[36] + D[39] + D[42];

				posMin = 6;
				mn = disteucl1[0];
				for (i = 0; i <= 7; i++)
				{
					if (mn>disteucl1[i])
					{
						mn = disteucl1[i];
						posMin = posicion[i];
					}

				}

				d_Pout[(Row * m + Col) * 3 + 0] = vectR[posMin];
				d_Pout[(Row * m + Col) * 3 + 1] = vectG[posMin];
				d_Pout[(Row * m + Col) * 3 + 2] = vectB[posMin];

				/*
				d_Pout[(Row * m + Col) * 3 + 0] = 255;
				d_Pout[(Row * m + Col) * 3 + 1] = 255;
				d_Pout[(Row * m + Col) * 3 + 2] = 255;
				*/

			}
			else {
				// si no es ruido la salida el el pixel central de la ventana
				d_Pout[(Row * m + Col) * 3 + 0] = vectR[12];
				d_Pout[(Row * m + Col) * 3 + 1] = vectG[12];
				d_Pout[(Row * m + Col) * 3 + 2] = vectB[12];

				/*
				d_Pout[(Row * m + Col) * 3 + 0] = 255;
				d_Pout[(Row * m + Col) * 3 + 1] = 255;
				d_Pout[(Row * m + Col) * 3 + 2] = 255;
				*/

			}
			/*
			d_Pout[(Row * m + Col) * 3 + 0] = r * 255;
			d_Pout[(Row * m + Col) * 3 + 1] = r * 255;
			d_Pout[(Row * m + Col) * 3 + 2] = r * 255;
			*/
		}
	}
}



float Big_Nabla(float nabla) {

	const float med1 = 60, var1=1000;

	if (nabla > med1 )		return 1;

	float aux = nabla - med1;

	return exp( -1 * (aux*aux)/ (2*var1));
}
float Small_Nabla(float nabla) {
	const float med2 = 10, var1 = 1000;

	if (nabla < med2)		return 1;

	float aux = nabla - med2;

	return exp(-1 * (aux*aux) / (2*var1));
}

float Small_Theta(float theta) {
	const float med1 = .1, var2 = .8;

	if (theta < med1)		return 1;

	float aux = theta - med1;

	return exp(-1 * (aux*aux) / (2*var2));
}

float Big_Theta(float theta) {
	const float med2 = .6, var2=.8;

	if (theta > med2)		return 1;

	float aux = theta - med2;

	return exp(-1 * (aux*aux) / (2*var2));
}

float Angle(float a, float b) {
	//Recordar que NormalizeAng	= 2*(255*255)
	float arriva = NormalizeAng + a*b;
	float abajo = sqrt(NormalizeAng + (a*a)) * sqrt(NormalizeAng + (b*b));
	float aux = arriva / abajo;
	if (aux >= 1)	aux = 1;
	float result = acos(aux);
	return result;
}

float Nabla(float a, float b) {
	return abs (a-b);
}


void swap(float *xp, float *yp)
{
	float temp = *xp; //nota aqui temp debe de ser float
	*xp = *yp;
	*yp = temp;
}

// Bubble Sort menor a mayor
void bubbleSort(float arr[], float rho[], int n)
{
	int i, j;
	for (i = 0; i < n - 1; i++)
		// Last i elements are already in place  
		for (j = 0; j < n - i - 1; j++)
			if (arr[j] > arr[j + 1]) {
				swap(&arr[j], &arr[j + 1]);
				swap(&rho[j], &rho[j + 1]);
			}
}



void FTSCF_CPU_Multi2(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads)
{
	int Col = 0, Row = 0, c = 0, d = 0, i = 0, j = 0;

	int x = 0, posicion[9], hold2 = 0, F = 0;
	float vectR[25], vectG[25], vectB[25];
	float vectR_Fil[9], vectG_Fil[9], vectB_Fil[9]; // falta agregar esto a la linea de openMP

	float disteucl = 0.0, disteucl1[9], hold;
	float mn, mx;
	int posMin = 0;

	float rs_R[8], rs_G[8], rs_B[8], r = 0, r_R = 0, r_G = 0, r_B = 0;
	float uBigR[3], uSmallR[3], uBigG[3], uSmallG[3], uBigB[3], uSmallB[3];

	float aux1R ,aux2R, aux1G, aux2G, aux1B , aux2B;

	float Ruido_R=0, Ruido_G = 0, Ruido_B = 0;

	float Pc_R, Pc_G, Pc_B, aux1=0,aux2=0;

	float rhoSumatoria = 0, rhoTotal = 0, sum = 0;

	//filtrado
	float rho[9], AuxArray[9];
	float sumR = 0;



	for (Row = 3; Row < n - 3; Row++) {
		//#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F ,x,disteucl1,disteucl,vectR,vectG,vectB,hold,hold2,posicion,posMin,mn,mnx) shared(d_Pout, d_Pin,n, m,channels,Row ) schedule(static)

		for (Col = 3; Col < m - 3; Col++) {
			//int tid = omp_get_thread_num();
			//hacer el arreglo
			F = 0;
			r = 0, r_R = 0, r_G = 0, r_B = 0;
			for (i = -2; i <= 2; i++) {
				for (j = -2; j <= 2; j++) {
					vectR[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
					vectG[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
					vectB[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
					F++;
				}
			}

			Pc_R = vectR[12];
			Pc_G = vectG[12];
			Pc_B = vectB[12];

			//N Nablas
			uBigR[0]	= Big_Nabla(Nabla(Pc_R, vectR[7]));
			uBigG[0]	= Big_Nabla(Nabla(Pc_G, vectG[7])); //Basico_1
			uBigB[0]	= Big_Nabla(Nabla(Pc_B, vectB[7]));

			uBigR[1]	= Big_Nabla(Nabla(Pc_R, vectR[11]));
			uBigG[1]	= Big_Nabla(Nabla(Pc_G, vectG[11])); //Basico_2
			uBigB[1]	= Small_Nabla(Nabla(Pc_B, vectB[11]));

			uBigR[2]	= Big_Nabla(Nabla(Pc_R, vectR[13]));
			uBigG[2]	= Big_Nabla(Nabla(Pc_G, vectG[13])); //Basico_3
			uBigB[2]	= Big_Nabla(Nabla(Pc_B, vectB[13]));

			uSmallR[0]	= Small_Nabla(Nabla(vectR[6], vectR[11]));
			uSmallG[0]	= Small_Nabla(Nabla(vectG[6], vectG[11])); //rel_1
			uSmallB[0]	= Small_Nabla(Nabla(vectB[6], vectB[11]));

			uSmallR[1]	= Small_Nabla(Nabla(vectR[8], vectR[13]));
			uSmallG[1]	= Small_Nabla(Nabla(vectG[8], vectG[13])); //rel_2
			uSmallB[1]	= Small_Nabla(Nabla(vectB[8], vectB[13]));
			
			aux1R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux1G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux1B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			//N Angles
			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[7]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[7])); //Basico_1
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[7]));

			uBigR[1] = Big_Theta(Angle(Pc_R, vectR[11]));
			uBigG[1] = Big_Theta(Angle(Pc_G, vectG[11])); //Basico_2
			uBigB[1] = Big_Theta(Angle(Pc_B, vectB[11]));

			uBigR[2] = Big_Theta(Angle(Pc_R, vectR[13]));
			uBigG[2] = Big_Theta(Angle(Pc_G, vectG[13])); //Basico_3
			uBigB[2] = Big_Theta(Angle(Pc_B, vectB[13]));

			uSmallR[0] = Small_Theta(Angle(vectR[6], vectR[11]));
			uSmallG[0] = Small_Theta(Angle(vectG[6], vectG[11])); //rel_1
			uSmallB[0] = Small_Theta(Angle(vectB[6], vectB[11]));

			uSmallR[1] = Small_Theta(Angle(vectR[8], vectR[13]));
			uSmallG[1] = Small_Theta(Angle(vectG[8], vectG[13])); //rel_2
			uSmallB[1] = Small_Theta(Angle(vectB[8], vectB[13])); /*aqui esta el error da nan*/


			aux2R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux2G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux2B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			rs_R[0] = std::min(aux1R, aux2R);
			rs_G[0] = std::min(aux1G, aux2G);
			rs_B[0] = std::min(aux1B, aux2B);

			//NE Nablas
			uBigR[0] = Big_Nabla(Nabla(Pc_R, vectR[8]));
			uBigG[0] = Big_Nabla(Nabla(Pc_G, vectG[8])); //Basico_1
			uBigB[0] = Big_Nabla(Nabla(Pc_B, vectB[8]));

			uBigR[1] = Big_Nabla(Nabla(Pc_R, vectR[6]));
			uBigG[1] = Big_Nabla(Nabla(Pc_G, vectG[6])); //Basico_2
			uBigB[1] = Big_Nabla(Nabla(Pc_B, vectB[6]));

			uBigR[2] = Big_Nabla(Nabla(Pc_R, vectR[18]));
			uBigG[2] = Big_Nabla(Nabla(Pc_G, vectG[18])); //Basico_3
			uBigB[2] = Big_Nabla(Nabla(Pc_B, vectB[18]));

			uSmallR[0] = Small_Nabla(Nabla(vectR[2], vectR[6]));
			uSmallG[0] = Small_Nabla(Nabla(vectG[2], vectG[6])); //rel_1
			uSmallB[0] = Small_Nabla(Nabla(vectB[2], vectB[6]));

			uSmallR[1] = Small_Nabla(Nabla(vectR[18], vectR[14]));
			uSmallG[1] = Small_Nabla(Nabla(vectG[18], vectG[14])); //rel_2
			uSmallB[1] = Small_Nabla(Nabla(vectB[18], vectB[14]));

			aux1R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux1G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux1B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			//NE Angles
			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[8]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[8])); //Basico_1
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[8]));

			uBigR[1] = Big_Theta(Angle(Pc_R, vectR[6]));
			uBigG[1] = Big_Theta(Angle(Pc_G, vectG[6])); //Basico_2
			uBigB[1] = Big_Theta(Angle(Pc_B, vectB[6]));

			uBigR[2] = Big_Theta(Angle(Pc_R, vectR[18]));
			uBigG[2] = Big_Theta(Angle(Pc_G, vectG[18])); //Basico_3
			uBigB[2] = Big_Theta(Angle(Pc_B, vectB[18]));

			uSmallR[0] = Small_Theta(Angle(vectR[2], vectR[6]));
			uSmallG[0] = Small_Theta(Angle(vectG[2], vectG[6])); //rel_1
			uSmallB[0] = Small_Theta(Angle(vectB[2], vectB[6]));

			uSmallR[1] = Small_Theta(Angle(vectR[18], vectR[14]));
			uSmallG[1] = Small_Theta(Angle(vectG[18], vectG[14])); //rel_2
			uSmallB[1] = Small_Theta(Angle(vectB[18], vectB[14]));
			
			aux2R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux2G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux2B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			rs_R[1] = std::min(aux1R, aux2R);
			rs_G[1] = std::min(aux1G, aux2G);
			rs_B[1] = std::min(aux1B, aux2B);

			//E
			uBigR[0] = Big_Nabla(Nabla(Pc_R, vectR[13]));
			uBigG[0] = Big_Nabla(Nabla(Pc_G, vectG[13])); //Basico_1
			uBigB[0] = Big_Nabla(Nabla(Pc_B, vectB[13]));

			uBigR[1] = Big_Nabla(Nabla(Pc_R, vectR[17]));
			uBigG[1] = Big_Nabla(Nabla(Pc_G, vectG[17])); //Basico_2
			uBigB[1] = Big_Nabla(Nabla(Pc_B, vectB[17]));

			uBigR[2] = Big_Nabla(Nabla(Pc_R, vectR[7]));
			uBigG[2] = Big_Nabla(Nabla(Pc_G, vectG[7])); //Basico_3
			uBigB[2] = Big_Nabla(Nabla(Pc_B, vectB[7]));

			uSmallR[0] = Small_Nabla(Nabla(vectR[7], vectR[8]));
			uSmallG[0] = Small_Nabla(Nabla(vectG[7], vectG[8])); //rel_1
			uSmallB[0] = Small_Nabla(Nabla(vectB[7], vectB[8]));

			uSmallR[1] = Small_Nabla(Nabla(vectR[17], vectR[18]));
			uSmallG[1] = Small_Nabla(Nabla(vectG[17], vectG[18])); //rel_2
			uSmallB[1] = Small_Nabla(Nabla(vectB[17], vectB[18]));

			aux1R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux1G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux1B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			//E Angles
			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[13]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[13])); //Basico_1
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[13]));

			uBigR[1] = Big_Theta(Angle(Pc_R, vectR[17]));
			uBigG[1] = Big_Theta(Angle(Pc_G, vectG[17])); //Basico_2
			uBigB[1] = Big_Theta(Angle(Pc_B, vectB[17]));

			uBigR[2] = Big_Theta(Angle(Pc_R, vectR[7]));
			uBigG[2] = Big_Theta(Angle(Pc_G, vectG[7])); //Basico_3
			uBigB[2] = Big_Theta(Angle(Pc_B, vectB[7]));

			uSmallR[0] = Small_Theta(Angle(vectR[7], vectR[8]));
			uSmallG[0] = Small_Theta(Angle(vectG[7], vectG[8])); //rel_1
			uSmallB[0] = Small_Theta(Angle(vectB[7], vectB[8]));

			uSmallR[1] = Small_Theta(Angle(vectR[17], vectR[18]));
			uSmallG[1] = Small_Theta(Angle(vectG[17], vectG[18])); //rel_2
			uSmallB[1] = Small_Theta(Angle(vectB[17], vectB[18]));

			aux2R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux2G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux2B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			rs_R[2] = std::min(aux1R, aux2R);
			rs_G[2] = std::min(aux1G, aux2G);
			rs_B[2] = std::min(aux1B, aux2B);


			//SE Nablas
			uBigR[0] = Big_Nabla(Nabla(Pc_R, vectR[18]));
			uBigG[0] = Big_Nabla(Nabla(Pc_G, vectG[18])); //Basico_1
			uBigB[0] = Big_Nabla(Nabla(Pc_B, vectB[18]));

			uBigR[1] = Big_Nabla(Nabla(Pc_R, vectR[8]));
			uBigG[1] = Big_Nabla(Nabla(Pc_G, vectG[8])); //Basico_2
			uBigB[1] = Big_Nabla(Nabla(Pc_B, vectB[8]));

			uBigR[2] = Big_Nabla(Nabla(Pc_R, vectR[16]));
			uBigG[2] = Big_Nabla(Nabla(Pc_G, vectG[16])); //Basico_3
			uBigB[2] = Big_Nabla(Nabla(Pc_B, vectB[16]));

			uSmallR[0] = Small_Nabla(Nabla(vectR[8], vectR[14]));
			uSmallG[0] = Small_Nabla(Nabla(vectG[8], vectG[14])); //rel_1
			uSmallB[0] = Small_Nabla(Nabla(vectB[8], vectB[14]));

			uSmallR[1] = Small_Nabla(Nabla(vectR[16], vectR[22]));
			uSmallG[1] = Small_Nabla(Nabla(vectG[16], vectG[22])); //rel_2
			uSmallB[1] = Small_Nabla(Nabla(vectB[16], vectB[22]));

			aux1R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux1G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux1B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			//SE Angles
			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[18]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[18])); //Basico_1
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[18]));

			uBigR[1] = Big_Theta(Angle(Pc_R, vectR[8]));
			uBigG[1] = Big_Theta(Angle(Pc_G, vectG[8])); //Basico_2
			uBigB[1] = Big_Theta(Angle(Pc_B, vectB[8]));

			uBigR[2] = Big_Theta(Angle(Pc_R, vectR[16]));
			uBigG[2] = Big_Theta(Angle(Pc_G, vectG[16])); //Basico_3
			uBigB[2] = Big_Theta(Angle(Pc_B, vectB[16]));

			uSmallR[0] = Small_Theta(Angle(vectR[8], vectR[14]));
			uSmallG[0] = Small_Theta(Angle(vectG[8], vectG[14])); //rel_1
			uSmallB[0] = Small_Theta(Angle(vectB[8], vectB[14]));

			uSmallR[1] = Small_Theta(Angle(vectR[16], vectR[22]));
			uSmallG[1] = Small_Theta(Angle(vectG[16], vectG[22])); //rel_2
			uSmallB[1] = Small_Theta(Angle(vectB[16], vectB[22]));

			aux2R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux2G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux2B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			rs_R[3] = std::min(aux1R, aux2R);
			rs_G[3] = std::min(aux1G, aux2G);
			rs_B[3] = std::min(aux1B, aux2B);

			//S Nabla
			uBigR[0] = Big_Nabla(Nabla(Pc_R, vectR[17]));
			uBigG[0] = Big_Nabla(Nabla(Pc_G, vectG[17])); //Basico_1
			uBigB[0] = Big_Nabla(Nabla(Pc_B, vectB[17]));

			uBigR[1] = Big_Nabla(Nabla(Pc_R, vectR[11]));
			uBigG[1] = Big_Nabla(Nabla(Pc_G, vectG[11])); //Basico_2
			uBigB[1] = Big_Nabla(Nabla(Pc_B, vectB[11]));

			uBigR[2] = Big_Nabla(Nabla(Pc_R, vectR[13]));
			uBigG[2] = Big_Nabla(Nabla(Pc_G, vectG[13])); //Basico_3
			uBigB[2] = Big_Nabla(Nabla(Pc_B, vectB[13]));

			uSmallR[0] = Small_Nabla(Nabla(vectR[11], vectR[16]));
			uSmallG[0] = Small_Nabla(Nabla(vectG[11], vectG[16])); //rel_1
			uSmallB[0] = Small_Nabla(Nabla(vectB[11], vectB[16]));

			uSmallR[1] = Small_Nabla(Nabla(vectR[13], vectR[18]));
			uSmallG[1] = Small_Nabla(Nabla(vectG[13], vectG[18])); //rel_2
			uSmallB[1] = Small_Nabla(Nabla(vectB[13], vectB[18]));

			aux1R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux1G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux1B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			//S Angles
			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[17]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[17])); //Basico_1
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[17]));

			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[11]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[11])); //Basico_2
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[11]));

			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[13]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[13])); //Basico_3
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[13]));

			uSmallR[0] = Small_Theta(Angle(vectR[11], vectR[16]));
			uSmallG[0] = Small_Theta(Angle(vectG[11], vectG[16])); //rel_1
			uSmallB[0] = Small_Theta(Angle(vectB[11], vectB[16]));

			uSmallR[1] = Small_Theta(Angle(vectR[13], vectR[18]));
			uSmallG[1] = Small_Theta(Angle(vectG[13], vectG[18])); //rel_2
			uSmallB[1] = Small_Theta(Angle(vectB[13], vectB[18]));

			aux2R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux2G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux2B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			rs_R[4] = std::min(aux1R, aux2R);
			rs_G[4] = std::min(aux1G, aux2G);
			rs_B[4] = std::min(aux1B, aux2B);

			//SW Nabla
			uBigR[0] = Big_Nabla(Nabla(Pc_R, vectR[16]));
			uBigG[0] = Big_Nabla(Nabla(Pc_G, vectG[16])); //Basico_1
			uBigB[0] = Big_Nabla(Nabla(Pc_B, vectB[16]));

			uBigR[1] = Big_Nabla(Nabla(Pc_R, vectR[18]));
			uBigG[1] = Big_Nabla(Nabla(Pc_G, vectG[18])); //Basico_2
			uBigB[1] = Big_Nabla(Nabla(Pc_B, vectB[18]));

			uBigR[2] = Big_Nabla(Nabla(Pc_R, vectR[6]));
			uBigG[2] = Big_Nabla(Nabla(Pc_G, vectG[6])); //Basico_3
			uBigB[2] = Big_Nabla(Nabla(Pc_B, vectB[6]));

			uSmallR[0] = Small_Nabla(Nabla(vectR[6], vectR[10])); //este da demasiado bajo (.00188)
			uSmallG[0] = Small_Nabla(Nabla(vectG[6], vectG[10])); //rel_1
			uSmallB[0] = Small_Nabla(Nabla(vectB[6], vectB[10]));

			uSmallR[1] = Small_Nabla(Nabla(vectR[22], vectR[18]));
			uSmallG[1] = Small_Nabla(Nabla(vectG[22], vectG[18])); //rel_2
			uSmallB[1] = Small_Nabla(Nabla(vectB[22], vectB[18]));

			aux1R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux1G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux1B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			//SW Angles
			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[16]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[16])); //Basico_1
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[16]));

			uBigR[1] = Big_Theta(Angle(Pc_R, vectR[18]));
			uBigG[1] = Big_Theta(Angle(Pc_G, vectG[18])); //Basico_2
			uBigB[1] = Big_Theta(Angle(Pc_B, vectB[18]));

			uBigR[2] = Big_Theta(Angle(Pc_R, vectR[6]));
			uBigG[2] = Big_Theta(Angle(Pc_G, vectG[6])); //Basico_3
			uBigB[2] = Big_Theta(Angle(Pc_B, vectB[6]));

			uSmallR[0] = Small_Theta(Angle(vectR[6], vectR[10]));
			uSmallG[0] = Small_Theta(Angle(vectG[6], vectG[10])); //rel_1
			uSmallB[0] = Small_Theta(Angle(vectB[6], vectB[10]));

			uSmallR[1] = Small_Theta(Angle(vectR[22], vectR[18]));
			uSmallG[1] = Small_Theta(Angle(vectG[22], vectG[18])); //rel_2
			uSmallB[1] = Small_Theta(Angle(vectB[22], vectB[18]));

			aux2R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux2G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux2B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			rs_R[5] = std::min(aux1R, aux2R);
			rs_G[5] = std::min(aux1G, aux2G);
			rs_B[5] = std::min(aux1B, aux2B);

			//W Nabla
			uBigR[0] = Big_Nabla(Nabla(Pc_R, vectR[11]));
			uBigG[0] = Big_Nabla(Nabla(Pc_G, vectG[11])); //Basico_1
			uBigB[0] = Big_Nabla(Nabla(Pc_B, vectB[11]));

			uBigR[1] = Big_Nabla(Nabla(Pc_R, vectR[7]));
			uBigG[1] = Big_Nabla(Nabla(Pc_G, vectG[7])); //Basico_2
			uBigB[1] = Big_Nabla(Nabla(Pc_B, vectB[7]));

			uBigR[2] = Big_Nabla(Nabla(Pc_R, vectR[17]));
			uBigG[2] = Big_Nabla(Nabla(Pc_G, vectG[17])); //Basico_3
			uBigB[2] = Big_Nabla(Nabla(Pc_B, vectB[17]));

			uSmallR[0] = Small_Nabla(Nabla(vectR[6], vectR[7]));
			uSmallG[0] = Small_Nabla(Nabla(vectG[6], vectG[7])); //rel_1
			uSmallB[0] = Small_Nabla(Nabla(vectB[6], vectB[7]));

			uSmallR[1] = Small_Nabla(Nabla(vectR[16], vectR[17]));
			uSmallG[1] = Small_Nabla(Nabla(vectG[16], vectG[17])); //rel_2
			uSmallB[1] = Small_Nabla(Nabla(vectB[16], vectB[17]));

			aux1R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux1G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux1B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			//W Angles
			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[11]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[11])); //Basico_1
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[11]));

			uBigR[1] = Big_Theta(Angle(Pc_R, vectR[7]));
			uBigG[1] = Big_Theta(Angle(Pc_G, vectG[7])); //Basico_2
			uBigB[1] = Big_Theta(Angle(Pc_B, vectB[7]));

			uBigR[2] = Big_Theta(Angle(Pc_R, vectR[17]));
			uBigG[2] = Big_Theta(Angle(Pc_G, vectG[17])); //Basico_3
			uBigB[2] = Big_Theta(Angle(Pc_B, vectB[17]));

			uSmallR[0] = Small_Theta(Angle(vectR[6], vectR[7]));
			uSmallG[0] = Small_Theta(Angle(vectG[6], vectG[7])); //rel_1
			uSmallB[0] = Small_Theta(Angle(vectB[6], vectB[7]));

			uSmallR[1] = Small_Theta(Angle(vectR[16], vectR[17]));
			uSmallG[1] = Small_Theta(Angle(vectG[16], vectG[17])); //rel_2
			uSmallB[1] = Small_Theta(Angle(vectB[16], vectB[17]));


			aux2R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux2G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux2B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			rs_R[6] = std::min(aux1R, aux2R);
			rs_G[6] = std::min(aux1G, aux2G);
			rs_B[6] = std::min(aux1B, aux2B);

			//NW Nabla
			uBigR[0] = Big_Nabla(Nabla(Pc_R, vectR[6]));
			uBigG[0] = Big_Nabla(Nabla(Pc_G, vectG[6])); //Basico_1
			uBigB[0] = Big_Nabla(Nabla(Pc_B, vectB[6]));

			uBigR[1] = Big_Nabla(Nabla(Pc_R, vectR[16]));
			uBigG[1] = Big_Nabla(Nabla(Pc_G, vectG[16])); //Basico_2
			uBigB[1] = Big_Nabla(Nabla(Pc_B, vectB[16]));

			uBigR[2] = Big_Nabla(Nabla(Pc_R, vectR[8]));
			uBigG[2] = Big_Nabla(Nabla(Pc_G, vectG[8])); //Basico_3
			uBigB[2] = Big_Nabla(Nabla(Pc_B, vectB[8]));

			uSmallR[0] = Small_Nabla(Nabla(vectR[2], vectR[8]));
			uSmallG[0] = Small_Nabla(Nabla(vectG[2], vectG[8])); //rel_1
			uSmallB[0] = Small_Nabla(Nabla(vectB[2], vectB[8]));

			uSmallR[1] = Small_Nabla(Nabla(vectR[16], vectR[10]));
			uSmallG[1] = Small_Nabla(Nabla(vectG[16], vectG[10])); //rel_2
			uSmallB[1] = Small_Nabla(Nabla(vectB[16], vectB[10]));

			aux1R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux1G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux1B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			//NW Angles
			uBigR[0] = Big_Theta(Angle(Pc_R, vectR[6]));
			uBigG[0] = Big_Theta(Angle(Pc_G, vectG[6])); //Basico_1
			uBigB[0] = Big_Theta(Angle(Pc_B, vectB[6]));

			uBigR[1] = Big_Theta(Angle(Pc_R, vectR[16]));
			uBigG[1] = Big_Theta(Angle(Pc_G, vectG[16])); //Basico_2
			uBigB[1] = Big_Theta(Angle(Pc_B, vectB[16]));

			uBigR[2] = Big_Theta(Angle(Pc_R, vectR[8]));
			uBigG[2] = Big_Theta(Angle(Pc_G, vectG[8])); //Basico_3
			uBigB[2] = Big_Theta(Angle(Pc_B, vectB[8]));

			uSmallR[0] = Small_Theta(Angle(vectR[2], vectR[8]));
			uSmallG[0] = Small_Theta(Angle(vectG[2], vectG[8])); //rel_1
			uSmallB[0] = Small_Theta(Angle(vectB[2], vectB[8]));

			uSmallR[1] = Small_Theta(Angle(vectR[16], vectR[10]));
			uSmallG[1] = Small_Theta(Angle(vectG[16], vectG[10])); //rel_2
			uSmallB[1] = Small_Theta(Angle(vectB[16], vectB[10]));

			aux2R = uBigR[0] * uSmallR[0] * uSmallR[1] * uBigR[1] * uBigR[2];
			aux2G = uBigG[0] * uSmallG[0] * uSmallG[1] * uBigG[1] * uBigG[2];
			aux2B = uBigB[0] * uSmallB[0] * uSmallB[1] * uBigB[1] * uBigB[2];

			rs_R[7] = std::min(aux1R, aux2R);
			rs_G[7] = std::min(aux1G, aux2G);
			rs_B[7] = std::min(aux1B, aux2B);

			for ( i = 0; i<8; i++)
			{
				if (rs_R[i] > r_R)
					r_R = rs_R[i];

				if (rs_G[i] > r_G)
					r_G = rs_G[i];

				if (rs_B[i] > r_B)
					r_B = rs_B[i];
			}

			if (r_R>=.3){
				Ruido_R = 255;
			}

			if (r_G >= .3) {
				Ruido_G = 255;
			}

			if (r_B >= .3) {
				Ruido_B = 255;
			}

			
			/////////////Filtrado/////////////

			
			if (Ruido_R ==255){
				int F;
				float weights[9], sum_weights = 0, hold2, suma = 0;
				
				for (j = 0; j <= 7; j++)
				{
					if (j == 4) {
						sum_weights += 0;//central pixel
					}
					else {
						sum_weights += (1 - rs_R[j]);
					}
				}
				sum_weights = (sum_weights + 3 * sqrt(1 - r_R)) / 2;
				weights[0] = (1 - rs_R[0]);
				weights[1] = (1 - rs_R[1]);
				weights[2] = (1 - rs_R[2]);
				weights[3] = (1 - rs_R[7]);
				weights[4] = 3 * sqrt(1 - r_R);
				weights[5] = (1 - rs_R[3]);
				weights[6] = (1 - rs_R[6]);
				weights[7] = (1 - rs_R[5]);
				weights[8] = (1 - rs_R[4]);


				vectR_Fil[0] = vectR[6];
				vectR_Fil[1] = vectR[7];
				vectR_Fil[2] = vectR[8];
				vectR_Fil[3] = vectR[11];
				vectR_Fil[4] = vectR[12];
				vectR_Fil[5] = vectR[13];
				vectR_Fil[6] = vectR[16];
				vectR_Fil[7] = vectR[17];
				vectR_Fil[8] = vectR[18];
					


				for (j = 0; j <= 8; j++)
				{
					for (x = 0; x <= 7; x++)
					{
						if (vectR_Fil[x] > vectR_Fil[x + 1])
						{
							hold = vectR_Fil[x];
							hold2 = weights[x];
							vectR_Fil[x] = vectR_Fil[x + 1];
							weights[x] = weights[x + 1];
							vectR_Fil[x + 1] = hold;
							weights[x + 1] = hold2;
						}
					}
				}

				for (j = 8; j >= 0; j--)
				{
					suma += weights[j];
					if (suma >= sum_weights)
					{
						if (j < 2)
						{
							sum_weights = sum_weights - (weights[0] + weights[1]);
							sum_weights = sum_weights / 2;
							suma = 0;
							for (F = 8; F >= 2; F--)
							{
								suma += weights[F];
								if (suma >= sum_weights)
								{
									d_Pout[(Row * m + Col) * channels + 0] = vectR_Fil[F];
									F = -1;
								}
							}
							j = -1;
						}
						else
						{
							d_Pout[(Row * m + Col) * channels + 0] = vectR_Fil[j];
							j = -1;
						}
						suma = -1;
					}
				}
			}
			else
			{
				d_Pout[(Row * m + Col) * channels + 0] = vectR[12];

			}


			
			if (Ruido_G == 255) {
				int F;
				float weights[9], sum_weights = 0, hold2, suma = 0;

				for (j = 0; j <= 7; j++)
				{
					if (j == 4) {
						sum_weights += 0;//central pixel
					}
					else {
						sum_weights += (1 - rs_G[j]);
					}
				}
				sum_weights = (sum_weights + 3 * sqrt(1 - r_G)) / 2;
				weights[0] = (1 - rs_G[0]);
				weights[1] = (1 - rs_G[1]);
				weights[2] = (1 - rs_G[2]);
				weights[3] = (1 - rs_G[7]);
				weights[4] = 3 * sqrt(1 - r_G);
				weights[5] = (1 - rs_G[3]);
				weights[6] = (1 - rs_G[6]);
				weights[7] = (1 - rs_G[5]);
				weights[8] = (1 - rs_G[4]);

				vectG_Fil[0] = vectG[6];
				vectG_Fil[1] = vectG[7];
				vectG_Fil[2] = vectG[8];
				vectG_Fil[3] = vectG[11];
				vectG_Fil[4] = vectG[12];
				vectG_Fil[5] = vectG[13];
				vectG_Fil[6] = vectG[16];
				vectG_Fil[7] = vectG[17];
				vectG_Fil[8] = vectG[18];

				for (j = 0; j <= 8; j++)
				{
					for (x = 0; x <= 7; x++)
					{
						if (vectG_Fil[x] > vectG_Fil[x + 1])
						{
							hold = vectG_Fil[x];
							hold2 = weights[x];
							vectG_Fil[x] = vectG_Fil[x + 1];
							weights[x] = weights[x + 1];
							vectG_Fil[x + 1] = hold;
							weights[x + 1] = hold2;
						}
					}
				}

				for (j = 8; j >= 0; j--)
				{
					suma += weights[j];
					if (suma >= sum_weights)
					{
						if (j < 2)
						{
							sum_weights = sum_weights - (weights[0] + weights[1]);
							sum_weights = sum_weights / 2;
							suma = 0;
							for (F = 8; F >= 2; F--)
							{
								suma += weights[F];
								if (suma >= sum_weights)
								{
									d_Pout[(Row * m + Col) * channels + 1] = vectG_Fil[F];
									F = -1;
								}
							}
							j = -1;
						}
						else
						{
							d_Pout[(Row * m + Col) * channels + 1] = vectG_Fil[j];
							j = -1;
						}
						suma = -1;
					}
				}
			}
			else
			{
				d_Pout[(Row * m + Col) * channels + 1] = vectG_Fil[12];

			}

			
			
			if (Ruido_B == 255) {
				int F;
				float weights[9], sum_weights = 0, hold2, suma = 0;

				for (j = 0; j <= 7; j++)
				{
					if (j == 4) {
						sum_weights += 0;//central pixel
					}
					else {
						sum_weights += (1 - rs_B[j]);
					}
				}
				sum_weights = (sum_weights + 3 * sqrt(1 - r_B)) / 2;
				weights[0] = (1 - rs_B[0]);
				weights[1] = (1 - rs_B[1]);
				weights[2] = (1 - rs_B[2]);
				weights[3] = (1 - rs_B[7]);
				weights[4] = 3 * sqrt(1 - r_B);
				weights[5] = (1 - rs_B[3]);
				weights[6] = (1 - rs_B[6]);
				weights[7] = (1 - rs_B[5]);
				weights[8] = (1 - rs_B[4]);

				vectB_Fil[0] = vectB[6];
				vectB_Fil[1] = vectB[7];
				vectB_Fil[2] = vectB[8];
				vectB_Fil[3] = vectB[11];
				vectB_Fil[4] = vectB[12];
				vectB_Fil[5] = vectB[13];
				vectB_Fil[6] = vectB[16];
				vectB_Fil[7] = vectB[17];
				vectB_Fil[8] = vectB[18];

				for (j = 0; j <= 8; j++)
				{
					for (x = 0; x <= 7; x++)
					{
						if (vectB_Fil[x] > vectB_Fil[x + 1])
						{
							hold = vectB_Fil[x];
							hold2 = weights[x];
							vectB_Fil[x] = vectB_Fil[x + 1];
							weights[x] = weights[x + 1];
							vectB_Fil[x + 1] = hold;
							weights[x + 1] = hold2;
						}
					}
				}

				for (j = 8; j >= 0; j--)
				{
					suma += weights[j];
					if (suma >= sum_weights)
					{
						if (j < 2)
						{
							sum_weights = sum_weights - (weights[0] + weights[1]);
							sum_weights = sum_weights / 2;
							suma = 0;
							for (F = 8; F >= 2; F--)
							{
								suma += weights[F];
								if (suma >= sum_weights)
								{
									d_Pout[(Row * m + Col) * channels + 2] = vectB_Fil[F];
									F = -1;
								}
							}
							j = -1;
						}
						else
						{
							d_Pout[(Row * m + Col) * channels + 2] = vectB_Fil[j];
							j = -1;
						}
						suma = -1;
					}
				}
			}
			else
			{
				d_Pout[(Row * m + Col) * channels + 2] = vectB[12];

			}
			


		}
	}
}

/*
for (j = 0; j <= 8; j++) {
for (i = 6; i <= 17; i++) {
if (vectR[i] >  vectR[i + 1]) {
hold = vectR[i];
hold2 = posicion[i];
vectR[i] = vectR[i + 1];
posicion[i] = posicion[i + 1];
disteucl1[i + 1] = hold;
posicion[i + 1] = hold2;
}
}
}
*/





#define min(a, b) ((a < b) ? a : b)
#define max(a, b) ((a > b) ? a : b) 


void FTSCF_CPU_Multi_Original(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads)
{

	int M = 0, j=0,x=0;
	int vectR[9], vectG[9], vectB[9], hold;
	int vect_R[25], vect_G[25], vect_B[25];
	float gam_small_1[18] = { 0 }, med_1, med_2, var_1, gam_big_1[18] = { 0 };
	float gam_small_2[18] = { 0 }, med1, med2, var1, gam_big_2[18] = { 0 };

	int array_R[25];
	int array_G[25];
	int array_B[25];

	int Row = 0, Col = 0;
	int F=0, i=0; 


	for (Row = 3; Row < n - 3; Row++) {
#pragma omp parallel for  num_threads(nThreads) private(Col ,i, j, F ,vectR,vectG,vectB,array_R,array_G,array_B,hold, x, M,gam_small_1,med_1, med_2, var_1, gam_big_1, gam_small_2, med1, med2, var1, gam_big_2) shared(d_Pout, d_Pin,n, m,channels,Row ) schedule(static)

		for (Col = 3; Col < m - 3; Col++) {
			//int tid = omp_get_thread_num();
			//hacer el arreglo
			F = 0;
			for (i = -2; i <= 2; i++) {
				for (j = -2; j <= 2; j++) {
					array_R[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 0];
					array_G[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 1];
					array_B[F] = d_Pin[((Row + i) * m + (Col + j)) * 3 + 2];
					F++;
				}
			}

			
			// se copia a continuacion solo los 8-vecinos
			M = 0;
			for (F = 6; F <= 8; F++){
				vectG[M] = (array_G[F]);
				vectR[M] = (array_R[F]);
				vectB[M] = (array_B[F]);
				M++;
			}
			for (F = 11; F <= 13; F++){
				vectG[M] = (array_G[F]);
				vectR[M] = (array_R[F]);
				vectB[M] = (array_B[F]);
				M++;
				}
			for (F = 16; F <= 18; F++){
				vectG[M] = (array_G[F]);
				vectR[M] = (array_R[F]);
				vectB[M] = (array_B[F]);
				M++;
			}
	

			float noreste_C_R, noreste_N1_R, noreste_N2_R, sur_C_R, sur_N1_R, sur_N2_R, noroeste_C_R, noroeste_N1_R, noroeste_N2_R;
			float este_C_R, este_N1_R, este_N2_R, oeste_C_R, oeste_N1_R, oeste_N2_R, sureste_C_R, sureste_N1_R, sureste_N2_R;
			float norte_C_R, norte_N1_R, norte_N2_R, suroeste_C_R, suroeste_N1_R, suroeste_N2_R;
			float suroeste_NW_R, suroeste_SE_R, sur_W_R, sur_E_R, sureste_SW_R, sureste_NE_R, este_S_R, este_N_R, noreste_SE_R, noreste_NW_R;
			float norte_W_R, norte_E_R, noroeste_NE_R, noroeste_SW_R, oeste_S_R, oeste_N_R;
			float noreste_C_G, noreste_N1_G, noreste_N2_G, sur_C_G, sur_N1_G, sur_N2_G, noroeste_C_G, noroeste_N1_G, noroeste_N2_G;
			float este_C_G, este_N1_G, este_N2_G, oeste_C_G, oeste_N1_G, oeste_N2_G, sureste_C_G, sureste_N1_G, sureste_N2_G;
			float norte_C_G, norte_N1_G, norte_N2_G, suroeste_C_G, suroeste_N1_G, suroeste_N2_G;
			float suroeste_NW_G, suroeste_SE_G, sur_W_G, sur_E_G, sureste_SW_G, sureste_NE_G, este_S_G, este_N_G, noreste_SE_G, noreste_NW_G;
			float norte_W_G, norte_E_G, noroeste_NE_G, noroeste_SW_G, oeste_S_G, oeste_N_G;
			float noreste_C_B, noreste_N1_B, noreste_N2_B, sur_C_B, sur_N1_B, sur_N2_B, noroeste_C_B, noroeste_N1_B, noroeste_N2_B;
			float este_C_B, este_N1_B, este_N2_B, oeste_C_B, oeste_N1_B, oeste_N2_B, sureste_C_B, sureste_N1_B, sureste_N2_B;
			float norte_C_B, norte_N1_B, norte_N2_B, suroeste_C_B, suroeste_N1_B, suroeste_N2_B;
			float suroeste_NW_B, suroeste_SE_B, sur_W_B, sur_E_B, sureste_SW_B, sureste_NE_B, este_S_B, este_N_B, noreste_SE_B, noreste_NW_B;
			float norte_W_B, norte_E_B, noroeste_NE_B, noroeste_SW_B, oeste_S_B, oeste_N_B;
			float largo[9], largo_1[9], largo_2[9], LARGO[9], LARGO_1[9], LARGO_2[9];
			float noise_R_R, noise_G_G, noise_B_B;
			int SW_C_B, SW_N1_B, SW_N2_B, SW_NW_B, SW_SE_B, S_C_B, S_N1_B, S_N2_B, S_W_B, S_E_B, SE_C_B, SE_N1_B, SE_N2_B, SE_SW_B, SE_NE_B;
			int E_C_B, E_N1_B, E_N2_B, E_S_B, E_N_B, NE_C_B, NE_N1_B, NE_N2_B, NE_SE_B, NE_NW_B, N_C_B, N_N1_B, N_N2_B, N_W_B, N_E_B;
			int NW_C_B, NW_N1_B, NW_N2_B, NW_NE_B, NW_SW_B, W_C_B, W_N1_B, W_N2_B, W_S_B, W_N_B;
			int SW_C_R, SW_N1_R, SW_N2_R, SW_NW_R, SW_SE_R, S_C_R, S_N1_R, S_N2_R, S_W_R, S_E_R, SE_C_R, SE_N1_R, SE_N2_R, SE_SW_R, SE_NE_R;
			int E_C_R, E_N1_R, E_N2_R, E_S_R, E_N_R, NE_C_R, NE_N1_R, NE_N2_R, NE_SE_R, NE_NW_R, N_C_R, N_N1_R, N_N2_R, N_W_R, N_E_R;
			int NW_C_R, NW_N1_R, NW_N2_R, NW_NE_R, NW_SW_R, W_C_R, W_N1_R, W_N2_R, W_S_R, W_N_R;
			int SW_C_G, SW_N1_G, SW_N2_G, SW_NW_G, SW_SE_G, S_C_G, S_N1_G, S_N2_G, S_W_G, S_E_G, SE_C_G, SE_N1_G, SE_N2_G, SE_SW_G, SE_NE_G;
			int E_C_G, E_N1_G, E_N2_G, E_S_G, E_N_G, NE_C_G, NE_N1_G, NE_N2_G, NE_SE_G, NE_NW_G, N_C_G, N_N1_G, N_N2_G, N_W_G, N_E_G;
			int NW_C_G, NW_N1_G, NW_N2_G, NW_NE_G, NW_SW_G, W_C_G, W_N1_G, W_N2_G, W_S_G, W_N_G;
			int AAA, BBB, CCC, cons1 = 255, cons2 = 255;

			/**************************************	AZUL	******************************************/
			SW_C_B = abs(array_B[6] - array_B[12]);
			SW_N1_B = abs(array_B[10] - array_B[16]);
			SW_N2_B = abs(array_B[2] - array_B[8]);
			SW_NW_B = abs(array_B[12] - array_B[16]);
			SW_SE_B = abs(array_B[12] - array_B[8]);
			S_C_B = abs(array_B[7] - array_B[12]);
			S_N1_B = abs(array_B[6] - array_B[11]);
			S_N2_B = abs(array_B[8] - array_B[13]);
			S_W_B = abs(array_B[12] - array_B[11]);
			S_E_B = abs(array_B[12] - array_B[13]);
			SE_C_B = abs(array_B[8] - array_B[12]);
			SE_N1_B = abs(array_B[2] - array_B[6]);
			SE_N2_B = abs(array_B[14] - array_B[18]);
			SE_SW_B = abs(array_B[12] - array_B[6]);
			SE_NE_B = abs(array_B[12] - array_B[18]);
			E_C_B = abs(array_B[13] - array_B[12]);
			E_N1_B = abs(array_B[8] - array_B[7]);
			E_N2_B = abs(array_B[18] - array_B[17]);
			E_S_B = abs(array_B[12] - array_B[7]);
			E_N_B = abs(array_B[12] - array_B[17]);
			NE_C_B = abs(array_B[18] - array_B[12]);
			NE_N1_B = abs(array_B[14] - array_B[8]);
			NE_N2_B = abs(array_B[22] - array_B[16]);
			NE_SE_B = abs(array_B[12] - array_B[8]);
			NE_NW_B = abs(array_B[12] - array_B[16]);
			N_C_B = abs(array_B[17] - array_B[12]);
			N_N1_B = abs(array_B[18] - array_B[13]);
			N_N2_B = abs(array_B[16] - array_B[11]);
			N_W_B = abs(array_B[12] - array_B[11]);
			N_E_B = abs(array_B[12] - array_B[13]);
			NW_C_B = abs(array_B[16] - array_B[12]);
			NW_N1_B = abs(array_B[22] - array_B[18]);
			NW_N2_B = abs(array_B[10] - array_B[6]);
			NW_NE_B = abs(array_B[12] - array_B[18]);
			NW_SW_B = abs(array_B[12] - array_B[6]);
			W_C_B = abs(array_B[11] - array_B[12]);
			W_N1_B = abs(array_B[16] - array_B[17]);
			W_N2_B = abs(array_B[6] - array_B[7]);
			W_S_B = abs(array_B[12] - array_B[7]);
			W_N_B = abs(array_B[12] - array_B[17]);

			SW_C_G = abs(array_G[6] - array_G[12]);
			SW_N1_G = abs(array_G[10] - array_G[16]);
			SW_N2_G = abs(array_G[2] - array_G[8]);
			SW_NW_G = abs(array_G[12] - array_G[16]);
			SW_SE_G = abs(array_G[12] - array_G[8]);
			S_C_G = abs(array_G[7] - array_G[12]);
			S_N1_G = abs(array_G[6] - array_G[11]);
			S_N2_G = abs(array_G[8] - array_G[13]);
			S_W_G = abs(array_G[12] - array_G[11]);
			S_E_G = abs(array_G[12] - array_G[13]);
			SE_C_G = abs(array_G[8] - array_G[12]);
			SE_N1_G = abs(array_G[2] - array_G[6]);
			SE_N2_G = abs(array_G[14] - array_G[18]);
			SE_SW_G = abs(array_G[12] - array_G[6]);
			SE_NE_G = abs(array_G[12] - array_G[18]);
			E_C_G = abs(array_G[13] - array_G[12]);
			E_N1_G = abs(array_G[8] - array_G[7]);
			E_N2_G = abs(array_G[18] - array_G[17]);
			E_S_G = abs(array_G[12] - array_G[7]);
			E_N_G = abs(array_G[12] - array_G[17]);
			NE_C_G = abs(array_G[18] - array_G[12]);
			NE_N1_G = abs(array_G[14] - array_G[8]);
			NE_N2_G = abs(array_G[22] - array_G[16]);
			NE_SE_G = abs(array_G[12] - array_G[8]);
			NE_NW_G = abs(array_G[12] - array_G[16]);
			N_C_G = abs(array_G[17] - array_G[12]);
			N_N1_G = abs(array_G[18] - array_G[13]);
			N_N2_G = abs(array_G[16] - array_G[11]);
			N_W_G = abs(array_G[12] - array_G[11]);
			N_E_G = abs(array_G[12] - array_G[13]);
			NW_C_G = abs(array_G[16] - array_G[12]);
			NW_N1_G = abs(array_G[22] - array_G[18]);
			NW_N2_G = abs(array_G[10] - array_G[6]);
			NW_NE_G = abs(array_G[12] - array_G[18]);
			NW_SW_G = abs(array_G[12] - array_G[6]);
			W_C_G = abs(array_G[11] - array_G[12]);
			W_N1_G = abs(array_G[16] - array_G[17]);
			W_N2_G = abs(array_G[6] - array_G[7]);
			W_S_G = abs(array_G[12] - array_G[7]);
			W_N_G = abs(array_G[12] - array_G[17]);

			SW_C_R = abs(array_R[6] - array_R[12]);
			SW_N1_R = abs(array_R[10] - array_R[16]);
			SW_N2_R = abs(array_R[2] - array_R[8]);
			SW_NW_R = abs(array_R[12] - array_R[16]);
			SW_SE_R = abs(array_R[12] - array_R[8]);
			S_C_R = abs(array_R[7] - array_R[12]);
			S_N1_R = abs(array_R[6] - array_R[11]);
			S_N2_R = abs(array_R[8] - array_R[13]);
			S_W_R = abs(array_R[12] - array_R[11]);
			S_E_R = abs(array_R[12] - array_R[13]);
			SE_C_R = abs(array_R[8] - array_R[12]);
			SE_N1_R = abs(array_R[2] - array_R[6]);
			SE_N2_R = abs(array_R[14] - array_R[18]);
			SE_SW_R = abs(array_R[12] - array_R[6]);
			SE_NE_R = abs(array_R[12] - array_R[18]);
			E_C_R = abs(array_R[13] - array_R[12]);
			E_N1_R = abs(array_R[8] - array_R[7]);
			E_N2_R = abs(array_R[18] - array_R[17]);
			E_S_R = abs(array_R[12] - array_R[7]);
			E_N_R = abs(array_R[12] - array_R[17]);
			NE_C_R = abs(array_R[18] - array_R[12]);
			NE_N1_R = abs(array_R[14] - array_R[8]);
			NE_N2_R = abs(array_R[22] - array_R[16]);
			NE_SE_R = abs(array_R[12] - array_R[8]);
			NE_NW_R = abs(array_R[12] - array_R[16]);
			N_C_R = abs(array_R[17] - array_R[12]);
			N_N1_R = abs(array_R[18] - array_R[13]);
			N_N2_R = abs(array_R[16] - array_R[11]);
			N_W_R = abs(array_R[12] - array_R[11]);
			N_E_R = abs(array_R[12] - array_R[13]);
			NW_C_R = abs(array_R[16] - array_R[12]);
			NW_N1_R = abs(array_R[22] - array_R[18]);
			NW_N2_R = abs(array_R[10] - array_R[6]);
			NW_NE_R = abs(array_R[12] - array_R[18]);
			NW_SW_R = abs(array_R[12] - array_R[6]);
			W_C_R = abs(array_R[11] - array_R[12]);
			W_N1_R = abs(array_R[16] - array_R[17]);
			W_N2_R = abs(array_R[6] - array_R[7]);
			W_S_R = abs(array_R[12] - array_R[7]);
			W_N_R = abs(array_R[12] - array_R[17]);

			if (((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[12])) == 0) suroeste_C_R = 0;
			else	suroeste_C_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[10] * array_R[16])) == 0) suroeste_N1_R = 0;
			else   suroeste_N1_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[10] * array_R[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[10], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[2] * array_R[8])) == 0) suroeste_N2_R = 0;
			else   suroeste_N2_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[2] * array_R[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[2], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[8], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[16])) == 0) suroeste_NW_R = 0;
			else	suroeste_NW_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[16], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[8])) == 0) suroeste_SE_R = 0;
			else	suroeste_SE_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[8], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[7] * array_R[12])) == 0) sur_C_R = 0;
			else	sur_C_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[7] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[7], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[11])) == 0) sur_N1_R = 0;
			else	sur_N1_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[11], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[8] * array_R[13])) == 0) sur_N2_R = 0;
			else   sur_N2_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[8] * array_R[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[13], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[11])) == 0) sur_W_R = 0;
			else	sur_W_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[11], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[13])) == 0) sur_E_R = 0;
			else	sur_E_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[13], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[8] * array_R[12])) == 0) sureste_C_R = 0;
			else	sureste_C_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[8] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[2])) == 0) sureste_N1_R = 0;
			else	sureste_N1_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[2])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[2], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[14] * array_R[18])) == 0) sureste_N2_R = 0;
			else	sureste_N2_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[14] * array_R[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[14], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[6])) == 0) sureste_SW_R = 0;
			else	sureste_SW_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[6])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[18])) == 0) sureste_NE_R = 0;
			else	sureste_NE_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[13] * array_R[12])) == 0) este_C_R = 0;
			else	este_C_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[13] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[8] * array_R[7])) == 0) este_N1_R = 0;
			else	este_N1_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[8] * array_R[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[8], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[7], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[18] * array_R[17])) == 0) este_N2_R = 0;
			else	este_N2_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[18] * array_R[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[18], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[7])) == 0) este_S_R = 0;
			else	este_S_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[7], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[17])) == 0) este_N_R = 0;
			else	este_N_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[18] * array_R[12])) == 0) noreste_C_R = 0;
			else	noreste_C_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[18] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[14] * array_R[8])) == 0) noreste_N1_R = 0;
			else	noreste_N1_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[14] * array_R[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[14], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[22] * array_R[16])) == 0) noreste_N2_R = 0;
			else	noreste_N2_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[22] * array_R[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[22], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[8])) == 0) noreste_SE_R = 0;
			else	noreste_SE_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[16])) == 0) noreste_NW_R = 0;
			else	noreste_NW_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[17] * array_R[12])) == 0) norte_C_R = 0;
			else	norte_C_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[17] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[18] * array_R[13])) == 0) norte_N1_R = 0;
			else	norte_N1_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[18] * array_R[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[18], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[16] * array_R[11])) == 0) norte_N2_R = 0;
			else	norte_N2_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[16] * array_R[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[16], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[13])) == 0) norte_E_R = 0;
			else	norte_E_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[11])) == 0) norte_W_R = 0;
			else	norte_W_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[16] * array_R[12])) == 0) noroeste_C_R = 0;
			else	noroeste_C_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[16] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[22] * array_R[18])) == 0) noroeste_N1_R = 0;
			else	noroeste_N1_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[22] * array_R[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[22], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[10])) == 0) noroeste_N2_R = 0;
			else	noroeste_N2_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[10])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[10], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[18])) == 0) noroeste_NE_R = 0;
			else	noroeste_NE_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[12])) == 0) noroeste_SW_R = 0;
			else	noroeste_SW_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[11] * array_R[12])) == 0) oeste_C_R = 0;
			else	oeste_C_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[11] * array_R[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[16] * array_R[17])) == 0) oeste_N1_R = 0;
			else	oeste_N1_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[16] * array_R[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[16], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[7])) == 0) oeste_N2_R = 0;
			else	oeste_N2_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[6] * array_R[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[7], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[17])) == 0) oeste_N_R = 0;
			else	oeste_N_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[17], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[7])) == 0) oeste_S_R = 0;
			else	oeste_S_R = acos(((cons1 + cons1) + (cons2*cons2) + (array_R[12] * array_R[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[7], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_R[12], 2))))));

			if (((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[12])) == 0) suroeste_C_G = 0;
			else	suroeste_C_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[10] * array_G[16])) == 0) suroeste_N1_G = 0;
			else   suroeste_N1_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[10] * array_G[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[10], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[2] * array_G[8])) == 0) suroeste_N2_G = 0;
			else   suroeste_N2_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[2] * array_G[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[2], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[8], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[16])) == 0) suroeste_NW_G = 0;
			else	suroeste_NW_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[16], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[8])) == 0) suroeste_SE_G = 0;
			else	suroeste_SE_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[8], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[7] * array_G[12])) == 0) sur_C_G = 0;
			else	sur_C_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[7] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[7], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[11])) == 0) sur_N1_G = 0;
			else	sur_N1_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[11], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[8] * array_G[13])) == 0) sur_N2_G = 0;
			else   sur_N2_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[8] * array_G[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[13], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[11])) == 0) sur_W_G = 0;
			else	sur_W_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[11], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[13])) == 0) sur_E_G = 0;
			else	sur_E_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[13], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[8] * array_G[12])) == 0) sureste_C_G = 0;
			else	sureste_C_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[8] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[2])) == 0) sureste_N1_G = 0;
			else	sureste_N1_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[2])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[2], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[14] * array_G[18])) == 0) sureste_N2_G = 0;
			else	sureste_N2_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[14] * array_G[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[14], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[6])) == 0) sureste_SW_G = 0;
			else	sureste_SW_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[6])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[18])) == 0) sureste_NE_G = 0;
			else	sureste_NE_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[13] * array_G[12])) == 0) este_C_G = 0;
			else	este_C_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[13] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[8] * array_G[7])) == 0) este_N1_G = 0;
			else	este_N1_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[8] * array_G[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[8], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[7], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[18] * array_G[17])) == 0) este_N2_G = 0;
			else	este_N2_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[18] * array_G[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[18], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[7])) == 0) este_S_G = 0;
			else	este_S_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[7], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[17])) == 0) este_N_G = 0;
			else	este_N_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[18] * array_G[12])) == 0) noreste_C_G = 0;
			else	noreste_C_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[18] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[14] * array_G[8])) == 0) noreste_N1_G = 0;
			else	noreste_N1_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[14] * array_G[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[14], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[22] * array_G[16])) == 0) noreste_N2_G = 0;
			else	noreste_N2_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[22] * array_G[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[22], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[8])) == 0) noreste_SE_G = 0;
			else	noreste_SE_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[16])) == 0) noreste_NW_G = 0;
			else	noreste_NW_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[17] * array_G[12])) == 0) norte_C_G = 0;
			else	norte_C_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[17] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[18] * array_G[13])) == 0) norte_N1_G = 0;
			else	norte_N1_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[18] * array_G[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[18], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[16] * array_G[11])) == 0) norte_N2_G = 0;
			else	norte_N2_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[16] * array_G[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[16], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[13])) == 0) norte_E_G = 0;
			else	norte_E_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[11])) == 0) norte_W_G = 0;
			else	norte_W_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[16] * array_G[12])) == 0) noroeste_C_G = 0;
			else	noroeste_C_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[16] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[22] * array_G[18])) == 0) noroeste_N1_G = 0;
			else	noroeste_N1_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[22] * array_G[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[22], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[10])) == 0) noroeste_N2_G = 0;
			else	noroeste_N2_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[10])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[10], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[18])) == 0) noroeste_NE_G = 0;
			else	noroeste_NE_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[12])) == 0) noroeste_SW_G = 0;
			else	noroeste_SW_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[11] * array_G[12])) == 0) oeste_C_G = 0;
			else	oeste_C_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[11] * array_G[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[16] * array_G[17])) == 0) oeste_N1_G = 0;
			else	oeste_N1_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[16] * array_G[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[16], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[7])) == 0) oeste_N2_G = 0;
			else	oeste_N2_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[6] * array_G[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[7], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[17])) == 0) oeste_N_G = 0;
			else	oeste_N_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[17], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[7])) == 0) oeste_S_G = 0;
			else	oeste_S_G = acos(((cons1 + cons1) + (cons2*cons2) + (array_G[12] * array_G[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[7], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_G[12], 2))))));

			if (((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[12])) == 0) suroeste_C_B = 0;
			else	suroeste_C_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[10] * array_B[16])) == 0) suroeste_N1_B = 0;
			else   suroeste_N1_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[10] * array_B[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[10], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[2] * array_B[8])) == 0) suroeste_N2_B = 0;
			else   suroeste_N2_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[2] * array_B[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[2], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[8], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[16])) == 0) suroeste_NW_B = 0;
			else	suroeste_NW_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[16], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[8])) == 0) suroeste_SE_B = 0;
			else	suroeste_SE_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[8], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[7] * array_B[12])) == 0) sur_C_B = 0;
			else	sur_C_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[7] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2)))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[7], 2)))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[11])) == 0) sur_N1_B = 0;
			else	sur_N1_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[11], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[8] * array_B[13])) == 0) sur_N2_B = 0;
			else   sur_N2_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[8] * array_B[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[13], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[11])) == 0) sur_W_B = 0;
			else	sur_W_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[11], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[13])) == 0) sur_E_B = 0;
			else	sur_E_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[13], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[8] * array_B[12])) == 0) sureste_C_B = 0;
			else	sureste_C_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[8] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[2])) == 0) sureste_N1_B = 0;
			else	sureste_N1_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[2])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[2], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[14] * array_B[18])) == 0) sureste_N2_B = 0;
			else	sureste_N2_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[14] * array_B[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[14], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[6])) == 0) sureste_SW_B = 0;
			else	sureste_SW_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[6])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[18])) == 0) sureste_NE_B = 0;
			else	sureste_NE_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[13] * array_B[12])) == 0) este_C_B = 0;
			else	este_C_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[13] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[8] * array_B[7])) == 0) este_N1_B = 0;
			else	este_N1_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[8] * array_B[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[8], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[7], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[18] * array_B[17])) == 0) este_N2_B = 0;
			else	este_N2_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[18] * array_B[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[18], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[7])) == 0) este_S_B = 0;
			else	este_S_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[7], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[17])) == 0) este_N_B = 0;
			else	este_N_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[18] * array_B[12])) == 0) noreste_C_B = 0;
			else	noreste_C_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[18] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[14] * array_B[8])) == 0) noreste_N1_B = 0;
			else	noreste_N1_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[14] * array_B[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[14], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[22] * array_B[16])) == 0) noreste_N2_B = 0;
			else	noreste_N2_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[22] * array_B[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[22], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[8])) == 0) noreste_SE_B = 0;
			else	noreste_SE_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[8])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[8], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[16])) == 0) noreste_NW_B = 0;
			else	noreste_NW_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[16])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[17] * array_B[12])) == 0) norte_C_B = 0;
			else	norte_C_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[17] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[18] * array_B[13])) == 0) norte_N1_B = 0;
			else	norte_N1_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[18] * array_B[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[18], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[16] * array_B[11])) == 0) norte_N2_B = 0;
			else	norte_N2_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[16] * array_B[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[16], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[13])) == 0) norte_E_B = 0;
			else	norte_E_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[13])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[13], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[11])) == 0) norte_W_B = 0;
			else	norte_W_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[11])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[16] * array_B[12])) == 0) noroeste_C_B = 0;
			else	noroeste_C_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[16] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[16], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[22] * array_B[18])) == 0) noroeste_N1_B = 0;
			else	noroeste_N1_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[22] * array_B[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[22], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[10])) == 0) noroeste_N2_B = 0;
			else	noroeste_N2_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[10])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[10], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[18])) == 0) noroeste_NE_B = 0;
			else	noroeste_NE_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[18])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[18], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[12])) == 0) noroeste_SW_B = 0;
			else	noroeste_SW_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[11] * array_B[12])) == 0) oeste_C_B = 0;
			else	oeste_C_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[11] * array_B[12])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[11], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[16] * array_B[17])) == 0) oeste_N1_B = 0;
			else	oeste_N1_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[16] * array_B[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[16], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[17], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[7])) == 0) oeste_N2_B = 0;
			else	oeste_N2_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[6] * array_B[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[7], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[6], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[17])) == 0) oeste_N_B = 0;
			else	oeste_N_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[17])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[17], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))))));
			if (((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[7])) == 0) oeste_S_B = 0;
			else	oeste_S_B = acos(((cons1 + cons1) + (cons2*cons2) + (array_B[12] * array_B[7])) / ((sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[7], 2))*(sqrt(pow(cons1, 2) + pow(cons1, 2) + pow(array_B[12], 2))))));
			/******	SUROESTE	******/

			med_1 = 1/*1*/, var_1 = 0.8/*0.1*/;
			med_2 = 0.1/*0.8*/;

			if (suroeste_C_R > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((suroeste_C_R)-med_1), 2) / (2 * var_1))));
			if (suroeste_N1_R < med_2) gam_small_1[0] = 1;
			else 	gam_small_1[0] = (exp(-(pow(((suroeste_N1_R)-med_2), 2) / (2 * var_1))));
			if (suroeste_N2_R < med_2) gam_small_1[1] = 1;
			else 	gam_small_1[1] = (exp(-(pow(((suroeste_N2_R)-med_2), 2) / (2 * var_1))));
			if (suroeste_NW_R > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((suroeste_NW_R)-med_1), 2) / (2 * var_1))));
			if (suroeste_SE_R > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((suroeste_SE_R)-med_1), 2) / (2 * var_1))));
			largo[0] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (sur_C_R > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((sur_C_R)-med_1), 2) / (2 * var_1))));
			if (sur_N1_R < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((sur_N1_R)-med_2), 2) / (2 * var_1))));
			if (sur_N2_R < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((sur_N2_R)-med_2), 2) / (2 * var_1))));
			if (sur_W_R > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((sur_W_R)-med_1), 2) / (2 * var_1))));
			if (sur_E_R > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((sur_E_R)-med_1), 2) / (2 * var_1))));
			largo[1] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (sureste_C_R > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((sureste_C_R)-med_1), 2) / (2 * var_1))));
			if (sureste_N1_R < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((sureste_N1_R)-med_2), 2) / (2 * var_1))));
			if (sureste_N2_R < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((sureste_N2_R)-med_2), 2) / (2 * var_1))));
			if (sureste_NE_R > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((sureste_NE_R)-med_1), 2) / (2 * var_1))));
			if (sureste_SW_R > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((sureste_SW_R)-med_1), 2) / (2 * var_1))));
			largo[2] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (este_C_R > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((este_C_R)-med_1), 2) / (2 * var_1))));
			if (este_N1_R < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((este_N1_R)-med_2), 2) / (2 * var_1))));
			if (este_N2_R < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((este_N2_R)-med_2), 2) / (2 * var_1))));
			if (este_N_R > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((este_N_R)-med_1), 2) / (2 * var_1))));
			if (este_S_R > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((este_S_R)-med_1), 2) / (2 * var_1))));
			largo[3] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (noreste_C_R > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((noreste_C_R)-med_1), 2) / (2 * var_1))));
			if (noreste_N1_R < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((noreste_N1_R)-med_2), 2) / (2 * var_1))));
			if (noreste_N2_R < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((noreste_N2_R)-med_2), 2) / (2 * var_1))));
			if (noreste_NW_R > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((noreste_NW_R)-med_1), 2) / (2 * var_1))));
			if (noreste_SE_R > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((noreste_SE_R)-med_1), 2) / (2 * var_1))));
			largo[4] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (norte_C_R > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((norte_C_R)-med_1), 2) / (2 * var_1))));
			if (norte_N1_R < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((norte_N1_R)-med_2), 2) / (2 * var_1))));
			if (norte_N2_R < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((norte_N2_R)-med_2), 2) / (2 * var_1))));
			if (norte_W_R > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((norte_W_R)-med_1), 2) / (2 * var_1))));
			if (norte_E_R > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((norte_E_R)-med_1), 2) / (2 * var_1))));
			largo[5] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (noroeste_C_R > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((noroeste_C_R)-med_1), 2) / (2 * var_1))));
			if (noroeste_N1_R < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((noroeste_N1_R)-med_2), 2) / (2 * var_1))));
			if (noroeste_N2_R < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((noroeste_N2_R)-med_2), 2) / (2 * var_1))));
			if (noroeste_NE_R > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((noroeste_NE_R)-med_1), 2) / (2 * var_1))));
			if (noroeste_SW_R > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((noroeste_SW_R)-med_1), 2) / (2 * var_1))));
			largo[6] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (oeste_C_R > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((oeste_C_R)-med_1), 2) / (2 * var_1))));
			if (oeste_N1_R < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((oeste_N1_R)-med_2), 2) / (2 * var_1))));
			if (oeste_N2_R < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((oeste_N2_R)-med_2), 2) / (2 * var_1))));
			if (oeste_N_R > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((oeste_N_R)-med_1), 2) / (2 * var_1))));
			if (oeste_S_R > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((oeste_S_R)-med_1), 2) / (2 * var_1))));
			largo[7] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (suroeste_C_G > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((suroeste_C_G)-med_1), 2) / (2 * var_1))));
			if (suroeste_N1_G < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((suroeste_N1_G)-med_2), 2) / (2 * var_1))));
			if (suroeste_N2_G < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((suroeste_N2_G)-med_2), 2) / (2 * var_1))));
			if (suroeste_NW_G > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((suroeste_NW_G)-med_1), 2) / (2 * var_1))));
			if (suroeste_SE_G > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((suroeste_SE_G)-med_1), 2) / (2 * var_1))));
			largo_1[0] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (sur_C_G > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((sur_C_G)-med_1), 2) / (2 * var_1))));
			if (sur_N1_G < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((sur_N1_G)-med_2), 2) / (2 * var_1))));
			if (sur_N2_G < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((sur_N2_G)-med_2), 2) / (2 * var_1))));
			if (sur_W_G > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((sur_W_G)-med_1), 2) / (2 * var_1))));
			if (sur_E_G > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((sur_E_G)-med_1), 2) / (2 * var_1))));
			largo_1[1] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (sureste_C_G > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((sureste_C_G)-med_1), 2) / (2 * var_1))));
			if (sureste_N1_G < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((sureste_N1_G)-med_2), 2) / (2 * var_1))));
			if (sureste_N2_G < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((sureste_N2_G)-med_2), 2) / (2 * var_1))));
			if (sureste_NE_G > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((sureste_NE_G)-med_1), 2) / (2 * var_1))));
			if (sureste_SW_G > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((sureste_SW_G)-med_1), 2) / (2 * var_1))));
			largo_1[2] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (este_C_G > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((este_C_G)-med_1), 2) / (2 * var_1))));
			if (este_N1_G < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((este_N1_G)-med_2), 2) / (2 * var_1))));
			if (este_N2_G < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((este_N2_G)-med_2), 2) / (2 * var_1))));
			if (este_N_G > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((este_N_G)-med_1), 2) / (2 * var_1))));
			if (este_S_G > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((este_S_G)-med_1), 2) / (2 * var_1))));
			largo_1[3] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (noreste_C_G > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((noreste_C_G)-med_1), 2) / (2 * var_1))));
			if (noreste_N1_G < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((noreste_N1_G)-med_2), 2) / (2 * var_1))));
			if (noreste_N2_G < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((noreste_N2_G)-med_2), 2) / (2 * var_1))));
			if (noreste_NW_G > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((noreste_NW_G)-med_1), 2) / (2 * var_1))));
			if (noreste_SE_G > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((noreste_SE_G)-med_1), 2) / (2 * var_1))));
			largo_1[4] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (norte_C_G > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((norte_C_G)-med_1), 2) / (2 * var_1))));
			if (norte_N1_G < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((norte_N1_G)-med_2), 2) / (2 * var_1))));
			if (norte_N2_G < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((norte_N2_G)-med_2), 2) / (2 * var_1))));
			if (norte_W_G > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((norte_W_G)-med_1), 2) / (2 * var_1))));
			if (norte_E_G > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((norte_E_G)-med_1), 2) / (2 * var_1))));
			largo_1[5] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (noroeste_C_G > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((noroeste_C_G)-med_1), 2) / (2 * var_1))));
			if (noroeste_N1_G < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((noroeste_N1_G)-med_2), 2) / (2 * var_1))));
			if (noroeste_N2_G < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((noroeste_N2_G)-med_2), 2) / (2 * var_1))));
			if (noroeste_NE_G > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((noroeste_NE_G)-med_1), 2) / (2 * var_1))));
			if (noroeste_SW_G > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((noroeste_SW_G)-med_1), 2) / (2 * var_1))));
			largo_1[6] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (oeste_C_G > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((oeste_C_G)-med_1), 2) / (2 * var_1))));
			if (oeste_N1_G < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((oeste_N1_G)-med_2), 2) / (2 * var_1))));
			if (oeste_N2_G < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((oeste_N2_G)-med_2), 2) / (2 * var_1))));
			if (oeste_N_G > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((oeste_N_G)-med_1), 2) / (2 * var_1))));
			if (oeste_S_G > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((oeste_S_G)-med_1), 2) / (2 * var_1))));
			largo_1[7] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (suroeste_C_B > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((suroeste_C_B)-med_1), 2) / (2 * var_1))));
			if (suroeste_N1_B < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((suroeste_N1_B)-med_2), 2) / (2 * var_1))));
			if (suroeste_N2_B < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((suroeste_N2_B)-med_2), 2) / (2 * var_1))));
			if (suroeste_NW_B > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((suroeste_NW_B)-med_1), 2) / (2 * var_1))));
			if (suroeste_SE_B > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((suroeste_SE_B)-med_1), 2) / (2 * var_1))));
			largo_2[0] = (gam_big_2[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_2[2]);
			if (sur_C_B > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((sur_C_B)-med_1), 2) / (2 * var_1))));
			if (sur_N1_B < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((sur_N1_B)-med_2), 2) / (2 * var_1))));
			if (sur_N2_B < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((sur_N2_B)-med_2), 2) / (2 * var_1))));
			if (sur_W_B > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((sur_W_B)-med_1), 2) / (2 * var_1))));
			if (sur_E_B > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((sur_E_B)-med_1), 2) / (2 * var_1))));
			largo_2[1] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (sureste_C_B > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((sureste_C_B)-med_1), 2) / (2 * var_1))));
			if (sureste_N1_B < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((sureste_N1_B)-med_2), 2) / (2 * var_1))));
			if (sureste_N2_B < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((sureste_N2_B)-med_2), 2) / (2 * var_1))));
			if (sureste_NE_B > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((sureste_NE_B)-med_1), 2) / (2 * var_1))));
			if (sureste_SW_B > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((sureste_SW_B)-med_1), 2) / (2 * var_1))));
			largo_2[2] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (este_C_B > med_1) gam_big_1[0] = 1;

			else	gam_big_1[0] = (exp(-(pow(((este_C_B)-med_1), 2) / (2 * var_1))));
			if (este_N1_B < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((este_N1_B)-med_2), 2) / (2 * var_1))));
			if (este_N2_B < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((este_N2_B)-med_2), 2) / (2 * var_1))));
			if (este_N_B > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((este_N_B)-med_1), 2) / (2 * var_1))));
			if (este_S_B > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((este_S_B)-med_1), 2) / (2 * var_1))));
			largo_2[3] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (noreste_C_B > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((noreste_C_B)-med_1), 2) / (2 * var_1))));
			if (noreste_N1_B < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((noreste_N1_B)-med_2), 2) / (2 * var_1))));
			if (noreste_N2_B < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((noreste_N2_B)-med_2), 2) / (2 * var_1))));
			if (noreste_NW_B > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((noreste_NW_B)-med_1), 2) / (2 * var_1))));
			if (noreste_SE_B > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((noreste_SE_B)-med_1), 2) / (2 * var_1))));
			largo_2[4] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (norte_C_B > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((norte_C_B)-med_1), 2) / (2 * var_1))));
			if (norte_N1_B < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((norte_N1_B)-med_2), 2) / (2 * var_1))));
			if (norte_N2_B < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((norte_N2_B)-med_2), 2) / (2 * var_1))));
			if (norte_W_B > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((norte_W_B)-med_1), 2) / (2 * var_1))));
			if (norte_E_B > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((norte_E_B)-med_1), 2) / (2 * var_1))));
			largo_2[5] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (noroeste_C_B > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((noroeste_C_B)-med_1), 2) / (2 * var_1))));
			if (noroeste_N1_B < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((noroeste_N1_B)-med_2), 2) / (2 * var_1))));
			if (noroeste_N2_B < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((noroeste_N2_B)-med_2), 2) / (2 * var_1))));
			if (noroeste_NE_B > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((noroeste_NE_B)-med_1), 2) / (2 * var_1))));
			if (noroeste_SW_B > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((noroeste_SW_B)-med_1), 2) / (2 * var_1))));
			largo_2[6] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);
			if (oeste_C_B > med_1) gam_big_1[0] = 1;
			else	gam_big_1[0] = (exp(-(pow(((oeste_C_B)-med_1), 2) / (2 * var_1))));
			if (oeste_N1_B < med_2) gam_small_1[0] = 1;
			else	gam_small_1[0] = (exp(-(pow(((oeste_N1_B)-med_2), 2) / (2 * var_1))));
			if (oeste_N2_B < med_2) gam_small_1[1] = 1;
			else	gam_small_1[1] = (exp(-(pow(((oeste_N2_B)-med_2), 2) / (2 * var_1))));
			if (oeste_N_B > med_1) gam_big_1[1] = 1;
			else	gam_big_1[1] = (exp(-(pow(((oeste_N_B)-med_1), 2) / (2 * var_1))));
			if (oeste_S_B > med_1) gam_big_1[2] = 1;
			else	gam_big_1[2] = (exp(-(pow(((oeste_S_B)-med_1), 2) / (2 * var_1))));
			largo_2[7] = (gam_big_1[0] * gam_small_1[0] * gam_small_1[1] * gam_big_1[1] * gam_big_1[2]);

			med1 = 60;
			med2 = 10/*20*/;
			var1 = 1000;
			if (SW_C_R > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((SW_C_R)-med1), 2) / (2 * var1))));
			if (SW_N1_R < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((SW_N1_R)-med2), 2) / (2 * var1))));
			if (SW_N2_R < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((SW_N2_R)-med2), 2) / (2 * var1))));
			if (SW_NW_R > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((SW_NW_R)-med1), 2) / (2 * var1))));
			if (SW_SE_R > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((SW_SE_R)-med1), 2) / (2 * var1))));
			LARGO[0] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (S_C_R > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((S_C_R)-med1), 2) / (2 * var1))));
			if (S_N1_R < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((S_N1_R)-med2), 2) / (2 * var1))));
			if (S_N2_R < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((S_N2_R)-med2), 2) / (2 * var1))));
			if (S_W_R > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((S_W_R)-med1), 2) / (2 * var1))));
			if (S_E_R > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((S_E_R)-med1), 2) / (2 * var1))));
			LARGO[1] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (SE_C_R > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((SE_C_R)-med1), 2) / (2 * var1))));
			if (SE_N1_R < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((SE_N1_R)-med2), 2) / (2 * var1))));
			if (SE_N2_R < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((SE_N2_R)-med2), 2) / (2 * var1))));
			if (SE_NE_R > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((SE_NE_R)-med1), 2) / (2 * var1))));
			if (SE_SW_R > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((SE_SW_R)-med1), 2) / (2 * var1))));
			LARGO[2] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (E_C_R > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((E_C_R)-med1), 2) / (2 * var1))));
			if (E_N1_R < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((E_N1_R)-med2), 2) / (2 * var1))));
			if (E_N2_R < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((E_N2_R)-med2), 2) / (2 * var1))));
			if (E_N_R > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((E_N_R)-med1), 2) / (2 * var1))));
			if (E_S_R > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((E_S_R)-med1), 2) / (2 * var1))));
			LARGO[3] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (NE_C_R > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((NE_C_R)-med1), 2) / (2 * var1))));
			if (NE_N1_R < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((NE_N1_R)-med2), 2) / (2 * var1))));
			if (NE_N2_R < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((NE_N2_R)-med2), 2) / (2 * var1))));
			if (NE_NW_R > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((NE_NW_R)-med1), 2) / (2 * var1))));
			if (NE_SE_R > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((NE_SE_R)-med1), 2) / (2 * var1))));
			LARGO[4] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (N_C_R > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((N_C_R)-med1), 2) / (2 * var1))));
			if (N_N1_R < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((N_N1_R)-med2), 2) / (2 * var1))));
			if (N_N2_R < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((N_N2_R)-med2), 2) / (2 * var1))));
			if (N_W_R > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((N_W_R)-med1), 2) / (2 * var1))));
			if (N_E_R > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((N_E_R)-med1), 2) / (2 * var1))));
			LARGO[5] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (NW_C_R > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((NW_C_R)-med1), 2) / (2 * var1))));
			if (NW_N1_R < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((NW_N1_R)-med2), 2) / (2 * var1))));
			if (NW_N2_R < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((NW_N2_R)-med2), 2) / (2 * var1))));
			if (NW_NE_R > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((NW_NE_R)-med1), 2) / (2 * var1))));
			if (NW_SW_R > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((NW_SW_R)-med1), 2) / (2 * var1))));
			LARGO[6] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (W_C_R > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((W_C_R)-med1), 2) / (2 * var1))));
			if (W_N1_R < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((W_N1_R)-med2), 2) / (2 * var1))));
			if (W_N2_R < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((W_N2_R)-med2), 2) / (2 * var1))));
			if (W_N_R > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((W_N_R)-med1), 2) / (2 * var1))));
			if (W_S_R > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((W_S_R)-med1), 2) / (2 * var1))));
			LARGO[7] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (SW_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((SW_C_G)-med1), 2) / (2 * var1))));
			if (SW_N1_G < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((SW_N1_G)-med2), 2) / (2 * var1))));
			if (SW_N2_G < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((SW_N2_G)-med2), 2) / (2 * var1))));
			if (SW_NW_G > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((SW_NW_G)-med1), 2) / (2 * var1))));
			if (SW_SE_G > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((SW_SE_G)-med1), 2) / (2 * var1))));
			LARGO_1[0] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (S_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((S_C_G)-med1), 2) / (2 * var1))));
			if (S_N1_G < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((S_N1_G)-med2), 2) / (2 * var1))));
			if (S_N2_G < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((S_N2_G)-med2), 2) / (2 * var1))));
			if (S_W_G > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((S_W_G)-med1), 2) / (2 * var1))));
			if (S_E_G > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((S_E_G)-med1), 2) / (2 * var1))));
			LARGO_1[1] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (SE_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((SE_C_G)-med1), 2) / (2 * var1))));
			if (SE_N1_G < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((SE_N1_G)-med2), 2) / (2 * var1))));
			if (SE_N2_G < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((SE_N2_G)-med2), 2) / (2 * var1))));
			if (SE_NE_G > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((SE_NE_G)-med1), 2) / (2 * var1))));
			if (SE_SW_G > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((SE_SW_G)-med1), 2) / (2 * var1))));
			LARGO_1[2] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (E_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((E_C_G)-med1), 2) / (2 * var1))));
			if (E_N1_G < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((E_N1_G)-med2), 2) / (2 * var1))));
			if (E_N2_G < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((E_N2_G)-med2), 2) / (2 * var1))));
			if (E_N_G > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((E_N_G)-med1), 2) / (2 * var1))));
			if (E_S_G > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((E_S_G)-med1), 2) / (2 * var1))));
			LARGO_1[3] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (NE_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((NE_C_G)-med1), 2) / (2 * var1))));
			if (NE_N1_G < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((NE_N1_G)-med2), 2) / (2 * var1))));
			if (NE_N2_G < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((NE_N2_G)-med2), 2) / (2 * var1))));
			if (NE_NW_G > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((NE_NW_G)-med1), 2) / (2 * var1))));
			if (NE_SE_G > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((NE_SE_G)-med1), 2) / (2 * var1))));
			LARGO_1[4] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (N_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((N_C_G)-med1), 2) / (2 * var1))));
			if (N_N1_G < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((N_N1_G)-med2), 2) / (2 * var1))));
			if (N_N2_G < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((N_N2_G)-med2), 2) / (2 * var1))));
			if (N_W_G > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((N_W_G)-med1), 2) / (2 * var1))));
			if (N_E_G > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((N_E_G)-med1), 2) / (2 * var1))));
			LARGO_1[5] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (NW_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((NW_C_G)-med1), 2) / (2 * var1))));
			if (NW_N1_G < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((NW_N1_G)-med2), 2) / (2 * var1))));
			if (NW_N2_G < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((NW_N2_G)-med2), 2) / (2 * var1))));
			if (NW_NE_G > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((NW_NE_G)-med1), 2) / (2 * var1))));
			if (NW_SW_G > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((NW_SW_G)-med1), 2) / (2 * var1))));
			LARGO_1[6] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (W_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((W_C_G)-med1), 2) / (2 * var1))));
			if (W_N1_G < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((W_N1_G)-med2), 2) / (2 * var1))));
			if (W_N2_G < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((W_N2_G)-med2), 2) / (2 * var1))));
			if (W_N_G > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((W_N_G)-med1), 2) / (2 * var1))));
			if (W_S_G > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((W_S_G)-med1), 2) / (2 * var1))));
			LARGO_1[7] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (SW_C_G > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((SW_C_B)-med1), 2) / (2 * var1))));
			if (SW_N1_B < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((SW_N1_B)-med2), 2) / (2 * var1))));
			if (SW_N2_B < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((SW_N2_B)-med2), 2) / (2 * var1))));
			if (SW_NW_B > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((SW_NW_B)-med1), 2) / (2 * var1))));
			if (SW_SE_B > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((SW_SE_B)-med1), 2) / (2 * var1))));
			LARGO_2[0] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (S_C_B > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((S_C_B)-med1), 2) / (2 * var1))));
			if (S_N1_B < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((S_N1_B)-med2), 2) / (2 * var1))));
			if (S_N2_B < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((S_N2_B)-med2), 2) / (2 * var1))));
			if (S_W_B > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((S_W_B)-med1), 2) / (2 * var1))));
			if (S_E_B > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((S_E_B)-med1), 2) / (2 * var1))));
			LARGO_2[1] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (SE_C_B > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((SE_C_B)-med1), 2) / (2 * var1))));
			if (SE_N1_B < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((SE_N1_B)-med2), 2) / (2 * var1))));
			if (SE_N2_B < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((SE_N2_B)-med2), 2) / (2 * var1))));
			if (SE_NE_B > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((SE_NE_B)-med1), 2) / (2 * var1))));
			if (SE_SW_B > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((SE_SW_B)-med1), 2) / (2 * var1))));
			LARGO_2[2] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (E_C_B > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((E_C_B)-med1), 2) / (2 * var1))));
			if (E_N1_B < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((E_N1_B)-med2), 2) / (2 * var1))));
			if (E_N2_B < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((E_N2_B)-med2), 2) / (2 * var1))));
			if (E_N_B > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((E_N_B)-med1), 2) / (2 * var1))));
			if (E_S_B > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((E_S_B)-med1), 2) / (2 * var1))));
			LARGO_2[3] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (NE_C_B > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((NE_C_B)-med1), 2) / (2 * var1))));
			if (NE_N1_B < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((NE_N1_B)-med2), 2) / (2 * var1))));
			if (NE_N2_B < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((NE_N2_B)-med2), 2) / (2 * var1))));
			if (NE_NW_B > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((NE_NW_B)-med1), 2) / (2 * var1))));
			if (NE_SE_B > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((NE_SE_B)-med1), 2) / (2 * var1))));
			LARGO_2[4] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (N_C_B > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((N_C_B)-med1), 2) / (2 * var1))));
			if (N_N1_B < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((N_N1_B)-med2), 2) / (2 * var1))));
			if (N_N2_B < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((N_N2_B)-med2), 2) / (2 * var1))));
			if (N_W_B > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((N_W_B)-med1), 2) / (2 * var1))));
			if (N_E_B > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((N_E_B)-med1), 2) / (2 * var1))));
			LARGO_2[5] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (NW_C_B > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((NW_C_B)-med1), 2) / (2 * var1))));
			if (NW_N1_B < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((NW_N1_B)-med2), 2) / (2 * var1))));
			if (NW_N2_B < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((NW_N2_B)-med2), 2) / (2 * var1))));
			if (NW_NE_B > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((NW_NE_B)-med1), 2) / (2 * var1))));
			if (NW_SW_B > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((NW_SW_B)-med1), 2) / (2 * var1))));
			LARGO_2[6] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);
			if (W_C_B > med1) gam_big_2[0] = 1;
			else	gam_big_2[0] = (exp(-(pow(((W_C_B)-med1), 2) / (2 * var1))));
			if (W_N1_B < med2) gam_small_2[0] = 1;
			else	gam_small_2[0] = (exp(-(pow(((W_N1_B)-med2), 2) / (2 * var1))));
			if (W_N2_B < med2) gam_small_2[1] = 1;
			else	gam_small_2[1] = (exp(-(pow(((W_N2_B)-med2), 2) / (2 * var1))));
			if (W_N_B > med1) gam_big_2[1] = 1;
			else	gam_big_2[1] = (exp(-(pow(((W_N_B)-med1), 2) / (2 * var1))));
			if (W_S_B > med1) gam_big_2[2] = 1;
			else	gam_big_2[2] = (exp(-(pow(((W_S_B)-med1), 2) / (2 * var1))));
			LARGO_2[7] = (gam_big_2[0] * gam_small_2[0] * gam_small_2[1] * gam_big_2[1] * gam_big_2[2]);

			float	mu_R_R[8], mu_G_G[8], mu_B_B[8];

			mu_R_R[0] = min(largo[0], LARGO[0]);
			mu_R_R[1] = min(largo[1], LARGO[1]);
			mu_R_R[2] = min(largo[2], LARGO[2]);
			mu_R_R[3] = min(largo[3], LARGO[3]);
			mu_R_R[4] = min(largo[4], LARGO[4]);
			mu_R_R[5] = min(largo[5], LARGO[5]);
			mu_R_R[6] = min(largo[6], LARGO[6]);
			mu_R_R[7] = min(largo[7], LARGO[7]);

			mu_G_G[0] = min(largo_1[0], LARGO_1[0]);
			mu_G_G[1] = min(largo_1[1], LARGO_1[1]);
			mu_G_G[2] = min(largo_1[2], LARGO_1[2]);
			mu_G_G[3] = min(largo_1[3], LARGO_1[3]);
			mu_G_G[4] = min(largo_1[4], LARGO_1[4]);
			mu_G_G[5] = min(largo_1[5], LARGO_1[5]);
			mu_G_G[6] = min(largo_1[6], LARGO_1[6]);
			mu_G_G[7] = min(largo_1[7], LARGO_1[7]);

			mu_B_B[0] = min(largo_2[0], LARGO_2[0]);
			mu_B_B[1] = min(largo_2[1], LARGO_2[1]);
			mu_B_B[2] = min(largo_2[2], LARGO_2[2]);
			mu_B_B[3] = min(largo_2[3], LARGO_2[3]);
			mu_B_B[4] = min(largo_2[4], LARGO_2[4]);
			mu_B_B[5] = min(largo_2[5], LARGO_2[5]);
			mu_B_B[6] = min(largo_2[6], LARGO_2[6]);
			mu_B_B[7] = min(largo_2[7], LARGO_2[7]);

			noise_R_R = max(max(max(max(max(max(max(mu_R_R[0], mu_R_R[1]), mu_R_R[2]), mu_R_R[3]), mu_R_R[4]), mu_R_R[5]), mu_R_R[6]), mu_R_R[7]);
			noise_G_G = max(max(max(max(max(max(max(mu_G_G[0], mu_G_G[1]), mu_G_G[2]), mu_G_G[3]), mu_G_G[4]), mu_G_G[5]), mu_G_G[6]), mu_G_G[7]);
			noise_B_B = max(max(max(max(max(max(max(mu_B_B[0], mu_B_B[1]), mu_B_B[2]), mu_B_B[3]), mu_B_B[4]), mu_B_B[5]), mu_B_B[6]), mu_B_B[7]);

			//printf( "%f",noise_B_B);

			if ((noise_B_B >= 0.3))
			{
			
				float weights[9], sum_weights = 0, hold2, suma = 0;
				for (j = 0; j <= 7; j++)
				{
					sum_weights += (1 - mu_B_B[j]);
				}
				sum_weights = (sum_weights + 3 * sqrt(1 - noise_B_B)) / 2;
				weights[0] = (1 - mu_B_B[0]);
				weights[1] = (1 - mu_B_B[1]);
				weights[2] = (1 - mu_B_B[2]);
				weights[3] = (1 - mu_B_B[7]);
				weights[4] = 3 * sqrt(1 - noise_B_B);
				weights[5] = (1 - mu_B_B[3]);
				weights[6] = (1 - mu_B_B[6]);
				weights[7] = (1 - mu_B_B[5]);
				weights[8] = (1 - mu_B_B[4]);

				for (j = 0; j <= 8; j++)
				{
					for (x = 0; x <= 7; x++)
					{
						if (vectB[x] > vectB[x + 1])
						{
							hold = vectB[x];
							hold2 = weights[x];
							vectB[x] = vectB[x + 1];
							weights[x] = weights[x + 1];
							vectB[x + 1] = hold;
							weights[x + 1] = hold2;
						}
					}
				}
				for (j = 8; j >= 0; j--)
				{
					suma += weights[j];
					if (suma >= sum_weights)
					{
						if (j < 2)
						{
							sum_weights = sum_weights - (weights[0] + weights[1]);
							sum_weights = sum_weights / 2;
							suma = 0;
							for (F = 8; F >= 2; F--)
							{
								suma += weights[F];
								if (suma > sum_weights)
								{
									d_Pout[(Row * m + Col) * channels + 2] = vectB[F];
									F = -1;
								}
							}
							j = -1;
						}
						else
						{
							d_Pout[(Row * m + Col) * channels + 2] = vectB[j];
							//d_Pout[(Row * m + Col) * channels + 0] = d_Pout[(Row * m + Col) * channels + 0];
							j = -1;
						}
						suma = -1;
					}
				}
				//		fwrite (&CCC, 1, 1, header_file);
			}
			else
			{
				d_Pout[(Row * m + Col) * channels + 2] = vectB[4];
				//d_Pout[(Row * m + Col) * channels + 0] = 0;

				//		fwrite (&CCC, 1, 1, header_file);
			}

			if (noise_G_G >= 0.3)
			{
	
				float weights[9], sum_weights = 0, hold2, suma = 0;
				for (j = 0; j <= 7; j++)
				{
					sum_weights += (1 - mu_G_G[j]);
				}
				sum_weights = (sum_weights + 3 * sqrt(1 - noise_G_G)) / 2;
				weights[0] = (1 - mu_G_G[0]);
				weights[1] = (1 - mu_G_G[1]);
				weights[2] = (1 - mu_G_G[2]);
				weights[3] = (1 - mu_G_G[7]);
				weights[4] = 3 * sqrt(1 - noise_G_G);
				weights[5] = (1 - mu_G_G[3]);
				weights[6] = (1 - mu_G_G[6]);
				weights[7] = (1 - mu_G_G[5]);
				weights[8] = (1 - mu_G_G[4]);
				for (j = 0; j <= 8; j++)
				{
					for (x = 0; x <= 7; x++)
					{
						if (vectG[x] > vectG[x + 1])
						{
							hold = vectG[x];
							hold2 = weights[x];
							vectG[x] = vectG[x + 1];
							weights[x] = weights[x + 1];
							vectG[x + 1] = hold;
							weights[x + 1] = hold2;
						}
					}
				}
				for (j = 8; j >= 0; j--)
				{
					suma += weights[j];
					if (suma >= sum_weights)
					{
						if (j < 2)
						{
							sum_weights = sum_weights - (weights[0] + weights[1]);
							sum_weights = sum_weights / 2;
							suma = 0;
							for (F = 8; F >= 2; F--)
							{
								suma += weights[F];
								if (suma >= sum_weights)
								{
									d_Pout[(Row * m + Col) * channels + 1] = vectG[F];
									F = -1;
								}
							}
							j = -1;
						}
						else
						{
							d_Pout[(Row * m + Col) * channels + 1] = vectG[j];
							j = -1;
						}
						suma = -1;
					}
				}
				//		fwrite (&BBB, 1, 1, header_file);
			}
			else
			{
				d_Pout[(Row * m + Col) * channels + 1] = vectG[4];
				//		fwrite (&BBB, 1, 1, header_file);
			}

			if (noise_R_R >= 0.3)
			{

				float weights[9], sum_weights = 0, hold2, suma = 0;
				for (j = 0; j <= 7; j++)
				{
					sum_weights += (1 - mu_R_R[j]);
				}
				sum_weights = (sum_weights + 3 * sqrt(1 - noise_R_R)) / 2;
				weights[0] = (1 - mu_R_R[0]);
				weights[1] = (1 - mu_R_R[1]);
				weights[2] = (1 - mu_R_R[2]);
				weights[3] = (1 - mu_R_R[7]);
				weights[4] = 3 * sqrt(1 - noise_R_R);
				weights[5] = (1 - mu_R_R[3]);
				weights[6] = (1 - mu_R_R[6]);
				weights[7] = (1 - mu_R_R[5]);
				weights[8] = (1 - mu_R_R[4]);
				for (j = 0; j <= 8; j++)
				{
					for (x = 0; x <= 7; x++)
					{
						if (vectR[x] > vectR[x + 1])
						{
							hold = vectR[x];
							hold2 = weights[x];
							vectR[x] = vectR[x + 1];
							weights[x] = weights[x + 1];
							vectR[x + 1] = hold;
							weights[x + 1] = hold2;
						}
					}
				}
				for (j = 8; j >= 0; j--)
				{
					suma += weights[j];
					if (suma >= sum_weights)
					{
						if (j < 2)
						{
							sum_weights = sum_weights - (weights[0] + weights[1]);
							sum_weights = sum_weights / 2;
							suma = 0;
							for (F = 8; F >= 2; F--)
							{
								suma += weights[F];
								if (suma > sum_weights)
								{
									d_Pout[(Row * m + Col) * channels + 0] = vectR[F];
									F = -1;
								}
							}
							j = -1;
						}
						else
						{
							d_Pout[(Row * m + Col) * channels + 0] = vectR[j];
							j = -1;
						}
						suma = -1;
					}
				}
				//      fwrite (&AAA, 1, 1, header_file);
			}
			else
			{
				d_Pout[(Row * m + Col) * channels + 0] = vectR[4];
				//		fwrite (&AAA, 1, 1, header_file);
			}

		}
	}
}