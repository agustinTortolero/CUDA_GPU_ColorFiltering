#pragma once
void FTSCF_CPU_Multi(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);
void FTSCF_CPU_Multi2(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);//Segundo intento
void FTSCF_CPU_Multi_Original(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);
