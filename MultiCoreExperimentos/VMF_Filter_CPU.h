#pragma once
void VMF_CPU_Single(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels);
void VMF_CPU_Single_Optimizado(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels);

void VMF_CPU_Multi(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);
void VMF_CPU_Multi_Optimizado(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);

void VMF_CPU_Multi_L1(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);
void VMF_CPU_Multi_L2(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);

void VMF_CPU_Multi_Reuso(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads);
void VMF_CPU_Multi_Reuso_Sort(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads);
void VMF_CPU_Multi_Sort(unsigned char* d_Pout, const unsigned char* d_Pin, int n, int m, int channels, int nThreads);