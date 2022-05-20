#pragma once
void PeerGroup_Multi(unsigned char* image_out, const unsigned char* image_in, unsigned char* image_Noise, int n, int m, int channels, int nThreads);
void Detection_FuzzyMetric(unsigned char* image_Noise, const unsigned char* image_in, int n, int m, int channels, int nThreads);
void AMF_Filtering(unsigned char* image_out, const unsigned char* image_in, unsigned char* image_Noise, int n, int m, int channels, int nThreads);
void VMF_Filtering(unsigned char* image_out, const unsigned char* image_in, unsigned char* image_Noise, int n, int m, int channels, int nThreads);
void AMF_Filtering_ExpWindow(unsigned char* image_out, const unsigned char* image_in, unsigned char* image_Noise, int n, int m, int channels, int nThreads);
void VMF_Filtering_Single(unsigned char* image_out, const unsigned char* image_in, unsigned char* image_Noise, int n, int m, int channels);
void FiltradoPropuesta(unsigned char* image_out, const unsigned char* image_in, unsigned char* image_Noise, int n, int m, int channels, int nThreads);
void MPGFMF(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);