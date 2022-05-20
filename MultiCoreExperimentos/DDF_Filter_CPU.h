#pragma once
void DDF_CPU_Single(unsigned char* image_out, unsigned char* image_in, int n, int m, int channels);
void DDF_CPU_Multi(unsigned char* image_out, const unsigned char* image_in, int n, int m, int channels, int nThreads);