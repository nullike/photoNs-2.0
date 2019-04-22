#ifndef PARTMESH_H
#define PARTMESH_H

/* 
 * The convolutions for particle-method is based on FFT method 
 * We employ the 2DCOMP lib for two dimensional FFT.
 * 
 * Reference : N. Li and S. Laizet. 2DECOMP&FFT – A highly scalable 2D
 *             decomposition library  and FFT interface. In Cray User 
 *             Group 2010 conference, Edinburgh, 2010 
 */

void partition_fft();
void partition_fft2();

void partmesh_thread();
void partmesh();

void powerspectrum(char powname[]) ;

#ifdef PMTHREAD
void* pm_thread(void *arg);
#endif

#define convolution convolution_
#define conv_pmonly conv_pmonly_
#define densitykspace densitykspace_


#endif
