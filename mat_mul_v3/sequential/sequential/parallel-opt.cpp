// yas-concurrent.cpp : Defines the entry point for the console application.


#include "stdafx.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
# include <omp.h>
#include <chrono>
#include <fstream>
#include <algorithm> 
#include <immintrin.h>

using namespace std::chrono;
using namespace std;

const int g_outputWidth = 9; // Set the output width for elements.
const int SAMPLE_SIZE = 6;	// Set sample size
const int g_outputPrecision = 15; // Output elements precision.
const int LEVEL_TWO_CACHE_SIZE_KB = 256; // L1cache size* no of physical cores
const int g_cacheBlockSize = floor(sqrt((LEVEL_TWO_CACHE_SIZE_KB * 1024) / 4));	// Block size for loop tilling/blocking.


//populate matrix with random values.
double** generateMatrix(int n){
	double x;
	double max = DBL_MAX;
	double min = DBL_MIN;
	double** matA = new double*[n];
	for (int i = 0; i < n; i++) {
		matA[i] = new double[n];
		for (int j = 0; j < n; j++) {
			double randVal = (double)rand() / RAND_MAX;
			matA[i][j] = min + randVal * (max - min);
		}
	}
	return matA;
}

//generate matrix for final result.
double** generateMatrixFinal(int n){
	double f;
	double** matA = new double*[n];
	for (int i = 0; i < n; i++) {
		matA[i] = new double[n];
		for (int j = 0; j < n; j++) {
			matA[i][j] = 0;
		}
	}
	return matA;
}

//matrix multiplication - Cache blocking (parameter),SIMD instruction (x8 float, x4 double) and OpenMP 'parallel for' on outermost loop.
double matrixMultiplicationOptimized_doubleSMID(double** A, double** B, double** C, int n) {
	int i, j, k, l;
	int limit0 = n; 			// Index i limit 
	int limit1 = n; 			// Index j limit
	int limit2 = n; 			// Index k limit
	int aux_i, aux_j, aux_k;
	int aux_limit_i; 	 			// Block index limit i
	int aux_limit_j; 	 			// Block index limit j
	int aux_limit_k; 	 			// Block index limit k
	int unroll_factor = g_cacheBlockSize;
	int unroll_limit; 	 			// Loop unroll index limit
	clock_t begin_time = clock();
	cout << "\nn: " << "clock started" << "\n";
#pragma omp parallel for private(i, j, k, aux_i, aux_j, aux_k, aux_limit_i, aux_limit_j, aux_limit_k, unroll_limit) shared (A,B,C)
	for (i = 0; i < limit0; i += g_cacheBlockSize) {
		// Blocking index i limit
		aux_limit_i = min((i + g_cacheBlockSize), limit0);

		for (j = 0; j < limit1; j += g_cacheBlockSize) {
			// Blocking index j limit
			aux_limit_j = min((j + g_cacheBlockSize), limit1);

			for (k = 0; k < limit2; k += g_cacheBlockSize) {
				// Blocking index k limit
				aux_limit_k = min((k + g_cacheBlockSize), limit2);

				unroll_limit = aux_limit_k - (unroll_factor - 1); // Unrolling by factor of 4

				for (aux_i = i; aux_i < aux_limit_i; ++aux_i) {
					for (aux_j = j; aux_j < aux_limit_j; ++aux_j) {

						double zero = 0;
						__m256d acc = _mm256_broadcast_sd(&zero);

						// Unrolling for k loop
						for (aux_k = k; aux_k < unroll_limit; aux_k += unroll_factor) {
							acc = _mm256_add_pd(acc,
								_mm256_mul_pd(_mm256_load_pd(&A[aux_i][aux_k]), _mm256_load_pd(&B[aux_j][aux_k])));
						}

						// Gather possible uncounted elements
						for (; aux_k < aux_limit_k; ++aux_k)
							C[aux_i][aux_j] += A[aux_i][aux_k] * B[aux_j][aux_k];


						// Sum up everything
						double acc_vet[4];

						_mm256_storeu_pd(acc_vet, acc);

						C[aux_i][aux_j] += acc_vet[0] + acc_vet[1] + acc_vet[2] + acc_vet[3];

					}
				}
			}
		}
	}
	//ofstream out("filename.txt", ios::out | ios::app);
	double t = float(clock() - begin_time);
	//out << "parallelSMID: " << t << "\n";
	//out.close();
	return t;
}

// transpose of a matrix
double** trn(double **b, double **t)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			t[j][i] = b[i][j];
		}
	}
	return t;
}



int _tmain(int argc, _TCHAR* argv[])
{

	ofstream out("filename.txt", ios::out | ios::app);
	out << "--------------STARTED--------------" << "\n";
	int start = 200, stop = 2000, step = 200;

	for (int n = start; n <= stop; n += step)
	{

		srand(time(NULL));
		cout << "\nn: " << n << "\n";
		//    initialize a Matrix of size 5
		double t3 = 0;

		int my_size = n;
		// Calculate time for n samples.
		for (int i = 0; i < SAMPLE_SIZE; i++){
			double **A = generateMatrix(my_size);
			cout << "\nn: " << "generated a" << "\n";
			double **B = generateMatrix(my_size);
			cout << "\nn: " << "generated b" << "\n";
			double **C = generateMatrixFinal(my_size);
			double **T = generateMatrixFinal(my_size);
			cout << "\nn: " << "generated c" << "\n";
			T = trn(B, T);
			t3 += matrixMultiplicationOptimized_doubleSMID(A, T, C, n);
			delete A;
			delete B;
			delete C;
			delete T;
				}
		double T3 = 0;
		T3 = t3 / SAMPLE_SIZE;
	
		out << "-----Mean Time----" << n << " - matrix size" << "\n";
		out << "parallelSMID: " << T3 << "\n";
		}
	out << "-----------FINISHED-----------------" << "\n";
	out.close();
	return 0;
}

