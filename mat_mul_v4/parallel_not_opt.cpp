// yas-concurrent.cpp : Defines the entry point for the console application.


#include "includes/stdafx.h"

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
#include <cfloat>
#include <limits>
#include <math.h>

using namespace std::chrono;
using namespace std;

const int SAMPLE_SIZE = 20;	// Set sample size

//populate matrix with random values.
double** generateMatrix(int n){ 
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
	double** matA = new double*[n];
	for (int i = 0; i < n; i++) {
		matA[i] = new double[n];
		for (int j = 0; j < n; j++) {
			matA[i][j] = 0;
		}
	}
	return matA;
}

//matrix multiplication - parallel
double matrixMultiplicationParallel(double** A, double** B, double** C, int n){
	int i, j, k;
	clock_t begin_time = clock();
	//	cout << "clock started" << "\n";
# pragma omp parallel shared ( A,B,C,n  ) // private ( i, j, k )
	{
# pragma omp for
		for (i = 0; i < n; i++) {
			//            cout<< i << ", " ;
			for (j = 0; j < n; j++) {
				for (k = 0; k < n; k++) {
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
	}
	//ofstream out("filename.txt", ios::out | ios::app);
	double t = float(clock() - begin_time);
	//out << "parallel: " << t << "\n";
	//out.close();
	return t;
}

int _tmain(int argc, _TCHAR* argv[])
{

	ofstream out("results_parallel_NOT_opt.txt", ios::out | ios::app);
	out << "Sample size = " << SAMPLE_SIZE << "\n";
	out << "--------------STARTED--------------" << "\n";
	int start = 200, stop = 2000, step = 200;

	for (int n = start; n <= stop; n += step)
	{

		srand(time(NULL));
		cout << "\nn: " << n << "\n";
		//    initialize a Matrix of size 5
		double t1 = 0;

		int my_size = n;
		// Calculate time for n samples.
		double times[SAMPLE_SIZE] = {};
		for (int i = 0; i < SAMPLE_SIZE; i++){
			//            cout << "sample #" << i << " : ";
			double **A = generateMatrix(my_size);
			//			cout << "\nn: " << "generated a" << "\n";
			double **B = generateMatrix(my_size);
			//			cout << "\nn: " << "generated b" << "\n";
			double **C = generateMatrixFinal(my_size);
			//			cout << "\nn: " << "generated c" << "\n";
			double single_sample_time = matrixMultiplicationParallel(A, B, C, n);
			//            cout << single_sample_time;
			times[i] = single_sample_time;

			t1 += single_sample_time;
			for (int i = 0; i < n; i++) {
				delete[] A[i];
				delete[] B[i];
				delete[] C[i];
			}
			delete[] A;
			delete[] B;
			delete[] C;
		}
		double T1 = 0;
		T1 = t1 / SAMPLE_SIZE;
		out << "-----Mean Time----" << n << " - matrix size" << "\n";
		cout << "mean = " << T1 << "\n";
		out << "Parallel mean: " << T1 << "\n";

		//calculating std deviation
		double sq = 0;
		for (int k = 0; k<SAMPLE_SIZE; k++){
			sq += (times[k] - T1)*(times[k] - T1);
		}
		double std_dev = sqrt(sq / SAMPLE_SIZE);
		cout << "standard deviation = " << std_dev << "\n";

		//calculating sample size
		double samplesize = ((196 * std_dev) / (5 * T1))* ((196 * std_dev) / (5 * T1));
		cout << "sample size = " << samplesize << "\n";
	}
	out << "-----------FINISHED-----------------" << "\n";
	out.close();
	return 0;
}



