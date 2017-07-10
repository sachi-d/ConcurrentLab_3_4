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
#include <cfloat>
#include <limits>
#include <math.h>

using namespace std::chrono;
using namespace std;

const int SAMPLE_SIZE = 20;	// Set sample size

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

//matrix multiplication - sequential
double matrixMultiplicationSequential(double** A, double** B, double** C, int n){
	int i, j, k, l;
	clock_t begin_time = clock();
	cout << "\nn: " << "clock started" << "\n";

	for (i = 0; i < n; i++) {
		for (j = 0; j <n; j++) {
			for (k = 0; k < n; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	//ofstream out("filename.txt", ios::out | ios::app);
	double t = float(clock() - begin_time);
	//out << "sequential: " << t << "\n";
	//out.close();
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
		double t0 = 0;
		int my_size = n;
		// Calculate time for n samples.
		double times[SAMPLE_SIZE]={};
		for (int i = 0; i < SAMPLE_SIZE; i++){
			double **A = generateMatrix(my_size);
			cout << "\nn: " << "generated a" << "\n";
			double **B = generateMatrix(my_size);
			cout << "\nn: " << "generated b" << "\n";
			double **C = generateMatrixFinal(my_size);
			double **T = generateMatrixFinal(my_size);
			cout << "\nn: " << "generated c" << "\n";
			double single_sample_time=matrixMultiplicationSequential(A, B, C, n);
			t0+=single_sample_time;
			delete A;
			delete B;
			delete C;
		}
		double T0 = 0;
		T0 = t0 / SAMPLE_SIZE;
		out << "-----Mean Time----" << n << " - matrix size" << "\n";
		out << "Sequential: " << T0 << "\n";

        //calculating std deviation
		double sq=0;
		for(int k=0;k<SAMPLE_SIZE;k++){
            sq+=(times[k]-T0)*(times[k]-T0);
		}
		double std_dev = sqrt(sq/SAMPLE_SIZE);
		out<< "standard deviation = " << std_dev << "\n";

		//calculating sample size
		double samplesize = ((196*std_dev)/(5*T0))* ((196*std_dev)/(5*T0)) ;
		out << "sample size = " << samplesize << "\n";
	}
	out << "-----------FINISHED-----------------" << "\n";
	out.close();
	return 0;
}

