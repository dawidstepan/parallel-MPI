#include <mpi.h>
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <iostream>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>

void fill_matrix(double* &,const int &); 
void calc_algo_row(double* &,const int &,const int &,const int &,const int &); 
void print_matrix(double* &,const int &,const int &,const int &);
void transpose(double* &,const int &);
void teststripes(double* &, double* &, int &, int &, int &, int &);
void calc_algo_serial(double* &, const int &, const int &);
void parallel(double* &, double* &, int &, int &, const int*, const int*, int &, int &, int &);
double check_sum(double* &A, const int &n) {

	int i, j; 
	double res = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			res += A[j + i*n];
		}		
	}
	return res; 
}

int main(int argc, char *argv[])
{
	int rank, numprocs, N, count, counts, remainder, stripsize, steps;
	double timer,timers, t_start, t_end, ts_start, ts_end = 0;
	double* matrix		= NULL;
	double* serial		= NULL;
	double* matrixPart	= NULL; 
	int* rscounts 		= 0;
	int* rsdispls 		= 0;
	int impl 		= 0;
	int proc		= 1; 

	srand(time(NULL));
	MPI_Init (&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     
	if (rank == 0) { 

		if(argc > 1) {
			N = atoi(argv[1]);
			steps = atoi(argv[2]);
			impl = atoi(argv[3]);
		   }
		   else {
			printf("Matrix dimension: "); 
			std::cin >> N;
			printf("Steps: ");
			std::cin >> steps;
		   } 

		if (N < 1 || steps < 1) return EXIT_FAILURE;
		
		matrix 	 = new double[N * N];
		serial 	 = new double[N * N];
		rscounts = new int[numprocs];
		rsdispls = new int[numprocs];

		fill_matrix(matrix,N);
		fill_matrix(serial,N);

		if(impl == 1) {
			ts_start = MPI_Wtime();
			calc_algo_serial(serial,N,steps);
			ts_end = MPI_Wtime();
			timers = ts_end - ts_start;
		}
		count = N / numprocs;
		remainder = N - count * numprocs;
		int sum = 0;
		for (int i = 0; i < numprocs; ++i) {
			rscounts[i] = (i < remainder) ? count + 1 : count;
			rsdispls[i] = sum * N;
			sum += rscounts[i];
			rscounts[i] = rscounts[i] * N;
		}		
	}
	
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (0 != rank) {
		count = N / numprocs;
		remainder = N - count * numprocs;
		
	}
	
	stripsize = rank < remainder ? count + 1 : count;
	counts = stripsize * N; 
	
	matrixPart = new double[stripsize * N];	
	MPI_Barrier(MPI_COMM_WORLD);
	t_start = MPI_Wtime();
		
	for (int i = 1; i <= steps; i++ ) {
		parallel(matrix, matrixPart, stripsize, N, rscounts, rsdispls, counts, rank, numprocs);		
	}

	t_end = MPI_Wtime();
	timer = t_end - t_start;
	MPI_Barrier(MPI_COMM_WORLD);

	delete[] matrixPart;
	if (0 != rank) {
		delete[] matrix;
		delete[] serial;
	}

	if (0 == rank)
	{
		double res = check_sum(matrix, N);
			
		if(impl == 1) {
			double ress = check_sum(serial, N);
			printf("Results: \n");  
			printf("Checksum parallel: %lf\n", res);
			printf("Checksum serial: %lf\n", ress);	
			printf("Processed time parallel: %lf\n", timer);
			printf("Processed time serial: %lf\n", timers);	
		}

		std::ofstream log;
  		log.open ("log.txt", std::ofstream::app);
  		log << N <<  "," << numprocs << "," << timer << "," << std::setprecision(20) << res << std::endl;
  		log.close();

		delete[] matrix; 
 		delete[] serial;
	}

	MPI_Finalize();

	return 0;
}

void print_matrix(double* &A,const int &col, const int &row, const int &n) { // Row major matrix
	
	int i, j; 
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
		printf("%f\t", A[j + i*n]);
		}
	printf("\n");

	}
}

void fill_matrix(double* &A, const int &n) {

	int i; 	
	for(i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
		  A[j + i*n] = double((i*2.0+j*3.0)/3.0);
		}
	}

	for (i = 0; i < n; i++) {
	      A[i]	     = 9.50;
	      A[i + (n-1)*n] = 99.50;
	      A[i*n]	     = 999.50;
	      A[i*n + (n-1)] = 9999.50;
	}
}

void transpose(double* &A, const int &n) {

	int i, j; 
	#pragma omp parallel for	
	for (i = 0; i < n; i++) {
			for (j = i+1; j < n; j++) {
				std::swap(A[j + i*n], A[i + j*n]);
			}
	}
}

void calc_algo_row(double* &A, const int &n, const int &row, const int &ProcRank, const int &ProcNum) {

	int i, j, start, end;
	
	if(ProcRank == 0) {
		start = 1; 
		end = row-1;		
	}
	else if(ProcRank == ProcNum - 1 ) {
		start = 0;
		end = row-2;
	}
	else{
		start = 0; 
		end = row-1;
	}

	// forward/downward
	#pragma omp parallel for	
	for (i = start; i <= end; i++) {
		for (j = 1; j <= n-2; j++) {
			A[i*n + j] = (A[i*n + j] + A[i*n + (j-1)]) / 2.0 +
				     (A[i*n + j] / 3.0 + A[i*n + (j-1)] / 3.0) / n - 
				     (A[i*n + j] / 5.0 + A[i*n + (j-1)] / 5.0) / n;
			}
		}
	
	// backward/unpward
	#pragma omp parallel for	
	for (i = start; i <= end; i++) {
		for (j = n-2; j >= 1; j--) {
			A[i*n + j] = (A[i*n + j] + A[i*n + (j+1)]) / 2.0 +
				     (A[i*n + j] / 3.0 + A[i*n + (j+1)] / 3.0) / n - 
				     (A[i*n + j] / 5.0 + A[i*n + (j+1)] / 5.0) / n;
			}
		}
}

void calc_algo_serial(double* &A, const int &n, const int &steps) {

	int s,i,j;
	for (s = 1; s <= steps; s++) {	
	      	// Phase 1: forward/backward along rows
	     	// forward
		for (i = 1; i <= n-2; i++) {
			for (j = 1; j <= n-2; j++) {
				A[i*n + j] = (A[i*n + j] + A[i*n + (j-1)]) / 2.0 +
					     (A[i*n + j] / 3.0 + A[i*n + (j-1)] / 3.0) / n - 
					     (A[i*n + j] / 5.0 + A[i*n + (j-1)] / 5.0) / n;
			}
		}

		// backward
		for (i = 1; i <= n-2; i++) {
			 for (j = n-2; j >= 1; j--) {
				A[i*n + j] = (A[i*n + j] + A[i*n + (j+1)]) / 2.0 +
					     (A[i*n + j] / 3.0 + A[i*n + (j+1)] / 3.0) / n - 
					     (A[i*n + j] / 5.0 + A[i*n + (j+1)] / 5.0) / n;
			}
		}
		
		// Phase 2: downward/upward along columns
		// downward
		for (i = 1; i <= n-2; i++) {
			for (j = 1; j <= n-2; j++) {
				A[i*n + j] = (A[i*n + j] + A[(i-1)*n + j]) / 2.0 +
					     (A[i*n + j] / 3.0 + A[(i-1)*n + j] / 3.0) / n - 
					     (A[i*n + j] / 5.0 + A[(i-1)*n + j] / 5.0) / n;
			} 
		}

		// upward
		for (i = n-2; i >= 1; i--) {
			for (j = 1; j <= n-2; j++) {
				A[i*n + j] = (A[i*n + j] + A[(i+1)*n + j]) / 2.0 +
					     (A[i*n + j] / 3.0 + A[(i+1)*n + j] / 3.0) / n - 
					     (A[i*n + j] / 5.0 + A[(i+1)*n + j] / 5.0) / n;
			}
		}
	}
}   

void teststripes(double* &A, double* &matrixPart, int &N, int &stripsize, int &rank, int &numprocs) { 

	if (rank == 0) { 
		 printf("Initial Matrix: \n"); 
		 print_matrix(A, N, N, N); 
	} 
	MPI_Barrier(MPI_COMM_WORLD); 
	for (int i=0; i<numprocs; i++) { 
	if (rank == i) { 
		printf("Teststripe:\n"); 
		std::cout << stripsize << std::endl;
		print_matrix(matrixPart, N, stripsize, N); 
		printf("\n"); 
		} 
	MPI_Barrier(MPI_COMM_WORLD);
	}  
} 

void parallel(double* &matrix, double* &matrixPart, int &stripsize, int &N, const int* rscounts, const int* rsdispls, int &counts, int &rank, int &numprocs) {
	
	MPI_Scatterv(matrix, rscounts, rsdispls, MPI_DOUBLE, matrixPart, stripsize * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	calc_algo_row(matrixPart, N, stripsize, rank, numprocs);
	//teststripes(matrix, matrixPart,N, stripsize, rank, numprocs);
	MPI_Gatherv(matrixPart,counts, MPI_DOUBLE, matrix, rscounts, rsdispls, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	if (0 == rank) {
		transpose(matrix,N);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Scatterv(matrix, rscounts, rsdispls, MPI_DOUBLE, matrixPart, stripsize * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	calc_algo_row(matrixPart, N, stripsize, rank, numprocs);
	//teststripes(matrix, matrixPart,N, stripsize, rank, numprocs);
	MPI_Gatherv(matrixPart,counts, MPI_DOUBLE, matrix, rscounts, rsdispls, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	if (0 == rank) {
		transpose(matrix,N);
	}
}






