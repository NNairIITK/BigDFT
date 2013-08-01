#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "s_gpu.h"

#include "type.h"

#define MASTER_TASK 0
#define MATRIX_W 320
// The height of the matrix must be a multiple of "number of MPI tasks - 1"
#define MATRIX_H 320
#define MAX_TASKS 128

int rank, nb_tasks;

// This function is implemented in add_kernel.cu
void callback_add ( void * );

// S_GPU trace functions
void callback_beg ( int i )
{
	printf ( ">> beg on card %d on rank %d\n", i, rank );
}
void callback_end ( int i )
{
	printf ( "<< end on card %d on rank %d\n", i, rank );
}
void callback_mg_beg ( int i )
{
	printf ( ">> beg mem gpu on card %d on rank %d\n", i, rank );
}
void callback_mg_end ( int i )
{
	printf ( "<< end mem gpu on card %d on rank %d\n", i, rank );
}
void callback_mc_beg ( int i )
{
	printf ( ">> beg memcpy on card %d on rank %d\n", i, rank );
}
void callback_mc_end ( int i )
{
	printf ( "<< end memcpy on card %d on rank %d\n", i, rank );
}


int main ( int narg, char *args[] )
{
	int i, j;

	// MPI initialization
	MPI_Init ( &narg, &args );
	MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
	MPI_Comm_size ( MPI_COMM_WORLD, &nb_tasks );

	if ( nb_tasks > MAX_TASKS )
	{
		printf ( "Too many MPI tasks: %d\n", nb_tasks );
		return 1;
	}
	if ( MATRIX_H % ( nb_tasks - 1 ) != 0 )
  {
    printf( "The height of the matrix must be a multiple of \"number of MPI tasks - 1\"\nmatrix height = %d\nnumber of MPI tasks = %d", MATRIX_H, nb_tasks);
  }

	MPI_Request ireq[MAX_TASKS];  // asynch request, assume size<128
	MPI_Status stats[MAX_TASKS];  // status of asynch communication
	MPI_Status stat;              // status of asynch communication

	// S_GPU initialization
	int gpushare, usegpu;
  // This initialization function must be called by all the MPI task
  // It must be called befaore any other s_gpu function
	sg_init ( &gpushare, &usegpu, rank );
  // tracer_callbacks_t is a structure which contains the pointers to the trace functions
	tracer_callbacks_t tracs;
  // This operation is crucial when using trace functions
  // If, for example, we decide to use only 4 of the 6 functions available for tracing,
  // undeclared function pointers must be set to NULL.
	memset ( &tracs, ( int ) NULL, sizeof ( tracer_callbacks_t ) );
	tracs.calc_beg = &callback_beg;
	tracs.calc_end = &callback_end;
	tracs.mem_gpu_beg = &callback_mg_beg;
	tracs.mem_gpu_end = &callback_mg_end;
//	tracs.memcpy_beg = &callback_mc_beg;
//	tracs.memcpy_end = &callback_mc_end;
//	add_tracer_callbacks ( &tracs );

  // For this example, we split the matrix in ( nb_tasks - 1 ) equal parts
  // Each mpi client task will compute MATRIX_H / ( nb_tasks - 1 ) rows of the matrix
	int nb_rows_by_task = MATRIX_H / ( nb_tasks - 1 );

	// Master node computation
	if ( rank == MASTER_TASK )
	{
		double *matrix_a, *matrix_b, *result_matrix;
		size_t matrix_size = MATRIX_H * MATRIX_W * sizeof ( double );
		matrix_a = ( double* ) malloc ( matrix_size );
		matrix_b = ( double* ) malloc ( matrix_size );
		result_matrix = ( double* ) malloc ( matrix_size );

    // Initialization of the matrices
    srand ( 0 );
		for ( i = 0; i < MATRIX_W; i++ )
		{
			for ( j = 0; j < MATRIX_H; j++ )
			{
//        matrix_a[i][j] = (double)(rand() / (RAND_MAX + 1.0));
//        matrix_b[i][j] = (double)(rand() / (RAND_MAX + 1.0));
				matrix_a[i *MATRIX_W + j] = ( double ) i;
				matrix_b[i *MATRIX_W + j] = ( double ) j;
			}
		}

		// Send intervals of matrices A and B to worker processes
		for ( i = 1; i < nb_tasks; i++ )
			MPI_Isend ( matrix_a + ( ( i - 1 ) *nb_rows_by_task*MATRIX_W ), nb_rows_by_task*MATRIX_W, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &ireq[i-1] );
		for ( i = 1; i < nb_tasks; i++ )
			MPI_Isend ( matrix_b + ( ( i - 1 ) *nb_rows_by_task*MATRIX_W ), nb_rows_by_task*MATRIX_W, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &ireq[nb_tasks -1 + i - 1] );
		MPI_Waitall ( nb_tasks * 2 - 2, ireq, stats );

    // Receive all the result intervals
		for ( i = 1; i < nb_tasks; i++ )
			MPI_Irecv ( result_matrix + ( ( i - 1 ) *nb_rows_by_task*MATRIX_W ), nb_rows_by_task*MATRIX_W, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &ireq[i-1] );
		MPI_Waitall ( nb_tasks - 1 , ireq, stats );

//		printf ( "=====\n" );
//		for ( i = 0; i < MATRIX_W; i++ )
//		{
//			for ( j = 0; j < MATRIX_H; j++ )
//			{
//				printf ( "%f ", result_matrix[i * MATRIX_W + j] );
//			}
//			printf ( "\n" );
//		}
//		printf ( "=====\n" );
	}
	else
	{
    // Client node
		double *_a_h, *_b_h, *_c_h; // Host pointers
		double *_a_p, *_b_p, *_c_p; // Pinned pointers
		double *_a_d, *_b_d, *_c_d; // Device pointers

    // Allocate memory to receive the interval on this MPI node
		_a_h = ( double* ) malloc ( nb_rows_by_task * MATRIX_W * sizeof ( double ) );
		_b_h = ( double* ) malloc ( nb_rows_by_task * MATRIX_W * sizeof ( double ) );
		_c_h = ( double* ) malloc ( nb_rows_by_task * MATRIX_W * sizeof ( double ) );

    // Receive the interval of the matrices A and B
		MPI_Recv ( _a_h, nb_rows_by_task*MATRIX_W, MPI_DOUBLE, MASTER_TASK, rank, MPI_COMM_WORLD, &stat );
		MPI_Recv ( _b_h, nb_rows_by_task*MATRIX_W, MPI_DOUBLE, MASTER_TASK, rank, MPI_COMM_WORLD, &stat );

    // In order to benefit the asynchroneous mem copy CPU -> GPU (and GPU -> CPU), the host memory must be a special memory
    // called pinned memory. Before any transfer to the GPU, the host memory is first copied to a pinned memory area.
    // Allocate pinned memory
    // This operation is immediate and blocking
		sg_cpu_malloc_pinned ( ( void** ) &_a_p, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2 );
		sg_cpu_malloc_pinned ( ( void** ) &_b_p, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2 );
		sg_cpu_malloc_pinned ( ( void** ) &_c_p, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2 );

    // Allocate memory on the GPU device
    // This operation is immediate and blocking
		sg_gpu_malloc ( ( void** ) &_a_d, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2 );
		sg_gpu_malloc ( ( void** ) &_b_d, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2 );
		sg_gpu_malloc ( ( void** ) &_c_d, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2 );

    // Create the stream
		sg_stream_ptr_t stream_ptr = sg_create_stream();
    // Send the second half of the part of the matrices A and B on the GPU for computation
    // This operation is added in the stream. The execution is delayed.
		sg_gpu_send_mem ( _a_d, &_a_h[nb_rows_by_task*MATRIX_W/2], _a_p, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2, stream_ptr );
		sg_gpu_send_mem ( _b_d, &_b_h[nb_rows_by_task*MATRIX_W/2], _b_p, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2, stream_ptr );

    // The parameter must be passed to the callback function through a structure.
    // This structure is defined in type.h
		param_t param;
		param.a_d = _a_d;
		param.b_d = _b_d;
		param.c_d = _c_d;
		param.N = nb_rows_by_task * MATRIX_W / 2;

    // Add a function call to the stream, with a struct as parameter.
    // This operation is added in the stream. The execution is delayed.
		sg_gpu_send_calc ( &callback_add, &param, sizeof ( param_t ), stream_ptr );

		// Retrieve results from device and store it in host array
    // This operation is added in the stream. The execution is delayed.
		sg_gpu_recv_mem ( &_c_h[nb_rows_by_task*MATRIX_W/2], _c_d, _c_p, nb_rows_by_task*MATRIX_W*sizeof ( double ) / 2, stream_ptr );

    // The operations added to the stream previously are now executed asynchroneously
    // If several streams have been created, they would be executed concurently
		sg_exec_all_streams();

    // In the same time, the MPI node compute the other half of the matrix addition
		for ( i = 0; i < nb_rows_by_task / 2; i++ )
		{
			for ( j = 0; j < MATRIX_W; j++ )
			{
				_c_h[i * MATRIX_W + j] = _a_h[i * MATRIX_W + j] + _b_h[i * MATRIX_W + j];
			}
		}

    // The MPI node send the result of its interval to the master node
		MPI_Send ( _c_h, nb_rows_by_task*MATRIX_W, MPI_DOUBLE, MASTER_TASK, rank, MPI_COMM_WORLD );
	}
	MPI_Finalize();
}

