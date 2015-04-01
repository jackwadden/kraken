#include "stdio.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess) 
        {
            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
        }
}

__global__ void vector_add(int *a, int *b, int length)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;

    if(id < length)
        a[id] = a[id] + b[id] + 5;
}

void kernel_wrapper(int *a, int *b, int length)
{

    int *a_d;
    int *b_d;

    int blockSize = 32;
    int numBlocks = (int)(length / blockSize);

    if(length % blockSize) {
        numBlocks++;
    }

    dim3 threads( blockSize, 1 );
    dim3 blocks( numBlocks, 1 );

    size_t byteLength = length * sizeof(int);

    gpuErrchk(cudaMalloc( (void **)&a_d, byteLength ));
    gpuErrchk(cudaMalloc( (void **)&b_d, byteLength ));

    gpuErrchk( cudaMemcpy( a_d, a, byteLength, cudaMemcpyHostToDevice ));
    gpuErrchk( cudaMemcpy( b_d, b, byteLength, cudaMemcpyHostToDevice ));

    vector_add<<< blocks, threads >>>( a_d, b_d , length);
    gpuErrchk( cudaPeekAtLastError() );

    gpuErrchk( cudaMemcpy( a, a_d, byteLength, cudaMemcpyDeviceToHost ));
    gpuErrchk( cudaMemcpy( b, b_d, byteLength, cudaMemcpyDeviceToHost ));

    gpuErrchk( cudaFree(a_d) );
    gpuErrchk( cudaFree(b_d) );

}
