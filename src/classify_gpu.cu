/*
 * -- Kraken-GPU -- 
 * Jack Wadden - 2015
 */

#include "stdio.h"
#include <stdlib.h>
#include <stdint.h>

// kraken
#include "seqreader.hpp"

using namespace std;
using namespace kraken;

DNASequence global;

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
        a[id] = a[id] + b[id] + 7;
}

__host__ void addRead(char * h_reads, 
                      uint32_t bytesPerRead, 
                      uint32_t readLen, 
                      const char * read, 
                      int readNum) {

    h_reads[bytesPerRead * readNum] = readLen;
    
    // fill remaining bits with bp from read
    // for every 4 bp
    for(int i = 0; i < bytesPerRead; i++) {
        char tmp = 0;
        for(int j = 0; j < 8/2; j++) {
            // shift in 00
            tmp = tmp << 2;
            // if there is a bp left, add it
            if((i*4 + j) < readLen) {
                switch(read[i*4 + j]) {
                case 0 : // A
                    tmp = tmp | 0;
                    break;
                case 1 : // T
                    tmp = tmp | 1;
                    break;
                case 2 : // G
                    tmp = tmp | 2;
                    break;
                case 3 : // C
                    tmp = tmp | 3;
                    break;
                }
            }
        }
        // store every 4 bp into a single char array after length
        h_reads[(bytesPerRead * readNum) + 1 + i] = tmp;
    }
}


/*
 * Kernel wrapper is called and given a ptr to the database file
 *  Database is passed to the GPU in the following format
 *  |read len (8)| read (512) |
 *
 */
void kernel_wrapper(int maxReadLen, vector<DNASequence> reads)
{

    // Convert vector of strings to GPU appropriate data structure
    int numReads = reads.size();
    const unsigned int bytesPerRead =  
        sizeof(char) + // spot for read len (8 bits for now) 
        ceil((maxReadLen*2)/8); // spot for read

    // Malloc appropriate size pinned array on host
    char *h_reads;
    gpuErrchk(cudaMallocHost( (void **)&h_reads, numReads * bytesPerRead));

    // Copy each read to compressed format on host
    for(int i = 0; i < numReads; i++) {
        // add read
        addRead(h_reads, 
                bytesPerRead, 
                reads[i].seq.length(), 
                reads[i].seq.c_str(), 
                i);
        printf("length: %c\n", h_reads[i*bytesPerRead]);        
    }



    /*
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
    */
    gpuErrchk( cudaFreeHost(h_reads) );
}
