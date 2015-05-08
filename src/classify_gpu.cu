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
    //test comment -Dor
    if(id < length)
        a[id] = a[id] + b[id] + 7;
}
__constant__ uint64_t INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;


// Code mostly from Jellyfish 1.6 source                                                                                             
__device__ uint64_t reverse_complement(uint64_t kmer, uint8_t n) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
}

// Code mostly from Jellyfish 1.6 source                                                                                             
__device__ uint64_t reverse_complement(uint64_t kmer) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - ((key_bits/2) << 1));
}


// Lexicographically smallest of k-mer and reverse comp. of k-mer                                                                    
__device__ uint64_t canonical_representation(uint64_t kmer, uint8_t n) {
  uint64_t revcom = reverse_complement(kmer, n);
  return kmer < revcom ? kmer : revcom;
}

__device__ uint64_t canonical_representation(uint64_t kmer) {
  uint64_t revcom = reverse_complement(kmer, key_bits/2);
  return kmer < revcom ? kmer : revcom;
}


//__device__ int get_kmer(int kmer_start, int *read)
__device__ uint64_t bin_key(uint64_t kmer, uint64_t idx_nt) {
  uint8_t nt = idx_nt;
  uint64_t xor_mask = INDEX2_XOR_MASK;
  uint64_t mask = 1 << (nt * 2);
  mask--;
  xor_mask &= mask;
  uint64_t min_bin_key = ~0;
  for (uint64_t i = 0; i < key_bits / 2 - nt + 1; i++) {
    uint64_t temp_bin_key = xor_mask ^ canonical_representation(kmer & mask, nt);
    if (temp_bin_key < min_bin_key)
      min_bin_key = temp_bin_key;
    kmer >>= 2;
  }
  return min_bin_key;
}

// Separate functions to avoid a conditional in the function                                                                         
// This probably isn't necessary...                                                                                                   
__device__ uint64_t bin_key(uint64_t kmer) {
  uint8_t nt = index_ptr->indexed_nt();
  uint8_t idx_type = index_ptr->index_type();
  uint64_t xor_mask = idx_type == 1 ? 0 : INDEX2_XOR_MASK;
  uint64_t mask = 1 << (nt * 2);
  mask--;
  xor_mask &= mask;
  uint64_t min_bin_key = ~0;
  for (uint64_t i = 0; i < key_bits / 2 - nt + 1; i++) {
    uint64_t temp_bin_key = xor_mask ^ canonical_representation(kmer & mask, nt);
    if (temp_bin_key < min_bin_key)
      min_bin_key = temp_bin_key;
    kmer >>= 2;
  }
  return min_bin_key;
}





__global__ void kernel(int *reads, int num_of_reads)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    int read_length = 0;
    int length_bits = 10;
    int ints_per_read = 16;
    int kmer_length = 31;

    if(blockIdx.x =< num_of_reads){
      if(threadIdx.x == 1){
	read_length = ( reads[blockIdx.x*int_per_read] >> sizeof(read[0])*8-length_bits) 
      }
          
      __syncthreads();
    
      __shared__ unsigned int kmers [read_length-kmer_length+1];
      
      if (threadIdx.x =< (read_length - kmer_length)){
	
	kmers[threadIdx.x] = get_kmer(threadIdx.x, reads)
	
	
	
	
      }
      
      
      

      
    }        
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
