/*
 * -- Kraken-GPU -- 
 * Jack Wadden - 2015
 */

#include "stdio.h"
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>

// kraken
#include "seqreader.hpp"
#include "krakendb.hpp"

using namespace std;
using namespace kraken;

__constant__ uint64_t INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess) 
        {
            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
        }
}


// Code mostly from Jellyfish 1.6 source                                      
__device__ uint64_t reverse_complement_gpu(uint64_t kmer, uint8_t n, uint64_t key_bits_d) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (n << 1));
}



// Code mostly from Jellyfish 1.6 source
inline __device__ uint64_t reverse_complement_gpu(uint64_t kmer, uint64_t key_bits_d) {
  kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
  kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
  kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
  kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
  kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
  return (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - ((key_bits_d/2) << 1));
}



// Lexicographically smallest of k-mer and reverse comp. of k-mer
inline __device__ uint64_t canonical_representation_gpu(uint64_t kmer, uint8_t n, uint64_t key_bits_d) {
    uint64_t revcom = reverse_complement_gpu(kmer, n, key_bits_d);
    return kmer < revcom ? kmer : revcom;
}



inline __device__ uint64_t canonical_representation_gpu(uint64_t kmer, uint64_t key_bits_d) {
    uint64_t revcom = reverse_complement_gpu(kmer, key_bits_d/2, key_bits_d);
  return kmer < revcom ? kmer : revcom;
}


/*
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
*/

// Separate functions to avoid a conditional in the function
// This probably isn't necessary...
__device__ uint64_t bin_key_gpu(uint64_t kmer, uint64_t key_bits_d, uint64_t nt) {
    /*
      uint8_t nt = index_ptr->indexed_nt();
      uint8_t idx_type = index_ptr->index_type();
      uint64_t xor_mask = idx_type == 1 ? 0 : INDEX2_XOR_MASK;
    */
    //ONLY SUPPORTS PREMADE DB WITH V2 IDX FOR NOW

    uint64_t xor_mask = INDEX2_XOR_MASK;
    uint64_t mask = 1 << ((uint8_t)nt * 2);
    mask--;
    xor_mask &= mask;
    uint64_t min_bin_key = ~0;
    for (uint64_t i = 0; i < key_bits_d / 2 - nt + 1; i++) {
        uint64_t temp_bin_key = xor_mask ^ canonical_representation_gpu(kmer & mask, nt, key_bits_d);
        if (temp_bin_key < min_bin_key)
            min_bin_key = temp_bin_key;
        kmer >>= 2;
    }

    return min_bin_key;
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
                case 'A' : // A
                    tmp = tmp | 0;
                    break;
                case 'C' : // T
                    tmp = tmp | 1;
                    break;
                case 'G' : // G
                    tmp = tmp | 2;
                    break;
                case 'T' : // C
                    tmp = tmp | 3;
                    break;
                }
                //cout << hex << read[i*4 + j] << endl;
            }
        }
        // store every 4 bp into a single char array after length
        //cout << hex << (int)tmp << endl;
        h_reads[(bytesPerRead * readNum) + 1 + i] = tmp;
    }
    //cout << endl;
}


__global__ void kernel(uint64_t * idx,
                       char * kdb,
                       uint32_t pair_sz,
                       uint32_t key_bits,
                       uint64_t key_len,
                       char *reads, 
                       int num_reads, 
                       int read_bytes,
                       int kmer_len,
                       uint64_t * output,
                       uint64_t key_bits_gpu,
                       uint64_t nt)
{

    int block_id = blockIdx.x + blockIdx.y * gridDim.x
        + gridDim.x * gridDim.y * blockIdx.z;
    int thread_id = block_id * blockDim.x + threadIdx.x;
    int work_size = 1;
    
    // FOR EVERY READ
    //for(int read_num = work_size * block_id;
    int read_num = work_size * block_id;

    //ACTIVATE THREADS
    int read_len;
    read_len = reads[read_num * read_bytes];
    if(threadIdx.x < read_len - kmer_len) {
        
        //EACH THREAD GETS A KMER FROM THE READ
        uint64_t kmer = 0;
        int cur_bp = 0;
        int limit_bp = threadIdx.x + kmer_len;
        for(int i = 0; i < read_bytes - 1; i++) {
            //char window = rd[i];
            char window = reads[read_num * read_bytes + 1 + i] ;
            for(int j = 6; j >= 0; j-=2) {
                if(cur_bp < limit_bp) {
                    // if there is a bp left, add it
                    if(cur_bp < read_len) {
                        // shift in 00
                        kmer <<= 2;
                        // grab two bits of window
                        kmer |= ((window >> j) & 3);
                    }
                }
                cur_bp++;
            }
        }
        uint64_t kmer_mask = ~0;
        kmer_mask >>= 2;
        kmer &= kmer_mask;
        
        //EACH THREAD CALCULATES A MINIMIZER FROM THE KMER
        uint64_t canon_kmer = canonical_representation_gpu(kmer, key_bits_gpu);
        uint64_t minimizer = bin_key_gpu(canon_kmer, key_bits_gpu, nt);
        //output[thread_id] = kmer;      
        //EACH THREAD GETS A HIGH AND LOW SEARCH SPACE
        int min = idx[minimizer];
        int max = idx[minimizer + 1] - 1;
        
        //EACH THREAD LOOKS FOR ITS TAXON VIA BINARY SEARCH
        // Binary search with large window
        uint64_t comp_kmer;
        int mid;
        uint32_t * taxon_ptr = 0;
        
        int found = 0;
        while (min + 15 <= max) {
            mid = min + (max - min) / 2;
            comp_kmer = 0x0000000000000000ULL;
            memcpy(&comp_kmer, kdb + pair_sz * mid, key_len);
            comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
            if (canon_kmer > comp_kmer)
                min = mid + 1;
            else if (canon_kmer < comp_kmer)
                max = mid - 1;
            else{
                taxon_ptr = (uint32_t *) (kdb + pair_sz * mid + key_len);
                found = 1;
                break;
            }        
        }
        
        // Linear search once window shrinks
        if(!found) {
            for (mid = min; mid <= max; mid++) {
                comp_kmer = 0x0000000000000000ULL;
                memcpy(&comp_kmer, kdb + pair_sz * mid, key_len);
                comp_kmer &= (1ull << key_bits) - 1;  // trim any excess
                if (canon_kmer == comp_kmer) {
                    taxon_ptr = (uint32_t *) (kdb + pair_sz * mid + key_len);
                }
            }
        }
               
        //output[thread_id * 2] = kmer;
        //output[thread_id * 2 + 1] = taxon_ptr ? *taxon_ptr : 0;
       
        //EACH THREAD UPDATES THE TAXON HISTOGRAM
        /*
          if(thread_id == 0) {
          for(int i = 0; i < 1000; i++) {
          output = kdb + pair_sz * i + key_len;
          }
          }
        */
    }
}



/*
 * Kernel wrapper is called and given a ptr to the database file
 *  Database is passed to the GPU in the following format
 *  |read len (8)| read (512) |
 *
 */
void kernel_wrapper(KrakenDB * database,
                    int maxReadLen, 
                    vector<DNASequence> reads, 
                    uint64_t key_bits_gpu, 
                    uint64_t nt)
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
        /*
        for(int j = 0; j < 63; j++)
            cout << j << ": " << dec <<  h_reads[i*bytesPerRead + j] << endl;        
        */
    }
    cout << hex << endl;
    
    printf("num reads: %d\n", numReads);

    // CUDA SETUP
    int blockSize = 256;
    int numBlocks = (int)(numReads);
    /*
    if(numReads % blockSize) {
        numBlocks++;
    }
    */

    dim3 threads( blockSize, 1, 1 );
    dim3 blocks( numBlocks, 1, 1 );

    // HOST ALLOCATE
    size_t byteLength = sizeof(char) * 2000000;
    uint64_t * h_output;
    gpuErrchk(cudaMallocHost( (void **)&h_output, byteLength ));

    // DEVICE ALLOCATE
    cudaEvent_t alloc_start, alloc_stop;
    cudaEventCreate(&alloc_start);
    cudaEventCreate(&alloc_stop);

    char * d_reads;
    uint64_t * d_output;
    uint64_t * d_idx;
    char * d_kdb;

    cudaEventRecord(alloc_start);
    gpuErrchk(cudaMalloc( (void **)&d_reads, numReads * bytesPerRead) );
    gpuErrchk(cudaMalloc( (void **)&d_output, byteLength) );
    // allocate for index array
    size_t idx_size = 536870928; //unclear how to get this number from the object yet
    gpuErrchk(cudaMalloc( (void **)&d_idx, idx_size) );
    // allocate for db
    size_t kdb_size = database->pair_size() * database->get_key_ct();
    gpuErrchk(cudaMalloc( (void **)&d_kdb, kdb_size) );

    cudaEventRecord(alloc_stop);
    cudaEventSynchronize(alloc_stop);
    float alloc_ms = 0;
    cudaEventElapsedTime(&alloc_ms, alloc_start, alloc_stop);
    cout<< "ALLOCATION TIME: " << alloc_ms << " ms" << endl;

    // MEM COPY TO
    cudaEvent_t copyto_start, copyto_stop;
    cudaEventCreate(&copyto_start);
    cudaEventCreate(&copyto_stop);

    cudaEventRecord(copyto_start);
    
    //
    gpuErrchk(cudaMemcpy(d_reads, h_reads, numReads * bytesPerRead, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_idx, database->get_index()->get_array(), idx_size, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_kdb, database->get_pair_ptr(), kdb_size, cudaMemcpyHostToDevice));

    
    cudaEventRecord(copyto_stop);
    cudaEventSynchronize(copyto_stop);
    float copyto_ms = 0;
    cudaEventElapsedTime(&copyto_ms, copyto_start, copyto_stop);
    cout<< "COPY TO TIME: " << copyto_ms << " ms" << endl;

    // KERNEL LAUNCH
    cudaEvent_t kernel_start, kernel_stop;
    cudaEventCreate(&kernel_start);
    cudaEventCreate(&kernel_stop);

    cout << "LAUNCHING KERNEL :: blocks(" 
         << dec 
         << blocks.x << "," 
         << blocks.y << "," 
         << blocks.z << ")" 
         << " threads("        
         << threads.x << "," 
         << threads.y << "," 
         << threads.z << ")" 
         << endl;

    cudaEventRecord(kernel_start);
    //
    kernel<<< blocks, threads>>>(d_idx,
                                                d_kdb,
                                                database->pair_size(),
                                                database->get_key_bits(),
                                                database->get_key_len(),
                                                d_reads,
                                                numReads,
                                                bytesPerRead,
                                                31,
                                                d_output,
                                                key_bits_gpu,
                                                nt);
    gpuErrchk( cudaPeekAtLastError() );

    cudaEventRecord(kernel_stop);
    cudaEventSynchronize(kernel_stop);
    float kernel_ms = 0;
    cudaEventElapsedTime(&kernel_ms, kernel_start, kernel_stop);
    cout<< "KERNEL TIME: " << kernel_ms << " ms" << endl;


    // MEM COPY FROM
    cudaEvent_t copyfrom_start, copyfrom_stop;
    cudaEventCreate(&copyfrom_start);
    cudaEventCreate(&copyfrom_stop);

    cudaEventRecord(copyfrom_start);

    //
    //gpuErrchk( cudaMemcpy( h_output, d_output, byteLength, cudaMemcpyDeviceToHost ));

    cudaEventRecord(copyfrom_stop);
    cudaEventSynchronize(copyfrom_stop);
    float copyfrom_ms = 0;
    cudaEventElapsedTime(&copyfrom_ms, copyfrom_start, copyfrom_stop);
    cout<< "COPY FROM TIME: " << copyfrom_ms << " ms" << endl;

    // PRINT OUTPUT
    for(int i = 0; i < 1024; i+=2) {
        //cout << "#" << dec << i/2 << ":gpu::kmer: " << hex << h_output[i] << " taxon: " << h_output[i+1] << endl;
    }
    /*
      bytes
    for(int i = 0; i < 500; i++) {
        printf("i:%d %02X\n",i, (int)h_output[i]);
        //cout << "i:" << i << setw(2) << setfill('0') << hex << (int) h_output[i] << endl ;

    }
    */
    // CLEANUP
    gpuErrchk( cudaFree(d_reads) );
    gpuErrchk( cudaFree(d_output) );    
    gpuErrchk( cudaFree(d_idx) );    
    gpuErrchk( cudaFree(d_kdb) );    

    gpuErrchk( cudaFreeHost(h_reads) );
    gpuErrchk( cudaFreeHost(h_output) );
}
