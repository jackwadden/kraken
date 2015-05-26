/*
 * -- Kraken-GPU -- 
 * Jack Wadden - 2015
 */

#include "stdio.h"
#include <stdlib.h>
#include <stdint.h>
#include <iomanip>

// kraken
#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include "seqreader.hpp"

using namespace std;
using namespace kraken;

__constant__ uint64_t INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;

extern string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig);
extern bool Print_classified;

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
                       uint64_t key_bits_gpu,
                       uint64_t key_len,
                       char *reads, 
                       int num_reads, 
                       int read_bytes,
                       int kmer_len,
                       uint32_t * output,
                       uint64_t nt)
{

    //SHARED DECLS
    __shared__ uint32_t taxons[256];
    __shared__ uint32_t hist[512];
    __shared__ uint32_t hit_counter;
    hit_counter = 0;
    int found = 0;
    uint32_t * taxon_ptr = 0;    

    // ID INITS
    int block_id = blockIdx.x + blockIdx.y * gridDim.x
        + gridDim.x * gridDim.y * blockIdx.z;
    int thread_id = block_id * blockDim.x + threadIdx.x;
    int work_size = 1;
    
    // FOR EVERY READ
    //for(int read_num = work_size * block_id;
    int read_num = work_size * block_id;


    int read_len;
    read_len = reads[read_num * read_bytes];
    output[threadIdx.x * 2] = 0;
    output[threadIdx.x * 2 + 1] = 0;
    taxons[threadIdx.x] = 0;

    //ACTIVATE THREADS
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
        
        while (min + 15 <= max) {
            mid = min + (max - min) / 2;
            comp_kmer = 0x0000000000000000ULL;
            memcpy(&comp_kmer, kdb + pair_sz * mid, key_len);
            comp_kmer &= (1ull << key_bits_gpu) - 1;  // trim any excess
            if (canon_kmer > comp_kmer)
                min = mid + 1;
            else if (canon_kmer < comp_kmer)
                max = mid - 1;
            else{
                taxon_ptr = (uint32_t *) (kdb + pair_sz * mid + key_len);
                taxons[threadIdx.x] = *taxon_ptr;
                atomicAdd(&hit_counter, 1);
                found = 1;
                break;
            }        
        }
        
        // Linear search once window shrinks
        if(!found) {
            for (mid = min; mid <= max; mid++) {
                comp_kmer = 0x0000000000000000ULL;
                memcpy(&comp_kmer, kdb + pair_sz * mid, key_len);
                comp_kmer &= (1ull << key_bits_gpu) - 1;  // trim any excess
                if (canon_kmer == comp_kmer) {
                    taxon_ptr = (uint32_t *) (kdb + pair_sz * mid + key_len);
                    taxons[threadIdx.x] = *taxon_ptr;
                    atomicAdd(&hit_counter, 1);
                }
            }
        }
    }
    
    __threadfence();

    /*
      if(threadIdx.x == 0){
      output[block_id * 2] = kmer;
      output[block_id * 2 + 1] = hit_count;
      }
    */

    //LEAD THREAD POPULATES HISTOGRAM (NOTE: NAIVE O(N^2) IMPL.)
    uint32_t uniq_counter = 0;
    uint32_t total_counter = 0;
    int hit = 0;
    if(threadIdx.x == 0) {
        // for all threads
        int thread_index = 0;
        while(total_counter < hit_counter && thread_index < blockDim.x) {
            //if the thread found a taxon
            uint32_t taxon = taxons[thread_index];
            if(taxon != 0) {
                total_counter++;
                //search in stack for taxon
                hit = 0;
                for(int j = 0; j < uniq_counter; j++){
                    if(hist[j * 2] == taxon){
                        hist[j * 2 + 1]++;
                        hit = 1;
                    }
                }
                if(!hit) {
                    hist[uniq_counter * 2] = taxon;
                    hist[uniq_counter * 2 + 1] = 1;
                    //inc top of stack
                    uniq_counter++;
                }
            }
            thread_index++;
        }
    }

    // THREADS COPY HISTOGRAM TO GLOBAL FOR EXPORT
    output[block_id * blockDim.x * 2 + threadIdx.x * 2] = hist[threadIdx.x * 2];
    output[block_id * blockDim.x * 2 + threadIdx.x * 2 + 1] = hist[threadIdx.x * 2 + 1];


    //EACH THREAD SUMS TAXON'S LTR PATH

    //LEAD THREAD FINDS MAX VIA REDUCTION

}

/*
 *  Read input is passed to the GPU in the following format
 *  |read len (8)| read (512) |
 *
 */
void process_file_gpu(char *filename, 
                      KrakenDB *database, 
                      map<uint32_t,uint32_t> &parent_map,
                      bool Fastq_input, 
                      size_t work_unit_size,
                      ostream *Kraken_output,
                      ostream *Classified_output,
                      ostream *Unclassified_output) {
    string file_str(filename);
    DNASequenceReader *reader;
    DNASequence dna;
    
    if (Fastq_input)
        reader = new FastqReader(file_str);
    else
        reader = new FastaReader(file_str);
    
    vector<DNASequence> reads;
    ostringstream koss, coss, uoss;
    bool Print_unclassified = false;
    //bool Print_classified = true;
    bool Print_kraken = false;
    bool Only_classified_kraken_output = true;
    bool Quick_mode = false;
    uint32_t Minimum_hit_count = 3;

    //INITIALIZE GPU
    uint32_t max_read_len = 251;
    // HOST ALLOCATE
    size_t byteLength = sizeof(uint32_t) * 2 * 256 * 2000000;
    uint32_t * h_output;
    gpuErrchk(cudaMallocHost( (void **)&h_output, byteLength ));

    // DEVICE ALLOCATE
    cudaEvent_t alloc_start, alloc_stop;
    cudaEventCreate(&alloc_start);
    cudaEventCreate(&alloc_stop);

    char *h_reads;
    char * d_reads;
    uint32_t * d_output;
    uint64_t * d_idx;
    char * d_kdb;

    cudaEventRecord(alloc_start);
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
    gpuErrchk(cudaMemcpy(d_idx, database->get_index()->get_array(), idx_size, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_kdb, database->get_pair_ptr(), kdb_size, cudaMemcpyHostToDevice));
    
    
    cudaEventRecord(copyto_stop);
    cudaEventSynchronize(copyto_stop);
    float copyto_ms = 0;
    cudaEventElapsedTime(&copyto_ms, copyto_start, copyto_stop);
    cout<< "COPY TO TIME: " << copyto_ms << " ms" << endl;

    //KERNEL LOOP
    while (reader->is_valid()) {
        reads.clear();
        size_t total_nt = 0;
        while (total_nt < work_unit_size) {
            dna = reader->next_sequence();
            if (! reader->is_valid())
                break;
            reads.push_back(dna);
            total_nt += dna.seq.size();
        }

        cout << "TOTAL NT: " << total_nt << endl;
        
        // Convert vector of strings to GPU appropriate data structure
        int numReads = reads.size();
        const unsigned int bytesPerRead =  
            sizeof(char) + // spot for read len (8 bits for now) 
            ceil((max_read_len * 2) / 8); // spot for read
        
        // Malloc appropriate size pinned array on host
        gpuErrchk(cudaMallocHost( (void **)&h_reads, numReads * bytesPerRead));
        gpuErrchk(cudaMalloc( (void **)&d_reads, numReads * bytesPerRead) );
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
        
        printf("PROCESSING #READS: %d\n", numReads);
        
        // CUDA SETUP
        int blockSize = 256;
        int numBlocks = (int)(numReads);
        
        dim3 threads( blockSize, 1, 1 );
        dim3 blocks( numBlocks, 1, 1 );
        
        //COPY WORK UNIT
        cudaEvent_t wucopyto_start, wucopyto_stop;
        cudaEventCreate(&wucopyto_start);
        cudaEventCreate(&wucopyto_stop);
        cudaEventRecord(wucopyto_start);
        //
        gpuErrchk(cudaMemcpy(d_reads, h_reads, numReads * bytesPerRead, cudaMemcpyHostToDevice));
        cudaEventRecord(wucopyto_stop);
        cudaEventSynchronize(wucopyto_stop);
        float wucopyto_ms = 0;
        cudaEventElapsedTime(&wucopyto_ms, wucopyto_start, wucopyto_stop);
        cout<< "WORK UNIT COPY TO TIME: " << wucopyto_ms << " ms" << endl;
        
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
                                     database->get_index()->indexed_nt());
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
        gpuErrchk( cudaMemcpy( h_output, d_output, byteLength, cudaMemcpyDeviceToHost ));
        
        cudaEventRecord(copyfrom_stop);
        cudaEventSynchronize(copyfrom_stop);
        float copyfrom_ms = 0;
        cudaEventElapsedTime(&copyfrom_ms, copyfrom_start, copyfrom_stop);
        cout<< "COPY FROM TIME: " << copyfrom_ms << " ms" << endl;
        
        /*
        for(int i = 0; i < 1024; i+=2) {
            cout << "#" << dec << i/2 << ":gpu::taxon: " << dec << h_output[i] << " hit_count: " << h_output[i+1] << endl;
        }
        */
        
        koss.str("");
        coss.str("");
        uoss.str("");    
        
        // RESOLVE TREE FOR EACH READ
        cout << "RESOLVING TREES FOR EACH READ..." << endl;
        uint32_t total_classified = 0;
        // TODO: Parameterize this
        size_t read_hist_size = 2 * 256; // two locations for each thread
        vector<uint32_t> taxa;
        vector<uint8_t> ambig_list;            
        map<uint32_t, uint32_t> hit_counts;

        for(int read = 0; read < numReads; read++) {
            //cout << read << endl;
            
            // create hit counts map
            uint32_t taxon = h_output[read * read_hist_size];
            uint32_t taxon_counter = 0;
            
            while(taxon != 0){
                taxa.push_back(taxon);
                ambig_list.push_back(0);
                hit_counts[taxon] = h_output[read * read_hist_size + taxon_counter * 2 + 1];
                taxon_counter++;
                taxon = h_output[read * read_hist_size + taxon_counter * 2];
            }

            // PRINT CLASSIFICATION
            uint32_t call = 0;
            uint32_t hits = 0;
            
            // quick mode not implemented now
            //cout << "starting resolve tree" << endl;
            if (Quick_mode)
                call = hits >= Minimum_hit_count ? taxon : 0;
            else
                call = resolve_tree(hit_counts, parent_map);

            //cout << "TAXA: " << taxa.size() << endl;
            //cout << "HITS: " << hit_counts.size() << endl;

            hit_counts.clear();
            ambig_list.clear();
            taxa.clear();


            //cout << "ending resolve tree" << endl;
            //cout << "CLASSIFICATION: " << call << endl;
            
            if (call)
#pragma omp atomic
                total_classified++;
            
            if (Print_unclassified || Print_classified) {
                ostringstream *oss_ptr = call ? &coss : &uoss;
                bool print = call ? Print_classified : Print_unclassified;
                if (print) {
                    if (Fastq_input) {
                        (*oss_ptr) << "@" << reads[read].header_line << endl
                                   << reads[read].seq << endl
                                   << "+" << endl
                                   << reads[read].quals << endl;
                    }
                    else {
                        (*oss_ptr) << ">" << reads[read].header_line << endl
                                   << reads[read].seq << endl;
                    }
                }
            }
            
            if (! Print_kraken)
                continue;
            
            if (call) {
                koss << "C\t";
            }
            else {
                if (Only_classified_kraken_output)
                    continue;
                koss << "U\t";
            }
            koss << reads[read].id << "\t" << call << "\t" << reads[read].seq.size() << "\t";
            
            if (Quick_mode) {
                koss << "Q:" << hits;
            }
            else {
                if (taxa.empty())
                    koss << "0:0";
                else
                    koss << hitlist_string(taxa, ambig_list);
            }
            
            koss << endl;
        
        }
        
        uint32_t total_sequences = 0;
        uint32_t total_bases = 0;
#pragma omp critical(write_output)
        {
            if (Print_kraken)
                (*Kraken_output) << koss.str();
            if (Print_classified)
                (*Classified_output) << coss.str();
            if (Print_unclassified)
                (*Unclassified_output) << uoss.str();
            total_sequences += work_unit_size;
            total_bases += total_nt;
            cerr << "\rProcessed " << total_sequences << " sequences (" << total_bases << " bp) ...";
            
        }
        
        
        if (total_nt == 0)
            break;
        
    }
            
            

    // PRINT OUTPUT
    /*
      for(int i = 0; i < 256; i+=2) {
      cout << "#" << dec << i/2 << ":gpu::kmer: " << hex << h_output[i] << " hit_count: " << h_output[i+1] << endl;
      }
    */
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
  
    delete reader;
}
