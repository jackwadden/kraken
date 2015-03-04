__global__ void vector_add(int *a, int *b, int length)
{
    int tx = threadIdx.x;

    if(tx < length) {
        a[tx] = a[tx] + b[tx];
    }
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

    int byteLength = length * sizeof(int);

    cudaMalloc( (void **)&a_d, byteLength );
    cudaMalloc( (void **)&b_d, byteLength );

    cudaMemcpy( a_d, a, byteLength, cudaMemcpyHostToDevice );
    cudaMemcpy( b_d, b, byteLength, cudaMemcpyHostToDevice );

    vector_add<<< blocks, threads >>>( a, b , length);

    cudaMemcpy( a, a_d, byteLength, cudaMemcpyDeviceToHost );
    cudaMemcpy( b, b_d, byteLength, cudaMemcpyDeviceToHost );

    cudaFree(a_d);
    cudaFree(b_d);
}
