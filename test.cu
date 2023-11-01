#include "stdio.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "curand_kernel.h"


#define BLOCKS 32


__global__ void setup_kernel(curandState *state, unsigned long seed)
{
    printf("hhh");
	// int tid = blockIdx.x *blockDim.x + threadIdx.x; //获取线程号0~blocks*THREAD_NUM-1  grid划分成1维，block划分为1维
	// curand_init(seed, tid, 0, &state[tid]);// initialize the state
    
}

__global__ void use(curandState *globalState)
{
    
	// unsigned int j;
	// int tid = blockIdx.x *blockDim.x + threadIdx.x; //获取线程号0~blocks*THREAD_NUM-1  grid划分成1维，block划分为1维
	// curandState localState = globalState[tid];
	// j = (curand(&localState));
	// printf("%u\n", j);
    printf("hhh");
}

int main()
{
	curandState* devStates;  //创建一个随机算法状态的对象
	cudaMalloc(&devStates, BLOCKS * THREAD_NUM * sizeof(curandState));
	srand(time(0));

	setup_kernel << <BLOCKS, THREAD_NUM >> > (devStates, rand()); // blocks number is POP. thread number is P
	use << < BLOCKS, THREAD_NUM >> > (devStates);
    
    return 0;
}

