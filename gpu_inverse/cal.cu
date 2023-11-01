#include<iostream>
#include<thread>
#include <cuda.h>
#include <curand.h>
//输入a，b

//求解\theta \phi

//矩阵求逆*初始变量（不要完全写出矩阵，无法储存）

//矩阵减去已经存在的变量

//对比前后的变化量，差距不大就输出


const int pre=16384;//矩阵大小
__constant__ double theta_gpu[pre+1];
__constant__ double phi_gpu[pre+1];
__constant__ double powb_gpu[pre];
double*partial_c_gpu;//用于累加暂时储存部分和
double*partial_c_cpu;
double*tmp_gpu;//用于储存向量中间计算结果
double*tmp_cpu;
int n=20;


const int threadsPerBlock = 256;
const int blocksPerGrid = min(32, (pre + threadsPerBlock - 1) / threadsPerBlock);


double b;
double a[pre];

void get_theta(){
    double theta_cpu[pre+1];
    theta_cpu[0]=1;
    theta_cpu[1]=a[1-1];
    double b2=b*b;
    for (int i = 2; i <= pre; i++)
    {
        theta_cpu[i]=a[i-1]*theta_cpu[i-1]-b2*theta_cpu[i-2];
    }
    cudaMemcpyToSymbol(theta_gpu, theta_cpu, pre * sizeof(double));
}

void get_phi(){
    double phi_cpu[pre+1];
    phi_cpu[pre]=1;
    phi_cpu[pre-1]=a[pre-1];
    double b2=b*b;
    for (int i = pre-2; i >= 0; i--)
    {
        phi_cpu[i]=a[i-1]*phi_cpu[i+1]-b2*phi_cpu[i+2];
    }
    cudaMemcpyToSymbol(phi_gpu, phi_cpu, pre * sizeof(double));
}

void get_powb(){
    double powb_cpu[pre];
    powb_cpu[0]=1;
    for (int i = 1; i < pre; i++)
    {powb_cpu[i]=b*powb_cpu[i-1];}
    cudaMemcpyToSymbol(powb_gpu, powb_cpu, pre * sizeof(double));
}




__global__ void mult_matrix_and_vert_gpu(double* in, double* out){
    int offset=blockIdx.x * blockDim.x + threadIdx.x;
    if(offset>=pre){
        return;
    }
    out[offset]=0;
    double tmp;
    for(int i=0;i<pre;i++){
        if(offset>i){
            tmp=powb_gpu[offset-i]*theta_gpu[i]*phi_gpu[offset]/theta_gpu[pre];
        }
        else{
            tmp=powb_gpu[offset-i]*theta_gpu[offset]*phi_gpu[i]/theta_gpu[pre];
        }
        out[offset]+=tmp*in[i];
    }
}
void mult_matrix_and_vert(double* in, double* out){
    mult_matrix_and_vert_gpu<<<blocksPerGrid,threadsPerBlock>>>(in,out);
}




__global__ void dot_gpu(double* a, double* b, double* c) {
    __shared__ double cache[threadsPerBlock];
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int cacheIndex = threadIdx.x;
    double temp = 0;
    while (tid < pre) {
        temp += a[tid] * b[tid];
        tid += blockDim.x * gridDim.x;
    }
    cache[cacheIndex] = temp;
    __syncthreads();
    int i = blockDim.x / 2;
    while (i != 0) {
        if (cacheIndex < i)
            cache[cacheIndex] += cache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }
    if (cacheIndex == 0)
        c[blockIdx.x] = cache[0];
}
double dot(double*a,double*b){
    dot_gpu<<<blocksPerGrid, threadsPerBlock>>>(a, b, partial_c_gpu);
    cudaMemcpy(partial_c_cpu, partial_c_gpu, blocksPerGrid * sizeof(double), cudaMemcpyDeviceToHost);
    double c=0;
    for (int i = 0; i < blocksPerGrid; i++){c += partial_c_cpu[i];}
    return c;
}


__global__ void sub_gpu(double*a,double*b,double*c){
    int offset=blockIdx.x * blockDim.x + threadIdx.x;
    if(offset>=pre){
        return;
    }
    a[offset]=b[offset]-c[offset];
}
void sub(double*a,double*b,double*c){
    sub_gpu<<<blocksPerGrid,threadsPerBlock>>>(a,b,c);
}
void sub(double*a,double*b){
    sub(a,a,b);
}



__global__ void multinum_gpu(double*a,double*b,double c){
    int offset=blockIdx.x * blockDim.x + threadIdx.x;
    if(offset>=pre){
        return;
    }
    a[offset]=b[offset]*c;
}
void multinum(double*a,double*b,double c){
    multinum_gpu<<<blocksPerGrid,threadsPerBlock>>>(a,b,c);
}
void multinum(double*a,double c){
    multinum_gpu<<<blocksPerGrid,threadsPerBlock>>>(a,a,c);
}


void cal_one_round(double* cal1,double* cal2,double* memery,int memery_n,bool& turn){
    double *in;
    double *out;
    if(turn){
        in=cal1;
        out=cal2;
    }
    else{
        in=cal2;
        out=cal1;
    }
    mult_matrix_and_vert(in,out);
    for(int i=0;i<memery_n;i++){
        double coff=dot(out,(memery+i*pre));
        multinum(tmp_gpu,(memery+i*pre),coff);
        sub(out,(memery+i*pre));
    }
    double normal=1/sqrt(dot(out,out));
    multinum(out,normal);
    turn=!turn;
}


__global__ void copyvec_gpu(double* a,double* b){
    int offset=blockIdx.x * blockDim.x + threadIdx.x;
    if(offset>=pre){
        return;
    }
    a[offset]=b[offset];
}
void copyvec(double* a,double* b){
    copyvec_gpu<<<blocksPerGrid,threadsPerBlock>>>(a,b);
}

__global__ void mycurand_gpu(double* a){
    int offset=blockIdx.x * blockDim.x + threadIdx.x;
    if(offset>=pre){
        return;
    }
    a[offset]=;
}
void mycurand(double* a){
    mycurand_gpu<<<blocksPerGrid,threadsPerBlock>>>(a);
}

void run(){
    std::ios::sync_with_stdio(false);
    std::thread thread_a([](){get_theta();});
    std::thread thread_b([](){get_phi();});
    std::thread thread_c([](){get_powb();});
    thread_a.detach();
    thread_b.detach();
    thread_c.detach();

    double* cal1;
    double* cal2;
    double* memery;
    cudaMalloc((void**)&cal1, pre * sizeof(double));
    cudaMalloc((void**)&cal2, pre * sizeof(double));
    cudaMalloc((void**)&memery, (pre*n) * sizeof(double));
    cudaMalloc((void**)&partial_c_gpu, blocksPerGrid * sizeof(double));
    partial_c_cpu = (double*)malloc(pre * sizeof(double));
    cudaMalloc((void**)&tmp_gpu, pre * sizeof(double));
    bool turn=true;
    for (int i = 0; i < n; i++)
    {
        //随机一下
        while (true)
        {

            mycurand(cal1);
            mycurand(cal2);
            cal_one_round(cal1,cal2,memery,i,turn);
            sub_gpu(tmp_gpu,cal1,cal2);
            dot(tmp_gpu,tmp_gpu);
            if(){
                break;
            }
        }
        copyvec(memery,cal1);
    }
    
}