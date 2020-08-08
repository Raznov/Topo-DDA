#include "definition.h"

__global__ void A2AsKernel(double *A, cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, int NxFFT, int NyFFT, int NzFFT){
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int index = k + j*NzFFT + i*NzFFT*NyFFT;
    
    if(i < NxFFT && j <NyFFT && k<NzFFT){
        A00[index].x = A[0+2*0+12*index];
        A00[index].y = A[1+2*0+12*index];

        A01[index].x = A[0+2*1+12*index];
        A01[index].y = A[1+2*1+12*index];

        A02[index].x = A[0+2*2+12*index];
        A02[index].y = A[1+2*2+12*index];

        A11[index].x = A[0+2*3+12*index];
        A11[index].y = A[1+2*3+12*index];

        A12[index].x = A[0+2*4+12*index];
        A12[index].y = A[1+2*4+12*index];

        A22[index].x = A[0+2*5+12*index];
        A22[index].y = A[1+2*5+12*index];
    }
    

}

void A2As(double *A, cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, int NxFFT, int NyFFT, int NzFFT){

    int tmpx = 10;
    int tmpy = 10;
    int tmpz = 10;
    dim3 dimBlock(tmpx, tmpy, tmpz);
    dim3 dimGrid(ceil((double)NxFFT/tmpx), ceil((double)NyFFT/tmpy), ceil((double)NzFFT/tmpz));   

    A2AsKernel<<<dimGrid, dimBlock>>>(A, A00, A01, A02, A11, A12, A22, NxFFT, NyFFT, NzFFT);

}

__global__ void B2BsKernel(double *bDev, cufftDoubleComplex *bxDev, cufftDoubleComplex *byDev, cufftDoubleComplex *bzDev, int NxFFT, int NyFFT, int NzFFT){
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int index = k + j*NzFFT + i*NzFFT*NyFFT;
    
    if(i < NxFFT && j <NyFFT && k<NzFFT){
        bxDev[index].x = bDev[0+2*0+6*index];
        bxDev[index].y = bDev[1+2*0+6*index];

        byDev[index].x = bDev[0+2*1+6*index];
        byDev[index].y = bDev[1+2*1+6*index];

        bzDev[index].x = bDev[0+2*2+6*index];
        bzDev[index].y = bDev[1+2*2+6*index];
    }
    

}

void B2Bs(double *bDev, cufftDoubleComplex *bxDev, cufftDoubleComplex *byDev, cufftDoubleComplex *bzDev, int NxFFT, int NyFFT, int NzFFT){
    int tmpx = 10;
    int tmpy = 10;
    int tmpz = 10;
    dim3 dimBlock(tmpx, tmpy, tmpz);
    dim3 dimGrid(ceil((double)NxFFT/tmpx), ceil((double)NyFFT/tmpy), ceil((double)NzFFT/tmpz));   

    B2BsKernel<<<dimGrid, dimBlock>>>(bDev, bxDev, byDev, bzDev, NxFFT, NyFFT, NzFFT);
}

__global__ void ConvKernel(cufftDoubleComplex *Convx, cufftDoubleComplex *Convy, cufftDoubleComplex *Convz, cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, cufftDoubleComplex *bxDev, cufftDoubleComplex *byDev, cufftDoubleComplex *bzDev, int NxFFT, int NyFFT, int NzFFT){
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int index = k + j*NzFFT + i*NzFFT*NyFFT;
    
    if(i < NxFFT && j <NyFFT && k<NzFFT){
        Convx[index] = cuCadd(cuCadd(cuCmul(A00[index], bxDev[index]), cuCmul(A01[index], byDev[index])), cuCmul(A02[index], bzDev[index]));
        Convy[index] = cuCadd(cuCadd(cuCmul(A01[index], bxDev[index]), cuCmul(A11[index], byDev[index])), cuCmul(A12[index], bzDev[index]));
        Convz[index] = cuCadd(cuCadd(cuCmul(A02[index], bxDev[index]), cuCmul(A12[index], byDev[index])), cuCmul(A22[index], bzDev[index]));
    }
}
    
void Conv(cufftDoubleComplex *Convx, cufftDoubleComplex *Convy, cufftDoubleComplex *Convz, cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, cufftDoubleComplex *bxDev, cufftDoubleComplex *byDev, cufftDoubleComplex *bzDev, int NxFFT, int NyFFT, int NzFFT){
    int tmpx = 10;
    int tmpy = 10;
    int tmpz = 10;
    dim3 dimBlock(tmpx, tmpy, tmpz);
    dim3 dimGrid(ceil((double)NxFFT/tmpx), ceil((double)NyFFT/tmpy), ceil((double)NzFFT/tmpz));   

    ConvKernel<<<dimGrid, dimBlock>>>(Convx, Convy, Convz, A00, A01, A02, A11, A12, A22, bxDev, byDev, bzDev, NxFFT, NyFFT, NzFFT);
}

__global__ void Conv2BKernel(cufftDoubleComplex *Convx, cufftDoubleComplex *Convy, cufftDoubleComplex *Convz, double *bDev, int NxFFT, int NyFFT, int NzFFT){
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int index = k + j*NzFFT + i*NzFFT*NyFFT;
    
    if(i < NxFFT && j <NyFFT && k<NzFFT){
        bDev[0+2*0+6*index] = Convx[index].x;
        bDev[1+2*0+6*index] = Convx[index].y;

        bDev[0+2*1+6*index] = Convy[index].x;
        bDev[1+2*1+6*index] = Convy[index].y;

        bDev[0+2*2+6*index] = Convz[index].x;
        bDev[1+2*2+6*index] = Convz[index].y;
    }
    

}

void Conv2B(cufftDoubleComplex *Convx, cufftDoubleComplex *Convy, cufftDoubleComplex *Convz, double *bDev, int NxFFT, int NyFFT, int NzFFT){
    int tmpx = 10;
    int tmpy = 10;
    int tmpz = 10;
    dim3 dimBlock(tmpx, tmpy, tmpz);
    dim3 dimGrid(ceil((double)NxFFT/tmpx), ceil((double)NyFFT/tmpy), ceil((double)NzFFT/tmpz));   

    Conv2BKernel<<<dimGrid, dimBlock>>>(Convx, Convy, Convz, bDev, NxFFT, NyFFT, NzFFT);
}

__global__ void APtoESumKernel(cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, 
    cufftDoubleComplex *PxDev, cufftDoubleComplex *PyDev, cufftDoubleComplex *PzDev,
    cufftDoubleComplex *ESumxDev, cufftDoubleComplex *ESumyDev, cufftDoubleComplex *ESumzDev, 
    int NxFFT, int NyFFT, int NzFFT,
    int NxA, int NyA, int NzA, 
    int index1, int index2, int index3, int deduction){
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    int j = threadIdx.y + blockIdx.y*blockDim.y;
    int k = threadIdx.z + blockIdx.z*blockDim.z;
    int index = k + j*NzFFT + i*NzFFT*NyFFT;
    int indexA = (k+index3) + (j+index2)*NzA +(i+index1)*NzA*NyA;
    if(i < NxFFT && j <NyFFT && k<NzFFT){
        if(deduction == 1){
            ESumxDev[index] = cuCadd(ESumxDev[index], cuCadd(cuCadd(cuCmul(A00[indexA], PxDev[0]), cuCmul(A01[indexA], PyDev[0])), cuCmul(A02[indexA], PzDev[0])));
            ESumyDev[index] = cuCadd(ESumyDev[index], cuCadd(cuCadd(cuCmul(A01[indexA], PxDev[0]), cuCmul(A11[indexA], PyDev[0])), cuCmul(A12[indexA], PzDev[0])));
            ESumzDev[index] = cuCadd(ESumzDev[index], cuCadd(cuCadd(cuCmul(A02[indexA], PxDev[0]), cuCmul(A12[indexA], PyDev[0])), cuCmul(A22[indexA], PzDev[0])));
        }
        if(deduction == -1){
            ESumxDev[index] = cuCsub(ESumxDev[index], cuCadd(cuCadd(cuCmul(A00[indexA], PxDev[0]), cuCmul(A01[indexA], PyDev[0])), cuCmul(A02[indexA], PzDev[0])));
            ESumyDev[index] = cuCsub(ESumyDev[index], cuCadd(cuCadd(cuCmul(A01[indexA], PxDev[0]), cuCmul(A11[indexA], PyDev[0])), cuCmul(A12[indexA], PzDev[0])));
            ESumzDev[index] = cuCsub(ESumzDev[index], cuCadd(cuCadd(cuCmul(A02[indexA], PxDev[0]), cuCmul(A12[indexA], PyDev[0])), cuCmul(A22[indexA], PzDev[0])));
        }
    }

}

void APtoESum(cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, 
    cufftDoubleComplex *PxDev, cufftDoubleComplex *PyDev, cufftDoubleComplex *PzDev,
    cufftDoubleComplex *ESumxDev, cufftDoubleComplex *ESumyDev, cufftDoubleComplex *ESumzDev, int NxFFT, int NyFFT, int NzFFT, 
    int NxA, int NyA, int NzA,
    int index1, int index2, int index3, int deduction){
        int tmpx = 10;
        int tmpy = 10;
        int tmpz = 10;
        dim3 dimBlock(tmpx, tmpy, tmpz);
        dim3 dimGrid(ceil((double)NxFFT/tmpx), ceil((double)NyFFT/tmpy), ceil((double)NzFFT/tmpz));   

        APtoESumKernel<<<dimGrid, dimBlock>>>(A00, A01, A02, A11, A12, A22, PxDev, PyDev, PzDev, ESumxDev, ESumyDev, ESumzDev, NxFFT, NyFFT, NzFFT, NxA, NyA, NzA, index1, index2, index3, deduction);

}