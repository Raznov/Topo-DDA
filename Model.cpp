#include "definition.h"
#define NUM_THREADS 6

int Model::get_N(){
    return N;
}
int Model::get_Nx(){
    return Nx;
}
int Model::get_Ny(){
    return Ny;
}
int Model::get_Nz(){
    return Nz;
}

VectorXcd* Model::get_P(){
    return &P;
}

VectorXi* Model::get_R(){
    return &R;
}

double Model::get_d(){
    return d;
}

Vector3d Model::get_nE0(){
    return n_E0;
}

Vector3d Model::get_nK(){
    return n_K;
}

double Model::get_E0(){
    return E0;
}

double Model::get_wl(){
    return lam;
}

VectorXcd* Model::get_Einternal(){
    return &Einternal;
}

Model::Model(Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_){
    

    time=0;
    ITERATION=0;
    Error=0.0;
    d=d_;
    lam=lam_;
    K=2*M_PI/lam;
    E0=E0_;
    n_K=n_K_;
    n_E0=n_E0_;
    material=material_;

    tie(Nx, Ny, Nz, N)=(*space_).get_Ns();
    list<Structure> *ln=(*space_).get_ln();
    R=VectorXi::Zero(3*N);

    RResultSwitch = false;
    
    VectorXd diel_tmp=VectorXd::Zero(3*N);
    // Why not directly use diel_old?
    
    list<Structure>::iterator it=(*ln).begin();
    int n1=0;
    for(int i=0;i<=int((*ln).size())-1;i++){
        int n2=((*((*it).get_geometry())).size());
        if((*it).get_para()==1){
            para_nums.push_back(n2);
            para_starts.push_back(n1);
        }
        if((*it).get_para()==2){
            para_dep_nums.push_back(n2);
            para_dep_starts.push_back(n1);
        }
        for(int j=0;j<=n2-1;j++){
            R(n1+j)=(*((*it).get_geometry()))(j);
            diel_tmp(n1+j)=(*((*it).get_diel()))(j);
        }
        it++;
        n1=n1+n2;
    }
    
    RResult = R;
    diel=VectorXcd::Zero(3*N);
    diel_old=VectorXd::Zero(3*N);

    diel_old_max = diel_old;
    diel_max = diel;

    for(int i=0;i<=3*N-1;i++){
        diel(i)=material(0)+diel_tmp(i)*(material(1)-material(0));
        diel_old(i)=diel_tmp(i);
    }
    
    P = VectorXcd::Zero(N*3);
    P_max = P;
    E = VectorXcd::Zero(N*3);
    Einternal = VectorXcd::Zero(N*3);
    EResult = VectorXcd::Zero(N*3);
    for (int i=0;i<N;i++) {
        E(3*i) = E0*n_E0(0)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);
        E(3*i+1) = E0*n_E0(1)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);
        E(3*i+2) = E0*n_E0(2)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);                                           
    }
    al = VectorXcd::Zero(N*3);
    for (int i=0;i<N*3;i++) {
        al(i)=1.0/Get_Alpha(lam,K,d,diel(i));
    }  
    al_max = al;
    verbose = true;

    
    //Allocation and part of initialization for FFT
    //NFFT:
    NxFFT = 2*Nx-1;
    NyFFT = 2*Ny-1;
    NzFFT = 2*Nz-1;
    NFFT = NxFFT*NyFFT*NzFFT;

    //Plan:
    if (cufftPlan3d(&Plan, NxFFT, NyFFT, NzFFT, CUFFT_Z2Z) != CUFFT_SUCCESS){
	    fprintf(stderr, "CUFFT error: Plan creation failed");
    }
    
    high_resolution_clock::time_point t_init = high_resolution_clock::now();
    //A_dic:
    AHos = new double[2*6*NFFT];
    for(int i=0; i<=Nx-1; i++){
        for(int j=0; j<=Ny-1; j++){
            for(int k=0; k<=Nz-1; k++){
                double x=d*i; double y=d*j; double z=d*k;
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*j; double z=d*k;
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }             
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*i; double y=d*(j-(2*Ny-1)); double z=d*k;
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }      
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*i; double y=d*j; double z=d*(k-(2*Nz-1));
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*(j-(2*Ny-1)); double z=d*k;
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*j; double z=d*(k-(2*Nz-1));
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*i; double y=d*(j-(2*Ny-1)); double z=d*(k-(2*Nz-1));
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*(j-(2*Ny-1)); double z=d*(k-(2*Nz-1));
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    
    cudaMalloc((void**)&ADev, sizeof(double)*2*6*NFFT);
    cudaMemcpy(ADev, AHos, sizeof(double)*2*6*NFFT, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&A00, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A01, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A02, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A11, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A12, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A22, sizeof(cufftDoubleComplex)*NFFT);
    A2As(ADev, A00, A01, A02, A11, A12, A22, NxFFT, NyFFT, NzFFT);
    high_resolution_clock::time_point t_init_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_init_end-t_init).count();
    cout<<"Cost "<<duration<<"ms to build and allocate A from host to device"<<endl;
    
    //FFT of As:
    if (cufftExecZ2Z(Plan, A00, A00, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A01, A01, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A02, A02, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A11, A11, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A12, A12, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A22, A22, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    
    

    //b:
    bHos = new double[2*3*NFFT];
    cudaMalloc((void**)&bDev, sizeof(double)*2*3*NFFT);
    cudaMalloc((void**)&bxDev, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&byDev, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&bzDev, sizeof(cufftDoubleComplex)*NFFT);
    
    

    //Conv:
    cudaMalloc((void**)&Convx, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&Convy, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&Convz, sizeof(cufftDoubleComplex)*NFFT);
    
    
}

Model::Model(Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_, VectorXi *RResult_){
    time=0;
    ITERATION=0;
    Error=0.0;
    d=d_;
    lam=lam_;
    K=2*M_PI/lam;
    E0=E0_;
    n_K=n_K_;
    n_E0=n_E0_;
    material=material_;

    tie(Nx, Ny, Nz, N)=(*space_).get_Ns();
    list<Structure> *ln=(*space_).get_ln();
    R=VectorXi::Zero(3*N);

    RResultSwitch = true;

    RResult = *RResult_;
    if(RResult.size()%3 != 0){
    
        cout<<"RResult does not have a size with 3*integer"<<endl;
        // This could be changed to throwing an exception.
    }
    
    VectorXd diel_tmp=VectorXd::Zero(3*N);
    // Why not directly use diel_old?
    
    list<Structure>::iterator it=(*ln).begin();
    int n1=0;
    for(int i=0;i<=int((*ln).size())-1;i++){
        int n2=((*((*it).get_geometry())).size());
        if((*it).get_para()==1){
            para_nums.push_back(n2);
            para_starts.push_back(n1);
        }
        if((*it).get_para()==2){
            para_dep_nums.push_back(n2);
            para_dep_starts.push_back(n1);
        }
        for(int j=0;j<=n2-1;j++){
            R(n1+j)=(*((*it).get_geometry()))(j);
            diel_tmp(n1+j)=(*((*it).get_diel()))(j);
        }
        it++;
        n1=n1+n2;
    }
    diel=VectorXcd::Zero(3*N);
    diel_old=VectorXd::Zero(3*N);
    for(int i=0;i<=3*N-1;i++){
        diel(i)=material(0)+diel_tmp(i)*(material(1)-material(0));
        diel_old(i)=diel_tmp(i);
    }
    
    diel_old_max = diel_old;
    diel_max = diel;

    P = VectorXcd::Zero(N*3);
    P_max = P;
    E = VectorXcd::Zero(N*3);
    Einternal = VectorXcd::Zero(N*3);
    EResult = VectorXcd::Zero(RResult.size());
    for (int i=0;i<N;i++) {
        E(3*i) = E0*n_E0(0)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);
        E(3*i+1) = E0*n_E0(1)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);
        E(3*i+2) = E0*n_E0(2)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);                                           
    }
    al = VectorXcd::Zero(N*3);
    for (int i=0;i<N*3;i++) {
        al(i)=1.0/Get_Alpha(lam,K,d,diel(i));
    }  
    al_max = al;
    verbose = true;
    

    //Allocation and part of initialization for FFT
    //NFFT:
    NxFFT = 2*Nx-1;
    NyFFT = 2*Ny-1;
    NzFFT = 2*Nz-1;
    NFFT = NxFFT*NyFFT*NzFFT;

    //Plan:
    if (cufftPlan3d(&Plan, NxFFT, NyFFT, NzFFT, CUFFT_Z2Z) != CUFFT_SUCCESS){
	    fprintf(stderr, "CUFFT error: Plan creation failed");
    }
    
    high_resolution_clock::time_point t_init = high_resolution_clock::now();
    //A_dic:
    AHos = new double[2*6*NFFT];
    for(int i=0; i<=Nx-1; i++){
        for(int j=0; j<=Ny-1; j++){
            for(int k=0; k<=Nz-1; k++){
                double x=d*i; double y=d*j; double z=d*k;
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*j; double z=d*k;
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }             
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*i; double y=d*(j-(2*Ny-1)); double z=d*k;
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }      
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*i; double y=d*j; double z=d*(k-(2*Nz-1));
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*(j-(2*Ny-1)); double z=d*k;
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*j; double z=d*(k-(2*Nz-1));
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*i; double y=d*(j-(2*Ny-1)); double z=d*(k-(2*Nz-1));
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*(j-(2*Ny-1)); double z=d*(k-(2*Nz-1));
                Matrix3cd Atmp = this->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    int index_imag = 1+2*l+12*(k+NzFFT*j+NzFFT*NyFFT*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                    
                }
            }
        }
    }
    cudaMalloc((void**)&ADev, sizeof(double)*2*6*NFFT);
    cudaMemcpy(ADev, AHos, sizeof(double)*2*6*NFFT, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&A00, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A01, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A02, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A11, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A12, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&A22, sizeof(cufftDoubleComplex)*NFFT);
    A2As(ADev, A00, A01, A02, A11, A12, A22, NxFFT, NyFFT, NzFFT);

    high_resolution_clock::time_point t_init_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_init_end-t_init).count();
    cout<<"Cost "<<duration<<"ms to build and allocate A from host to device"<<endl;
    
    //FFT of As:
    if (cufftExecZ2Z(Plan, A00, A00, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A01, A01, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A02, A02, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A11, A11, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A12, A12, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A22, A22, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    
    //b:
    bHos = new double[2*3*NFFT];
    cudaMalloc((void**)&bDev, sizeof(double)*2*3*NFFT);
    cudaMalloc((void**)&bxDev, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&byDev, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&bzDev, sizeof(cufftDoubleComplex)*NFFT);

    //Conv:
    cudaMalloc((void**)&Convx, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&Convy, sizeof(cufftDoubleComplex)*NFFT);
    cudaMalloc((void**)&Convz, sizeof(cufftDoubleComplex)*NFFT);
         
    
}

Model::~Model(){
        //cout<<"fuck"<<endl;
        delete [] AHos;
        cudaFree(ADev); 
        //cout<<"fuck"<<endl;
        cudaFree(A00);
        cudaFree(A01);
        cudaFree(A02);
        cudaFree(A11);
        cudaFree(A12);
        cudaFree(A22);
        delete [] bHos;
        cudaFree(bDev);
        cudaFree(bxDev);
        cudaFree(byDev);
        cudaFree(bzDev);
        cudaFree(Convx);
        cudaFree(Convy);
        cudaFree(Convz);
        cufftDestroy(Plan);
        
}

Matrix3cd Model::A_dic_generator(double x,double y,double z){
    Matrix3cd result(3,3);
    double xsquare = x*x; double ysquare = y*y; double zsquare = z*z;
    double rnorm = sqrt(xsquare+ysquare+zsquare);

    if (rnorm==0.0) {
        result(0,0) = 0.0;
        result(1,1) = 0.0;
        result(2,2) = 0.0;
        return result;
    }
    else {
        double xy = x*y; double yz = y*z; double zx = z*x;
        double rsquare = rnorm*rnorm;
        double rcubic = rnorm*rnorm*rnorm;
        std::complex<double> const1 = 1i*K*rnorm;
        const1 = exp(const1)/rcubic;
        std::complex<double> const2 = (1.0-K*rnorm*1i);
        const2 = const2/rsquare;
        std::complex<double> const3 = K*K-3.0*const2;

        result(0,0) = const2*(rsquare-3*xsquare)-K*K*(ysquare+zsquare);
        result(0,1) = xy*const3;
        result(0,2) = zx*const3;
        result(1,1) = const2*(rsquare-3*ysquare)-K*K*(zsquare+xsquare);
        result(1,2) = yz*const3;
        result(2,2) = const2*(rsquare-3*zsquare)-K*K*(xsquare+ysquare);

        result(1,0) = result(0,1);
        result(2,0) = result(0,2);
        result(2,1) = result(1,2);

        result = const1*result;
        return result;
    }
}

VectorXcd Model::Aproduct(VectorXcd &b){
    //! in geometry use np.meshgrid with 'ij'  
    for(int i=0; i<=2*3*NFFT-1; i++){
        bHos[i] = 0.0;
    }
    for(int i=0;i<=N-1;i++){
        int x=round(R(3*i));
        int y=round(R(3*i+1));
        int z=round(R(3*i+2));
        for(int j=0; j<=2; j++){
            int index_real = 0+2*j+2*3*(z+NzFFT*(y+NyFFT*x));
            int index_imag = 1+2*j+2*3*(z+NzFFT*(y+NyFFT*x));
            bHos[index_real] = b(3*i+j).real();
            bHos[index_imag] = b(3*i+j).imag();
        }
    }
    cudaMemcpy(bDev, bHos, sizeof(double)*2*3*NFFT, cudaMemcpyHostToDevice);

    B2Bs(bDev, bxDev, byDev, bzDev, NxFFT, NyFFT, NzFFT);
    if (cufftExecZ2Z(Plan, bxDev, bxDev, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, byDev, byDev, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, bzDev, bzDev, CUFFT_FORWARD) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }

    Conv(Convx, Convy, Convz, A00, A01, A02, A11, A12, A22, bxDev, byDev, bzDev, NxFFT, NyFFT, NzFFT);

    if (cufftExecZ2Z(Plan, Convx, Convx, CUFFT_INVERSE) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Inverse failed");	
    }
    if (cufftExecZ2Z(Plan, Convy, Convy, CUFFT_INVERSE) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Inverse failed");	
    }
    if (cufftExecZ2Z(Plan, Convz, Convz, CUFFT_INVERSE) != CUFFT_SUCCESS){
	        fprintf(stderr, "CUFFT error: ExecZ2Z Inverse failed");	
    }
    
    Conv2B(Convx, Convy, Convz, bDev, NxFFT, NyFFT, NzFFT);
    cudaMemcpy(bHos, bDev, sizeof(double)*2*3*NFFT, cudaMemcpyDeviceToHost);
    
   
    VectorXcd result(3*N);
    for(int i=0;i<=N-1;i++){
        int x=round(R(3*i));
        int y=round(R(3*i+1));
        int z=round(R(3*i+2));
        int position=z+NzFFT*(y+NyFFT*x);
        for(int j=0 ;j<=2; j++){
            result(3*i+j)=(bHos[0+2*j+2*3*position] + bHos[1+2*j+2*3*position]*1.0i)/double(NFFT) + b(3*i+j)*al(3*i+j);
        }
    
    }
    return result;

}


void Model::bicgstab(int MAX_ITERATION,double MAX_ERROR){
    if (verbose) {
        cout << "--------------Calculation start. Iterative method used: BICGSTAB---------------" << endl;
        cout << endl;
    }

    
    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    //fftw_init_threads();                                                       //////////Initialize the multi-thread
    //fftw_plan_with_nthreads(NUM_THREADS);
    cout<<"Threads"<<NUM_THREADS<<endl;;
    VectorXcd p=VectorXcd::Zero(N*3); 
    VectorXcd t = VectorXcd::Zero(N*3);
    VectorXcd w = VectorXcd::Zero(N*3);
    VectorXcd r = VectorXcd::Zero(N*3);
    VectorXcd r0 = VectorXcd::Zero(N*3);
    VectorXcd rl = VectorXcd::Zero(N*3);
    VectorXcd y = VectorXcd::Zero(N*3);
    VectorXcd u = VectorXcd::Zero(N*3);
    VectorXcd z = VectorXcd::Zero(N*3);
    VectorXcd x = VectorXcd::Zero(N*3);
    std::complex<double> alpha;
    std::complex<double> beta;
    std::complex<double> eta;
    std::complex<double> zeta;
    
    




    VectorXcd Ax0 = this->Aproduct(P);
    
    r = E-Ax0;
    r0 = r;
    p = r;
    VectorXcd Ap0 = this->Aproduct(p);


    alpha = r0.dot(r)/r0.dot(Ap0);
    t = r-alpha*Ap0;
    VectorXcd At0 = this->Aproduct(t);


    zeta = At0.dot(t)/At0.dot(At0);
    u = zeta*Ap0;
    z = zeta*r-alpha*u;
    P = P+alpha*p+z;                                    //this will directly change P in this.
    rl = r;
    r = t-zeta*At0;

    VectorXcd Ap = VectorXcd::Zero(N*3);
    VectorXcd At = VectorXcd::Zero(N*3);
    for (int it=0;it<=MAX_ITERATION-1;it++) {
        if (verbose && (it+1)%100==0) {
            cout << "                Iter " << it+1 << ", error=" << r.norm()/E.norm() << " MAX_ERROR="<<MAX_ERROR<<endl;
        }
        beta = (alpha/zeta)*r0.dot(r)/r0.dot(rl);
        p = r+beta*(p-u);
        Ap = this->Aproduct(p);
        alpha = r0.dot(r)/r0.dot(Ap);
        t = r-alpha*Ap;
        At = this->Aproduct(t);

        zeta = At.dot(t)/At.dot(At);
        u = zeta*Ap;
        z = zeta*r-alpha*u;
        P = P+alpha*p+z;
        rl = r;
        r = t-zeta*At;

        if (r.norm()/E.norm()<=MAX_ERROR) {
            Error = r.norm()/E.norm();
            ITERATION = it+1;
            if (verbose) {
                high_resolution_clock::time_point t_end = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(t_end-t_start).count();
                time = duration_cast<milliseconds>(t_end-t_start).count();
                cout << "--------------Calculation finished. Duration: " << duration/1000.0 << "s.-------------------" << endl;
                ofstream fout;
                fout.open("DDATime.txt", fstream::app);
                fout<<N<<" "<<duration/1000.0<<endl;
                fout.close();

                cout << "              Error: "<<Error<<endl;
                cout << "              Iteration: "<<ITERATION<<endl;
                cout << endl;
            }
            return;
        }
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    time = duration_cast<milliseconds>(t_end-t_start).count();
    cout<<"                ERROR:does not converge in "<<MAX_ITERATION<<" iterations"<<endl;
    return;
}

void Model::change_E(VectorXcd E_){
    E=E_;
}

void Model::reset_E(){
    for (int i=0;i<N;i++) {
        E(3*i) = E0*n_E0(0)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);
        E(3*i+1) = E0*n_E0(1)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);
        E(3*i+2) = E0*n_E0(2)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);                                            
    }
}

double Model::get_step_length(VectorXd gradients, double epsilon){
    int para_size = para_nums.size();
    int para_dep_size = para_dep_nums.size();
    double cur_step = epsilon;
    int position = 0;
    cout << para_size << ' ' << para_dep_size << endl;
    if(para_dep_size!=0){
        if(para_dep_size!=para_size){
            cout<<"ERROR: para_dep_size not equal para_size"<<endl;
        }
        list<int>::iterator it1=para_nums.begin();
        list<int>::iterator it2=para_starts.begin();
        for(int i=0;i<=para_size-1;i++){
            int para_begin=round((*it2)/3);
            int para_number=round((*it1)/3);
            for(int j=0;j<=para_number-1;j++){
                int position1=(j+para_begin);
                if (gradients(position)>0){
                    cout << "1" << endl;
                    cur_step = min(cur_step, 0.1/gradients(position));
                    if (diel_old(3*position1)<0.95){
                        cur_step = min(epsilon,1.0*(1.0-diel_old(3*position1))/gradients(position));
                    }    
                }
                else if (gradients(position)<0){
                    cout << "2" << endl;
                    cur_step = min(cur_step, -0.1/gradients(position));
                    if (diel_old(3*position1)>0.05){
                        cur_step = min(epsilon,1.0*(-diel_old(3*position1))/gradients(position));
                    }
                }
                position=position+1;
            }
            it1++;
            it2++;
        }  
    }
    else{
        list<int>::iterator it1=para_nums.begin();
        list<int>::iterator it2=para_starts.begin();
        for(int i=0;i<=para_size-1;i++){
            int para_begin=round((*it2)/3);
            int para_number=round((*it1)/3);
            for(int j=0;j<=para_number-1;j++){
                int position1=(j+para_begin);
                if (gradients(position)>0){
                    cur_step = min(cur_step, 0.1/gradients(position));
                    if (diel_old(3*position1)<0.95){
                        cur_step = min(cur_step,1.0*(1.0-diel_old(3*position1))/gradients(position));
                    }    
                }
                else if (gradients(position)<0){
                    cur_step = min(cur_step, -0.1/gradients(position));
                    if (diel_old(3*position1)>0.05){
                        cur_step = min(cur_step,1.0*(-diel_old(3*position1))/gradients(position));
                    }
                }
                position=position+1;
            }
            it1++;
            it2++;
        }
    }
    return cur_step;
}


void Model::change_para_diel(VectorXd step){
    
    int para_size=para_nums.size();
    int para_dep_size=para_dep_nums.size();
    cout << para_size << ' ' << para_dep_size << endl;
    if(para_dep_size!=0){
        if(para_dep_size!=para_size){
            cout<<"ERROR: para_dep_size not equal para_size"<<endl;
        }
        list<int>::iterator it1=para_nums.begin();
        list<int>::iterator it2=para_starts.begin();
        list<int>::iterator it3=para_dep_nums.begin();
        list<int>::iterator it4=para_dep_starts.begin();
        int position=0;
        for(int i=0;i<=para_size-1;i++){
            int times=round((*it3)/(*it1));
            int para_begin=round((*it2)/3);
            int para_number=round((*it1)/3);
            int para_dep_begin=round((*it4)/3);
            for(int j=0;j<=para_number-1;j++){
                int position1=(j+para_begin);
                diel_old(3*position1)+=step(position);
                if(diel_old(3*position1)>=1){
                    diel_old(3*position1)=1;
                }
                if(diel_old(3*position1)<=0){
                    diel_old(3*position1)=0;
                }

                diel_old(3*position1+1)=diel_old(3*position1);
                diel_old(3*position1+2)=diel_old(3*position1);
                

                diel(3*position1)=material(0)+diel_old(3*position1)*(material(1)-material(0));
                diel(3*position1+1)=diel(3*position1);
                diel(3*position1+2)=diel(3*position1);
                
                al(3*position1)=1.0/Get_Alpha(lam,K,d,diel(3*position1));
                al(3*position1+1)=1.0/Get_Alpha(lam,K,d,diel(3*position1+1));
                al(3*position1+2)=1.0/Get_Alpha(lam,K,d,diel(3*position1+2));

                for(int k=0;k<=times-1;k++){
                    int position2=j+para_dep_begin+k*para_number;
                    diel_old(3*position2)=diel_old(3*position1);
                    diel_old(3*position2+1)=diel_old(3*position1);
                    diel_old(3*position2+2)=diel_old(3*position1);

                    diel(3*position2)=material(0)+diel_old(3*position2)*(material(1)-material(0));
                    diel(3*position2+1)=diel(3*position2);
                    diel(3*position2+2)=diel(3*position2);

                    al(3*position2)=1.0/Get_Alpha(lam,K,d,diel(3*position2));
                    al(3*position2+1)=1.0/Get_Alpha(lam,K,d,diel(3*position2+1));
                    al(3*position2+2)=1.0/Get_Alpha(lam,K,d,diel(3*position2+2));
                    

                }
                position=position+1;
            }
            it1++;
            it2++;
            it3++;
            it4++;
        }
    }
    else{
        
        list<int>::iterator it1=para_nums.begin();
        list<int>::iterator it2=para_starts.begin();
        int position=0;
        for(int i=0;i<=para_size-1;i++){
            int para_begin=round((*it2)/3);
            int para_number=round((*it1)/3);

            for(int j=0;j<=para_number-1;j++){
                int position1=(j+para_begin);
                diel_old(3*position1)+=step(position);
                if(diel_old(3*position1)>=1){
                    diel_old(3*position1)=1;
                }
                if(diel_old(3*position1)<=0){
                    diel_old(3*position1)=0;
                }

                diel_old(3*position1+1)=diel_old(3*position1);
                diel_old(3*position1+2)=diel_old(3*position1);
                

                diel(3*position1)=material(0)+diel_old(3*position1)*(material(1)-material(0));
                diel(3*position1+1)=diel(3*position1);
                diel(3*position1+2)=diel(3*position1);
                
                al(3*position1)=1.0/Get_Alpha(lam,K,d,diel(3*position1));
                al(3*position1+1)=1.0/Get_Alpha(lam,K,d,diel(3*position1+1));
                al(3*position1+2)=1.0/Get_Alpha(lam,K,d,diel(3*position1+2));

                position=position+1;

            }
            it1++;
            it2++;
        }
        
    }

}

tuple<list<int>, list<int>, list<int>, list<int>> Model::get_para_info(){
    return make_tuple(para_nums, para_starts, para_dep_nums, para_dep_starts);
}

void Model::solve_E(){
    if(RResultSwitch == true){
        for (int j = 0;j<int(round(RResult.size()/3));j++){
            double x = d*RResult(3*j);
            double y = d*RResult(3*j+1);
            double z = d*RResult(3*j+2);
            Vector3cd sum=Vector3cd::Zero();
            Vector3cd E_ext=Vector3cd::Zero();
            E_ext(0) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            E_ext(1) = E0*n_E0(1)*(cos(K*(n_K(0)*y+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            E_ext(2) = E0*n_E0(2)*(cos(K*(n_K(0)*z+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            for (int i=0;i<N;i++){
                double rx=x-d*R(3*i);                  //R has no d in it, so needs to time d
                double ry=y-d*R(3*i+1);
                double rz=z-d*R(3*i+2);
                Matrix3cd A=this->A_dic_generator(rx,ry,rz);
                sum(0)+=(A(0,0)*P(3*i)+A(0,1)*P(3*i+1)+A(0,2)*P(3*i+2));
                sum(1)+=(A(1,0)*P(3*i)+A(1,1)*P(3*i+1)+A(1,2)*P(3*i+2));
                sum(2)+=(A(2,0)*P(3*i)+A(2,1)*P(3*i+1)+A(2,2)*P(3*i+2));
            }
            EResult(3*j) = E_ext(0)-sum(0);
            EResult(3*j+1) = E_ext(1)-sum(1);
            EResult(3*j+2) = E_ext(2)-sum(2);
        }
    }
    else{
        EResult = Einternal;
    }

    
}

void Model::update_E_in_structure(){
    for(int i=0;i<=3*N-1;i++){

        Einternal(i)=al(i)*P(i);
    }
}

void Model::output_to_file(){
    
    ofstream fout("Model_results.txt");
    fout<<Nx<<endl<<Ny<<endl<<Nz<<endl<<N<<endl;
    fout<<R<<endl;
    for(int i=0;i<=P.size()-1;i++){
        if(P(i).imag()<0){
            fout<<P(i).real()<<P(i).imag()<<"j"<<endl;
        }
        else{
            fout<<P(i).real()<<"+"<<P(i).imag()<<"j"<<endl;
        }
        
    }
    fout<<diel_old<<endl;
    fout<<d<<endl;
    fout<<lam<<endl;
    fout<<n_K<<endl;
    fout<<n_E0<<endl;
    fout<<EResult.size()<<endl;
    fout<<RResult<<endl;
    for(int i=0;i<=EResult.size()-1;i++){
        if(EResult(i).imag()<0){
            fout<<EResult(i).real()<<EResult(i).imag()<<"j"<<endl;
        }
        else{
            fout<<EResult(i).real()<<"+"<<EResult(i).imag()<<"j"<<endl;
        }
        
    }
    fout.close();
}

void Model::output_to_file(string save_position, int iteration){
    
    string name;
    name=save_position+"Model_results"+to_string(iteration)+".txt";
    ofstream fout(name);
    fout<<Nx<<endl<<Ny<<endl<<Nz<<endl<<N<<endl;
    fout<<R<<endl;
    for(int i=0;i<=P.size()-1;i++){
        if(P(i).imag()<0){
            fout<<P(i).real()<<P(i).imag()<<"j"<<endl;
        }
        else{
            fout<<P(i).real()<<"+"<<P(i).imag()<<"j"<<endl;
        }
        
    }
    fout<<diel_old<<endl;
    fout<<d<<endl;
    fout<<lam<<endl;
    fout<<n_K<<endl;
    fout<<n_E0<<endl;
    fout<<EResult.size()<<endl;
    fout<<RResult<<endl;
    for(int i=0;i<=EResult.size()-1;i++){
        if(EResult(i).imag()<0){
            fout<<EResult(i).real()<<EResult(i).imag()<<"j"<<endl;
        }
        else{
            fout<<EResult(i).real()<<"+"<<EResult(i).imag()<<"j"<<endl;
        }
        
    }
    fout.close();
}
