#include "definition.h"
#define NUM_THREADS 6



AProductCore::AProductCore(Space* space_, double d_, double lam_, Vector2cd material_){
    
    space = space_;
    d=d_;
    lam=lam_;
    K=2*M_PI/lam;
    material = material_;

    cout << "(d=" << d << ") " << "(lam=" << lam << ") " << "(K=" << K << ") " << endl;

    tie(Nx, Ny, Nz, N)=(*space_).get_Ns();
    list<Structure> *ln=(*space_).get_ln();
    R=VectorXi::Zero(3*N);
    RDep = VectorXi::Zero(3 * N);
    VectorXd diel_tmp = VectorXd::Zero(3 * N);

    //-----------------------------------------------------------------Input strs-------------------------------------------------------------
    list<Structure>::iterator it = (*ln).begin();

    int PositionParaN = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = ((*((*it).get_geometry())).size());
        if ((*it).get_para() == 1) {
            PositionParaN += int(round(n2 / 3));
        }
        it++;
    }
    PositionPara = VectorXi::Zero(PositionParaN);

    //Build R, RDep
    int n1 = 0;
    it = (*ln).begin();
    int PositionParaPos = 0;
    VectorXi RPara = VectorXi::Zero(3 * N);
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = ((*((*it).get_geometry())).size());
        if ((*it).get_para() == 1) {
            para_nums.push_back(n2);
            para_starts.push_back(n1);
            for (int i = int(round(n1 / 3)); i <= int(round(n1 / 3)) + int(round(n2 / 3)) - 1; i++) {
                PositionPara(PositionParaPos) = i;
                PositionParaPos += 1;
            }
        }
        if ((*it).get_para() == 2) {
            para_dep_nums.push_back(n2);
            para_dep_starts.push_back(n1);
        }
        for (int j = 0; j <= n2 - 1; j++) {
            R(n1 + j) = (*((*it).get_geometry()))(j);
            diel_tmp(n1 + j) = (*((*it).get_diel()))(j);
            RDep(n1 + j) = (*((*it).get_geometry_dep()))(j);
            RPara(n1 + j) = (*it).get_para();
        }
        it++;
        n1 = n1 + n2;
    }


    for (int i = 0; i <= PositionPara.size() - 1; i++) {
        list<int> tmpPos;
        int Parax = R(3 * PositionPara(i));
        int Paray = R(3 * PositionPara(i) + 1);
        int Paraz = R(3 * PositionPara(i) + 2);
        for (int j = 0; j <= N - 1; j++) {
            int jx = R(3 * j);
            int jy = R(3 * j + 1);
            int jz = R(3 * j + 2);
            int Depx = RDep(3 * j);
            int Depy = RDep(3 * j + 1);
            int Depz = RDep(3 * j + 2);
            if (RPara(3 * j) == 2) {
                if (Depx == Parax && Depy == Paray && Depz == Paraz) {
                    tmpPos.push_back(j);
                }
            }
        }
        PositionDep.push_back(tmpPos);
    }

    int PositionDepN = 0;
    list<list<int>>::iterator it1PositionDep = PositionDep.begin();
    for (int i = 0; i <= int(PositionDep.size()) - 1; i++) {
        PositionDepN += (*it1PositionDep).size();
        it1PositionDep++;
    }
    if (PositionDepN + PositionParaN != N) {
        cout << "PositionDepN = " << PositionDepN << endl;
        cout << "PositionParaN = " << PositionParaN << endl;
        cout << "In AProductCore::AProductCore(Space* space_, double d_, double lam_, Vector2cd material_) : PositionDepN + PositionParaN! = N" << endl;
    }

    //---------------------------------------------------initial diel------------------------------------
    diel = VectorXcd::Zero(3 * N);
    diel_old = VectorXd::Zero(3 * N);

    diel_old_max = diel_old;
    diel_max = diel;

    for (int i = 0; i <= 3 * N - 1; i++) {
        diel(i) = material(0) + diel_tmp(i) * (material(1) - material(0));
        diel_old(i) = diel_tmp(i);
    }
    
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

AProductCore::AProductCore(Space* space_, double d_, double lam_, Vector2cd material_, int MAXm_, int MAXn_, double Lm_, double Ln_) {

    space = space_;
    d = d_;
    lam = lam_;
    K = 2 * M_PI / lam;
    material = material_;

    cout << "(d=" << d << ") " << "(lam=" << lam << ") " << "(K=" << K << ") " << endl;

    tie(Nx, Ny, Nz, N) = (*space_).get_Ns();
    list<Structure>* ln = (*space_).get_ln();
    R = VectorXi::Zero(3 * N);
    RDep = VectorXi::Zero(3 * N);
    VectorXd diel_tmp = VectorXd::Zero(3 * N);

    //-----------------------------------------------------------------Input strs-------------------------------------------------------------
    list<Structure>::iterator it = (*ln).begin();

    int PositionParaN = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = ((*((*it).get_geometry())).size());
        if ((*it).get_para() == 1) {
            PositionParaN += int(round(n2 / 3));
        }
        it++;
    }
    PositionPara = VectorXi::Zero(PositionParaN);

    //Build R, RDep
    int n1 = 0;
    it = (*ln).begin();
    int PositionParaPos = 0;
    VectorXi RPara = VectorXi::Zero(3 * N);
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = ((*((*it).get_geometry())).size());
        if ((*it).get_para() == 1) {
            para_nums.push_back(n2);
            para_starts.push_back(n1);
            for (int i = int(round(n1 / 3)); i <= int(round(n1 / 3)) + int(round(n2 / 3)) - 1; i++) {
                PositionPara(PositionParaPos) = i;
                PositionParaPos += 1;
            }
        }
        if ((*it).get_para() == 2) {
            para_dep_nums.push_back(n2);
            para_dep_starts.push_back(n1);
        }
        for (int j = 0; j <= n2 - 1; j++) {
            R(n1 + j) = (*((*it).get_geometry()))(j);
            diel_tmp(n1 + j) = (*((*it).get_diel()))(j);
            RDep(n1 + j) = (*((*it).get_geometry_dep()))(j);
            RPara(n1 + j) = (*it).get_para();
        }
        it++;
        n1 = n1 + n2;
    }


    for (int i = 0; i <= PositionPara.size() - 1; i++) {
        list<int> tmpPos;
        int Parax = R(3 * PositionPara(i));
        int Paray = R(3 * PositionPara(i) + 1);
        int Paraz = R(3 * PositionPara(i) + 2);
        for (int j = 0; j <= N - 1; j++) {
            int jx = R(3 * j);
            int jy = R(3 * j + 1);
            int jz = R(3 * j + 2);
            int Depx = RDep(3 * j);
            int Depy = RDep(3 * j + 1);
            int Depz = RDep(3 * j + 2);
            if (RPara(3 * j) == 2) {
                if (Depx == Parax && Depy == Paray && Depz == Paraz) {
                    tmpPos.push_back(j);
                }
            }
        }
        PositionDep.push_back(tmpPos);
    }

    int PositionDepN = 0;
    list<list<int>>::iterator it1PositionDep = PositionDep.begin();
    for (int i = 0; i <= int(PositionDep.size()) - 1; i++) {
        PositionDepN += (*it1PositionDep).size();
        it1PositionDep++;
    }
    if (PositionDepN + PositionParaN != N) {
        cout << "PositionDepN = " << PositionDepN << endl;
        cout << "PositionParaN = " << PositionParaN << endl;
        cout << "In AProductCore::AProductCore(Space* space_, double d_, double lam_, Vector2cd material_) : PositionDepN + PositionParaN! = N" << endl;
    }

    //---------------------------------------------------initial diel------------------------------------
    diel = VectorXcd::Zero(3 * N);
    diel_old = VectorXd::Zero(3 * N);

    diel_old_max = diel_old;
    diel_max = diel;

    for (int i = 0; i <= 3 * N - 1; i++) {
        diel(i) = material(0) + diel_tmp(i) * (material(1) - material(0));
        diel_old(i) = diel_tmp(i);
    }

    //Allocation and part of initialization for FFT
    //NFFT:
    NxFFT = 2 * Nx - 1;
    NyFFT = 2 * Ny - 1;
    NzFFT = 2 * Nz - 1;
    NFFT = NxFFT * NyFFT * NzFFT;

    //Plan:
    if (cufftPlan3d(&Plan, NxFFT, NyFFT, NzFFT, CUFFT_Z2Z) != CUFFT_SUCCESS) {
        fprintf(stderr, "CUFFT error: Plan creation failed");
    }

    high_resolution_clock::time_point t_init = high_resolution_clock::now();
    //A_dic:
    AHos = new double[2 * 6 * NFFT];
    for (int i = 0; i <= Nx - 1; i++) {
        for (int j = 0; j <= Ny - 1; j++) {
            for (int k = 0; k <= Nz - 1; k++) {
                double x = d * i; double y = d * j; double z = d * k;
                Matrix3cd Atmp = Matrix3cd::Zero();
                for (int m = -MAXm; m <= MAXm; m++) {
                    for (int n = -MAXn; n <= MAXn; n++) {
                        Atmp = Atmp + this->A_dic_generator(x, y, z, m, n);
                    }
                }

                int first[6] = { 0, 0, 0, 1, 1, 2 };
                int second[6] = { 0, 1, 2, 1, 2, 2 };
                for (int l = 0; l <= 5; l++) {
                    int index_real = 0 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    int index_imag = 1 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();

                }
            }
        }
    }
    for (int i = Nx; i <= 2 * Nx - 2; i++) {
        for (int j = 0; j <= Ny - 1; j++) {
            for (int k = 0; k <= Nz - 1; k++) {
                double x = d * (i - (2 * Nx - 1)); double y = d * j; double z = d * k;
                Matrix3cd Atmp = Matrix3cd::Zero();
                for (int m = -MAXm; m <= MAXm; m++) {
                    for (int n = -MAXn; n <= MAXn; n++) {
                        Atmp = Atmp + this->A_dic_generator(x, y, z, m, n);
                    }
                }
                int first[6] = { 0, 0, 0, 1, 1, 2 };
                int second[6] = { 0, 1, 2, 1, 2, 2 };
                for (int l = 0; l <= 5; l++) {
                    int index_real = 0 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    int index_imag = 1 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();

                }
            }
        }
    }
    for (int i = 0; i <= Nx - 1; i++) {
        for (int j = Ny; j <= 2 * Ny - 2; j++) {
            for (int k = 0; k <= Nz - 1; k++) {
                double x = d * i; double y = d * (j - (2 * Ny - 1)); double z = d * k;
                Matrix3cd Atmp = Matrix3cd::Zero();
                for (int m = -MAXm; m <= MAXm; m++) {
                    for (int n = -MAXn; n <= MAXn; n++) {
                        Atmp = Atmp + this->A_dic_generator(x, y, z, m, n);
                    }
                }
                int first[6] = { 0, 0, 0, 1, 1, 2 };
                int second[6] = { 0, 1, 2, 1, 2, 2 };
                for (int l = 0; l <= 5; l++) {
                    int index_real = 0 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    int index_imag = 1 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();

                }
            }
        }
    }
    for (int i = 0; i <= Nx - 1; i++) {
        for (int j = 0; j <= Ny - 1; j++) {
            for (int k = Nz; k <= 2 * Nz - 2; k++) {
                double x = d * i; double y = d * j; double z = d * (k - (2 * Nz - 1));
                Matrix3cd Atmp = Matrix3cd::Zero();
                for (int m = -MAXm; m <= MAXm; m++) {
                    for (int n = -MAXn; n <= MAXn; n++) {
                        Atmp = Atmp + this->A_dic_generator(x, y, z, m, n);
                    }
                }
                int first[6] = { 0, 0, 0, 1, 1, 2 };
                int second[6] = { 0, 1, 2, 1, 2, 2 };
                for (int l = 0; l <= 5; l++) {
                    int index_real = 0 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    int index_imag = 1 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();

                }
            }
        }
    }
    for (int i = Nx; i <= 2 * Nx - 2; i++) {
        for (int j = Ny; j <= 2 * Ny - 2; j++) {
            for (int k = 0; k <= Nz - 1; k++) {
                double x = d * (i - (2 * Nx - 1)); double y = d * (j - (2 * Ny - 1)); double z = d * k;
                Matrix3cd Atmp = Matrix3cd::Zero();
                for (int m = -MAXm; m <= MAXm; m++) {
                    for (int n = -MAXn; n <= MAXn; n++) {
                        Atmp = Atmp + this->A_dic_generator(x, y, z, m, n);
                    }
                }
                int first[6] = { 0, 0, 0, 1, 1, 2 };
                int second[6] = { 0, 1, 2, 1, 2, 2 };
                for (int l = 0; l <= 5; l++) {
                    int index_real = 0 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    int index_imag = 1 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();

                }
            }
        }
    }
    for (int i = Nx; i <= 2 * Nx - 2; i++) {
        for (int j = 0; j <= Ny - 1; j++) {
            for (int k = Nz; k <= 2 * Nz - 2; k++) {
                double x = d * (i - (2 * Nx - 1)); double y = d * j; double z = d * (k - (2 * Nz - 1));
                Matrix3cd Atmp = Matrix3cd::Zero();
                for (int m = -MAXm; m <= MAXm; m++) {
                    for (int n = -MAXn; n <= MAXn; n++) {
                        Atmp = Atmp + this->A_dic_generator(x, y, z, m, n);
                    }
                }
                int first[6] = { 0, 0, 0, 1, 1, 2 };
                int second[6] = { 0, 1, 2, 1, 2, 2 };
                for (int l = 0; l <= 5; l++) {
                    int index_real = 0 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    int index_imag = 1 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();

                }
            }
        }
    }
    for (int i = 0; i <= Nx - 1; i++) {
        for (int j = Ny; j <= 2 * Ny - 2; j++) {
            for (int k = Nz; k <= 2 * Nz - 2; k++) {
                double x = d * i; double y = d * (j - (2 * Ny - 1)); double z = d * (k - (2 * Nz - 1));
                Matrix3cd Atmp = Matrix3cd::Zero();
                for (int m = -MAXm; m <= MAXm; m++) {
                    for (int n = -MAXn; n <= MAXn; n++) {
                        Atmp = Atmp + this->A_dic_generator(x, y, z, m, n);
                    }
                }
                int first[6] = { 0, 0, 0, 1, 1, 2 };
                int second[6] = { 0, 1, 2, 1, 2, 2 };
                for (int l = 0; l <= 5; l++) {
                    int index_real = 0 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    int index_imag = 1 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();

                }
            }
        }
    }
    for (int i = Nx; i <= 2 * Nx - 2; i++) {
        for (int j = Ny; j <= 2 * Ny - 2; j++) {
            for (int k = Nz; k <= 2 * Nz - 2; k++) {
                double x = d * (i - (2 * Nx - 1)); double y = d * (j - (2 * Ny - 1)); double z = d * (k - (2 * Nz - 1));
                Matrix3cd Atmp = Matrix3cd::Zero();
                for (int m = -MAXm; m <= MAXm; m++) {
                    for (int n = -MAXn; n <= MAXn; n++) {
                        Atmp = Atmp + this->A_dic_generator(x, y, z, m, n);
                    }
                }
                int first[6] = { 0, 0, 0, 1, 1, 2 };
                int second[6] = { 0, 1, 2, 1, 2, 2 };
                for (int l = 0; l <= 5; l++) {
                    int index_real = 0 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    int index_imag = 1 + 2 * l + 12 * (k + NzFFT * j + NzFFT * NyFFT * i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();

                }
            }
        }
    }

    cudaMalloc((void**)&ADev, sizeof(double) * 2 * 6 * NFFT);
    cudaMemcpy(ADev, AHos, sizeof(double) * 2 * 6 * NFFT, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&A00, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&A01, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&A02, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&A11, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&A12, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&A22, sizeof(cufftDoubleComplex) * NFFT);
    A2As(ADev, A00, A01, A02, A11, A12, A22, NxFFT, NyFFT, NzFFT);
    high_resolution_clock::time_point t_init_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_init_end - t_init).count();
    cout << "Cost " << duration << "ms to build and allocate A from host to device" << endl;

    //FFT of As:
    if (cufftExecZ2Z(Plan, A00, A00, CUFFT_FORWARD) != CUFFT_SUCCESS) {
        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A01, A01, CUFFT_FORWARD) != CUFFT_SUCCESS) {
        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A02, A02, CUFFT_FORWARD) != CUFFT_SUCCESS) {
        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A11, A11, CUFFT_FORWARD) != CUFFT_SUCCESS) {
        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A12, A12, CUFFT_FORWARD) != CUFFT_SUCCESS) {
        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }
    if (cufftExecZ2Z(Plan, A22, A22, CUFFT_FORWARD) != CUFFT_SUCCESS) {
        fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed");
    }



    //b:
    bHos = new double[2 * 3 * NFFT];
    cudaMalloc((void**)&bDev, sizeof(double) * 2 * 3 * NFFT);
    cudaMalloc((void**)&bxDev, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&byDev, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&bzDev, sizeof(cufftDoubleComplex) * NFFT);



    //Conv:
    cudaMalloc((void**)&Convx, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&Convy, sizeof(cufftDoubleComplex) * NFFT);
    cudaMalloc((void**)&Convz, sizeof(cufftDoubleComplex) * NFFT);


}

AProductCore::~AProductCore(){
        delete [] AHos;
        cudaFree(ADev);
        cudaFree(ADev); 
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

Matrix3cd AProductCore::A_dic_generator(double x,double y,double z){
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

Matrix3cd AProductCore::A_dic_generator(double x, double y, double z, int m, int n) {
    x = x + m * Lm;
    y = y + n * Ln;
    return this->A_dic_generator(x, y, z);
}

VectorXcd AProductCore::Aproduct(VectorXcd &b){
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
            result(3*i+j)=(bHos[0+2*j+2*3*position] + bHos[1+2*j+2*3*position]*1.0i)/double(NFFT); // + b(3*i+j)*al(3*i+j)
        }
    
    }
    return result;

}

void AProductCore::UpdateStr(VectorXd step) {

    cout << "step in UpdateStr" << step.mean() << endl;

    int para_size = para_nums.size();
    int para_dep_size = para_dep_nums.size();
    //cout << para_size << ' ' << para_dep_size << endl;
    if (para_dep_size != 0) {
        if (PositionPara.size() != PositionDep.size()) {
            cout << "In Model::change_para_diel(VectorXd step) : PositionPara.size() != PositionDep.size()" << endl;
        }

        list<list<int>>::iterator it1 = PositionDep.begin();
        for (int i = 0; i <= PositionPara.size() - 1; i++) {
            int position1 = PositionPara(i);
            diel_old(3 * position1) += step(i);
            if (diel_old(3 * position1) >= 1) {
                diel_old(3 * position1) = 1;
            }
            if (diel_old(3 * position1) <= 0) {
                diel_old(3 * position1) = 0;
            }

            diel_old(3 * position1 + 1) = diel_old(3 * position1);
            diel_old(3 * position1 + 2) = diel_old(3 * position1);


            diel(3 * position1) = material(0) + diel_old(3 * position1) * (material(1) - material(0));
            diel(3 * position1 + 1) = diel(3 * position1);
            diel(3 * position1 + 2) = diel(3 * position1);

            list<int>::iterator it2 = (*it1).begin();
            for (int j = 0; j <= (*it1).size() - 1; j++) {
                int position2 = (*it2);
                diel_old(3 * position2) = diel_old(3 * position1);
                diel_old(3 * position2 + 1) = diel_old(3 * position1);
                diel_old(3 * position2 + 2) = diel_old(3 * position1);

                diel(3 * position2) = material(0) + diel_old(3 * position2) * (material(1) - material(0));
                diel(3 * position2 + 1) = diel(3 * position2);
                diel(3 * position2 + 2) = diel(3 * position2);

                it2++;
            }
            it1++;
        }
    }
    else {

        list<int>::iterator it1 = para_nums.begin();
        list<int>::iterator it2 = para_starts.begin();
        int position = 0;
        for (int i = 0; i <= para_size - 1; i++) {
            int para_begin = round((*it2) / 3);
            int para_number = round((*it1) / 3);

            for (int j = 0; j <= para_number - 1; j++) {
                int position1 = (j + para_begin);
                diel_old(3 * position1) += step(position);
                if (diel_old(3 * position1) >= 1) {
                    diel_old(3 * position1) = 1;
                }
                if (diel_old(3 * position1) <= 0) {
                    diel_old(3 * position1) = 0;
                }

                diel_old(3 * position1 + 1) = diel_old(3 * position1);
                diel_old(3 * position1 + 2) = diel_old(3 * position1);


                diel(3 * position1) = material(0) + diel_old(3 * position1) * (material(1) - material(0));
                diel(3 * position1 + 1) = diel(3 * position1);
                diel(3 * position1 + 2) = diel(3 * position1);

                position = position + 1;

            }
            it1++;
            it2++;
        }

    }

}

void AProductCore::output_to_file() {

    ofstream fout("AProductCore.txt");
    fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
    fout << R << endl;
    fout << diel_old << endl;
    fout << d << endl;
    fout << lam << endl;
    fout.close();
}

void AProductCore::output_to_file(string save_position, int iteration) {

    string name;
    name=save_position+"AProductCoreit" + to_string(iteration) + ".txt";
    ofstream fout(name);
    fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
    fout << R << endl;
    fout << diel_old << endl;
    fout << d << endl;
    fout << lam << endl;
    fout.close();
}

int AProductCore::get_N() {
    return N;
}
int AProductCore::get_Nx() {
    return Nx;
}
int AProductCore::get_Ny() {
    return Ny;
}
int AProductCore::get_Nz() {
    return Nz;
}
tuple<list<int>, list<int>, list<int>, list<int>> AProductCore::get_para_info() {
    return make_tuple(para_nums, para_starts, para_dep_nums, para_dep_starts);
}
VectorXi* AProductCore::get_R() {
    return &R;
}
double AProductCore::get_d() {
    return d;
}
Space* AProductCore::get_space() {
    return space;
}
double AProductCore::get_lam() {
    return lam;
}
VectorXcd* AProductCore::get_diel() {
    return &diel;
}
list<list<int>>* AProductCore::get_PositionDep() {
    return &PositionDep;
}
VectorXi* AProductCore::get_PositionPara() {
    return &PositionPara;
}
list<int>* AProductCore::get_para_nums() {
    return &para_nums;
}
list<int>* AProductCore::get_para_starts() {
    return &para_starts;
}
list<int>* AProductCore::get_para_dep_nums() {
    return &para_dep_nums;
}
list<int>* AProductCore::get_para_dep_starts() {
    return &para_dep_starts;
}
VectorXd* AProductCore::get_diel_old() {
    return &diel_old;
}
Vector2cd* AProductCore::get_material() {
    return &material;
}
VectorXcd* AProductCore::get_diel_max() {
    return &diel_max;
}
VectorXd* AProductCore::get_diel_old_max() {
    return &diel_old_max;
}
