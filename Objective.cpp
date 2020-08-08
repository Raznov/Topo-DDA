#include "definition.h"



ObjectivePointE::ObjectivePointE(list<double> parameters, EvoModel *model_, bool HavePenalty_){
    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        PointEParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    x=PointEParameters(0);
    y=PointEParameters(1);
    z=PointEParameters(2);
    Have_Devx = false;
    model = model_;
    d = model->get_d();
    N = model->get_N();
    P = model->get_P();
    R = model->get_R();    
    Vector3d n_E0 = model->get_nE0();
    Vector3d n_K = model->get_nK();
    double E0 = model->get_E0();
    double lam = model->get_wl();
    double K = 2*M_PI/lam;
    E_sum = Vector3cd::Zero();
    E_ext = Vector3cd::Zero();
    E_ext(0) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(1) = E0*n_E0(1)*(cos(K*(n_K(0)*y+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(2) = E0*n_E0(2)*(cos(K*(n_K(0)*z+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);   
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectivePointE::SingleResponse(int idx, bool deduction){
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx=x-d*(*R)(3*idx);                  //R has no d in it, so needs to time d
    double ry=y-d*(*R)(3*idx+1);
    double rz=z-d*(*R)(3*idx+2);
    //cout << rx << "," << ry << "," << rz << idx << endl;
    Matrix3cd A=model->A_dic_generator(rx,ry,rz);
    if (deduction == false){
        E_sum(0)-=(A(0,0)*(*P)(3*idx)+A(0,1)*(*P)(3*idx+1)+A(0,2)*(*P)(3*idx+2));
        E_sum(1)-=(A(1,0)*(*P)(3*idx)+A(1,1)*(*P)(3*idx+1)+A(1,2)*(*P)(3*idx+2));
        E_sum(2)-=(A(2,0)*(*P)(3*idx)+A(2,1)*(*P)(3*idx+1)+A(2,2)*(*P)(3*idx+2));
    }
    else{
        E_sum(0)+=(A(0,0)*(*P)(3*idx)+A(0,1)*(*P)(3*idx+1)+A(0,2)*(*P)(3*idx+2));
        E_sum(1)+=(A(1,0)*(*P)(3*idx)+A(1,1)*(*P)(3*idx+1)+A(1,2)*(*P)(3*idx+2));
        E_sum(2)+=(A(2,0)*(*P)(3*idx)+A(2,1)*(*P)(3*idx+1)+A(2,2)*(*P)(3*idx+2));    
    }
}

double ObjectivePointE::GroupResponse(){
    if (Have_Penalty){
        return (E_sum).norm()-model->L1Norm();
    }
    else{
        return (E_sum).norm();
    }
    
}

double ObjectivePointE::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectivePointE::Reset(){
    E_sum(0) = E_ext(0);
    E_sum(1) = E_ext(1);
    E_sum(2) = E_ext(2);
}


ObjectiveSurfaceEExp::ObjectiveSurfaceEExp(list<double> parameters, EvoModel* model_, bool HavePenalty_){
    VectorXd SurfaceEExpParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        SurfaceEExpParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    z=SurfaceEExpParameters(0);                 //The z where the surfaceEexp is calculated(must be inside the entire calculated space)
    exponent = SurfaceEExpParameters(1);
    Have_Devx = false;
    model = model_;
    d = model->get_d();
    N = model->get_N();
    P = model->get_P();
    R = model->get_R();                 //The entire calculated space   
    Einternal = model->get_Einternal();                   //E of points in R

    //--------------------------------------------                 Extract the obj surface from the entire calculated space
    
    nz = int(round(z/d));
    Nobj = 0;
    for(int i=0; i<=N-1; i++){
        if((*R)(3*i+2) == nz){
            Nobj += 1;
        }
    }
    Robj = VectorXi::Zero(Nobj*3);
    E_sum = VectorXcd::Zero(Nobj*3);
    int Posobj = 0;
    for(int i=0; i<=N-1; i++){
        if((*R)(3*i+2) == nz){
            Robj(3*Posobj) = (*R)(3*i);
            Robj(3*Posobj+1) = (*R)(3*i+1);
            Robj(3*Posobj+2) = (*R)(3*i+2);
            Posobj += 1;
        }

    }
    
    //-----------------------------------------------



    Vector3d n_E0 = model->get_nE0();
    Vector3d n_K = model->get_nK();
    double E0 = model->get_E0();
    double lam = model->get_wl();
    double K = 2*M_PI/lam;
    


}

void ObjectiveSurfaceEExp::SingleResponse(int idx, bool deduction){         //Must always run one GetVal before any SingleResponse because Reset give the original value of partial direvative
    double rz=z-d*(*R)(3*idx+2);            //R has no d in it, so needs to time d
    for(int i=0; i<=Nobj-1; i++){
        double x = Robj(3*i);
        double y = Robj(3*i+1);
        double rx=x-d*(*R)(3*idx);                  
        double ry=y-d*(*R)(3*idx+1);

        Matrix3cd A=model->A_dic_generator(rx,ry,rz);
        if (deduction == false){
            E_sum(3*i)-=(A(0,0)*(*P)(3*idx)+A(0,1)*(*P)(3*idx+1)+A(0,2)*(*P)(3*idx+2));
            E_sum(3*i+1)-=(A(1,0)*(*P)(3*idx)+A(1,1)*(*P)(3*idx+1)+A(1,2)*(*P)(3*idx+2));
            E_sum(3*i+2)-=(A(2,0)*(*P)(3*idx)+A(2,1)*(*P)(3*idx+1)+A(2,2)*(*P)(3*idx+2));
        }
        else{
            E_sum(3*i)+=(A(0,0)*(*P)(3*idx)+A(0,1)*(*P)(3*idx+1)+A(0,2)*(*P)(3*idx+2));
            E_sum(3*i+1)+=(A(1,0)*(*P)(3*idx)+A(1,1)*(*P)(3*idx+1)+A(1,2)*(*P)(3*idx+2));
            E_sum(3*i+2)+=(A(2,0)*(*P)(3*idx)+A(2,1)*(*P)(3*idx+1)+A(2,2)*(*P)(3*idx+2));    
        }
    }


}

double ObjectiveSurfaceEExp::GroupResponse(){
    if (Have_Penalty){
        return Average(&E_sum, Nobj, exponent)-model->L1Norm();
    }
    else{
        return Average(&E_sum, Nobj, exponent);
    }
}

double ObjectiveSurfaceEExp::GetVal(){
    Reset();
    return GroupResponse();
}

void ObjectiveSurfaceEExp::Reset(){
    Einternal = model->get_Einternal();
    int Posobj = 0;
    for(int i=0; i<=N-1; i++){
        if((*R)(3*i+2) == nz){
            E_sum(3*Posobj) = (*Einternal)(3*i);
            E_sum(3*Posobj+1) = (*Einternal)(3*i+1);
            E_sum(3*Posobj+2) = (*Einternal)(3*i+2);
            Posobj += 1;
        }
    }
}


ObjectiveExtSurfaceEExp_CPU::ObjectiveExtSurfaceEExp_CPU(list<double> parameters, EvoModel *model_, bool HavePenalty_){
    VectorXd ExtSurfaceEExpParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        ExtSurfaceEExpParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    ExtSurfaceEExpRz = ExtSurfaceEExpParameters(0);                   //(Focus, form the bottom of str to the obj plane(in xy plane))
    exponent = ExtSurfaceEExpParameters(1);
    ratio = int(round(ExtSurfaceEExpParameters(2)));

    Have_Devx = false;
    model = model_;
    d = model->get_d();
    N = model->get_N();
    Nx = model->get_Nx();
    Ny = model->get_Ny();
    Nz = model->get_Nz();
    P = model->get_P();
    R = model->get_R();    
    Vector3d n_E0 = model->get_nE0();
    Vector3d n_K = model->get_nK();
    double E0 = model->get_E0();
    double lam = model->get_wl();
    double K = 2*M_PI/lam;
    
    int Nx_obj = int(floor(Nx / ratio));
    int Ny_obj = int(floor(Ny / ratio));
    Nobj = Nx_obj * Ny_obj;
    Robj = VectorXi::Zero(Nobj * 3);
    E_sum = VectorXcd::Zero(Nobj * 3);
    E_ext = VectorXcd::Zero(Nobj * 3);

    int Posobj = 0;

    for (int i = 0; i <= Nx_obj - 1; i++) {
        for (int j = 0; j <= Ny_obj - 1; j++) {
            Robj(3 * Posobj) = i * ratio;
            Robj(3 * Posobj + 1) = j * ratio;
            Robj(3 * Posobj + 2) = int(round(ExtSurfaceEExpRz / d));
            double x = d * Robj(3 * Posobj);
            double y = d * Robj(3 * Posobj + 1);
            double z = d * Robj(3 * Posobj + 2);
            E_ext(3 * Posobj) = E0 * n_E0(0) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
            E_ext(3 * Posobj + 1) = E0 * n_E0(1) * (cos(K * (n_K(0) * y + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
            E_ext(3 * Posobj + 2) = E0 * n_E0(2) * (cos(K * (n_K(0) * z + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);

            Posobj += 1;
        }
    }
    

    distance0 = int(round(ExtSurfaceEExpRz/d)) - (Nz-1);
    vector<vector<vector<Matrix3cd>>> tmp(2*Nx-1,vector<vector<Matrix3cd>>(2*Ny-1,vector<Matrix3cd>(Nz,Matrix3cd::Zero())));
    A_dic = tmp;
    for(int i=0; i<=2*Nx-2; i++){
        for(int j=0; j<=2*Ny-2; j++){
            for(int k=0; k<=Nz-1; k++){
                double x=d*(i-(Nx-1)); double y=d*(j-(Ny-1)); double z=d*(k + distance0);
                A_dic[i][j][k] = model->A_dic_generator(x,y,z);
            }
        }
    }


}

void ObjectiveExtSurfaceEExp_CPU::SingleResponse(int idx, bool deduction){
    for(int i=0; i<=Nobj-1; i++){
        
        int rnx=Robj(3*i)-(*R)(3*idx)+(Nx-1);                  //R has no d in it, so needs to time d
        int rny=Robj(3*i+1)-(*R)(3*idx+1)+(Ny-1);
        int rnz=Robj(3*i+2)-(*R)(3*idx+2)-distance0;

        //cout<<"distance0"<<distance0<<endl;
        //cout<<"Robjx"<<Robj(3*i)<<endl;
        //cout<<"Robjy"<<Robj(3*i+1)<<endl;
        //cout<<"Robjz"<<Robj(3*i+2)<<endl;
        //cout<<"Rx"<<(*R)(3*idx)<<endl;
        //cout<<"Ry"<<(*R)(3*idx+1)<<endl;
        //cout<<"Rz"<<(*R)(3*idx+2)<<endl;
        //cout<<"rnx "<<rnx<<endl;
        //cout<<"rny "<<rny<<endl;
        //cout<<"rnz "<<rnz<<endl;

        Matrix3cd A=A_dic[rnx][rny][rnz];
        if (deduction == false){
            E_sum(3*i)-=(A(0,0)*(*P)(3*idx)+A(0,1)*(*P)(3*idx+1)+A(0,2)*(*P)(3*idx+2));
            E_sum(3*i+1)-=(A(1,0)*(*P)(3*idx)+A(1,1)*(*P)(3*idx+1)+A(1,2)*(*P)(3*idx+2));
            E_sum(3*i+2)-=(A(2,0)*(*P)(3*idx)+A(2,1)*(*P)(3*idx+1)+A(2,2)*(*P)(3*idx+2));
        }
        else{
            E_sum(3*i)+=(A(0,0)*(*P)(3*idx)+A(0,1)*(*P)(3*idx+1)+A(0,2)*(*P)(3*idx+2));
            E_sum(3*i+1)+=(A(1,0)*(*P)(3*idx)+A(1,1)*(*P)(3*idx+1)+A(1,2)*(*P)(3*idx+2));
            E_sum(3*i+2)+=(A(2,0)*(*P)(3*idx)+A(2,1)*(*P)(3*idx+1)+A(2,2)*(*P)(3*idx+2));    
        }

    }
}

double ObjectiveExtSurfaceEExp_CPU::GroupResponse(){
    if (Have_Penalty){
        return Average(&E_sum, Nobj, exponent)-model->L1Norm();
    }
    else{
        return Average(&E_sum, Nobj, exponent);
    }  
}

double ObjectiveExtSurfaceEExp_CPU::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
    }
    
    /*
    //----bugfix
    VectorXd x = VectorXd::Zero(Nobj);
    for(int i=0; i<=Nobj-1; i++){
        Vector3cd tmp = Vector3cd::Zero();
        tmp(0) = (E_sum)(3*i);
        tmp(1) = (E_sum)(3*i+1);
        tmp(2) = (E_sum)(3*i+2);
        x(i) = pow(tmp.norm(), exponent);
        
    }
    cout<<"Max x"<<x.maxCoeff();
    */
    return GroupResponse();
}

void ObjectiveExtSurfaceEExp_CPU::Reset(){
    for(int i=0; i<=Nobj-1; i++){
        E_sum(3*i) = E_ext(3*i);
        E_sum(3*i+1) = E_ext(3*i+1);
        E_sum(3*i+2) = E_ext(3*i+2);
    }
}


ObjectiveExtSurfaceEExp::ObjectiveExtSurfaceEExp(list<double> parameters, EvoModel *model_, bool HavePenalty_){
    VectorXd ExtSurfaceEExpParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        ExtSurfaceEExpParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    ExtSurfaceEExpRz = ExtSurfaceEExpParameters(0);                   //(Focus, form the bottom of str to the obj plane(in xy plane))
    exponent = ExtSurfaceEExpParameters(1);

    Have_Devx = false;
    model = model_;
    d = model->get_d();
    N = model->get_N();
    Nx = model->get_Nx();
    Ny = model->get_Ny();
    Nz = model->get_Nz();
    P = model->get_P();
    R = model->get_R();    
    Vector3d n_E0 = model->get_nE0();
    Vector3d n_K = model->get_nK();
    double E0 = model->get_E0();
    double lam = model->get_wl();
    double K = 2*M_PI/lam;
    
    Nobj = 0;
    for(int i=0; i<=N-1; i++){
        if((*R)(3*i+2) == 0){
            Nobj += 1;
        }
    }
    Robj = VectorXi::Zero(Nobj*3);
    E_sum = VectorXcd::Zero(Nobj*3);
    E_ext = VectorXcd::Zero(Nobj*3);

    int Posobj = 0;
    for(int i=0; i<=N-1; i++){
        if((*R)(3*i+2) == 0){
            Robj(3*Posobj) = (*R)(3*i);
            Robj(3*Posobj+1) = (*R)(3*i+1);
            Robj(3*Posobj+2) = int(round(ExtSurfaceEExpRz/d));
            double x = d*Robj(3*Posobj);
            double y = d*Robj(3*Posobj+1);
            double z = d*Robj(3*Posobj+2);
            E_ext(3*Posobj) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            E_ext(3*Posobj+1) = E0*n_E0(1)*(cos(K*(n_K(0)*y+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            E_ext(3*Posobj+2) = E0*n_E0(2)*(cos(K*(n_K(0)*z+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i); 

            Posobj += 1;
        }
    }
    

    distance0 = int(round(ExtSurfaceEExpRz/d)) - (Nz-1);
    NxA = 2*Nx-1;
    NyA = 2*Ny-1;
    NzA = Nz;
    NA = NxA*NyA*NzA;

    AHos = new double[2*6*NA];
    for(int i=0; i<=NxA-1; i++){
        for(int j=0; j<=NyA-1; j++){
            for(int k=0; k<=NzA-1; k++){
                double x=d*(i-(Nx-1)); double y=d*(j-(Ny-1)); double z=d*(k + distance0);
                Matrix3cd Atmp = model->A_dic_generator(x,y,z);
                int first[6] = {0, 0, 0, 1, 1, 2};
                int second[6] = {0, 1, 2, 1, 2, 2};
                for(int l=0; l<=5; l++){
                    int index_real = 0+2*l+12*(k+NzA*j+NzA*NyA*i);
                    int index_imag = 1+2*l+12*(k+NzA*j+NzA*NyA*i);
                    AHos[index_real] = Atmp(first[l], second[l]).real();
                    AHos[index_imag] = Atmp(first[l], second[l]).imag();
                }
            }
        }
    }

    cudaMalloc((void**)&ADev, sizeof(double)*2*6*NA);
    cudaMemcpy(ADev, AHos, sizeof(double)*2*6*NA, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&A00, sizeof(cufftDoubleComplex)*NA);
    cudaMalloc((void**)&A01, sizeof(cufftDoubleComplex)*NA);
    cudaMalloc((void**)&A02, sizeof(cufftDoubleComplex)*NA);
    cudaMalloc((void**)&A11, sizeof(cufftDoubleComplex)*NA);
    cudaMalloc((void**)&A12, sizeof(cufftDoubleComplex)*NA);
    cudaMalloc((void**)&A22, sizeof(cufftDoubleComplex)*NA);
    A2As(ADev, A00, A01, A02, A11, A12, A22, NxA, NyA, NzA);
    
    PHos = new double[2*3];
    cudaMalloc((void**)&PDev, sizeof(double)*2*3);
    cudaMalloc((void**)&PxDev, sizeof(cufftDoubleComplex));
    cudaMalloc((void**)&PyDev, sizeof(cufftDoubleComplex));
    cudaMalloc((void**)&PzDev, sizeof(cufftDoubleComplex));

    ESumHos = new double[2*3*Nobj];
    cudaMalloc((void**)&ESumDev, sizeof(double)*2*3*Nobj);
    cudaMalloc((void**)&ESumxDev, sizeof(cufftDoubleComplex)*Nobj);
    cudaMalloc((void**)&ESumyDev, sizeof(cufftDoubleComplex)*Nobj);
    cudaMalloc((void**)&ESumzDev, sizeof(cufftDoubleComplex)*Nobj);
    
    for(int i=0; i<=2*3*Nobj-1; i++){
        ESumHos[i] = 0.0;
    }

}

void ObjectiveExtSurfaceEExp::SingleResponse(int idx, bool deduction){

    PHos[0] = (*P)(3*idx).real();
    PHos[1] = (*P)(3*idx).imag();
    PHos[2] = (*P)(3*idx+1).real();
    PHos[3] = (*P)(3*idx+1).imag();
    PHos[4] = (*P)(3*idx+2).real();
    PHos[5] = (*P)(3*idx+2).imag();
    cudaMemcpy(PDev, PHos, sizeof(double)*2*3, cudaMemcpyHostToDevice);
    cudaMemcpy(ESumDev, ESumHos, sizeof(double)*2*3*Nobj, cudaMemcpyHostToDevice);
    B2Bs(PDev, PxDev, PyDev, PzDev, 1, 1, 1);
    B2Bs(ESumDev, ESumxDev, ESumyDev, ESumzDev, Nx, Ny, 1);
    
    //int rnx=Robj(3*i)-(*R)(3*idx)+(Nx-1);                  //R has no d in it, so needs to time d
    //int rny=Robj(3*i+1)-(*R)(3*idx+1)+(Ny-1);
    //int rnz=Robj(3*i+2)-(*R)(3*idx+2)-distance0;
    
    int index1 = -(*R)(3*idx)+(Nx-1);
    int index2 = -(*R)(3*idx+1)+(Ny-1);
    int index3 = -(*R)(3*idx+2)-distance0;

    if(deduction == false){
        APtoESum(A00, A01, A02, A11, A12, A22, PxDev, PyDev, PzDev, ESumxDev, ESumyDev, ESumzDev, Nx, Ny, 1, NxA, NyA, NzA, index1, index2, index3, -1);
    }
    else{
        APtoESum(A00, A01, A02, A11, A12, A22, PxDev, PyDev, PzDev, ESumxDev, ESumyDev, ESumzDev, Nx, Ny, 1, NxA, NyA, NzA, index1, index2, index3, 1);
    }
    
    Conv2B(ESumxDev, ESumyDev, ESumzDev, ESumDev, Nx, Ny, 1);
    cudaMemcpy(ESumHos, ESumDev, sizeof(double)*2*3*Nobj, cudaMemcpyDeviceToHost);
}

double ObjectiveExtSurfaceEExp::GroupResponse(){
    if (Have_Penalty){
        for(int i=0; i<=Nobj-1; i++){
            E_sum(3*i) = ESumHos[0+2*0+2*3*i] + ESumHos[1+2*0+2*3*i]*1.0i;
            E_sum(3*i+1) = ESumHos[0+2*1+2*3*i] + ESumHos[1+2*1+2*3*i]*1.0i;
            E_sum(3*i+2) = ESumHos[0+2*2+2*3*i] + ESumHos[1+2*2+2*3*i]*1.0i;
        }
        return Average(&E_sum, Nobj, exponent)-model->L1Norm();
    }
    else{
        for(int i=0; i<=Nobj-1; i++){
            E_sum(3*i) = ESumHos[0+2*0+2*3*i] + ESumHos[1+2*0+2*3*i]*1.0i;
            E_sum(3*i+1) = ESumHos[0+2*1+2*3*i] + ESumHos[1+2*1+2*3*i]*1.0i;
            E_sum(3*i+2) = ESumHos[0+2*2+2*3*i] + ESumHos[1+2*2+2*3*i]*1.0i;
        }
        return Average(&E_sum, Nobj, exponent);
    }  
}

double ObjectiveExtSurfaceEExp::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
    }
    
    /*
    //----bugfix
    VectorXd x = VectorXd::Zero(Nobj);
    for(int i=0; i<=Nobj-1; i++){
        Vector3cd tmp = Vector3cd::Zero();
        tmp(0) = (E_sum)(3*i);
        tmp(1) = (E_sum)(3*i+1);
        tmp(2) = (E_sum)(3*i+2);
        x(i) = pow(tmp.norm(), exponent);
        
    }
    cout<<"Max x"<<x.maxCoeff();
    */
    return GroupResponse();
}

void ObjectiveExtSurfaceEExp::Reset(){
    for(int i=0; i<=Nobj-1; i++){
        E_sum(3*i) = E_ext(3*i);
        E_sum(3*i+1) = E_ext(3*i+1);
        E_sum(3*i+2) = E_ext(3*i+2);
        ESumHos[0+2*0+2*3*i] = E_sum(3*i).real();
        ESumHos[1+2*0+2*3*i] = E_sum(3*i).imag();
        ESumHos[0+2*1+2*3*i] = E_sum(3*i+1).real();
        ESumHos[1+2*1+2*3*i] = E_sum(3*i+1).imag();
        ESumHos[0+2*2+2*3*i] = E_sum(3*i+2).real();
        ESumHos[1+2*2+2*3*i] = E_sum(3*i+2).imag();
    }
}

ObjectiveExtSurfaceEExp::~ObjectiveExtSurfaceEExp(){
        cout<<"fuck"<<endl;
        //cout<<"fuck"<<endl;
        delete [] AHos;
        delete [] ADev; 
        //cout<<"fuck"<<endl;
        cudaFree(A00);
        cudaFree(A01);
        cudaFree(A02);
        cudaFree(A11);
        cudaFree(A12);
        cudaFree(A22);
        delete [] PHos;
        delete [] PDev;
        cudaFree(PxDev);
        cudaFree(PyDev);
        cudaFree(PzDev);
        delete [] ESumHos;
        delete [] ESumDev;
        cudaFree(ESumxDev);
        cudaFree(ESumyDev);
        cudaFree(ESumzDev);
        cout<<"fuck"<<endl;
}










double Average(VectorXcd* E, int N, double exponent){
    VectorXd x = VectorXd::Zero(N);
    for(int i = 0; i<=N-1; i++){
        Vector3cd tmp = Vector3cd::Zero();
        tmp(0) = (*E)(3*i);
        tmp(1) = (*E)(3*i+1);
        tmp(2) = (*E)(3*i+2);
        x(i) = exp(pow(tmp.norm(), exponent));
    }
    double avg = x.sum()/N;
    return avg;
}