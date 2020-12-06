#include "definition.h"



ObjectivePointEDDAModel::ObjectivePointEDDAModel(list<double> parameters, DDAModel *model_, EvoDDAModel* evomodel_, bool HavePenalty_){
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
    Have_Devx = true;
    model = model_;
    evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();
    P = (*model).get_P();
    R = (*Core).get_R();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    double E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    double K = 2*M_PI/lam;
    E_sum = Vector3cd::Zero();
    E_ext = Vector3cd::Zero();
    E_ext(0) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(1) = E0*n_E0(1)*(cos(K*(n_K(0)*y+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(2) = E0*n_E0(2)*(cos(K*(n_K(0)*z+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);   
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectivePointEDDAModel::SingleResponse(int idx, bool deduction){
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx=x-d*(*R)(3*idx);                  //R has no d in it, so needs to time d
    double ry=y-d*(*R)(3*idx+1);
    double rz=z-d*(*R)(3*idx+2);
    //cout << rx << "," << ry << "," << rz << idx << endl;
    AProductCore* Core = (*model).get_Core();
    Matrix3cd A=(*Core).A_dic_generator(rx,ry,rz);
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

double ObjectivePointEDDAModel::GroupResponse(){
    if (Have_Penalty){
        return (E_sum).norm()-(*evomodel).L1Norm();
    }
    else{
        return (E_sum).norm();
    }
    
}

double ObjectivePointEDDAModel::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectivePointEDDAModel::Reset(){
    E_sum(0) = E_ext(0);
    E_sum(1) = E_ext(1);
    E_sum(2) = E_ext(2);
}




ObjectiveIntegratedEDDAModel::ObjectiveIntegratedEDDAModel(list<double> parameters, DDAModel *model_, EvoDDAModel* evomodel_, bool HavePenalty_){
    Have_Penalty = HavePenalty_;
    Have_Devx = true;
    model = model_;
    evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    N = (*Core).get_N();
    Nx = (*Core).get_Nx();
    Ny = (*Core).get_Ny();
    Nz = (*Core).get_Nz();
    d = (*Core).get_d();
    al = (*model).get_al();
    P = (*model).get_P();
    E = VectorXcd::Zero(N*3);
    E_int = 0;
}

void ObjectiveIntegratedEDDAModel::SingleResponse(int idx, bool deduction){
    return;
}

double ObjectiveIntegratedEDDAModel::GroupResponse(){
    for(int i=0;i<N;i++){
        E(i*3) = (*al)(i*3)*(*P)(i*3);
        E(i*3+1) = (*al)(i*3+1)*(*P)(i*3+1);
        E(i*3+2) = (*al)(i*3+2)*(*P)(i*3+2);
    }
    E_int = 0;
    double diel_sum = 0;
    for(int i=0;i<N;i++){
        if((i/(Nx*Ny))>=9){
            double E_sum_temp = 0;
            for(int j=0;j<3;j++){
                E_sum_temp += pow(abs(E(3*i+j)),2);
            }
            diel_sum += (*((*model).get_Core()->get_diel_old()))(i*3);
            E_int += pow(E_sum_temp,2)*((*((*model).get_Core()->get_diel_old()))(i*3)+0.0001)/4.0; //prevent nan result for devp calculation.
        }
    }  
    // E_int /= diel_sum;
    if (Have_Penalty){
        return E_int-(*evomodel).L1Norm();
    }
    else{
        return E_int;
    }
    
}

double ObjectiveIntegratedEDDAModel::GetVal(){
    Reset();
    for(int i=0;i<N;i++){
        E(i*3) = (*al)(i*3)*(*P)(i*3);
        E(i*3+1) = (*al)(i*3+1)*(*P)(i*3+1);
        E(i*3+2) = (*al)(i*3+2)*(*P)(i*3+2);
    }
    return GroupResponse();
}

void ObjectiveIntegratedEDDAModel::Reset(){
    E_int = 0;
    E = VectorXcd::Zero(N*3);
}





ObjectiveMDipoleDDAModel::ObjectiveMDipoleDDAModel(list<double> parameters, DDAModel *model_, EvoDDAModel* evomodel_, bool HavePenalty_){
    VectorXd MDipoleParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        MDipoleParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    x=MDipoleParameters(0);
    y=MDipoleParameters(1);
    z=MDipoleParameters(2);
    Have_Devx = true;
    model = model_;
    evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();
    P = (*model).get_P();
    R = (*Core).get_R();
    al = (*model).get_al();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    double E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    double K = 2*M_PI/lam;
    Dipole_sum = Vector3cd::Zero();
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectiveMDipoleDDAModel::SingleResponse(int idx, bool deduction){
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx=x-d*(*R)(3*idx);                  //R has no d in it, so needs to time d
    double ry=y-d*(*R)(3*idx+1);
    double rz=z-d*(*R)(3*idx+2);
    Vector3cd D=Vector3cd::Zero();
    //cout << rx << "," << ry << "," << rz << idx << endl;
    D(0) = ((*al)(idx*3)+1.0)*(*P)(idx*3)/5.0;
    D(1) = ((*al)(idx*3+1)+1.0)*(*P)(idx*3+1)/5.0;
    D(2) = ((*al)(idx*3+2)+1.0)*(*P)(idx*3+2)/5.0;
    if (deduction == false){
        Dipole_sum(0) -= ry*D(2)-rz*D(1);
        Dipole_sum(1) -= rz*D(0)-rx*D(2);
        Dipole_sum(2) -= rx*D(1)-ry*D(0);
    }
    else{
        Dipole_sum(0) += ry*D(2)-rz*D(1);
        Dipole_sum(1) += rz*D(0)-rx*D(2);
        Dipole_sum(2) += rx*D(1)-ry*D(0); 
    }
}

double ObjectiveMDipoleDDAModel::GroupResponse(){
    if (Have_Penalty){
        return pow((Dipole_sum).norm(),2)-(*evomodel).L1Norm();
    }
    else{
        return pow((Dipole_sum).norm(),2);
    }
    
}

double ObjectiveMDipoleDDAModel::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectiveMDipoleDDAModel::Reset(){
    Dipole_sum = Vector3cd::Zero();
}




ObjectiveTDipoleDDAModel::ObjectiveTDipoleDDAModel(list<double> parameters, DDAModel *model_, EvoDDAModel* evomodel_, bool HavePenalty_){
    VectorXd TDipoleParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        TDipoleParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    x=TDipoleParameters(0);
    y=TDipoleParameters(1);
    z=TDipoleParameters(2);
    Have_Devx = true;
    model = model_;
    evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();
    P = (*model).get_P();
    R = (*Core).get_R();
    al = (*model).get_al();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    double E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    double K = 2*M_PI/lam;
    Dipole_sum = Vector3cd::Zero();
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectiveTDipoleDDAModel::SingleResponse(int idx, bool deduction){
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx=x-d*(*R)(3*idx);                  //R has no d in it, so needs to time d
    double ry=y-d*(*R)(3*idx+1);
    double rz=z-d*(*R)(3*idx+2);
    double rnorm2 = rx*rx+ry*ry+rz*rz;
    Vector3cd D=Vector3cd::Zero();
    //cout << rx << "," << ry << "," << rz << idx << endl;
    D(0) = ((*al)(idx*3)+1.0)*(*P)(idx*3);
    D(1) = ((*al)(idx*3+1)+1.0)*(*P)(idx*3+1);
    D(2) = ((*al)(idx*3+2)+1.0)*(*P)(idx*3+2);
    if (deduction == false){
        Dipole_sum(0) -= (rx*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(0))/5.0;
        Dipole_sum(1) -= (ry*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(1))/5.0;
        Dipole_sum(2) -= (rz*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(2))/5.0;
    }
    else{
        Dipole_sum(0) += (rx*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(0))/5.0;
        Dipole_sum(1) += (ry*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(1))/5.0;
        Dipole_sum(2) += (rz*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(2))/5.0;
    }
}

double ObjectiveTDipoleDDAModel::GroupResponse(){
    if (Have_Penalty){
        return pow((Dipole_sum).norm(),2)-(*evomodel).L1Norm();
    }
    else{
        return pow((Dipole_sum).norm(),2);
    }
    
}

double ObjectiveTDipoleDDAModel::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectiveTDipoleDDAModel::Reset(){
    Dipole_sum = Vector3cd::Zero();
}


ObjectiveEDipoleDDAModel::ObjectiveEDipoleDDAModel(list<double> parameters, DDAModel *model_, EvoDDAModel* evomodel_, bool HavePenalty_){
    VectorXd EDipoleParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        EDipoleParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    x=EDipoleParameters(0);
    y=EDipoleParameters(1);
    z=EDipoleParameters(2);
    Have_Devx = true;
    model = model_;
    evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();
    P = (*model).get_P();
    R = (*Core).get_R();
    al = (*model).get_al();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    double E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    double K = 2*M_PI/lam;
    Dipole_sum = Vector3cd::Zero();
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectiveEDipoleDDAModel::SingleResponse(int idx, bool deduction){
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx=x-d*(*R)(3*idx);                  //R has no d in it, so needs to time d
    double ry=y-d*(*R)(3*idx+1);
    double rz=z-d*(*R)(3*idx+2);
    double rnorm2 = rx*rx+ry*ry+rz*rz;
    Vector3cd D=Vector3cd::Zero();
    //cout << rx << "," << ry << "," << rz << idx << endl;
    D(0) = ((*al)(idx*3)+1.0)*(*P)(idx*3);
    D(1) = ((*al)(idx*3+1)+1.0)*(*P)(idx*3+1);
    D(2) = ((*al)(idx*3+2)+1.0)*(*P)(idx*3+2);
    if (deduction == false){
        Dipole_sum(0) -= D(0);
        Dipole_sum(1) -= D(1);
        Dipole_sum(2) -= D(2);
    }
    else{
        Dipole_sum(0) += D(0);
        Dipole_sum(1) += D(1);
        Dipole_sum(2) += D(2);
    }
}

double ObjectiveEDipoleDDAModel::GroupResponse(){
    if (Have_Penalty){
        return pow((Dipole_sum).norm(),2)-(*evomodel).L1Norm();
    }
    else{
        return pow((Dipole_sum).norm(),2);
    }
    
}

double ObjectiveEDipoleDDAModel::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectiveEDipoleDDAModel::Reset(){
    Dipole_sum = Vector3cd::Zero();
}



ObjectiveMDipole_TDipoleDDAModel::ObjectiveMDipole_TDipoleDDAModel(list<double> parameters, DDAModel *model_, EvoDDAModel* evomodel_, bool HavePenalty_){
    VectorXd MDipole_TDipoleParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        MDipole_TDipoleParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    x=MDipole_TDipoleParameters(0);
    y=MDipole_TDipoleParameters(1);
    z=MDipole_TDipoleParameters(2);
    Have_Devx = true;
    model = model_;
    evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();
    P = (*model).get_P();
    R = (*Core).get_R();
    al = (*model).get_al();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    double E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    K = 2*M_PI/lam;
    MDipole_sum = Vector3cd::Zero();
    TDipole_sum = Vector3cd::Zero();
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectiveMDipole_TDipoleDDAModel::SingleResponse(int idx, bool deduction){
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx=x-d*(*R)(3*idx);                  //R has no d in it, so needs to time d
    double ry=y-d*(*R)(3*idx+1);
    double rz=z-d*(*R)(3*idx+2);
    double rnorm2 = rx*rx+ry*ry+rz*rz;
    Vector3cd D=Vector3cd::Zero();
    //cout << rx << "," << ry << "," << rz << idx << endl;
    D(0) = ((*al)(idx*3)+1.0)*(*P)(idx*3)/5.0;
    D(1) = ((*al)(idx*3+1)+1.0)*(*P)(idx*3+1)/5.0;
    D(2) = ((*al)(idx*3+2)+1.0)*(*P)(idx*3+2)/5.0;
    if (deduction == false){
        MDipole_sum(0) -= ry*D(2)-rz*D(1);
        MDipole_sum(1) -= rz*D(0)-rx*D(2);
        MDipole_sum(2) -= rx*D(1)-ry*D(0);
        TDipole_sum(0) -= (rx*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(0))/5.0;
        TDipole_sum(1) -= (ry*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(1))/5.0;
        TDipole_sum(2) -= (rz*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(2))/5.0;
    }
    else{
        MDipole_sum(0) += ry*D(2)-rz*D(1);
        MDipole_sum(1) += rz*D(0)-rx*D(2);
        MDipole_sum(2) += rx*D(1)-ry*D(0);
        TDipole_sum(0) += (rx*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(0))/5.0;
        TDipole_sum(1) += (ry*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(1))/5.0;
        TDipole_sum(2) += (rz*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(2))/5.0;
    }
}

double ObjectiveMDipole_TDipoleDDAModel::GroupResponse(){
    if (Have_Penalty){
        return pow((MDipole_sum).norm(),2)-pow((TDipole_sum).norm(),2)*pow(K,2)-(*evomodel).L1Norm();
    }
    else{
        return pow((MDipole_sum).norm(),2)-pow((TDipole_sum).norm(),2)*pow(K,2);
    }
    
}

double ObjectiveMDipole_TDipoleDDAModel::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectiveMDipole_TDipoleDDAModel::Reset(){
    MDipole_sum = Vector3cd::Zero();
    TDipole_sum = Vector3cd::Zero();
}


ObjectiveEDipole_TDipoleDDAModel::ObjectiveEDipole_TDipoleDDAModel(list<double> parameters, DDAModel *model_, EvoDDAModel* evomodel_, bool HavePenalty_){
    VectorXd EDipole_TDipoleParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it=(parameters).begin();
    for(int i=0;i<=int((parameters).size()-1);i++){
        EDipole_TDipoleParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    x=EDipole_TDipoleParameters(0);
    y=EDipole_TDipoleParameters(1);
    z=EDipole_TDipoleParameters(2);
    Have_Devx = true;
    model = model_;
    evomodel = evomodel_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();
    P = (*model).get_P();
    R = (*Core).get_R();
    al = (*model).get_al();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    double E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    K = 2*M_PI/lam;
    EDipole_sum = Vector3cd::Zero();
    TDipole_sum = Vector3cd::Zero();
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectiveEDipole_TDipoleDDAModel::SingleResponse(int idx, bool deduction){
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx=x-d*(*R)(3*idx);                  //R has no d in it, so needs to time d
    double ry=y-d*(*R)(3*idx+1);
    double rz=z-d*(*R)(3*idx+2);
    double rnorm2 = rx*rx+ry*ry+rz*rz;
    Vector3cd D=Vector3cd::Zero();
    //cout << rx << "," << ry << "," << rz << idx << endl;
    D(0) = ((*al)(idx*3)+1.0)*(*P)(idx*3)/5.0;
    D(1) = ((*al)(idx*3+1)+1.0)*(*P)(idx*3+1)/5.0;
    D(2) = ((*al)(idx*3+2)+1.0)*(*P)(idx*3+2)/5.0;
    if (deduction == false){
        EDipole_sum(0) -= D(0);
        EDipole_sum(1) -= D(1);
        EDipole_sum(2) -= D(2);
        TDipole_sum(0) -= (rx*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(0))/5.0;
        TDipole_sum(1) -= (ry*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(1))/5.0;
        TDipole_sum(2) -= (rz*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(2))/5.0;
    }
    else{
        EDipole_sum(0) += D(0);
        EDipole_sum(1) += D(1);
        EDipole_sum(2) += D(2);
        TDipole_sum(0) += (rx*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(0))/5.0;
        TDipole_sum(1) += (ry*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(1))/5.0;
        TDipole_sum(2) += (rz*(rx*D(0)+ry*D(1)+rz*D(2))-2*(rnorm2)*D(2))/5.0;
    }
}

double ObjectiveEDipole_TDipoleDDAModel::GroupResponse(){
    if (Have_Penalty){
        return pow((EDipole_sum).norm(),2)-pow((TDipole_sum).norm(),2)*pow(K,2)-(*evomodel).L1Norm();
    }
    else{
        return pow((EDipole_sum).norm(),2)-pow((TDipole_sum).norm(),2)*pow(K,2);
    }
    
}

double ObjectiveEDipole_TDipoleDDAModel::GetVal(){
    Reset();
    for (int idx=0;idx<N;idx++){
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectiveEDipole_TDipoleDDAModel::Reset(){
    EDipole_sum = Vector3cd::Zero();
    TDipole_sum = Vector3cd::Zero();
}



