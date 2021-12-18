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
    Have_Devx = false;
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
    cout << "lam" << lam << endl;
    double K = 2*M_PI/lam;
    E_sum = Vector3cd::Zero();                                                                             //�ǲ���E_sum���˼�E_ext�ˣ� It is actually in Rest.
    E_ext = Vector3cd::Zero();
    E_ext(0) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(1) = E0*n_E0(1)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(2) = E0*n_E0(2)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);   
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




ObjectivePointListEDDAModel::ObjectivePointListEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it = (parameters).begin();
    if ((parameters.size() % 3) != 0) {
        cout << "ERROR: (parameters.size() % 3) != 0" << endl;
    }
    PNum = int(round(parameters.size() / 3));
    x = VectorXd::Zero(PNum);
    y = VectorXd::Zero(PNum);
    z = VectorXd::Zero(PNum);
    E_sum.resize(PNum, 3);
    E_ext.resize(PNum, 3);
    for (int i = 0; i <= PNum - 1; i++) {
        x(i) = (*it);
        it++;
        y(i) = (*it);
        it++;
        z(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    
    Have_Devx = false;
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
    cout << "lam" << lam << endl;
    double K = 2 * M_PI / lam;

    for (int i = 0; i <= PNum - 1; i++) {
        E_ext(i, 0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) + sin(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) * 1i);
        E_ext(i, 1) = E0 * n_E0(1) * (cos(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) + sin(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) * 1i);
        E_ext(i, 2) = E0 * n_E0(2) * (cos(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) + sin(K * (n_K(0) * x(i) + n_K(1) * y(i) + n_K(2) * z(i))) * 1i);
    }
    
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectivePointListEDDAModel::SingleResponse(int idx, bool deduction) {
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    AProductCore* Core = (*model).get_Core();
    for (int i = 0; i <= PNum - 1; i++) {
        double rx = x(i) - d * (*R)(3 * idx);                  //R has no d in it, so needs to time d
        double ry = y(i) - d * (*R)(3 * idx + 1);
        double rz = z(i) - d * (*R)(3 * idx + 2);
        Matrix3cd A = (*Core).A_dic_generator(rx, ry, rz);
        if (deduction == false) {
            E_sum(i, 0) -= (A(0, 0) * (*P)(3 * idx) + A(0, 1) * (*P)(3 * idx + 1) + A(0, 2) * (*P)(3 * idx + 2));
            E_sum(i, 1) -= (A(1, 0) * (*P)(3 * idx) + A(1, 1) * (*P)(3 * idx + 1) + A(1, 2) * (*P)(3 * idx + 2));
            E_sum(i, 2) -= (A(2, 0) * (*P)(3 * idx) + A(2, 1) * (*P)(3 * idx + 1) + A(2, 2) * (*P)(3 * idx + 2));
        }
        else {
            E_sum(i, 0) += (A(0, 0) * (*P)(3 * idx) + A(0, 1) * (*P)(3 * idx + 1) + A(0, 2) * (*P)(3 * idx + 2));
            E_sum(i, 1) += (A(1, 0) * (*P)(3 * idx) + A(1, 1) * (*P)(3 * idx + 1) + A(1, 2) * (*P)(3 * idx + 2));
            E_sum(i, 2) += (A(2, 0) * (*P)(3 * idx) + A(2, 1) * (*P)(3 * idx + 1) + A(2, 2) * (*P)(3 * idx + 2));
        }
    }

    
    //cout << rx << "," << ry << "," << rz << idx << endl;
    
    
    
}

double ObjectivePointListEDDAModel::GroupResponse() {
    double result = 0;
    if (Have_Penalty) {
        for (int i = 0; i <= PNum - 1; i++) {
            Vector3cd tmp;
            tmp(0) = E_sum(i, 0);
            tmp(1) = E_sum(i, 1);
            tmp(2) = E_sum(i, 2);
            //result = result+ sqrt(pow(tmp(0).real(), 2) + pow(tmp(0).imag(), 2)+ pow(tmp(1).real(), 2) + pow(tmp(1).imag(), 2)+ pow(tmp(2).real(), 2) + pow(tmp(2).imag(), 2));
            result += tmp.norm();
        }
        
        return result/PNum - (*evomodel).L1Norm();
    }
    else {
        for (int i = 0; i <= PNum - 1; i++) {
            Vector3cd tmp;
            tmp(0) = E_sum(i, 0);
            tmp(1) = E_sum(i, 1);
            tmp(2) = E_sum(i, 2);
            //result = result + sqrt(pow(tmp(0).real(), 2) + pow(tmp(0).imag(), 2) + pow(tmp(1).real(), 2) + pow(tmp(1).imag(), 2) + pow(tmp(2).real(), 2) + pow(tmp(2).imag(), 2));
            result += tmp.norm();
        }
        return result / PNum;
    }

}

double ObjectivePointListEDDAModel::GetVal() {
    Reset();
    for (int idx = 0; idx < N; idx++) {
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectivePointListEDDAModel::Reset() {
    for (int i = 0; i <= PNum - 1; i++) {
        E_sum(i, 0) = E_ext(i, 0);
        E_sum(i, 1) = E_ext(i, 1);
        E_sum(i, 2) = E_ext(i, 2);
    }
    
}




ObjectivePointIDDAModel::ObjectivePointIDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it = (parameters).begin();
    for (int i = 0; i <= int((parameters).size() - 1); i++) {
        PointEParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    x = PointEParameters(0);
    y = PointEParameters(1);
    z = PointEParameters(2);
    Have_Devx = false;
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
    cout << "lam" << lam << endl;
    double K = 2 * M_PI / lam;
    E_sum = Vector3cd::Zero();
    E_ext = Vector3cd::Zero();
    E_ext(0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
    E_ext(1) = E0 * n_E0(1) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
    E_ext(2) = E0 * n_E0(2) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
    // cout << R(3*5444+2) << "this" << endl;
}

void ObjectivePointIDDAModel::SingleResponse(int idx, bool deduction) {
    //VectorXcd P = model->get_P();
    //VectorXi R = model->get_R();
    double rx = x - d * (*R)(3 * idx);                  //R has no d in it, so needs to time d
    double ry = y - d * (*R)(3 * idx + 1);
    double rz = z - d * (*R)(3 * idx + 2);
    //cout << rx << "," << ry << "," << rz << idx << endl;
    AProductCore* Core = (*model).get_Core();
    Matrix3cd A = (*Core).A_dic_generator(rx, ry, rz);
    if (deduction == false) {
        E_sum(0) -= (A(0, 0) * (*P)(3 * idx) + A(0, 1) * (*P)(3 * idx + 1) + A(0, 2) * (*P)(3 * idx + 2));
        E_sum(1) -= (A(1, 0) * (*P)(3 * idx) + A(1, 1) * (*P)(3 * idx + 1) + A(1, 2) * (*P)(3 * idx + 2));
        E_sum(2) -= (A(2, 0) * (*P)(3 * idx) + A(2, 1) * (*P)(3 * idx + 1) + A(2, 2) * (*P)(3 * idx + 2));
    }
    else {
        E_sum(0) += (A(0, 0) * (*P)(3 * idx) + A(0, 1) * (*P)(3 * idx + 1) + A(0, 2) * (*P)(3 * idx + 2));
        E_sum(1) += (A(1, 0) * (*P)(3 * idx) + A(1, 1) * (*P)(3 * idx + 1) + A(1, 2) * (*P)(3 * idx + 2));
        E_sum(2) += (A(2, 0) * (*P)(3 * idx) + A(2, 1) * (*P)(3 * idx + 1) + A(2, 2) * (*P)(3 * idx + 2));
    }
}

double ObjectivePointIDDAModel::GroupResponse() {
    if (Have_Penalty) {
        return pow((E_sum).norm(),2) - (*evomodel).L1Norm();
    }
    else {
        return pow((E_sum).norm(),2);
    }

}

double ObjectivePointIDDAModel::GetVal() {
    Reset();
    for (int idx = 0; idx < N; idx++) {
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void ObjectivePointIDDAModel::Reset() {
    E_sum(0) = E_ext(0);
    E_sum(1) = E_ext(1);
    E_sum(2) = E_ext(2);
}



ObjectiveIntegratedEDDAModel::ObjectiveIntegratedEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
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
    E = VectorXcd::Zero(N * 3);
    R = (*Core).get_R();
    E_int = 0;
}

void ObjectiveIntegratedEDDAModel::SingleResponse(int idx, bool deduction) {
    return;
}

double ObjectiveIntegratedEDDAModel::GroupResponse() {
    for (int i = 0; i < N; i++) {
        E(i * 3) = (*al)(i * 3) * (*P)(i * 3);
        E(i * 3 + 1) = (*al)(i * 3 + 1) * (*P)(i * 3 + 1);
        E(i * 3 + 2) = (*al)(i * 3 + 2) * (*P)(i * 3 + 2);
    }
    E_int = 0;
    double diel_sum = 0;
    
    for (int i = 0; i < N; i++) {
        if ((*R)(3*i+2) >= 8) {
            double E_sum_temp = 0;
            for (int j = 0; j < 3; j++) {
                E_sum_temp += pow(abs(E(3 * i + j)), 2);
            }
            diel_sum += (*((*model).get_Core()->get_diel_old()))(i * 3);
            E_int += pow(E_sum_temp, 2) * ((*((*model).get_Core()->get_diel_old()))(i * 3) + 0.0001) / 4.0; //prevent nan result for devp calculation.
        }
    }
    
    E_int = log(E_int);


    // E_int /= diel_sum;
    if (Have_Penalty) {
        return E_int - (*evomodel).L1Norm();
    }
    else {
        return E_int;
    }

}

double ObjectiveIntegratedEDDAModel::GetVal() {
    Reset();
    for (int i = 0; i < N; i++) {
        E(i * 3) = (*al)(i * 3) * (*P)(i * 3);
        E(i * 3 + 1) = (*al)(i * 3 + 1) * (*P)(i * 3 + 1);
        E(i * 3 + 2) = (*al)(i * 3 + 2) * (*P)(i * 3 + 2);
    }
    return GroupResponse();
}

void ObjectiveIntegratedEDDAModel::Reset() {
    E_int = 0;
    E = VectorXcd::Zero(N * 3);
}



ObjectiveMidAvgEDDAModel::ObjectiveMidAvgEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it = (parameters).begin();
    for (int i = 0; i <= int((parameters).size() - 1); i++) {
        PointEParameters(i) = (*it);
        it++;
    }
    r = PointEParameters(0);
    
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
    E = VectorXcd::Zero(N * 3);
    R = (*Core).get_R();

}

void ObjectiveMidAvgEDDAModel::SingleResponse(int idx, bool deduction) {
    return;
}

double ObjectiveMidAvgEDDAModel::GroupResponse() {
    for (int i = 0; i < N; i++) {
        E(i * 3) = (*al)(i * 3) * (*P)(i * 3);
        E(i * 3 + 1) = (*al)(i * 3 + 1) * (*P)(i * 3 + 1);
        E(i * 3 + 2) = (*al)(i * 3 + 2) * (*P)(i * 3 + 2);
    }

    double avg_num = 0;

    for (int i = 0; i < N; i++) {
        if ((48<=(*R)(3 * i) * d<=118 && 48 <= (*R)(3 * i+1) * d <= 118) || ((212 <= (*R)(3 * i) * d <= 282 && 212 <= (*R)(3 * i + 1) * d <= 282))) {
            double E_sum_temp = 0;
            
            for (int j = 0; j < 3; j++) {
                E_sum_temp += pow(abs(E(3 * i + j)), 2);
            }
            E_sum_temp = sqrt(E_sum_temp);

            E_avg += E_sum_temp * (pow((*((*model).get_Core()->get_diel_old()))(i * 3), 1) + 0.0001);
            avg_num += 1;
        }
    }

    E_avg = E_avg;


    // E_int /= diel_sum;
    if (Have_Penalty) {
        return E_avg - (*evomodel).L1Norm();
    }
    else {
        return E_avg;
    }

}

double ObjectiveMidAvgEDDAModel::GetVal() {
    Reset();
    for (int i = 0; i < N; i++) {
        E(i * 3) = (*al)(i * 3) * (*P)(i * 3);
        E(i * 3 + 1) = (*al)(i * 3 + 1) * (*P)(i * 3 + 1);
        E(i * 3 + 2) = (*al)(i * 3 + 2) * (*P)(i * 3 + 2);
    }
    return GroupResponse();
}

void ObjectiveMidAvgEDDAModel::Reset() {
    E_avg = 0;
    E = VectorXcd::Zero(N * 3);
}

Objectivescattering0D::Objectivescattering0D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
    Paralength = (parameters).size();
    VectorXd FOMParameters = VectorXd::Zero(Paralength);
    list<double>::iterator it = (parameters).begin();
    for (int i = 0; i <= int(Paralength - 1); i++) {
        FOMParameters(i) = (*it);
        it++;
    }
    if (Paralength % 3 != 0) {
        cout << "FOMscattering2D ERROR: parameter must be times of 3." << endl;
    }

    Have_Penalty = HavePenalty_;
    Have_Devx = false;
    model = model_;
    evomodel = evomodel_;

    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();                   //Number of dipoles
    P = (*model).get_P();
    R = (*Core).get_R();
    Vector3d n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    cout << "lam" << lam << endl;
    K = 2 * M_PI / lam;
    

    for (int i = 0; i <= int(round(Paralength / 3) - 1); i++) {
        Vector3d n_K_tmp;
        n_K_tmp(0) = FOMParameters(3 * i);
        n_K_tmp(1) = FOMParameters(3 * i + 1);
        n_K_tmp(2) = FOMParameters(3 * i + 2);
        n_K_s_l.push_back(n_K_tmp);

        Matrix3d FconstM;
        double nkx = n_K_tmp(0);
        double nky = n_K_tmp(1);
        double nkz = n_K_tmp(2);
        double K3 = pow(K, 3);
        FconstM(0, 0) = K3 * (1 - nkx * nkx);
        FconstM(0, 1) = -K3 * nkx * nky;
        FconstM(0, 2) = -K3 * nkx * nkz;
        FconstM(1, 1) = K3 * (1 - nky * nky);
        FconstM(1, 2) = -K3 * nky * nkz;
        FconstM(2, 2) = K3 * (1 - nkz * nkz);
        FconstM(1, 0) = FconstM(0, 1);
        FconstM(2, 0) = FconstM(0, 2);
        FconstM(2, 1) = FconstM(1, 2);
        FconstM_l.push_back(FconstM);

        Vector3cd PSum_tmp;
        PSum_tmp = Vector3cd::Zero();
        PSum_l.push_back(PSum_tmp);
        
    }

    



}

void Objectivescattering0D::SingleResponse(int idx, bool deduction) {
    //list<Matrix3d>::iterator it1 = (FconstM_l).begin();
    list<Vector3d>::iterator it2 = (n_K_s_l).begin();
    list<Vector3cd>::iterator it3 = (PSum_l).begin();
    for (int i = 0; i <= int(round(Paralength / 3) - 1); i++) {
        //Matrix3d FconstM = (*it1);
        Vector3d n_K_s = (*it2);
        double nkx = n_K_s(0);
        double nky = n_K_s(1);
        double nkz = n_K_s(2);
        double phaseterm = -d * K * (nkx * ((*R)(3 * idx) - 32.5) + nky * ((*R)(3 * idx + 1) - 32.5) + nkz * ((*R)(3 * idx + 2) - 32.5));  //From equation 17. Time term will be eliminated
        complex<double> phase = cos(phaseterm) + sin(phaseterm) * 1i;
        if (deduction == false) {
            (*it3)(0) += (*P)(3 * idx) * phase;
            (*it3)(1) += (*P)(3 * idx + 1) * phase;
            (*it3)(2) += (*P)(3 * idx + 2) * phase;

        }
        else {
            (*it3)(0) -= (*P)(3 * idx) * phase;
            (*it3)(1) -= (*P)(3 * idx + 1) * phase;
            (*it3)(2) -= (*P)(3 * idx + 2) * phase;
        }
        it2++;
        it3++;
    }
}

double Objectivescattering0D::GroupResponse() {
    if (Have_Penalty) {
        return log10(this->FTUCnsquare() / (pow(K, 2) * pow(E0, 2))) - (*evomodel).L1Norm();
    }
    else {
        return log10(this->FTUCnsquare() / (pow(K, 2) * pow(E0, 2)));                           //K does not depend on scattering angle, so it is fine to divide it here.
    }

}

double Objectivescattering0D::GetVal() {
    Reset();
    for (int idx = 0; idx < N; idx++) {
        SingleResponse(idx, false);
    }
    return GroupResponse();
}

void Objectivescattering0D::Reset() {
    for (list<Vector3cd>::iterator it = PSum_l.begin(); it != PSum_l.end(); it++) {
        (*it) = Vector3cd::Zero();
    }
}


double Objectivescattering0D::FTUCnsquare() {

    list<Matrix3d>::iterator it1 = (FconstM_l).begin();
    list<Vector3cd>::iterator it2 = (PSum_l).begin();
    double result = 0.0;
    for (int i = 0; i <= int(round(Paralength / 3) - 1); i++) {
        Vector3cd FTUC;
        FTUC(0) = (*it1)(0, 0) * (*it2)(0) + (*it1)(0, 1) * (*it2)(1) + (*it1)(0, 2) * (*it2)(2);
        FTUC(1) = (*it1)(1, 0) * (*it2)(0) + (*it1)(1, 1) * (*it2)(1) + (*it1)(1, 2) * (*it2)(2);
        FTUC(2) = (*it1)(2, 0) * (*it2)(0) + (*it1)(2, 1) * (*it2)(1) + (*it1)(2, 2) * (*it2)(2);


        it1++;
        it2++;
        result += norm(FTUC(0)) + norm(FTUC(1)) + norm(FTUC(2));                                 //In C++, norm is the square of magnitude.
        

    }

    return result/double(round(Paralength / 3));
    
}

/*
Objectivescattering2D::Objectivescattering2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_) {
    VectorXd PointEParameters = VectorXd::Zero((parameters).size());
    list<double>::iterator it = (parameters).begin();
    for (int i = 0; i <= int((parameters).size() - 1); i++) {
        PointEParameters(i) = (*it);
        it++;
    }
    Have_Penalty = HavePenalty_;
    x = PointEParameters(0);
    y = PointEParameters(1);
    z = PointEParameters(2);
    Have_Devx = false;
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
    cout << "lam" << lam << endl;
    double K = 2 * M_PI / lam;
    E_sum = Vector3cd::Zero();
    E_ext = Vector3cd::Zero();
    E_ext(0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
    E_ext(1) = E0 * n_E0(1) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
    E_ext(2) = E0 * n_E0(2) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
    // cout << R(3*5444+2) << "this" << endl;
}

void Objectivescattering2D::SingleResponse(int idx, bool deduction) {
}

double Objectivescattering2D::GroupResponse() {
}

double Objectivescattering2D::GetVal() {
    Reset();
    for (int idx = 0; idx < N; idx++) {
        SingleResponse(idx, false);
        // cout << E_sum(0) << endl;
    }
    return GroupResponse();
}

void Objectivescattering2D::Reset() {
    E_sum(0) = E_ext(0);
    E_sum(1) = E_ext(1);
    E_sum(2) = E_ext(2);
}

double Objectivescattering2D::FTUCnsquare() {
    double 
}
*/









