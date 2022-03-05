FOMscattering2D::FOMscattering2D(list<double> parameters, DDAModel* model_) {
    Paralength = (parameters).size();
    FOMParameters = VectorXd::Zero(Paralength);
    list<double>::iterator it = (parameters).begin();
    for (int i = 0; i <= int(Paralength - 1); i++) {
        FOMParameters(i) = (*it);
        it++;
    }
    if (Paralength % 2 != 0){
        cout << "FOMscattering2D ERROR: parameter must be times of 2." << endl;
    }

    
    

    model = model_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();                   //Number of dipoles
    P = (*model).get_P();
    R = (*Core).get_R();
    n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    cout << "lam" << lam << endl;
    K = (*Core).get_K();

    double Lm = (*Core).get_Lm();
    double Ln = (*Core).get_Ln();
    ATUC = Lm * Ln;
    
    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
        Vector3d n_K_tmp;

        n_K_tmp(0) = 2 * M_PI * FOMParameters(2 * i) / (Lm * K) + n_K(0);
        n_K_tmp(1) = 2 * M_PI * FOMParameters(2 * i + 1) / (Ln * K) + n_K(1);
        n_K_tmp(2) = sqrt(1 - pow(n_K_tmp(0), 2) - pow(n_K_tmp(1), 2));
        n_K_s_l.push_back(n_K_tmp);
    }
}

list<double> FOMscattering2D::GetVal() {
    list<double> result;
    int listlength = (n_K_s_l).size();
    list<Vector3d>::iterator it = n_K_s_l.begin();
    
    for (int i = 0; i <= listlength - 1; i++) {
        if (FOMParameters(2 * i) == 0 && FOMParameters(2 * i + 1) == 0) {
            cout << "There is a normal output" << endl;
            Vector3cd tmp = this->FTUC((*it));
            double ksz = (*it)(2);                                 //x in the paper is actually z in our case.
            //cout << tmp << endl;
            tmp = tmp * (2.0 * M_PI * 1.0i) / (K * K * ATUC * abs(ksz));
            //cout << tmp << endl;
            //tmp = tmp + E0 * n_E0;
            //cout << tmp << endl;
            double tmpresult = (norm(tmp(0)) + norm(tmp(1)) + norm(tmp(2))) / pow(E0, 2);  //00 order transmission
            //cout << tmpresult << endl;
            result.push_back(tmpresult);
            cout << tmpresult << endl;
        }  
        else {
            Vector3cd tmp = this->FTUC((*it));
            double ksz = (*it)(2);                                 //x in the paper is actually z in our case.
            tmp = tmp * (2.0 * M_PI * 1.0i) / (K * K * ATUC * abs(ksz));
            double tmpresult = (norm(tmp(0)) + norm(tmp(1)) + norm(tmp(2))) / pow(E0, 2);  //00 order transmission
            result.push_back(tmpresult);
            
        } 
        it++;
    }

    return result;
}

Vector3cd FOMscattering2D::FTUC(Vector3d n_K_s) {
    Matrix3d FconstM;
    double nkx = n_K_s(0);
    double nky = n_K_s(1);
    double nkz = n_K_s(2);
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

    Vector3cd PSum;
    PSum = Vector3cd::Zero();
    for (int i = 0; i <= N-1; i++) {
        double phaseterm = -d * K * (nkx * (*R)(3 * i) + nky * (*R)(3 * i + 1) + nkz * (*R)(3 * i + 2));  //From equation 17. Time term will be eliminated
        complex<double> phase = cos(phaseterm) + sin(phaseterm) * 1i;
        PSum(0) += (*P)(3 * i) * phase;
        PSum(1) += (*P)(3 * i + 1) * phase;
        PSum(2) += (*P)(3 * i + 2) * phase;
    }

    Vector3cd FTUC;
    FTUC(0) = FconstM(0, 0) * PSum(0) + FconstM(0, 1) * PSum(1) + FconstM(0, 2) * PSum(2);
    FTUC(1) = FconstM(1, 0) * PSum(0) + FconstM(1, 1) * PSum(1) + FconstM(1, 2) * PSum(2);
    FTUC(2) = FconstM(2, 0) * PSum(0) + FconstM(2, 1) * PSum(1) + FconstM(2, 2) * PSum(2);

    //double result = norm(FTUC(0)) + norm(FTUC(1)) + norm(FTUC(2));                                 //In C++, norm is the square of magnitude.
    return FTUC;
}


FOMreflect2D::FOMreflect2D(list<double> parameters, DDAModel* model_) {
    Paralength = (parameters).size();
    FOMParameters = VectorXd::Zero(Paralength);
    list<double>::iterator it = (parameters).begin();
    for (int i = 0; i <= int(Paralength - 1); i++) {
        FOMParameters(i) = (*it);
        it++;
    }
    if (Paralength % 2 != 0) {
        cout << "FOMscattering2D ERROR: parameter must be times of 2." << endl;
    }




    model = model_;
    AProductCore* Core = (*model).get_Core();
    d = (*Core).get_d();
    N = (*Core).get_N();                   //Number of dipoles
    P = (*model).get_P();
    R = (*Core).get_R();
    n_E0 = (*model).get_nE0();
    Vector3d n_K = (*model).get_nK();
    E0 = (*model).get_E0();
    double lam = (*Core).get_lam();
    cout << "lam" << lam << endl;
    K = (*Core).get_K();

    double Lm = (*Core).get_Lm();
    double Ln = (*Core).get_Ln();
    ATUC = Lm * Ln;

    for (int i = 0; i <= int(round(Paralength / 2) - 1); i++) {
        Vector3d n_K_tmp;

        n_K_tmp(0) = 2 * M_PI * FOMParameters(2 * i) / (Lm * K) + n_K(0);
        n_K_tmp(1) = 2 * M_PI * FOMParameters(2 * i + 1) / (Ln * K) + n_K(1);
        n_K_tmp(2) = -sqrt(1 - pow(n_K_tmp(0), 2) - pow(n_K_tmp(1), 2));
        n_K_s_l.push_back(n_K_tmp);
    }
}

list<double> FOMreflect2D::GetVal() {
    list<double> result;
    int listlength = (n_K_s_l).size();
    list<Vector3d>::iterator it = n_K_s_l.begin();

    for (int i = 0; i <= listlength - 1; i++) {
        if (FOMParameters(2 * i) == 0 && FOMParameters(2 * i + 1) == 0) {
            cout << "There is a normal output" << endl;
            Vector3cd tmp = this->FTUC((*it));
            double ksz = (*it)(2);                                 //x in the paper is actually z in our case.
            //cout << tmp << endl;
            tmp = tmp * (2.0 * M_PI * 1.0i) / (K * K * ATUC * abs(ksz));
            tmp = tmp + E0 * n_E0;
            //cout << tmp << endl;
            //cout << tmp << endl;
            double tmpresult = (norm(tmp(0)) + norm(tmp(1)) + norm(tmp(2))) / pow(E0, 2);  //00 order transmission
            //cout << tmpresult << endl;
            result.push_back(tmpresult);
            cout << tmpresult << endl;
        }
        else {
            Vector3cd tmp = this->FTUC((*it));
            double ksz = (*it)(2);                                 //x in the paper is actually z in our case.
            tmp = tmp * (2.0 * M_PI * 1.0i) / (K * K * ATUC * abs(ksz));
            double tmpresult = (norm(tmp(0)) + norm(tmp(1)) + norm(tmp(2))) / pow(E0, 2);  //00 order transmission
            result.push_back(tmpresult);

        }
        it++;
    }

    return result;
}

Vector3cd FOMreflect2D::FTUC(Vector3d n_K_s) {
    Matrix3d FconstM;
    double nkx = n_K_s(0);
    double nky = n_K_s(1);
    double nkz = n_K_s(2);
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

    Vector3cd PSum;
    PSum = Vector3cd::Zero();
    for (int i = 0; i <= N - 1; i++) {
        double phaseterm = -d * K * (nkx * (*R)(3 * i) + nky * (*R)(3 * i + 1) + nkz * (*R)(3 * i + 2));  //From equation 17. Time term will be eliminated
        complex<double> phase = cos(phaseterm) + sin(phaseterm) * 1i;
        PSum(0) += (*P)(3 * i) * phase;
        PSum(1) += (*P)(3 * i + 1) * phase;
        PSum(2) += (*P)(3 * i + 2) * phase;
    }

    Vector3cd FTUC;
    FTUC(0) = FconstM(0, 0) * PSum(0) + FconstM(0, 1) * PSum(1) + FconstM(0, 2) * PSum(2);
    FTUC(1) = FconstM(1, 0) * PSum(0) + FconstM(1, 1) * PSum(1) + FconstM(1, 2) * PSum(2);
    FTUC(2) = FconstM(2, 0) * PSum(0) + FconstM(2, 1) * PSum(1) + FconstM(2, 2) * PSum(2);

    //double result = norm(FTUC(0)) + norm(FTUC(1)) + norm(FTUC(2));                                 //In C++, norm is the square of magnitude.
    return FTUC;
}



FOMscattering0D::FOMscattering0D(list<double> parameters, DDAModel* model_) {
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


    for (int i = 0; i <= int(round(Paralength / 3) - 1); i++) {
        Vector3d n_K_tmp;
        n_K_tmp(0) = FOMParameters(3 * i);
        n_K_tmp(1) = FOMParameters(3 * i + 1);
        n_K_tmp(2) = FOMParameters(3 * i + 2);
        n_K_s_l.push_back(n_K_tmp);
    }

    model = model_;
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
    K = (*Core).get_K();

}

list<double> FOMscattering0D::GetVal() {
    list<double> result;
    int listlength = (n_K_s_l).size();
    list<Vector3d>::iterator it = n_K_s_l.begin();
    for (int i = 0; i <= listlength - 1; i++) {
        double tmp = this->FTUCnsquare((*it));
        double tmpresult = tmp / (pow(K, 2) * pow(E0, 2));
        result.push_back(tmpresult);
        //cout << "(" << (*it)(0) << " " << (*it)(1) << " " << (*it)(2) << ") " << tmp<<" "<<tmpresult<<endl;

        it++;


    }

    return result;
}

double FOMscattering0D::FTUCnsquare(Vector3d n_K_s) {
    Matrix3d FconstM;
    double nkx = n_K_s(0);
    double nky = n_K_s(1);
    double nkz = n_K_s(2);
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

    Vector3cd PSum;
    PSum = Vector3cd::Zero();
    for (int i = 0; i <= N-1; i++) {
        double phaseterm = -d * K * (nkx * ((*R)(3 * i) - 32.5) + nky * ((*R)(3 * i + 1) - 32.5) + nkz * ((*R)(3 * i + 2) - 32.5));  //From equation 17. Time term will be eliminated
        complex<double> phase = cos(phaseterm) + sin(phaseterm) * 1i;
        PSum(0) += (*P)(3 * i) * phase;
        PSum(1) += (*P)(3 * i + 1) * phase;
        PSum(2) += (*P)(3 * i + 2) * phase;

        
    }

    Vector3cd FTUC;
    FTUC(0) = FconstM(0, 0) * PSum(0) + FconstM(0, 1) * PSum(1) + FconstM(0, 2) * PSum(2);
    FTUC(1) = FconstM(1, 0) * PSum(0) + FconstM(1, 1) * PSum(1) + FconstM(1, 2) * PSum(2);
    FTUC(2) = FconstM(2, 0) * PSum(0) + FconstM(2, 1) * PSum(1) + FconstM(2, 2) * PSum(2);

    double result = norm(FTUC(0)) + norm(FTUC(1)) + norm(FTUC(2));                                 //In C++, norm is the square of magnitude.
    
    cout << "compare:" << FconstM << endl;
    cout << PSum << endl;
    return result;
}