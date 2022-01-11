#include "definition.h"
#define NUM_THREADS 6



DDAModel::DDAModel(AProductCore* AProductCore_, Vector3d n_K_, double E0_, Vector3d n_E0_) {
    
    Core = AProductCore_;
    time=0;
    ITERATION=0;
    Error=0.0;
    E0=E0_;
    n_K=n_K_;
    n_E0=n_E0_;

    //cout << "lam" << (*Core).get_lam() << endl;

    //cout << "E0=" << E0 << endl;
    //cout << "n_K" << endl << n_K << endl;
    //cout << "n_E0" << endl << n_E0 << endl;
    int N = (*Core).get_N();
    int Nx = (*Core).get_Nx();
    int Ny = (*Core).get_Ny();
    int Nz = (*Core).get_Nz();
    double lam = (*Core).get_lam();
    double K = 2 * M_PI / lam;
    double d = (*Core).get_d();
    VectorXi* R = (*Core).get_R();
    VectorXd* diel_old = (*((*Core).get_CStr())).get_diel_old();
    VectorXcd* material = (*Core).get_material();

    RResultSwitch = false;
    RResult = *R;
    
    P = VectorXcd::Zero(N*3);
    P_max = P;
    E = VectorXcd::Zero(N*3);
    Einternal = VectorXcd::Zero(N*3);
    EResult = VectorXcd::Zero(N*3);
    cout<<"fuck"<<endl;
    for (int i=0;i<N;i++) {
        E(3 * i) = E0 * n_E0(0) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
        E(3 * i + 1) = E0 * n_E0(1) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
        E(3 * i + 2) = E0 * n_E0(2) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
    }
    al = VectorXcd::Zero(N*3);
    diel = VectorXcd::Zero(N * 3);
    cout<<"fuck"<<endl;
    cout<<"al size: "<<al.size()<<endl;
    cout<<"diel size: "<<diel.size()<<endl;
    cout<<"diel_old size: "<<diel_old.size()<<endl;
    for (int i = 0; i < N * 3; i++) {
        cout<<i<<endl;
        int labelfloor = int(floor((*diel_old)(i)));
        std::complex<double> diel_tmp = (*material)(labelfloor) + ((*diel_old)(i)-double(labelfloor)) * ((*material)(labelfloor + 1) - (*material)(labelfloor));
        diel(i) = diel_tmp;
        al(i) = 1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K);
    }  
    cout<<"fuck"<<endl;
    //cout << "al" << al(0) << endl;

    al_max = al;
    verbose = true;
}

DDAModel::DDAModel(AProductCore* AProductCore_, Vector3d n_K_, double E0_, Vector3d n_E0_, VectorXi* RResult_) {
    
    Core = AProductCore_;
    time = 0;
    ITERATION = 0;
    Error = 0.0;
    E0 = E0_;
    n_K = n_K_;
    n_E0 = n_E0_;

    cout << "E0=" << E0 << endl;
    cout << "n_K" << n_K << endl;
    cout << "n_E0" << n_E0 << endl;
    int N = (*Core).get_N();
    int Nx = (*Core).get_Nx();
    int Ny = (*Core).get_Ny();
    int Nz = (*Core).get_Nz();
    double lam = (*Core).get_lam();
    double K = 2 * M_PI / lam;
    double d = (*Core).get_d();
    VectorXi* R = (*Core).get_R();
    VectorXd* diel_old = (*((*Core).get_CStr())).get_diel_old();
    VectorXcd* material = (*Core).get_material();
    RResultSwitch = true;
    RResult = *RResult_;
    if(RResult.size()%3 != 0){
    
        cout<<"RResult does not have a size with 3*integer"<<endl;
        // This could be changed to throwing an exception.
    }
    
    

    P = VectorXcd::Zero(N*3);
    P_max = P;
    E = VectorXcd::Zero(N*3);
    Einternal = VectorXcd::Zero(N*3);
    EResult = VectorXcd::Zero(RResult.size());
    for (int i = 0; i < N; i++) {
        E(3 * i) = E0 * n_E0(0) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
        E(3 * i + 1) = E0 * n_E0(1) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
        E(3 * i + 2) = E0 * n_E0(2) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
    }
    al = VectorXcd::Zero(N * 3);
    diel = VectorXcd::Zero(N * 3);
    for (int i = 0; i < N * 3; i++) {
        int labelfloor = int(floor((*diel_old)(i)));
        std::complex<double> diel_tmp = (*material)(labelfloor) + ((*diel_old)(i) - double(labelfloor)) * ((*material)(labelfloor + 1) - (*material)(labelfloor));
        diel(i) = diel_tmp;
        al(i) = 1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K);
    }
    al_max = al;
    verbose = true;
    
}

void DDAModel::bicgstab(int MAX_ITERATION,double MAX_ERROR){
    if (verbose) {
        cout << "--------------Calculation start. Iterative method used: BICGSTAB---------------" << endl;
        cout << endl;
    }
    int N = (*Core).get_N();
    
    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    //fftw_init_threads();                                                       //////////Initialize the multi-thread
    //fftw_plan_with_nthreads(NUM_THREADS);
    //cout<<"Threads"<<NUM_THREADS<<endl;;
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
    
    
    //Always starts with P=0 to avoid strange behaviour
    //P = VectorXcd::Zero(N * 3);



    VectorXcd Ax0 = Aproductwithalb(P);
    
    r = E-Ax0;
    r0 = r;
    p = r;
    VectorXcd Ap0 = Aproductwithalb(p);


    alpha = r0.dot(r)/r0.dot(Ap0);
    t = r-alpha*Ap0;
    VectorXcd At0 = Aproductwithalb(t);


    zeta = At0.dot(t)/At0.dot(At0);
    u = zeta*Ap0;
    z = zeta*r-alpha*u;
    P = P+alpha*p+z;                                    //this will directly change P in this.
    rl = r;
    r = t-zeta*At0;

    VectorXcd Ap = VectorXcd::Zero(N*3);
    VectorXcd At = VectorXcd::Zero(N*3);
    for (int it=0;it<=MAX_ITERATION-1;it++) {
        if (verbose && (it+1)%1000==0) {
            //cout << "r.norm(): " << r.norm() << endl;
            //cout << "E.norm(): " << E.norm() << endl;
            cout << "                Iter " << it+1 << ", error=" << r.norm()/E.norm() << " MAX_ERROR="<<MAX_ERROR<<endl;
        }

        complex<double> r0dotrl = r0.dot(rl);
        if (r0dotrl.real() == 0 && r0dotrl.imag() == 0) {
            Error = r.norm() / E.norm();
            //cout << "r.norm(): " << r.norm() << endl;
            //cout << "E.norm(): " << E.norm() << endl;
            ITERATION = it + 1;
            if (verbose) {
                high_resolution_clock::time_point t_end = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(t_end - t_start).count();
                time = duration_cast<milliseconds>(t_end - t_start).count();
                cout << "----------------------------r0dotrl==(0,0), nan is going to occur so stop now------------------------" << endl;
                cout << "--------------Calculation finished. Duration: " << duration / 1000.0 << "s.-------------------" << endl;
                //ofstream fout;
                //fout.open("DDATime.txt", fstream::app);
                //fout<<N<<" "<<duration/1000.0<<endl;
                //fout << ITERATION << " " << duration / 1000.0 / ITERATION << endl;
                //fout.close();

                cout << "              Error: " << Error << endl;
                cout << "              Iteration: " << ITERATION << endl;
                cout << endl;
            }
            return;
        }

        beta = (alpha/zeta)*r0.dot(r)/r0.dot(rl);
        p = r+beta*(p-u);
        Ap = Aproductwithalb(p);
        alpha = r0.dot(r)/r0.dot(Ap);
        t = r-alpha*Ap;
        At = Aproductwithalb(t);

        zeta = At.dot(t)/At.dot(At);
        u = zeta*Ap;
        z = zeta*r-alpha*u;
        P = P+alpha*p+z;
        rl = r;
        r = t-zeta*At;

        if (r.norm()/E.norm()<=MAX_ERROR) {
            Error = r.norm()/E.norm();
            cout << "r.norm(): " << r.norm() << endl;
            cout << "E.norm(): " << E.norm() << endl;
            ITERATION = it+1;
            if (verbose) {
                high_resolution_clock::time_point t_end = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(t_end-t_start).count();
                time = duration_cast<milliseconds>(t_end-t_start).count();
                cout << "--------------Calculation finished. Duration: " << duration/1000.0 << "s.-------------------" << endl;
                //ofstream fout;
                //fout.open("DDATime.txt", fstream::app);
                //fout<<N<<" "<<duration/1000.0<<endl;
                //fout << ITERATION << " " << duration / 1000.0 / ITERATION << endl;
                //fout.close();

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

void DDAModel::bicgstab(int MAX_ITERATION, double MAX_ERROR, int EVOITERATION) {
    if (verbose) {
        cout << "--------------Calculation start. Iterative method used: BICGSTAB---------------" << endl;
        cout << endl;
    }
    int N = (*Core).get_N();

    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    //fftw_init_threads();                                                       //////////Initialize the multi-thread
    //fftw_plan_with_nthreads(NUM_THREADS);
    //cout<<"Threads"<<NUM_THREADS<<endl;;
    VectorXcd p = VectorXcd::Zero(N * 3);
    VectorXcd t = VectorXcd::Zero(N * 3);
    VectorXcd w = VectorXcd::Zero(N * 3);
    VectorXcd r = VectorXcd::Zero(N * 3);
    VectorXcd r0 = VectorXcd::Zero(N * 3);
    VectorXcd rl = VectorXcd::Zero(N * 3);
    VectorXcd y = VectorXcd::Zero(N * 3);
    VectorXcd u = VectorXcd::Zero(N * 3);
    VectorXcd z = VectorXcd::Zero(N * 3);
    VectorXcd x = VectorXcd::Zero(N * 3);
    std::complex<double> alpha;
    std::complex<double> beta;
    std::complex<double> eta;
    std::complex<double> zeta;


    //Always starts with P=0 to avoid strange behaviour
    //P = VectorXcd::Zero(N * 3);

    ofstream foutnew(".\\p330-lam542-beta8-TiO2-InE-circle-fordebug\\BUGINFO.txt");

    

    VectorXcd Ax0 = Aproductwithalb(P);



    r = E - Ax0;



    r0 = r;
    p = r;
    VectorXcd Ap0 = Aproductwithalb(p);



    alpha = r0.dot(r) / r0.dot(Ap0);



    t = r - alpha * Ap0;
    VectorXcd At0 = Aproductwithalb(t);



    zeta = At0.dot(t) / At0.dot(At0);



    u = zeta * Ap0;
    z = zeta * r - alpha * u;
    P = P + alpha * p + z;                                    //this will directly change P in this.
    rl = r;
    r = t - zeta * At0;

    



    VectorXcd Ap = VectorXcd::Zero(N * 3);
    VectorXcd At = VectorXcd::Zero(N * 3);
    for (int it = 0; it <= MAX_ITERATION - 1; it++) {
        if (verbose && (it + 1) % 1000 == 0) {
            cout << "r.norm(): " << r.norm() << endl;
            cout << "E.norm(): " << E.norm() << endl;
            cout << "                Iter " << it + 1 << ", error=" << r.norm() / E.norm() << " MAX_ERROR=" << MAX_ERROR << endl;
        }
        beta = (alpha / zeta) * r0.dot(r) / r0.dot(rl);

        /*
        if ((EVOITERATION == 88) && (it == 4849)) {
            //cout << "WE FIND 88" << endl;
            foutnew << "-----------------beta at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << beta << endl;

            //foutnew << "-----------------p at " << "it" + to_string(it) << "---------------" << endl;
            //foutnew << p << endl;

            //foutnew << "-----------------Ap at " << "it" + to_string(it) << "---------------" << endl;
            //foutnew << Ap << endl;

            foutnew << "-----------------alpha at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << alpha << endl;

            //foutnew << "-----------------t at " << "it" + to_string(it) << "---------------" << endl;
            //foutnew << t << endl;

            //foutnew << "-----------------At at " << "it" + to_string(it) << "---------------" << endl;
            //foutnew << At << endl;

            foutnew << "-----------------zeta at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << zeta << endl;

            foutnew << "-----------------alpha/zeta at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << alpha/zeta << endl;

            foutnew << "-----------------r1 at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << rl << endl;

            foutnew << "-----------------r0 at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << r0 << endl;

            //foutnew << "-----------------u at " << "it" + to_string(it) << "---------------" << endl;
            //foutnew << u << endl;

            //foutnew << "-----------------z at " << "it" + to_string(it) << "---------------" << endl;
            //foutnew << z << endl;

            //foutnew << "-----------------P at " << "it" + to_string(it) << "---------------" << endl;
            //foutnew << P << endl;

            foutnew << "-----------------r at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << r << endl;

            foutnew << "-----------------r0dotr at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << r0.dot(r) << endl;

            foutnew << "-----------------r0dotrl at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << r0.dot(rl) << endl;

            foutnew << "-----------------r0dotr/r0dotrl at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << r0.dot(r) / r0.dot(rl) << endl;
        }
        */
        p = r + beta * (p - u);
        Ap = Aproductwithalb(p);
        alpha = r0.dot(r) / r0.dot(Ap);
        t = r - alpha * Ap;
        At = Aproductwithalb(t);

        zeta = At.dot(t) / At.dot(At);
        u = zeta * Ap;
        z = zeta * r - alpha * u;
        P = P + alpha * p + z;
        rl = r;
        r = t - zeta * At;
        

        if (EVOITERATION == 88) {
            //cout << "WE FIND 88" << endl;
            foutnew << "-----------------r at " << "it" + to_string(it) << "---------------" << endl;
            foutnew << r.norm() << endl;
        }

        if (r.norm() / E.norm() <= MAX_ERROR) {
            Error = r.norm() / E.norm();
            cout << "r.norm(): " << r.norm() << endl;
            cout << "E.norm(): " << E.norm() << endl;
            ITERATION = it + 1;
            if (verbose) {
                high_resolution_clock::time_point t_end = high_resolution_clock::now();
                auto duration = duration_cast<milliseconds>(t_end - t_start).count();
                time = duration_cast<milliseconds>(t_end - t_start).count();
                cout << "--------------Calculation finished. Duration: " << duration / 1000.0 << "s.-------------------" << endl;
                //ofstream fout;
                //fout.open("DDATime.txt", fstream::app);
                //fout<<N<<" "<<duration/1000.0<<endl;
                //fout << ITERATION << " " << duration / 1000.0 / ITERATION << endl;
                //fout.close();

                cout << "              Error: " << Error << endl;
                cout << "              Iteration: " << ITERATION << endl;
                cout << endl;
            }
            return;
        }
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    time = duration_cast<milliseconds>(t_end - t_start).count();
    cout << "                ERROR:does not converge in " << MAX_ITERATION << " iterations" << endl;

    foutnew.close();


    return;
}

void DDAModel::change_E(VectorXcd E_){
    E=E_;
}

void DDAModel::reset_E(){
    int N = (*Core).get_N();
    double lam = (*Core).get_lam();
    double K = 2 * M_PI / lam;
    double d = (*Core).get_d();
    VectorXi* R = (*Core).get_R();
    for (int i=0;i<N;i++) {
        E(3 * i) = E0 * n_E0(0) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
        E(3 * i + 1) = E0 * n_E0(1) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
        E(3 * i + 2) = E0 * n_E0(2) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
    }
}

void DDAModel::UpdateAlpha() {

    VectorXd* diel_old = (*((*Core).get_CStr())).get_diel_old();
    VectorXcd* material = (*Core).get_material();
    double lam = (*Core).get_lam();
    double K = 2 * M_PI / lam;
    double d = (*Core).get_d();
    for (int i = 0; i <= (*diel_old).size() - 1; i++) {
        int labelfloor = int(floor((*diel_old)(i)));
        std::complex<double> diel_tmp = (*material)(labelfloor) + ((*diel_old)(i) - double(labelfloor)) * ((*material)(labelfloor + 1) - (*material)(labelfloor));
        diel(i) = diel_tmp;
        al(i) = 1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K);
    }
}

void DDAModel::UpdateAlphaSingle(int idx) {

    VectorXd* diel_old = (*((*Core).get_CStr())).get_diel_old();
    VectorXcd* material = (*Core).get_material();
    double lam = (*Core).get_lam();
    double K = 2 * M_PI / lam;
    double d = (*Core).get_d();

    /*
    CoreStructure* CStr = (*Core).get_CStr();
    SpacePara* spacepara = ((*CStr).get_spacepara());
    vector<list<int>>* Paratogeometry = (*spacepara).get_Paratogeometry();

    list<int>::iterator it = (*Paratogeometry)[idx].begin();
    for (int j = 0; j <= (*Paratogeometry)[idx].size() - 1; j++) {
        int position = *it;
        std::complex<double> diel_tmp = (*material)(0) + (*diel_old)(3 * position) * ((*material)(1) - (*material)(0));
        diel(3*position) = diel_tmp;
        diel(3 * position + 1) = diel_tmp;
        diel(3 * position + 2) = diel_tmp;
        al(3*position) = 1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K);
        al(3 * position + 1) = al(3 * position);
        al(3 * position + 2) = al(3 * position);
        it++;
    }
    */
    int labelfloor = int(floor((*diel_old)(3 * idx)));
    std::complex<double> diel_tmp = (*material)(labelfloor) + ((*diel_old)(3 * idx) - double(labelfloor)) * ((*material)(labelfloor + 1) - (*material)(labelfloor));
    diel(3 * idx) = diel_tmp;
    diel(3 * idx + 1) = diel_tmp;
    diel(3 * idx + 2) = diel_tmp;
    al(3 * idx) = 1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K);
    al(3 * idx + 1) = al(3 * idx);
    al(3 * idx + 2) = al(3 * idx);

}

void DDAModel::solve_E(){
    if(RResultSwitch == true){
        int N = (*Core).get_N();
        double d = (*Core).get_d();
        double lam = (*Core).get_lam();
        double K = 2 * M_PI / lam;
        VectorXi* R = (*Core).get_R();
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
                double rx = x - d * (*R)(3 * i);                  //R has no d in it, so needs to time d
                double ry = y - d * (*R)(3 * i + 1);
                double rz = z - d * (*R)(3 * i + 2);
                Matrix3cd A=(*Core).A_dic_generator(rx,ry,rz);
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

void DDAModel::update_E_in_structure(){
    int N = (*Core).get_N();
    
    for(int i=0;i<=3*N-1;i++){
        std::complex<double> lorentzfactor = 2.0 + diel(i);
        lorentzfactor = lorentzfactor / 3.0;
        Einternal(i) = al(i) * P(i)/lorentzfactor;
    }
}

VectorXcd DDAModel::Aproductwithalb(VectorXcd& b) {
    if (b.size() != al.size()) {
        cout << "In Aproductwithalb, bsize not equal to alsize" << endl;
    }
    VectorXcd result = VectorXcd::Zero(b.size());
    for (int i = 0; i <= al.size() - 1; i++) {
        result(i) = b(i) * al(i);
    }
    return (*Core).Aproduct(b) + result;
}


void DDAModel::output_to_file(){
    
    (*((*Core).get_CStr())).output_to_file();
    ofstream fout("Model_results.txt");
    for(int i=0;i<=P.size()-1;i++){
        if(P(i).imag()<0){
            fout<<P(i).real()<<P(i).imag()<<"j"<<endl;
        }
        else{
            fout<<P(i).real()<<"+"<<P(i).imag()<<"j"<<endl;
        }
        
    }
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

void DDAModel::output_to_file(string save_position, int iteration, int ModelLabel){
    
    string name;

    name = save_position + "Model_output\\" + "Model_results" + to_string(ModelLabel) + "it" + to_string(iteration) + ".txt";
    ofstream fout(name);
    for(int i=0;i<=P.size()-1;i++){
        if(P(i).imag()<0){
            fout<<P(i).real()<<P(i).imag()<<"j"<<endl;
        }
        else{
            fout<<P(i).real()<<"+"<<P(i).imag()<<"j"<<endl;
        }
        
    }
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

void DDAModel::output_to_file(string save_position, int iteration) {

    string name;
    name = save_position + "Model_results" + "it" + to_string(iteration) + ".txt";
    //name = save_position + "Model_output_verify\\" + "Model_results" + "it" + to_string(iteration) + ".txt";
    ofstream fout(name);
    /*
    for (int i = 0; i <= P.size() - 1; i++) {
        if (P(i).imag() < 0) {
            fout << P(i).real() << P(i).imag() << "j" << endl;
        }
        else {
            fout << P(i).real() << "+" << P(i).imag() << "j" << endl;
        }

    }
    */
    for (int i = 0; i <= EResult.size() - 1; i++) {
        if (EResult(i).imag() < 0) {
            fout << EResult(i).real() << EResult(i).imag() << "j" << endl;
        }
        else {
            fout << EResult(i).real() << "+" << EResult(i).imag() << "j" << endl;
        }

    }
    
    for (int i = 0; i <= P.size() - 1; i++) {
        if (P(i).imag() < 0) {
            fout << P(i).real() << P(i).imag() << "j" << endl;
        }
        else {
            fout << P(i).real() << "+" << P(i).imag() << "j" << endl;
        }

    }
    
    fout.close();
}

void DDAModel::output_to_file(string save_position, double wavelength, int iteration) {

    string name;

    name = save_position + "Model_output" + to_string(int(wavelength)) + "\\Model_results" + "it" + to_string(iteration) + ".txt";
    ofstream fout(name);
    /*
    for (int i = 0; i <= P.size() - 1; i++) {
        if (P(i).imag() < 0) {
            fout << P(i).real() << P(i).imag() << "j" << endl;
        }
        else {
            fout << P(i).real() << "+" << P(i).imag() << "j" << endl;
        }

    }
    */
    for (int i = 0; i <= EResult.size() - 1; i++) {
        if (EResult(i).imag() < 0) {
            fout << EResult(i).real() << EResult(i).imag() << "j" << endl;
        }
        else {
            fout << EResult(i).real() << "+" << EResult(i).imag() << "j" << endl;
        }

    }
    fout.close();
}

void DDAModel::InitializeP(VectorXcd& Initializer) {
    P = Initializer;
}

VectorXcd* DDAModel::get_P() {
    return &P;
}

Vector3d DDAModel::get_nE0() {
    return n_E0;
}

Vector3d DDAModel::get_nK() {
    return n_K;
}

double DDAModel::get_E0() {
    return E0;
}

VectorXcd* DDAModel::get_Einternal() {
    return &Einternal;
}

AProductCore* DDAModel::get_Core() {
    return Core;
}

VectorXcd* DDAModel::get_al() {
    return &al;
}

VectorXcd* DDAModel::get_al_max() {
    return &al_max;
}

VectorXcd* DDAModel::get_P_max() {
    return &P_max;
}

int DDAModel::get_ITERATION() {
    return ITERATION;
}

//----------------------------------------heritage from AProductCore------------------------------------

int DDAModel::get_N() {
    return (*Core).get_N();
}

int DDAModel::get_Nx() {
    return (*Core).get_Nx();
}

int DDAModel::get_Ny(){
    return (*Core).get_Ny();
}

int DDAModel::get_Nz() {
    return (*Core).get_Nz();
}

VectorXi* DDAModel::get_R() {
    return (*Core).get_R();
}

double DDAModel::get_d() {
    return (*Core).get_d();
}


double DDAModel::get_lam() {
    return (*Core).get_lam();
}

VectorXd* DDAModel::get_diel_old() {
    return (*Core).get_diel_old();
}

VectorXcd* DDAModel::get_material() {
    return (*Core).get_material();
}

VectorXd* DDAModel::get_diel_old_max() {
    return (*Core).get_diel_old_max();
}

SpacePara* DDAModel::get_spacepara() {
    return (*((*Core).get_CStr())).get_spacepara();
}
