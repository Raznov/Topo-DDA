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

    cout << "E0=" << E0 << endl;
    cout << "n_K" << endl << n_K << endl;
    cout << "n_E0" << endl << n_E0 << endl;
    int N = (*Core).get_N();
    int Nx = (*Core).get_Nx();
    int Ny = (*Core).get_Ny();
    int Nz = (*Core).get_Nz();
    double lam = (*Core).get_lam();
    double K = 2 * M_PI / lam;
    double d = (*Core).get_d();
    VectorXi* R = (*Core).get_R();
    VectorXcd* diel = (*Core).get_diel();

    RResultSwitch = false;
    RResult = *R;
    
    P = VectorXcd::Zero(N*3);
    P_store = P;
    P_max = P;
    E = VectorXcd::Zero(N*3);
    Einternal = VectorXcd::Zero(N*3);
    EResult = VectorXcd::Zero(N*3);
    for (int i=0;i<N;i++) {
        E(3 * i) = E0 * n_E0(0) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
        E(3 * i + 1) = E0 * n_E0(1) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
        E(3 * i + 2) = E0 * n_E0(2) * (cos(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) + sin(K * d * (n_K(0) * (*R)(3 * i) + n_K(1) * (*R)(3 * i + 1) + n_K(2) * (*R)(3 * i + 2))) * 1i);
    }
    al = VectorXcd::Zero(N*3);
    for (int i=0;i<N*3;i++) {
        al(i) = 1.0 / Get_Alpha(lam, K, d, (*diel)(i), n_E0, n_K);
    }  

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
    VectorXcd* diel = (*Core).get_diel();
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
    for (int i = 0; i < N * 3; i++) {
        al(i) = 1.0 / Get_Alpha(lam, K, d, (*diel)(i), n_E0, n_K);
    }
    al_max = al;
    verbose = true;
    
}

int DDAModel::bicgstab(int MAX_ITERATION,double MAX_ERROR){
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
    P = VectorXcd::Zero(N * 3);



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
        if (verbose && (it+1)%100==0) {
            cout << "                Iter " << it+1 << ", error=" << r.norm()/E.norm() << " MAX_ERROR="<<MAX_ERROR<<endl;
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
            return ITERATION;
        }
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    time = duration_cast<milliseconds>(t_end-t_start).count();
    cout<<"                ERROR:does not converge in "<<MAX_ITERATION<<" iterations"<<endl;
    return -1;;
}

void DDAModel::save_P(){
    P_store = P;
    int N = (*Core).get_N();
    P = VectorXcd::Zero(N * 3);
}

void DDAModel::set_P(){
    P = P_store;
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
    double lam = (*Core).get_lam();
    double K = 2 * M_PI / lam;
    double d = (*Core).get_d();
    VectorXcd* diel = (*Core).get_diel();
    for (int i = 0; i <= (*diel).size() - 1; i++) {
        al(i) = 1.0 / Get_Alpha(lam, K, d, (*diel)(i), n_E0, n_K);
    }
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

        Einternal(i)=al(i)*P(i);
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
    
    (*Core).output_to_file();
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

    name = save_position + "Model_results" + to_string(ModelLabel) + "it" + to_string(iteration) + ".txt";
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