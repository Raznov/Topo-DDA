#include "definition.h"

pair<VectorXi, VectorXd> InputInitial(string open_position, string model_label) {
    string name1 = open_position + "Commondata.txt";
    string name2 = open_position + "CoreStructure\\CoreStructure" + model_label + ".txt";
    ofstream TotalTime;
    TotalTime.open(open_position + "TotalTime.txt");

    ifstream fin1(name1), fin2(name2);
    int Nxtmp, Nytmp, Nztmp;
    int Ntmp;
    fin1 >> Nxtmp;
    fin1 >> Nytmp;
    fin1 >> Nztmp;
    fin1 >> Ntmp;
    cout << "Input geometry size: " << Ntmp << endl;
    VectorXi geometrytmp = VectorXi::Zero(3 * Ntmp);
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin1 >> geometrytmp(3 * i);
        fin1 >> geometrytmp(3 * i + 1);
        fin1 >> geometrytmp(3 * i + 2);
    }
    VectorXd dielinput = VectorXd::Zero(3 * Ntmp);
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin2 >> dielinput(3 * i);
        fin2 >> dielinput(3 * i + 1);
        fin2 >> dielinput(3 * i + 2);
    }
    fin1.close();
    fin2.close();
    return pair<VectorXi, VectorXd>(geometrytmp, dielinput);
}

//One Evo data generation
list<double> generatefocus(int min_num, int max_num, Vector3d lower_bound, Vector3d upper_bound, bool sym, double d) {
    //Center
    Vector3d center = (lower_bound + upper_bound) / 2;
    //Number of focus
    int num = round(((double)rand() / RAND_MAX) * (max_num - min_num)) + min_num;
    list<double> result;
    if (sym) {
        //For sym=ture, 2*num is generated. Symmetric distribution.
        for (int i = 0; i <= num - 1; i++) {
            double point1[3];
            double point2[3];
            for (int j = 0; j <= 2; j++) {
                point1[j] = ((double)rand() / RAND_MAX) * (upper_bound(j) - lower_bound(j)) + lower_bound(j);
                point2[j] = 2 * center(j) - point1[j];
                //Point2 wont get out of boundary.
                if (point2[j] > upper_bound[j]) {
                    point2[j] = upper_bound[j];
                }
                if (point2[j] < lower_bound[j]) {
                    point2[j] = lower_bound[j];
                }
            }
            //Push the points into the list.
            for (int j = 0; j <= 2; j++) {
                result.push_back(point1[j] * d);
            }

            for (int j = 0; j <= 2; j++) {
                result.push_back(point2[j] * d);
            }
        }

    }
    else {
        //When sym=false, 1*num focus is generated.
        for (int i = 0; i <= num - 1; i++) {
            for (int j = 0; j <= 2; j++) {
                double tmp = ((double)rand() / RAND_MAX) * (upper_bound(j) - lower_bound(j)) + lower_bound(j);
                result.push_back(tmp * d);
            }
        }
    }

    /*
    list<double>::iterator it = result.begin();
    cout << endl;
    for (int i = 0; i <= result.size() - 1; i++) {
        cout << (*it) << ",";
        it++;
    }
    cout << endl;
    */
    return result;
}

void Evo_Focus(SpacePara* spacepara_tmp, CoreStructure* CStr, DDAModel* TestModel, string save_position, int start_num, int max_evo,
    int min_num, int max_num, Vector3d lower_bound, Vector3d upper_bound, bool sym  //Parameters for focus generation
) {

    (*CStr).UpdateStr(spacepara_tmp);
    //(*CStr).output_to_file(save_position, start_num, "Simple");
    (*TestModel).UpdateAlpha();
    double d = (*CStr).get_d();

    double epsilon = 100;
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<DDAModel*> ModelpointerList;
    ModelpointerList.push_back(TestModel);

    list<string> ObjectFunctionNames{ "PointEList" };
    list<double> ObjectParameters1 = generatefocus(min_num, max_num, lower_bound, upper_bound, sym, d);
    list<list<double>*> ObjectParameters{ &ObjectParameters1 };

    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, CStr, ModelpointerList);
    
    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = max_evo;

    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam", start_num);

    return;
    
}

/*
void Evo_single(string save_position, Vector3i bind, Vector3d l, int MAX_ITERATION_EVO, Vector3d move_focus) {
    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    Vector3d center;
    center << l(0) / 2, l(1) / 2, l(2) / 2;
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);
    cout << center << endl;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Structure s1(S.get_total_space(), l, center);
    S = S + s1;
    SpacePara spacepara(&S, bind, "ONES");
    double d = 25;
    double E0 = 1.0;
    double epsilon = 10;
    double focus = (l(2) + 2) * d;   //nm       
    //double focus = (l(2) - 6) * d;
    cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    
    list<string> ObjectFunctionNames{ "PointE" };
    double exponent = 2;
    double ratio = 4;
    //list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };
    list<double> ObjectParameter{ center(0) * d + move_focus(0) * d,center(1) * d + move_focus(1) * d,focus + move_focus(2) * d };
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    Vector3d n_K;
    Vector3d n_E0;
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "2.5", lam(0), "nm");
    double nback = 1.0;
    AProductCore Core1(&CStr, lam(0), material, nback, "LDR");
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;
    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {
            for (int j = 0; j <= phi_num - 1; j++) {
                if (theta(i) != 0) {
                    double theta_tmp = theta(i) * M_PI / 180;
                    double phi_tmp = phi(j) * M_PI / 180;
                    n_K << sin(theta_tmp) * cos(phi_tmp), sin(theta_tmp)* sin(phi_tmp), cos(theta_tmp);
                    n_E0 = nEPerpinXZ(theta_tmp, phi_tmp);
                    if (CheckPerp(n_E0, n_K) == false) {
                        cout << "----------------------------------------theta" << theta[i] << "phi" << phi[j] << "Not perpendicular---------------------------------------" << endl;
                    }
                    if (k == 0) {
                        AngleInfo << theta[i] << endl;
                        AngleInfo << phi[j] << endl;
                        nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
                    }
                    DDAModel Model(Core, n_K, E0, n_E0);
                    ModelList.push_back(Model);
                }
            }
        }
        double theta_tmp = 0 * M_PI / 180;
        double phi_tmp = 0 * M_PI / 180;
        n_K << sin(theta_tmp) * cos(phi_tmp), sin(theta_tmp)* sin(phi_tmp), cos(theta_tmp);
        n_E0 = nEPerpinXZ(theta_tmp, phi_tmp);
        if (CheckPerp(n_E0, n_K) == false) {
            cout << "----------------------------------------theta" << 0 << "phi" << 0 << "Not perpendicular---------------------------------------" << endl;
        }
        if (k == 0) {
            AngleInfo << 0.0 << endl;
            AngleInfo << 0.0 << endl;
            nEInfo << n_E0(0) << " " << n_E0(1) << " " << n_E0(2) << endl;
        }
        DDAModel Model(Core, n_K, E0, n_E0);
        ModelList.push_back(Model);

        it++;
    }
    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;
    list<DDAModel>::iterator it1 = ModelList.begin();
    for (int i = 0; i <= ModelList.size() - 1; i++) {
        ModelpointerList.push_back(&(*it1));
        it1++;
    }
    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, HaveOriginHeritage, HaveAdjointHeritage, PenaltyFactor, save_position, &CStr, ModelpointerList);
    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return;
}
*/

void eval_FOM(string name, DDAModel* TestModel, list<double> theta, list<double> phi) {
    ofstream fout(name);

    list<double> ObjectParameter;
    list<double>::iterator itheta = theta.begin();
    
    for (int i = 0; i <= theta.size() - 1; i++) {
        double ttheta = (*itheta);
        itheta++;
        list<double>::iterator iphi = phi.begin();
        for (int j = 0; j <= phi.size() - 1; j++) {
            double tphi = (*iphi);
            iphi++;
            ObjectParameter.push_back(sin(ttheta) * cos(tphi));
            ObjectParameter.push_back(sin(ttheta) * sin(tphi));
            ObjectParameter.push_back(cos(ttheta));
            //cout << "(" << sin(ttheta) * cos(tphi) << " " << sin(ttheta) * sin(tphi) << " " << cos(ttheta) << ")" << endl;
        }
    }

    FOMscattering0D FOMcal(ObjectParameter, TestModel);
    list<double> FOMresults = FOMcal.GetVal();

    itheta = theta.begin();
    list<double>::iterator iresult = FOMresults.begin();
    for (int i = 0; i <= theta.size() - 1; i++) {
        double ttheta = (*itheta);
        itheta++;
        list<double>::iterator iphi = phi.begin();
        for (int j = 0; j <= phi.size() - 1; j++) {
            double tphi = (*iphi);
            double tresult = (*iresult);
            iphi++;
            iresult++;
            fout << (ttheta * 360) / (2 * M_PI) << " " << (tphi * 360) / (2 * M_PI) << " " << tresult << endl;
        }
    }

    fout.close();

}

void eval_FOM_2Dperiod(string name, DDAModel* TestModel, list<double> theta, list<double> phi) {
    ofstream fout(name);

    list<double> ObjectParameter;
    list<double>::iterator itheta = theta.begin();

    for (int i = 0; i <= theta.size() - 1; i++) {
        double ttheta = (*itheta);
        itheta++;
        list<double>::iterator iphi = phi.begin();
        for (int j = 0; j <= phi.size() - 1; j++) {
            double tphi = (*iphi);
            iphi++;
            ObjectParameter.push_back(sin(ttheta) * cos(tphi));
            ObjectParameter.push_back(sin(ttheta) * sin(tphi));
            ObjectParameter.push_back(cos(ttheta));
            //cout << "(" << sin(ttheta) * cos(tphi) << " " << sin(ttheta) * sin(tphi) << " " << cos(ttheta) << ")" << endl;
        }
    }

    FOMscattering2D FOMcal(ObjectParameter, TestModel);
    list<double> FOMresults = FOMcal.GetVal();

    itheta = theta.begin();
    list<double>::iterator iresult = FOMresults.begin();
    for (int i = 0; i <= theta.size() - 1; i++) {
        double ttheta = (*itheta);
        itheta++;
        list<double>::iterator iphi = phi.begin();
        for (int j = 0; j <= phi.size() - 1; j++) {
            double tphi = (*iphi);
            double tresult = (*iresult);
            iphi++;
            iresult++;
            fout << (ttheta * 360) / (2 * M_PI) << " " << (tphi * 360) / (2 * M_PI) << " " << tresult << endl;
        }
    }

    fout.close();

}
//Find the Max and Min of the input geometry in each direction
MatrixXi find_scope_3_dim(VectorXi* x) {
    int N = round(int((*x).size()) / 3);
    MatrixXi result(3, 2);
    result(0, 0) = (*x)(0);
    result(0, 1) = (*x)(0);
    result(1, 0) = (*x)(1);
    result(1, 1) = (*x)(1);
    result(2, 0) = (*x)(2);
    result(2, 1) = (*x)(2);

    for (int i = 0; i <= N - 1; i++) {
        if (result(0, 0) >= (*x)(3 * i)) {
            result(0, 0) = (*x)(3 * i);
        }
        if (result(0, 1) <= (*x)(3 * i)) {
            result(0, 1) = (*x)(3 * i);
        }
        if (result(1, 0) >= (*x)(3 * i + 1)) {
            result(1, 0) = (*x)(3 * i + 1);
        }
        if (result(1, 1) <= (*x)(3 * i + 1)) {
            result(1, 1) = (*x)(3 * i + 1);
        }
        if (result(2, 0) >= (*x)(3 * i + 2)) {
            result(2, 0) = (*x)(3 * i + 2);
        }
        if (result(2, 1) <= (*x)(3 * i + 2)) {
            result(2, 1) = (*x)(3 * i + 2);
        }
    }
    return result;
}

VectorXd initial_diel_func(string initial_diel, int N) {
    VectorXd diel;
    if (initial_diel.compare("ZEROS") == 0) {
        diel = VectorXd::Zero(N);
    }
    else if (initial_diel.compare("0.5") == 0) {
        diel = VectorXd::Ones(N);
        diel = 0.5 * diel;
    }
    else if (initial_diel.compare("ONES") == 0) {
        diel = VectorXd::Ones(N);
    }
    else if (initial_diel.compare("RANDOM") == 0) {
        diel = VectorXd::Zero(N);
        srand(time(0));
        for (int i = 0; i <= N - 1; i++) {
            double r = ((double)rand() / (RAND_MAX));
            diel(i) = r;
        }

    }
    else {
        diel = VectorXd::Zero(N);
        cout << "The initial type given does not match any of the built in method" << endl;
    }
    return diel;
}

double initial_diel_func(string initial_diel) {
    VectorXd diel;
    if (initial_diel.compare("ZEROS") == 0) {
        return 0.0;
    }
    else if (initial_diel.compare("0.5") == 0) {
        return 0.5;
    }
    else if (initial_diel.compare("ONES") == 0) {
        return 1.0;
    }
    else if (initial_diel.compare("RANDOM") == 0) {
        return ((double)rand() / (RAND_MAX));

    }
    else {
        cout << "The initial type given does not match any of the built in method" << endl;
        return 0.0;
    }
}

VectorXi build_a_bulk(int Nx, int Ny, int Nz){
    VectorXi result=VectorXi::Zero(3*Nx*Ny*Nz);
    for(int x=0;x<=Nx-1;x++){
        for(int y=0;y<=Ny-1;y++){
            for(int z=0;z<=Nz-1;z++){
                int position=z+Nz*(y+Ny*x);
                result(3*position)=x;
                result(3*position+1)=y;
                result(3*position+2)=z;
            }
        }
    }   
    return result; 
}


complex<double> Get_material(string mat, double wl, string unit){
    map<string,double> unit_dic;
    unit_dic.insert(pair<string,double>("nm",1.0e9));
    unit_dic.insert(pair<string,double>("um",1.0e6));
    unit_dic.insert(pair<string,double>("m",1.0));
    map<string,string> diel_dic;
    diel_dic.insert(pair<string,string>("Ag","diel/Ag (Silver) - CRC (raw)"));
    diel_dic.insert(pair<string,string>("Al","diel/Al (Aluminium) - Palik (raw)"));
    diel_dic.insert(pair<string,string>("Au","diel/Au (Gold) - Johnson and Christy (raw)"));
    diel_dic.insert(pair<string,string>("Si","diel/Si (Silicon) - Palik (raw)"));
    diel_dic.insert(pair<string,string>("SiO2","diel/SiO2 (Glass) - Palik (raw)"));
    diel_dic.insert(pair<string,string>("TiO2", "diel/TiO2_ALD (raw)"));
    diel_dic.insert(pair<string,string>("Air","diel/Air"));
    diel_dic.insert(pair<string,string>("1.5","diel/Diel1.5"));
    diel_dic.insert(pair<string,string>("2.0","diel/Diel2.0"));
    diel_dic.insert(pair<string,string>("2.5","diel/Diel2.5"));
    diel_dic.insert(pair<string,string>("3.0","diel/Diel3.0"));
    diel_dic.insert(pair<string,string>("3.5","diel/Diel3.5"));
    diel_dic.insert(pair<string,string>("4.0","diel/Diel4.0"));
    diel_dic.insert(pair<string,string>("4.5","diel/Diel4.5"));
    diel_dic.insert(pair<string,string>("5.0","diel/Diel5.0"));
    diel_dic.insert(pair<string, string>("H2O", "diel/H2O"));
    diel_dic.insert(pair<string, string>("TiN", "diel/TiN-Pfluger"));
    diel_dic.insert(pair<string, string>("Ti", "diel/Ti-JC"));
    wl=wl/unit_dic[unit];
    string real="Re_eps.txt";
    string imag="Im_eps.txt";
    string mat_real_name=diel_dic[mat]+real;
    string mat_imag_name=diel_dic[mat]+imag; 
    double mat_real,mat_imag,a,b,up,down,up_value,down_value;
    
    up=1.0;
    down=0.0;

    ifstream mat_real_file;
    mat_real_file.open(mat_real_name);
    while(mat_real_file>>a>>b){
        if(a<=wl&&a>down){
            down=a;
            down_value=b;
        }
        else if(a>wl&&a<up){
            up=a;
            up_value=b;
        }

    }
    mat_real_file.close();
    if(up==1.0||down==0.0){
        cout<<"ERROR, wavelength out of range of the diel file provided."<<endl;
    }
    else{
        mat_real=down_value+(wl-down)*(up_value-down_value)/(up-down); 
        cout << "real:" << mat_real << endl;
    }
    up=1.0;
    down=0.0;

    ifstream mat_imag_file;
    mat_imag_file.open(mat_imag_name);
    while(mat_imag_file>>a>>b){
        if(a<=wl&&a>down){
            down=a;
            down_value=b;
        }
        else if(a>wl&&a<up){
            up=a;
            up_value=b;
        }

    }
    mat_imag_file.close();
    if(up==1.0||down==0.0){
        cout<<"ERROR, wavelength out of range of the diel file provided."<<endl;
    }
    else{
        mat_imag=down_value+(wl-down)*(up_value-down_value)/(up-down); 
        cout << "imag:" << mat_imag << endl;
    }
    
    complex<double> result=mat_real+1.0i*mat_imag;
    return result;
}

Vector2cd Get_2_material(string sub, string mat, double wl, string unit){
    Vector2cd result;
    result(0)=Get_material(sub,wl,unit);
    result(1)=Get_material(mat,wl,unit);
    return result;
}

VectorXcd Get_X_material(list<string> mat_l, double wl, string unit) {
    list<string>::iterator it = mat_l.begin();
    VectorXcd result = VectorXcd::Zero(mat_l.size());
    int i = 0;
    while (it != mat_l.end()) {
        result(i) = Get_material(*it, wl, unit);
        i++;
        it++;
    }
    return result;
}

complex<double> Get_Alpha(double lam, double K, double d, complex<double> diel, Vector3d n_E0, Vector3d n_K){
    double b1 = -1.891531;
    double b2 = 0.1648469;
    double b3 = -1.7700004;
    //cout<<"lam"<<lam<<"K"<<K<<"d"<<d<<endl;
    std::complex<double> a1=(3*pow(d,3)/(4*M_PI))*(diel-1.0)/(diel+2.0);
    //cout<<a1<<endl;
    std::complex<double> result=0.0+(2.0/3.0)*pow(K*d,3)*1.0i;
    double S = pow(n_E0(0) * n_K(0), 2) + pow(n_E0(1) * n_K(1), 2) + pow(n_E0(2) * n_K(2), 2);
    result = 1.0 + (a1 / pow(d, 3)) * ((b1 + diel * b2 + diel * b3 * S) * pow(K * d, 2) - result);
    //cout<<result<<endl;
    result=a1/result;
    return result;
}

complex<double> Get_Alpha_FCD(double lam, double K, double d, complex<double> diel) {
    complex<double> M = (4.0 / 3.0) * pow((K * d), 2) + (2.0 / 3.0) * (1.0i + (1 / M_PI) * log((M_PI - K * d) / (M_PI + K * d))) * pow(K * d, 3);
    // complex<double> kappa = (diel - 1.0) / (4 * M_PI);
    // complex<double> result = pow(d, 3) * kappa / (1.0 + (4.0 * M_PI / 3.0 - M) * kappa);
    complex<double> a1 = (3 * pow(d, 3) / (4 * M_PI)) * (diel - 1.0) / (diel + 2.0);
    complex<double> result = 1.0 - M * a1 / pow(d, 3);
    result = a1 / result;
    return result;
}

/*
ArrayXcd FFT(int nx,int ny,int nz,ArrayXcd in,int _direction){
    //0 for forward, 1 for backward
    //FFTW_FORWARD and FFTW_BACKWARD are FFTW defined integer
    //normalize is needed because backward needs to be normalized to be the same as original

    int direction[2]={FFTW_FORWARD,FFTW_BACKWARD};
    int N=nx*ny*nz;
    int normalize[2]={1,N};
    fftw_complex* in1;
    fftw_complex* out;
    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for(int i=0;i<=N-1;i++){
        in1[i][0]=in(i).real();
        in1[i][1]=in(i).imag();
    }
    fftw_plan p1;
    p1=fftw_plan_dft_3d(nx,ny,nz,in1,out,direction[_direction],FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);
    ArrayXcd result(N);
    for(int i=0;i<=N-1;i++){
        result(i)=out[i][0]+out[i][1]*1.0i;
    }
    result=result/normalize[_direction];
    fftw_free(in1);
    fftw_free(out);

    return result;

}
*/

bool CheckPerp(Vector3d v1, Vector3d v2) {
    if (abs(v1.dot(v2)) <= 0.0001) {
        //cout << v1.dot(v2) << endl;
        return true;
    }
    else {
        return false;
    }
}

Vector3d nEPerpinXZ(double theta, double phi) {
    double tmp = cos(theta) / sqrt(pow(sin(theta), 2) * pow(cos(phi), 2) + pow(cos(theta), 2));
    double x[2] = { -tmp,tmp };
    double z[2] = { -sqrt(1 - tmp * tmp),sqrt(1 - tmp * tmp) };
    Vector3d k, nE;
    k << sin(theta) * cos(phi), sin(theta)* sin(phi), cos(theta);
    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= 1; j++) {
            nE << x[i], 0.0, z[j];
            if (CheckPerp(k, nE) == true && x[i] >= 0.0) {
                return nE;
            }
        }
    }

    cout << "ERROR : perp nE not found in Vector3d nEPerpinXZ(double theta, double phi)" << endl;

    return nE;
}

int makedirect(string name) {
    const char* tmp = name.c_str();
    return _mkdir(tmp);
}

list<double> makelist(double start, double end, double interval) {
    list<double> result;
    int number = floor((end - start) / interval + 1);
    for (int i = 0; i <= number - 1; i++) {
        result.push_back(start + i * interval);
    }
    return result;
}

list<double> makelist(double start, double end, int number) {
    list<double> result;
    double interval = (end - start) / (double(number) - 1.0);
    for (int i = 0; i <= number - 1; i++) {
        //cout << start + i * interval << endl;
        result.push_back(start + i * interval);
    }
    return result;
}

double exp_update(const double x, const double x_max, const double y_min, const double y_max) {
    return y_min + (y_max - y_min) * exp((x / x_max - 1.0));
}

double piecewise_update(const double x, const double x_max, const double y_min, const double y_max) {
    if (x <= 0.5 * x_max) {
        return y_min;
    }
    else if(0.5 * x_max < x&& x <= 0.7 * x_max) {
        return y_min+(y_max-y_min)/20;
    }
    else if (0.7 * x_max < x && x <= 0.8 * x_max) {
        return y_min + (y_max - y_min) / 10;
    }
    else {
        return y_max;
    }

}

double linear_update(const double x, const double x_max, const double y_min, const double y_max) {
    return y_min + (y_max - y_min) * x / x_max;
}