
//##########################################################DDA CALCULATION TEMPLATES###################################################
//DDAModel verification; updated after adding FCD and binding; updated on 2021-6-20
int main() {

    int Nx, Ny, Nz;
    Nx = 56; Ny = 56; Nz = 56;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r;
    Vector3d center;

    d = 2;
    r = 55 / d;
    center << Nx / 2, Ny / 2, Nz / 2;

    Structure s1(S.get_total_space(), r, center);



    S = S + s1;

    double ratioESItoG = 1 / (2.998 * pow(10, 4));
    double rationmtocm = 1 / (pow(10, 7));
    double ratioESItoGnm = (1 / (2.998 * pow(10, 4))) / (pow(10, 7));



    double lam = 525;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");

    double epsilon = 100;

    double focus = 350;   //nm       


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 100;

    //d = rationmtocm * d;
    //lam = rationmtocm * lam;
    //E0 = ratioESItoGnm * E0;
    //focus = rationmtocm * focus;


    list<string> ObjectFunctionNames{ "PointEratio" };

    double exponent = 2;
    double ratio = 4;

    //list<double> ObjectParameter1{ focus, exponent, ratio };
    list<double> ObjectParameter2{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter2 };
    string save_position = "";

    Vector3i bind(1, 1, 1);
    //SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);

    SpacePara spacepara(&S, bind, "ONES");

    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();
    TestModel.output_to_file();

    double x, y, z;
    y = center(1) * d;
    z = center(2) * d;


    N = Core.get_N();
    VectorXcd* P = TestModel.get_P();
    VectorXi* R = Core.get_R();
    double K = 2 * M_PI / lam;
    Vector3cd E_sum = Vector3cd::Zero();
    Vector3cd E_ext = Vector3cd::Zero();

    string name = "DDAModelEfieldscand=" + to_string(d) + ".txt";
    ofstream fout(name);
    for (x = center(0) * d; x <= center(0) * d + 350; x += 2) {
        cout << "x" << x - center(0) * d << endl;
        E_ext(0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_ext(1) = E0 * n_E0(1) * (cos(K * (n_K(0) * y + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_ext(2) = E0 * n_E0(2) * (cos(K * (n_K(0) * z + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_sum(0) = E_ext(0);
        E_sum(1) = E_ext(1);
        E_sum(2) = E_ext(2);
        for (int i = 0; i <= N - 1; i++) {
            double rx = x - d * (*R)(3 * i);                  //R has no d in it, so needs to time d
            double ry = y - d * (*R)(3 * i + 1);
            double rz = z - d * (*R)(3 * i + 2);
            Matrix3cd A = Core.A_dic_generator(rx, ry, rz);
            E_sum(0) -= (A(0, 0) * (*P)(3 * i) + A(0, 1) * (*P)(3 * i + 1) + A(0, 2) * (*P)(3 * i + 2));
            E_sum(1) -= (A(1, 0) * (*P)(3 * i) + A(1, 1) * (*P)(3 * i + 1) + A(1, 2) * (*P)(3 * i + 2));
            E_sum(2) -= (A(2, 0) * (*P)(3 * i) + A(2, 1) * (*P)(3 * i + 1) + A(2, 2) * (*P)(3 * i + 2));
        }
        double normE = E_sum.norm();
        cout << "E" << normE << endl;
        fout << x - center(0) * d << " " << normE << endl;
    }
    fout.close();

    //EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material, AMatrixMethod);

    //TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}
//2021-6-20 feeding heterodimer back d=10
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();


    double ld = 375;
    double tTiO2 = 180;
    double d = 10;
    double disp = 155;

    Vector3d l;
    Vector3d center;
    l << ld / d, ld / d, tTiO2 / d;
    //l << 40.0, 40.0, 8.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 1); Ny = round(l(1) + 1); Nz = round(l(2) + 1);
    cout << center << endl;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    Structure s1(S.get_total_space(), l, center);



    S = S + s1;

    Vector2d center1, center2;
    center1 << round(115 / d), round(115 / d);
    center2 << round(center1(0) + disp * cos(70 * PI / 180) / d), round(center1(1) + disp * sin(70 * PI / 180) / d);
    list<double> r_list{ 30 / d, 40 / d };
    list<Vector2d> center_list{ center1, center2 };
    list<string> initials_list{ "ZEROS", "ZEROS" };
    double r = 150 / d;

    Vector3i bind(1, 1, round(l(2)));
    //SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);


    SpacePara spacepara(&S, bind, "ONES", &initials_list, &r_list, &center_list);


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    //double focus = (l(2) + 2) * d;
    //cout << focus << endl;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 500;

    list<string> ObjectFunctionNames{ "IntegratedE" };


    Vector3d n_K;
    Vector3d n_E0;

    n_K << 0.0, 0.0, -1.0;
    n_E0 << 1.0, 0.0, 0.0;


    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    double lam_min = 630;
    double lam_max = 660;
    double lam_interval = 1;
    int lam_num = round((lam_max - lam_min) / lam_interval) + 1;

    CoreStructure CStr(&spacepara, d);
    string save_position = ".\\dimer-wavelengthscan-d=10\\";
    string name = save_position + "IntegratedEwavedepend=" + to_string(d) + ".txt";

    CStr.output_to_file(save_position, 0);
    //S.show_something_about_Structures();
    ofstream fout(name);

    for (int i = 0; i <= lam_num; i++) {
        double lam = lam_min + i * lam_interval;
        Vector2cd material;
        material = Get_2_material("Air", "TiO2", lam, "nm");
        int m, n;
        double Lm, Ln;
        m = 50;
        n = 50;
        Lm = (Nx + 1) * d;
        Ln = (Ny + 1) * d;
        AProductCore Core(&CStr, lam, material, m, n, Lm, Ln, "FCD");
        DDAModel TestModel(&Core, n_K, E0, n_E0);
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position, 0, lam);

        VectorXcd* E_internal = TestModel.get_Einternal();
        VectorXi* R = CStr.get_R();
        N = CStr.get_N();
        double E_int = 0;
        for (int j = 0; j < N; j++) {
            if ((*R)(3 * j + 2) >= 0) {
                double E_sum_temp = 0;
                for (int k = 0; k < 3; k++) {
                    E_sum_temp += pow(abs((*E_internal)(3 * j + k)), 2);
                }

                E_int += pow(E_sum_temp, 2) * ((*(CStr.get_diel_old()))(j * 3) + 0.0001) / 4.0;  //prevent nan result for devp calculation.
                //cout << "j: " << j << " " << E_int << endl;

            }
        }

        E_int = log(E_int);
        cout << "E_int" << E_int << endl;
        fout << lam << " " << E_int << endl;
    }
    fout.close();


    //AProductCore Core1(&CStr, lam(0), material, "LDR")






    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}
//2021-6-21 import structure from binary files converted from COMSOL and do DDA calculation
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();


    string open_position = ".\\170-air-d=10\\";
    string name1 = open_position + "CoreStructre-fromComsol-170air.txt";

    ifstream fin1(name1);
    int Nx, Ny, Nz;
    int Ntmp;
    fin1 >> Nx;
    fin1 >> Ny;
    fin1 >> Nz;
    fin1 >> Ntmp;
    VectorXi geometry = VectorXi::Zero(3 * Ntmp);
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin1 >> geometry(3 * i);
        fin1 >> geometry(3 * i + 1);
        fin1 >> geometry(3 * i + 2);
    }
    double d;
    fin1 >> d;
    fin1.close();


    string name2, name3, name4;
    name2 = open_position + "geometryPara-fromComsol-170air.txt";
    name3 = open_position + "Para-fromComsol-170air.txt";
    name4 = open_position + "bind-fromComsol-170air.txt";
    ifstream fin2(name2), fin3(name3), fin4(name4);

    VectorXi geometryPara = VectorXi::Zero(Ntmp);
    int parasize;
    fin3 >> parasize;
    VectorXd Para = VectorXd::Zero(parasize);
    Vector3i bind = Vector3i::Zero();
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin2 >> geometryPara(i);
    }
    for (int i = 0; i <= parasize - 1; i++) {
        fin3 >> Para(i);
    }
    for (int i = 0; i <= 2; i++) {
        fin4 >> bind(i);
    }
    fin2.close();
    fin3.close();
    fin4.close();


    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Structure s1(S.get_total_space(), &geometry);

    S = S + s1;


    SpacePara spacepara(&S, bind, &geometryPara, &Para);


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    //double focus = (l(2) + 2) * d;
    //cout << focus << endl;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 500;

    list<string> ObjectFunctionNames{ "IntegratedE" };


    Vector3d n_K;
    Vector3d n_E0;

    n_K << 0.0, 0.0, -1.0;
    n_E0 << 1.0, 0.0, 0.0;


    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    double lam_min = 600;
    double lam_max = 680;
    double lam_interval = 1;
    int lam_num = round((lam_max - lam_min) / lam_interval) + 1;

    CoreStructure CStr(&spacepara, d);
    string save_position = ".\\170-air-d=10\\";
    string name = save_position + "IntegratedEwavedepend=" + to_string(d) + ".txt";

    CStr.output_to_file(save_position, 0);
    //S.show_something_about_Structures();
    ofstream fout(name);

    for (int i = 0; i <= lam_num; i++) {
        double lam = lam_min + i * lam_interval;
        Vector2cd material;
        material = Get_2_material("Air", "TiO2", lam, "nm");
        int m, n;
        double Lm, Ln;
        m = 50;
        n = 50;
        Lm = (Nx + 1) * d;
        Ln = (Ny + 1) * d;
        AProductCore Core(&CStr, lam, material, m, n, Lm, Ln, "FCD");
        DDAModel TestModel(&Core, n_K, E0, n_E0);
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position, 0, lam);

        VectorXcd* E_internal = TestModel.get_Einternal();
        VectorXi* R = CStr.get_R();
        N = CStr.get_N();
        double E_int = 0.0;
        for (int j = 0; j < N; j++) {
            if ((*R)(3 * j + 2) >= 0) {
                double E_sum_temp = 0;
                for (int k = 0; k < 3; k++) {
                    E_sum_temp += pow(abs((*E_internal)(3 * j + k)), 2);
                }

                E_int += pow(E_sum_temp, 3 / 2) * ((*(CStr.get_diel_old()))(j * 3));  //prevent nan result for devp calculation.


            }
        }

        E_int = log(E_int * pow(d, 3));   //unit: nm^3
        cout << "E_int" << E_int << endl;
        fout << lam << " " << E_int << endl;
    }
    fout.close();


    //AProductCore Core1(&CStr, lam(0), material, "LDR")






    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();

    return 0;

}
//verify nn results
int main() {

    list<int> itchecklist{ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130 ,140, 150, 160, 170, 180, 190, 200 };

    for (list<int>::const_iterator itcheckit = itchecklist.begin(); itcheckit != itchecklist.end(); itcheckit++) {
        int itcheck = *itcheckit;
        //string open_position = "E:\\EM-Field-Surrogate-Solver-based-on-DDA\\32328-opt-nn-200-10-10-learn1\\";
        string open_position = "E:\\EM-Field-Surrogate-Solver-based-on-DDA\\Local files\\Optimization results\\4layer-lr0.0005\\32328-opt-nn-100-1-10-0.7-learn0.5-initialONE-weightdecay0.0-pointEmax\\";
        string name1 = open_position + "Commondata.txt";

        ifstream fin1(name1);
        int Nx, Ny, Nz;
        int Ntmp;
        fin1 >> Nx;
        fin1 >> Ny;
        fin1 >> Nz;
        fin1 >> Ntmp;
        VectorXi geometry = VectorXi::Zero(3 * Ntmp);
        for (int i = 0; i <= Ntmp - 1; i++) {
            fin1 >> geometry(3 * i);
            fin1 >> geometry(3 * i + 1);
            fin1 >> geometry(3 * i + 2);
        }
        double d;
        fin1 >> d;
        fin1.close();


        string name2, name3;
        name2 = open_position + "CoreStructure\\CoreStructure" + to_string(itcheck) + ".txt";
        name3 = open_position + "Model_output\\Model_resultsit" + to_string(itcheck) + ".txt";
        ifstream fin2(name2), fin3(name3);

        VectorXi geometryPara = VectorXi::Zero(Ntmp);
        int parasize = Ntmp;
        VectorXd Para = VectorXd::Zero(parasize);
        Vector3i bind(1, 1, 1);

        for (int i = 0; i <= Ntmp - 1; i++) {
            double x = geometry(3 * i);
            double y = geometry(3 * i + 1);
            double z = geometry(3 * i + 2);
            int pos = z + Nz * (y + Ny * x);
            geometryPara(i) = pos;
        }

        for (int i = 0; i <= parasize - 1; i++) {
            double fake;
            fin2 >> Para(i);
            fin2 >> fake;
            fin2 >> fake;
        }
        cout << parasize << endl;
        VectorXi FreeparatoPara = VectorXi::Zero(parasize);
        for (int i = 0; i <= parasize - 1; i++) {
            FreeparatoPara(i) = i;
        }
        fin2.close();
        fin3.close();


        int N = 0;
        VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
        list<Structure> ln;
        Space S(&total_space, Nx, Ny, Nz, N, &ln);

        Structure s1(S.get_total_space(), &geometry);

        S = S + s1;


        SpacePara spacepara(&S, bind, &geometryPara, &Para, &FreeparatoPara);


        double E0 = 1.0;


        //double epsilon = 1;

        //double focus = (l(2) + 2) * d;   //nm       
        //double focus = (l(2) + 2) * d;
        //cout << focus << endl;

        int MAX_ITERATION_DDA = 10000;
        double MAX_ERROR = 0.00001;


        Vector3d n_K;
        Vector3d n_E0;

        n_K << 0.0, 0.0, 1.0;
        n_E0 << 1.0, 0.0, 0.0;



        CoreStructure CStr(&spacepara, d);

        string save_position1 = open_position + "CoreStructure_verify\\";
        string save_position2 = open_position + "Model_output_verify\\";

        CStr.output_to_file(save_position1, itcheck, "simplify");
        double lam = 700;
        Vector2cd material;
        material = Get_2_material("Air", "SiO2", lam, "nm");
        AProductCore Core(&CStr, lam, material, "LDR");
        DDAModel TestModel(&Core, n_K, E0, n_E0);
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position2, itcheck);
    }

    return 0;

}
//SiO2 sphere far field scattering
int main() {

    int Nx, Ny, Nz;
    Nx = 65; Ny = 65; Nz = 65;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r;
    Vector3d center;

    d = 2;
    r = 60 / d;
    center << Nx / 2, Ny / 2, Nz / 2;

    Structure s1(S.get_total_space(), r, center);
    Vector3d l;
    l << 120 / d, 120 / d, 120 / d;
    //Structure s1(S.get_total_space(), r, center);


    S = S + s1;



    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES");
    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();



    string save_position = ".\\SiO2-sphere-scatter-r60-d2-lam500\\";
    CStr.output_to_file(save_position + "CoreStructure\\", 0, "simple");
    TestModel.output_to_file(save_position + "Model_output\\", 0);


    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    string name = save_position + "theta90phi0to360.txt";
    list<double> theta{ M_PI / 2 };
    double start = 0;
    double end = 2 * M_PI;
    int number = 200;
    list<double> phi = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta, phi);

    name = save_position + "theta-9090phi0.txt";
    list<double> phi1{ 0.0 };
    number = 100;
    start = 0.0;
    end = M_PI;
    list<double> theta1 = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta1, phi1);

    name = save_position + "theta-9090phi180.txt";
    list<double> phi2{ M_PI };
    start = M_PI;
    end = 0.0;
    list<double> theta2 = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta2, phi2);











    return 0;

}
//SiO2 block far field scattering
int main() {

    int Nx, Ny, Nz;
    Nx = 42; Ny = 42; Nz = 42;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r;
    Vector3d center;

    d = 25;
    //r = 60 / d;
    center << Nx / 2, Ny / 2, Nz / 2;


    Vector3d l;
    l << 39, 39, 39;
    //Structure s1(S.get_total_space(), r, center);
    Structure s1(S.get_total_space(), l, center);

    S = S + s1;



    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES");
    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();



    string save_position = ".\\1um1um1um-SiO2-scatter-d25-lam500\\";
    CStr.output_to_file(save_position + "CoreStructure\\", 0, "simple");
    TestModel.output_to_file(save_position + "Model_output\\", 0);


    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    string name = save_position + "theta90phi0to360.txt";
    list<double> theta{ M_PI / 2 };
    double start = 0;
    double end = 2 * M_PI;
    int number = 200;
    list<double> phi = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta, phi);

    name = save_position + "theta-9090phi0.txt";
    list<double> phi1{ 0.0 };
    number = 100;
    start = 0.0;
    end = M_PI;
    list<double> theta1 = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta1, phi1);

    name = save_position + "theta-9090phi180.txt";
    list<double> phi2{ M_PI };
    start = M_PI;
    end = 0.0;
    list<double> theta2 = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta2, phi2);











    return 0;

}
//dimer
int main() {

    int Nx, Ny, Nz;
    Nx = 125; Ny = 65; Nz = 65;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r1, r2;
    double gap;
    Vector3d center1, center2;
    d = 2;
    gap = 40 / d;
    r2 = 30 / d;
    r1 = 60 / d;
    center1 << (20 + r2 + gap + r1), Ny / 2, Nz / 2;

    Structure s1(S.get_total_space(), r1, center1);
    //Vector3d l;
    //l << 120 / d, 120 / d, 120 / d;

    center2 << 20, Ny / 2, Nz / 2;

    Structure s2(S.get_total_space(), r2, center2);


    S = S + s1;
    S = S + s2;



    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES");
    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();



    string save_position = ".\\SiO2-dimer-scatter-r3060-d2-lam500\\";
    CStr.output_to_file(save_position + "CoreStructure\\", 0, "simple");
    TestModel.output_to_file(save_position + "Model_output\\", 0);


    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    string name = save_position + "theta90phi0to360.txt";
    list<double> theta{ M_PI / 2 };
    double start = 0;
    double end = 2 * M_PI;
    int number = 200;
    list<double> phi = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta, phi);

    name = save_position + "theta-9090phi0.txt";
    list<double> phi1{ 0.0 };
    number = 100;
    start = 0.0;
    end = M_PI;
    list<double> theta1 = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta1, phi1);

    name = save_position + "theta-9090phi180.txt";
    list<double> phi2{ M_PI };
    start = M_PI;
    end = 0.0;
    list<double> theta2 = makelist(start, end, number);
    eval_FOM(name, &TestModel, theta2, phi2);











    return 0;

}
//gold sphere far field scattering spectrum
int main() {

    int Nx, Ny, Nz;
    Nx = 65; Ny = 65; Nz = 65;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r;
    Vector3d center;

    d = 2;
    r = 60 / d;
    center << Nx / 2, Ny / 2, Nz / 2;

    Structure s1(S.get_total_space(), r, center);
    Vector3d l;
    l << 120 / d, 120 / d, 120 / d;
    //Structure s1(S.get_total_space(), r, center);


    S = S + s1;


    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES");
    CoreStructure CStr(&spacepara, d);

    string save_position = ".\\Au-sphere-scatter-r60-d2-lamscan\\";
    list<double> wave_result;
    ofstream fout_wave;
    list<double> theta_wave{ M_PI / 2 };
    list<double> phi_wave{ M_PI / 2 };
    for (double lam = 200; lam <= 650; lam += 50) {
        Vector2cd material = Get_2_material("Air", "Au", lam, "nm");
        int iteration = round((lam - 200) / 50);
        AProductCore Core(&CStr, lam, material, "FCD");
        DDAModel TestModel(&Core, n_K, E0, n_E0);

        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();

        CStr.output_to_file(save_position + "CoreStructure\\", iteration, "simple");
        TestModel.output_to_file(save_position + "Model_output\\", iteration);

        ofstream Common;
        Common.open(save_position + "Commondata.txt");
        Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
        Common << (spacepara.get_geometry()) << endl;
        Common << d << endl;
        Common << n_E0 << endl;
        Common << n_K << endl;

        string name = save_position + to_string(iteration) + "theta90phi.txt";
        list<double> theta{ M_PI / 2 };
        double start = 0;
        double end = 2 * M_PI;
        int number = 200;
        list<double> phi = makelist(start, end, number);
        eval_FOM(name, &TestModel, theta, phi);

        name = save_position + to_string(iteration) + "theta-9090phi0.txt";
        list<double> phi1{ 0.0 };
        number = 100;
        start = 0.0;
        end = M_PI;
        list<double> theta1 = makelist(start, end, number);
        eval_FOM(name, &TestModel, theta1, phi1);

        name = save_position + to_string(iteration) + "theta-9090phi180.txt";
        list<double> phi2{ M_PI };
        start = M_PI;
        end = 0.0;
        list<double> theta2 = makelist(start, end, number);
        eval_FOM(name, &TestModel, theta2, phi2);



        list<double> ObjectParameter{ -0.5,0.8660254038,0 };

        FOMscattering0D FOMcal(ObjectParameter, &TestModel);
        list<double> FOMresults = FOMcal.GetVal();
        list<double>::iterator it_wave = FOMresults.begin();
        wave_result.push_back(*it_wave);

    }
    string name = save_position + "wavescan.txt";
    fout_wave.open(name);
    list<double>::iterator it_wave = wave_result.begin();
    for (double lam = 200; lam <= 650; lam += 50) {
        fout_wave << lam << " " << (*it_wave) << endl;
        it_wave++;
    }

    return 0;

}
//2D periodic scattering
int main() {

    int Nx, Ny, Nz;
    Nx = 65; Ny = 65; Nz = 65;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r;
    Vector3d center;

    d = 2;
    r = 60 / d;
    center << Nx / 2, Ny / 2, Nz / 2;

    Structure s1(S.get_total_space(), r, center);
    Vector3d l;
    l << 120 / d, 120 / d, 120 / d;
    //Structure s1(S.get_total_space(), r, center);
    double period = 800;
    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = period;
    Ln = period;

    S = S + s1;



    double lam = 400;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES");
    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, m, n, Lm, Ln, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();



    string save_position = ".\\SiO2-sphere-scatter-r60-d2-lam400-m800n800\\";
    CStr.output_to_file(save_position + "CoreStructure\\", 0, "simple");
    TestModel.output_to_file(save_position + "Model_output\\", 0);


    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    string name = save_position + "theta90phi0to360.txt";
    list<double> theta{ M_PI / 2 };
    double start = 0;
    double end = 2 * M_PI;
    int number = 200;
    list<double> phi = makelist(start, end, number);
    eval_FOM_2Dperiod(name, &TestModel, theta, phi);

    name = save_position + "theta0to180phi0.txt";
    list<double> phi1{ 0.0 };
    number = 200;
    start = 0.0;
    end = M_PI;
    list<double> theta1 = makelist(start, end, number);
    eval_FOM_2Dperiod(name, &TestModel, theta1, phi1);

    name = save_position + "theta180to0phi180.txt";
    list<double> phi2{ M_PI };
    start = M_PI;
    end = 0.0;
    list<double> theta2 = makelist(start, end, number);
    eval_FOM_2Dperiod(name, &TestModel, theta2, phi2);











    return 0;

}
//2D periodic scattering wave scan
int main() {

    int Nx, Ny, Nz;
    Nx = 40; Ny = 40; Nz = 40;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r;
    Vector3d center;

    d = 20;
    r = 300 / d;
    center << Nx / 2, Ny / 2, Nz / 2;

    Structure s1(S.get_total_space(), r, center);
    Vector3d l;
    l << 120 / d, 120 / d, 120 / d;
    //Structure s1(S.get_total_space(), r, center);
    double period = 800;
    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = period;
    Ln = period;

    S = S + s1;



    //double lam = 400;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES");
    CoreStructure CStr(&spacepara, d);
    string save_position = ".\\SiO2-sphere-scatter-r300-d20-lamscan-m800n800\\";
    CStr.output_to_file(save_position + "CoreStructure\\", 0, "simple");
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;
    int iteration = 0;

    list<double> ObjectParameter{ 0,0,1,0,0,1,1,1 };
    list<string> name_l;
    //list<ofstream> fout_l;
    list<list<double>> FOMresults_l;
    int out_file_num = int(round(ObjectParameter.size() / 2));
    list<double>::iterator ipara = ObjectParameter.begin();
    for (int i = 0; i <= out_file_num - 1; i++) {
        string name = "T";
        name = name + to_string(int(*ipara));
        ipara++;
        name = name + to_string(int(*ipara));
        ipara++;
        name = name + ".txt";
        name = save_position + name;
        name_l.push_back(name);
        list<double> tmp;
        FOMresults_l.push_back(tmp);
    }

    list<double> lam_l;
    double lam_min = 400;
    double lam_max = 800;
    double lam_step = 25;
    for (double lam = lam_min; lam <= lam_max; lam += lam_step) {
        lam_l.push_back(lam);
    }

    for (list<double>::iterator itlam = lam_l.begin(); itlam != lam_l.end(); itlam++) {
        double lam = *itlam;
        Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
        AProductCore Core(&CStr, lam, material, m, n, Lm, Ln, "LDR");
        DDAModel TestModel(&Core, n_K, E0, n_E0);
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position + "Model_output\\", iteration);

        iteration += 1;

        FOMscattering2D FOMcal(ObjectParameter, &TestModel);
        list<double> FOMresults = FOMcal.GetVal();
        list<double>::iterator iresult = FOMresults.begin();
        list<list<double>>::iterator iresult_l = FOMresults_l.begin();
        for (int i = 0; i <= out_file_num - 1; i++) {
            (*iresult_l).push_back(*iresult);
            iresult++;
            iresult_l++;
        }


    }

    list<string>::iterator iname = name_l.begin();
    list<list<double>>::iterator iresult_l = FOMresults_l.begin();
    for (int i = 0; i <= out_file_num - 1; i++) {
        ofstream fout(*iname);
        list<double>::iterator iresult = (*iresult_l).begin();
        for (list<double>::iterator itlam = lam_l.begin(); itlam != lam_l.end(); itlam++) {
            fout << (*itlam) << " " << (*iresult) << endl;
            iresult++;
        }
        iresult_l++;
        iname++;
    }
















    return 0;

}
//periodic 2D transmission and reflectance
int main() {




    Vector3d l;
    double d;
    double r = 10;
    Vector3d center;

    d = 10;
    l << 19.0, 19.0, 19.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Structure s1(S.get_total_space(), l, center);
    //Structure s1(S.get_total_space(), r, center);

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 20 * d;
    Ln = 20 * d;

    S = S + s1;



    //double lam = 400;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;

    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    Vector3i bind(1, 1, 1);
    SpacePara spacepara(&S, bind, "ONES", "ZEROS", r, "CYLINDER");
    CoreStructure CStr(&spacepara, d);
    string save_position = ".\\Air-TiN-cylinder-scatter-r100-d10-lamscan-m200n200\\";
    CStr.output_to_file(save_position + "CoreStructure\\", 0, "simple");
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;
    int iteration = 0;

    list<double> ObjectParameter{ 0,0 };
    list<string> name_l;
    list<string> name_l_R;
    //list<ofstream> fout_l;
    list<list<double>> FOMresults_l;
    list<list<double>> FOMresults_l_R;
    int out_file_num = int(round(ObjectParameter.size() / 2));
    list<double>::iterator ipara = ObjectParameter.begin();
    for (int i = 0; i <= out_file_num - 1; i++) {
        string name = "T";
        name = name + to_string(int(*ipara));
        ipara++;
        name = name + to_string(int(*ipara));
        ipara++;
        name = name + ".txt";
        name = save_position + name;
        name_l.push_back(name);
        list<double> tmp;
        FOMresults_l.push_back(tmp);
    }
    for (int i = 0; i <= out_file_num - 1; i++) {
        string name = "R";
        name = name + to_string(int(*ipara));
        ipara++;
        name = name + to_string(int(*ipara));
        ipara++;
        name = name + ".txt";
        name = save_position + name;
        name_l_R.push_back(name);
        list<double> tmp;
        FOMresults_l_R.push_back(tmp);
    }

    list<double> lam_l;
    double lam_min = 400;
    double lam_max = 800;
    double lam_step = 25;
    for (double lam = lam_min; lam <= lam_max; lam += lam_step) {
        lam_l.push_back(lam);
    }

    for (list<double>::iterator itlam = lam_l.begin(); itlam != lam_l.end(); itlam++) {
        double lam = *itlam;
        Vector2cd material = Get_2_material("Air", "TiN", lam, "nm");
        AProductCore Core(&CStr, lam, material, m, n, Lm, Ln, "FCD");
        DDAModel TestModel(&Core, n_K, E0, n_E0);
        TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel.update_E_in_structure();
        TestModel.solve_E();
        TestModel.output_to_file(save_position + "Model_output\\", iteration);

        iteration += 1;

        FOMscattering2D FOMcal(ObjectParameter, &TestModel);
        list<double> FOMresults = FOMcal.GetVal();
        list<double>::iterator iresult = FOMresults.begin();
        list<list<double>>::iterator iresult_l = FOMresults_l.begin();
        for (int i = 0; i <= out_file_num - 1; i++) {
            (*iresult_l).push_back(*iresult);
            iresult++;
            iresult_l++;
        }

        FOMreflect2D FOMcal_R(ObjectParameter, &TestModel);
        list<double> FOMresults_R = FOMcal_R.GetVal();
        list<double>::iterator iresult_R = FOMresults_R.begin();
        list<list<double>>::iterator iresult_l_R = FOMresults_l_R.begin();
        for (int i = 0; i <= out_file_num - 1; i++) {
            (*iresult_l_R).push_back(*iresult_R);
            iresult_R++;
            iresult_l_R++;
        }


    }

    list<string>::iterator iname = name_l.begin();
    list<list<double>>::iterator iresult_l = FOMresults_l.begin();
    for (int i = 0; i <= out_file_num - 1; i++) {
        ofstream fout(*iname);
        list<double>::iterator iresult = (*iresult_l).begin();
        for (list<double>::iterator itlam = lam_l.begin(); itlam != lam_l.end(); itlam++) {
            fout << (*itlam) << " " << (*iresult) << endl;
            iresult++;
        }
        iresult_l++;
        iname++;
    }

    list<string>::iterator iname_R = name_l_R.begin();
    list<list<double>>::iterator iresult_l_R = FOMresults_l_R.begin();
    for (int i = 0; i <= out_file_num - 1; i++) {
        ofstream fout(*iname_R);
        list<double>::iterator iresult_R = (*iresult_l_R).begin();
        for (list<double>::iterator itlam_R = lam_l.begin(); itlam_R != lam_l.end(); itlam_R++) {
            fout << (*itlam_R) << " " << (*iresult_R) << endl;
            iresult_R++;
        }
        iresult_l++;
        iname++;
    }
















    return 0;

}
//##########################################################INVERSE DESIGN TEMPLATES###################################################
//Optimization of focal metasurfaces: size 2um 2um 400nm at 500nm for diel2d5, nonperiodic
int main() {

    string save_position = ".\\thick400-diel2d5-phi0theta0-lam500-size2000-focus-m150\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 16.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.
    cout << center << endl;
    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(1, 1, 1);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 25;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 10;                                                //Fixed learning rate of the optimization.
    double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 200;                                        //Number of topology optimization.
    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };  //Focal spot postition.
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 0;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    Vector2cd material;
    material = Get_2_material("Air", "2.5", lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;
    Structure s1(S.get_total_space(), l, center);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    S = S + s1;                                                         //Add the geometry into the space.
    SpacePara spacepara(&S, bind, "ONES");                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "PointE" };                       //Name of the object function.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    AProductCore Core1(&CStr, lam(0), material, "LDR");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//Optimization of focal metasurfaces: size 2um 2um 400nm at 500nm for diel2d5, nonperiodic New 
int main() {

    string save_position = ".\\thick400-diel2d5-phi0theta0-lam500-size2000-focus-m150-new\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 16.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.
    cout << center << endl;
    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(1, 1, 1);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 25;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 10;                                                //Fixed learning rate of the optimization.
    double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 20;                                        //Number of topology optimization.
    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };  //Focal spot postition.
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 0;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    VectorXcd material;
    list<string> mat_l{ "Air", "2.5" };
    material = Get_X_material(mat_l, lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;
    Structure s1(S.get_total_space(), l, center);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    S = S + s1;                                                         //Add the geometry into the space.
    VectorXi* s1geometry = s1.get_geometry();
    list<VectorXi*> FPGeometryl{ s1geometry };
    list<VectorXi*> BPGeometryl{};
    list<double> BParal{};
    SpacePara spacepara(&S, bind, "ONES", FPGeometryl, BPGeometryl, BParal);                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "PointE" };                       //Name of the object function.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    AProductCore Core1(&CStr, lam(0), material, "LDR");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//Optimization of focal metasurfaces: size 2um 2um 400nm at 500nm for diel2d5 2D, nonperiodic filtered and binarized
int main() {

    string save_position = ".\\thick400-diel2d5-phi0theta0-lam500-size2000-focus450-2D\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 16.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.
    cout << center << endl;
    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(1, 1, 17);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 25;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 10;                                                //Fixed learning rate of the optimization.
    double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 200;                                        //Number of topology optimization.
    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };  //Focal spot postition.
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 0;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    VectorXcd material;
    list<string> mat_l{ "Air", "2.5" };
    material = Get_X_material(mat_l, lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;
    Structure s1(S.get_total_space(), l, center);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    S = S + s1;                                                         //Add the geometry into the space.
    VectorXi* s1geometry = s1.get_geometry();
    list<VectorXi*> FPGeometryl{ s1geometry };
    list<VectorXi*> BPGeometryl{};
    list<double> BParal{};
    bool Filter = true;
    double r_f = 2.0;
    FilterOption filteropt(0.0, 80.0, 0.5, "piecewise", r_f);
    SpacePara spacepara(&S, bind, "ONES", FPGeometryl, BPGeometryl, BParal, Filter, &filteropt);                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "PointE" };                       //Name of the object function.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    AProductCore Core1(&CStr, lam(0), material, "LDR");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//Optimization of periodic metasurfaces: 2021-4-8, USE this to regenerate the double arrow structure. Set integration object function for z>=8 because it is 2layer.
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 21.0, 21.0, 9.0;
    //l << 40.0, 40.0, 8.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 1); Ny = round(l(1) + 1); Nz = round(l(2) + 1);
    cout << center << endl;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    Structure s1(S.get_total_space(), l, center);



    S = S + s1;
    double d = 15;
    double r = 150 / d;

    Vector3i bind(1, 1, 10);
    SpacePara spacepara(&S, bind, "ONES", "ZEROS", r);

    //SpacePara spacepara(&S, bind, "RANDOM");


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 88;

    list<string> ObjectFunctionNames{ "IntegratedE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 70.0 };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = ".\\p330-lam542-beta7-TiO2-InE-2layer-circle\\";

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
    lam << 542;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "TiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 22 * d;
    Ln = 22 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");
    //AProductCore Core1(&CStr, lam(0), material, "FCD");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

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


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

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

    return 0;

}
//Optimization of periodic metasurfaces: with filter and binarization
int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 21.0, 21.0, 9.0;
    //l << 40.0, 40.0, 8.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 1); Ny = round(l(1) + 1); Nz = round(l(2) + 1);
    cout << center << endl;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    Structure s1(S.get_total_space(), l, center);



    S = S + s1;
    double d = 15;
    double r = 150 / d;




    Vector3i bind(1, 1, 10);
    VectorXi* s1geometry = s1.get_geometry();
    list<VectorXi*> FPGeometryl{ s1geometry };
    list<VectorXi*> BPGeometryl{};
    list<double> BParal{};
    bool Filter = true;
    double r_f = 2.0;
    FilterOption filteropt(0.0, 80.0, 0.5, "piecewise", r_f);
    SpacePara spacepara(&S, bind, "RANDOM", FPGeometryl, BPGeometryl, BParal, Filter, &filteropt);

    //SpacePara spacepara(&S, bind, "RANDOM");


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 100;

    list<string> ObjectFunctionNames{ "IntegratedE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 70.0 };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = ".\\p330-lam542-TiO2-InE-2layer-rand-filter-beta80\\";

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
    lam << 542;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "TiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 22 * d;
    Ln = 22 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");

    CorePointList.push_back(&Core1);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

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


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;


    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;


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

    return 0;

}
//This code for level up optimizing
int main() {



    string open_position = ".\\thick100-phi0theta0-lam500-size4000\\";
    string name1 = open_position + "CoreStructure\\CoreStructure99.txt";
    ofstream TotalTime;
    TotalTime.open(open_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    ifstream fin1(name1);
    int Nxtmp, Nytmp, Nztmp;
    int Ntmp;
    fin1 >> Nxtmp;
    fin1 >> Nytmp;
    fin1 >> Nztmp;
    fin1 >> Ntmp;
    cout << Ntmp << endl;
    VectorXi geometrytmp = VectorXi::Zero(3 * Ntmp);
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin1 >> geometrytmp(3 * i);
        fin1 >> geometrytmp(3 * i + 1);
        fin1 >> geometrytmp(3 * i + 2);
    }
    VectorXd dielinput = VectorXd::Zero(3 * Ntmp);
    for (int i = 0; i <= Ntmp - 1; i++) {
        fin1 >> dielinput(3 * i);
        fin1 >> dielinput(3 * i + 1);
        fin1 >> dielinput(3 * i + 2);
    }
    fin1.close();

    Vector3d l;
    Vector3d center;

    l << 180.0, 180.0, 4.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    Nx = round(l(0) + 1); Ny = round(l(1) + 1); Nz = round(l(2) + 1);
    //cout << center << endl;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    Structure s1(S.get_total_space(), l, center);



    S = S + s1;
    double d = 25;
    double r = 150 / d;

    Vector3i bind(1, 1, 1);

    SpacePara spacepara(&S, "ZEROS", geometrytmp, dielinput);


    double E0 = 1.0;


    double epsilon = 100;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 100;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = ".\\thick100-phi0theta0-lam500-size4500\\";

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
    material = Get_2_material("Air", "SiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    //int m, n;
    //double Lm, Ln;
    //m = 50;
    //n = 50;
    //Lm = 22 * d;
    //Ln = 22 * d;
    AProductCore Core1(&CStr, lam(0), material, "LDR");
    //AProductCore Core1(&CStr, lam(0), material, "FCD");
    //material = Get_2_material("Air", "2.5", lam(1), "nm");
    //AProductCore Core2(&CStr, lam(1), material, "LDR");
    //material = Get_2_material("Air", "2.5", lam(2), "nm");
    //AProductCore Core3(&CStr, lam(2), material, "LDR");
    CorePointList.push_back(&Core1);
    //CorePointList.push_back(&Core2);
    //CorePointList.push_back(&Core3);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

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


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

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

    return 0;

}
//Example for multi-task iteration
int main() {

    string save_position;
    Vector3i bind;
    Vector3d l;
    int MAX_ITERATION_EVO = 200;

    //-------------------------------------------bind2------------------------------------------
    //thick100

    save_position = ".\\thick100-diel2d5-phi0theta0-lam500-size1000-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 40.0, 40.0, 4.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick100-diel2d5-phi0theta0-lam500-size1500-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 60.0, 60.0, 4.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick100-diel2d5-phi0theta0-lam500-size2000-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 80.0, 80.0, 4.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick100-diel2d5-phi0theta0-lam500-size2500-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 100.0, 100.0, 4.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);
    //thick200
    save_position = ".\\thick200-diel2d5-phi0theta0-lam500-size1000-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 40.0, 40.0, 8.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick200-diel2d5-phi0theta0-lam500-size1500-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 60.0, 60.0, 8.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick200-diel2d5-phi0theta0-lam500-size2000-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 80.0, 80.0, 8.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick200-diel2d5-phi0theta0-lam500-size2500-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 100.0, 100.0, 8.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);
    //thick300
    save_position = ".\\thick300-diel2d5-phi0theta0-lam500-size1000-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 40.0, 40.0, 12.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick300-diel2d5-phi0theta0-lam500-size1500-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 60.0, 60.0, 12.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick300-diel2d5-phi0theta0-lam500-size2000-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 80.0, 80.0, 12.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick300-diel2d5-phi0theta0-lam500-size2500-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 100.0, 100.0, 12.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);
    //thick400
    save_position = ".\\thick400-diel2d5-phi0theta0-lam500-size1000-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 40.0, 40.0, 16.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick400-diel2d5-phi0theta0-lam500-size1500-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 60.0, 60.0, 16.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = ".\\thick400-diel2d5-phi0theta0-lam500-size2500-focus50-bind2\\";
    bind << 2, 2, 2;
    l << 100.0, 100.0, 16.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);


    save_position = ".\\thick400-diel2d5-phi0theta0-lam500-size2000-focus50-bind2-797915\\";
    bind << 2, 2, 2;
    l << 79.0, 79.0, 15.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);


    return 0;

}
//0D far field design
int main() {

    string save_position = ".\\1um1um1um-SiO2-phi0theta0-lam500-theta90\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 39.0, 39.0, 39.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.
    cout << center << endl;
    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(1, 1, 1);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 25;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 10;                                                //Fixed learning rate of the optimization.
    //double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    //cout << focus << endl;
    double stheta, sphi;
    stheta = M_PI / 2;
    sphi = 0.0;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 20;                                        //Number of topology optimization.
    list<double> ObjectParameter{ sin(stheta) * cos(sphi),sin(stheta) * sin(sphi),cos(stheta) };  //Focal spot postition.
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 0;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    Vector2cd material;
    material = Get_2_material("Air", "SiO2", lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;
    Structure s1(S.get_total_space(), l, center);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    S = S + s1;                                                         //Add the geometry into the space.
    SpacePara spacepara(&S, bind, "ONES");                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "scattering0D" };                       //Name of the object function.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    AProductCore Core1(&CStr, lam(0), material, "LDR");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//0D far field design y-dir 2D extrude
int main() {

    string save_position = ".\\1um1um1um-SiO2-phi0theta0-lam500-theta90-2Dextrusion-yinput\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 39.0, 39.0, 39.0;    //Size of the initialization block. 81*81*17 pixels in total.

    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 3);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    cout << center << endl;
    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(1, 1, 40);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 25;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 10;                                                //Fixed learning rate of the optimization.
    //double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    //cout << focus << endl;
    double stheta, sphi;
    stheta = M_PI / 2;
    sphi = 0.0;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 100;                                        //Number of topology optimization.
    list<double> ObjectParameter{ sin(stheta) * cos(sphi),sin(stheta) * sin(sphi),cos(stheta) };  //Focal spot postition.
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = true;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 90;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 90;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    Vector2cd material;
    material = Get_2_material("Air", "SiO2", lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;
    Structure s1(S.get_total_space(), l, center);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    S = S + s1;                                                         //Add the geometry into the space.
    SpacePara spacepara(&S, bind, "ONES");                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "scattering0D" };                       //Name of the object function.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    AProductCore Core1(&CStr, lam(0), material, "LDR");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//Abs design
int main() {

    string save_position = ".\\200nm3-gold-water-abs-lam600\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 19.0, 19.0, 19.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.

    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(1, 1, 20);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 10;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 10;                                                //Fixed learning rate of the optimization.
    //double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    //cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 20;                                        //Number of topology optimization.
    list<double> ObjectParameter{ 0.0 };  //Focal spot postition.
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 0;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 600;
    Vector2cd material;
    material = Get_2_material("H2O", "Au", lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;
    Structure s1(S.get_total_space(), l, center);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    S = S + s1;                                                         //Add the geometry into the space.
    SpacePara spacepara(&S, bind, "RANDOM");                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "Abs" };                       //Name of the object function.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    AProductCore Core1(&CStr, lam(0), material, "FCD");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//periodic Abs design
int main() {

    string save_position = ".\\200nm3-Water-Gold-abs-lam500-period-random\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 19.0, 19.0, 19.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.

    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(4, 4, 20);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 10;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 10;                                                //Fixed learning rate of the optimization.
    //double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    //cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 20;                                        //Number of topology optimization.
    list<double> ObjectParameter{ 0.0 };  //Focal spot postition.
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 0;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    Vector2cd material;
    material = Get_2_material("H20", "Au", lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;
    Structure s1(S.get_total_space(), l, center);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    S = S + s1;                                                         //Add the geometry into the space.
    SpacePara spacepara(&S, bind, "RANDOM");                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "Abs" };                       //Name of the object function.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 20 * d;
    Ln = 20 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//periodic Abs design with substrate
int main() {

    string save_position = ".\\100nm2-200nm-Water-TiN-abs-lam500-period-random\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 9.0, 9.0, 19.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.

    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(1, 1, 12);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 10;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 10;                                                //Fixed learning rate of the optimization.
    //double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    //cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 20;                                        //Number of topology optimization.
    list<double> ObjectParameter{ 0.0 };  //Focal spot postition.
    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 0;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    Vector2cd material;
    material = Get_2_material("H2O", "TiN", lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;

    Vector3d l1, l2, center1, center2;
    l1 << 9.0, 9.0, 11.0;
    center1 << l1(0) / 2, l1(1) / 2, 13.5;
    l2 << 9.0, 9.0, 7;
    center2 << l2(0) / 2, l2(1) / 2, 3.5;

    Structure s1(S.get_total_space(), l1, center1);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    Structure s2(S.get_total_space(), l2, center2);
    S = S + s1;                                                         //Add the geometry into the space.
    S = S + s2;

    SpacePara spacepara(&S, bind, "RANDOM", "ONES");                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "Abs" };                       //Name of the object function.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 20 * d;
    Ln = 20 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//periodic PointE design with substrate
int main() {

    string save_position = ".\\200nm3-Water-TiN-PointE7-lam500-period-ones-2D\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 19.0, 19.0, 19.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.

    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.
    Vector3i bind(1, 1, 20);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.
                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 10;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 100;                                                //Fixed learning rate of the optimization.
    //double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    //cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;                                     //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 100;                                        //Number of topology optimization.

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 180;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    Vector2cd material;
    material = Get_2_material("H2O", "TiN", lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;

    Vector3d l1, l2, center1, center2;
    l1 << 19.0, 19.0, 11.0;
    center1 << l1(0) / 2, l1(1) / 2, 13.5;
    l2 << 19.0, 19.0, 7;
    center2 << l2(0) / 2, l2(1) / 2, 3.5;

    Structure s1(S.get_total_space(), l1, center1);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    Structure s2(S.get_total_space(), l2, center2);
    S = S + s1;                                                         //Add the geometry into the space.
    S = S + s2;

    SpacePara spacepara(&S, bind, "ONES", "ONES");                              //Initialize with material index of 1, which is permittivity=2.5 in this case.
                                                                        //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "PointE" };                       //Name of the object function.
    list<double> ObjectParameter{ center2(0) * d,center2(1) * d,7.0 * d };  //Focal spot postition.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 20 * d;
    Ln = 20 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
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

    return 0;

}
//periodic PointE design wo 3layer and filter and binarization
int main() {

    string save_position = ".\\200nm2-sub100nm-Air-TiN-Ti-PointE7-lam500-period-ones-3D\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 19.0, 19.0, 29.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.

    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.

                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 10;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 100;                                                //Fixed learning rate of the optimization.
    //double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    //cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;
    //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 100;                                        //Number of topology optimization.

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 180;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    VectorXcd material;
    list<string> mat_l{ "Air", "TiN", "Ti" };
    material = Get_X_material(mat_l, lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;

    Vector3d l1, l2, l3, center1, center2, center3;
    l1 << 19.0, 19.0, 11.0;
    center1 << l1(0) / 2, l1(1) / 2, 23.5;
    l2 << 19.0, 19.0, 7.0;
    center2 << l2(0) / 2, l2(1) / 2, 13.5;
    l3 << 19.0, 19.0, 9.0;
    center3 << l3(0) / 2, l3(1) / 2, 4.5;

    Structure s1(S.get_total_space(), l1, center1);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    Structure s2(S.get_total_space(), l2, center2);
    Structure s3(S.get_total_space(), l3, center3);
    S = S + s1;                                                         //Add the geometry into the space.
    S = S + s2;
    S = S + s3;

    Vector3i bind(1, 1, 1);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.


    VectorXi* s1geometry = s1.get_geometry();
    VectorXi* s2geometry = s2.get_geometry();
    VectorXi* s3geometry = s3.get_geometry();
    list<VectorXi*> FPGeometryl{ s1geometry };
    list<VectorXi*> BPGeometryl{ s2geometry, s3geometry };
    list<double> BParal{ 1.0, 2.0 };
    bool Filter = false;
    double r_f = 2.0;
    FilterOption filteropt(0.0, 80.0, 0.5, "piecewise", r_f);

    SpacePara spacepara(&S, bind, "ONES", FPGeometryl, BPGeometryl, BParal, Filter);
    //SpacePara is where the parameter<->geometry link is established.

    list<string> ObjectFunctionNames{ "PointE" };                       //Name of the object function.
    list<double> ObjectParameter{ center2(0) * d,center2(1) * d,17.0 * d };  //Focal spot postition.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 20 * d;
    Ln = 20 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
            for (int j = 0; j <= phi_num - 1; j++) {
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

    return 0;

}
//periodic PointE design with 3layer and filter and binarization
int main() {

    string save_position = ".\\200nm2-sub100nm-Air-TiN-Ti-PointE7-lam500-period-ones-r_filter2-beta80\\";       //output file
    Vector3d l;
    Vector3d center;
    l << 19.0, 19.0, 29.0;    //Size of the initialization block. 81*81*17 pixels in total.
    center << l(0) / 2, l(1) / 2, l(2) / 2;      //Center of the block.
    int Nx, Ny, Nz;
    Nx = round(l(0) + 3); Ny = round(l(1) + 3); Nz = round(l(2) + 1);   //Size of the design space. Notice that this sets the limits for coordinates, 
                                                                        //but actual pixels outside the 81*81*7 are not included in simulation. 
                                                                        //However, if geometry has pixels outside this space, they will be cut off.

    int N = 0;                                                          //N counts the number of pixels in the geometries simulated. Initail is one.

                                                                        //This can be used to control the finest feature size of the designed structure.
    double d = 10;                                                      //Size of pixel. Here 25nm.   
    double E0 = 1.0;                                                    //Input field amplitude. 1V/m.
    double epsilon = 100;                                                //Fixed learning rate of the optimization.
    //double focus = (l(2) + 2) * d;   //nm                               //Focal spot is 50nm higher than the upper boundary of the intialization block.
    //cout << focus << endl;
    int MAX_ITERATION_DDA = 100000;
    //Number of maximum DDA iterations.
    double MAX_ERROR = 0.00001;                                         //Maximum error of DDA.
    int MAX_ITERATION_EVO = 100;                                        //Number of topology optimization.

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 0.0001;
    Vector3d n_K;                                                       //Wave vector direction.
    Vector3d n_E0;                                                      //Electric field polarization direction.
    int theta_num = 1;                                                  //Number of theta
    VectorXd theta(theta_num);
    theta << 180;                                                         //Theta values.
    int phi_num = 1;                                                    //Number of phi
    VectorXd phi(phi_num);
    phi << 0;                                                           //Phi values. For details of how theta and phi are defined, read tutorial.
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 500;
    VectorXcd material;
    list<string> mat_l{ "Air", "TiN", "Ti" };
    material = Get_X_material(mat_l, lam(0), "nm");              //Air as substrate. material with permittivity of 2.5 as design material.


    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);                    //Total space is a block.
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Vector3i direction;

    Vector3d l1, l2, l3, center1, center2, center3;
    l1 << 19.0, 19.0, 11.0;
    center1 << l1(0) / 2, l1(1) / 2, 23.5;
    l2 << 19.0, 19.0, 7.0;
    center2 << l2(0) / 2, l2(1) / 2, 13.5;
    l3 << 19.0, 19.0, 9.0;
    center3 << l3(0) / 2, l3(1) / 2, 4.5;

    Structure s1(S.get_total_space(), l1, center1);                       //Initialize the block which is 80*80*16 in terms of intervals or 81*81*17 in terms of pixels.
    Structure s2(S.get_total_space(), l2, center2);
    Structure s3(S.get_total_space(), l3, center3);
    S = S + s1;                                                         //Add the geometry into the space.
    S = S + s2;
    S = S + s3;

    Vector3i bind(1, 1, 12);                                             //binding in x,y,z. 2 means every 2 pixels will have the same material index. 3 means every 3.


    VectorXi* s1geometry = s1.get_geometry();
    VectorXi* s2geometry = s2.get_geometry();
    VectorXi* s3geometry = s3.get_geometry();
    list<VectorXi*> FPGeometryl{ s1geometry };
    list<VectorXi*> BPGeometryl{ s2geometry, s3geometry };
    list<double> BParal{ 1.0, 2.0 };
    bool Filter = true;
    double r_f = 2.0;
    FilterOption filteropt(0.0, 80.0, 0.5, "piecewise", r_f);
    SpacePara spacepara(&S, bind, "ONES", FPGeometryl, BPGeometryl, BParal, Filter, &filteropt);
    //SpacePara is where the parameter<->geometry link is established.
    list<string> ObjectFunctionNames{ "PointE" };                       //Name of the object function.
    list<double> ObjectParameter{ center2(0) * d,center2(1) * d,17.0 * d };  //Focal spot postition.
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;
    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");
    CoreStructure CStr(&spacepara, d);                                 //Used as a interface to update parameters       
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 20 * d;
    Ln = 20 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");                //Matrix vector product is carried out in AProductCore class. So in this step, wavelength and actual permittivity comes in.
    CorePointList.push_back(&Core1);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {                          //Each wavelength needs one new AProductCore. But if a structure working for multiple wavelengths
                                                                      //needs to be designed, these AProductCore should share the same CoreStructure so that the parameters
                                                                      //can be updated based on a gradient avergaed for all wavlengths. The majority of memory consumption is 
                                                                      //for the A matrix which is wavlength dependent, so having many wavelengts can use huge memory. 3 wavlengths 
                                                                      //for current structure should take up to 2 GB.
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {                    //However, as you can see. For single wavelength, different angles and polarizations can share the same AProductCore.
                                                                      //And one AProductCore can be shared among several DDAModel, such that the memory consumption for multiple-angle optimization
                                                                      //should be similar to single-angle case. Of course, this is optimization for performance for a range of angles averaged. 
            for (int j = 0; j <= phi_num - 1; j++) {
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

    return 0;

}

//##########################################################NN DATA SET GENERATION###################################################
//radnom 3D SiO2 several wavelength
int main() {

    srand((unsigned)(time(0)));

    int Nx, Ny, Nz;
    Nx = 32; Ny = 32; Nz = 8;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;

    Vector3d center;
    Vector3d l;

    d = 20;

    center << double(Nx - 1) / 2, double(Ny - 1) / 2, double(Nz - 1) / 2;
    l << Nx - 1, Ny - 1, Nz - 1;

    Structure s1(S.get_total_space(), l, center);



    S = S + s1;

    double lam1 = 650;
    double lam2 = 700;
    double lam3 = 750;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material1 = Get_2_material("Air", "SiO2", lam1, "nm");
    Vector2cd material2 = Get_2_material("Air", "SiO2", lam2, "nm");
    Vector2cd material3 = Get_2_material("Air", "SiO2", lam3, "nm");




    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;


    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    Vector3i bind(1, 1, 1);
    int maxnumber = 15;
    int minnumber = 6;
    int number = 0;
    double limitx1 = 2;
    double limitx2 = 5;
    double limity1 = 2;
    double limity2 = 5;
    double limitz1 = 2;
    double limitz2 = 5;

    SpacePara spacepara(&S, bind, maxnumber, limitx1, limitx2, limity1, limity2, limitz1, limitz2);


    CoreStructure CStr(&spacepara, d);
    AProductCore Core1(&CStr, lam1, material1, "LDR");
    AProductCore Core2(&CStr, lam2, material2, "LDR");
    AProductCore Core3(&CStr, lam3, material3, "LDR");
    DDAModel TestModel1(&Core1, n_K, E0, n_E0);
    DDAModel TestModel2(&Core2, n_K, E0, n_E0);
    DDAModel TestModel3(&Core3, n_K, E0, n_E0);

    string save_position = ".\\SiO2-rects-3D-650~750\\";

    TestModel1.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel1.update_E_in_structure();
    TestModel1.solve_E();
    TestModel2.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel2.update_E_in_structure();
    TestModel2.solve_E();
    TestModel3.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel3.update_E_in_structure();
    TestModel3.solve_E();
    //TestModel.output_to_file(save_position, 0, 0);
    //CStr.output_to_file(save_position, 0);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    int num_model = 10;
    int start_num = 0;

    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    for (int i = 0; i <= num_model - 1; i++) {
        number = round(((double)rand() / RAND_MAX) * (maxnumber - minnumber + 1)) + minnumber;
        SpacePara spacepara_tmp(&S, bind, number, limitx1, limitx2, limity1, limity2, limitz1, limitz2, spacepara.get_geometryPara());
        CStr.UpdateStr(&spacepara_tmp);
        CStr.output_to_file(save_position, start_num + i + 1, "Simple");
        TestModel1.UpdateAlpha();
        TestModel1.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel1.update_E_in_structure();
        TestModel1.solve_E();
        TestModel1.output_to_file(save_position, lam1, start_num + i + 1);

        TestModel2.UpdateAlpha();
        TestModel2.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel2.update_E_in_structure();
        TestModel2.solve_E();
        TestModel2.output_to_file(save_position, lam2, start_num + i + 1);

        TestModel3.UpdateAlpha();
        TestModel3.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
        TestModel3.update_E_in_structure();
        TestModel3.solve_E();
        TestModel3.output_to_file(save_position, lam3, start_num + i + 1);
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();





    return 0;

}
//random generation for  focal based str
int main() {

    srand((unsigned)(time(0)));

    int Nx, Ny, Nz;
    Nx = 32; Ny = 32; Nz = 8;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;

    Vector3d center;
    Vector3d l;

    d = 20;

    center << double(Nx - 1) / 2, double(Ny - 1) / 2, double(Nz - 1) / 2;
    l << Nx - 1, Ny - 1, Nz - 1;

    Structure s1(S.get_total_space(), l, center);



    S = S + s1;

    double lam = 700;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");




    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;


    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    Vector3i bind(1, 1, 1);

    SpacePara spacepara(&S, bind, "ONES");

    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    string save_position = ".\\32328-topogenerated-dataset\\";

    //TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    //TestModel.update_E_in_structure();
    //TestModel.solve_E();
    //TestModel.output_to_file(save_position, 0, 0);
    //CStr.output_to_file(save_position, 0);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    int num_var = 500;
    int min_num = 1;       //Smallest number of focus
    int max_num = 4;       //Biggest number of focus
    bool sym = false;
    int start_num = 0;
    int max_evo = 10;
    Vector3d lower_bound, upper_bound;

    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    lower_bound << center(0), center(1), 4;
    upper_bound << center(0), center(1), 4;
    Evo_Focus(&spacepara, &CStr, &TestModel, save_position, start_num, max_evo, 1, 1, lower_bound, upper_bound, sym);


    lower_bound << 1, 1, 1;
    upper_bound << Nx - 2, Ny - 2, Nz - 2;
    for (int i = 0; i <= num_var - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, "ONES", spacepara.get_geometryPara());
        start_num += max_evo;
        Evo_Focus(&spacepara_tmp, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);
    }

    for (int i = 0; i <= num_var - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, "RANDOM", spacepara.get_geometryPara());
        start_num += max_evo;
        Evo_Focus(&spacepara_tmp, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);
    }

    sym = true;
    max_num = 2;
    for (int i = 0; i <= num_var - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, "ONES", spacepara.get_geometryPara());
        start_num += max_evo;
        Evo_Focus(&spacepara_tmp, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);
    }

    for (int i = 0; i <= num_var - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, "RANDOM", spacepara.get_geometryPara());
        start_num += max_evo;
        Evo_Focus(&spacepara_tmp, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);
    }


    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();





    return 0;

}


//############################################################TMP###########################################################
//2021-11-5
int main() {

    srand((unsigned)(time(0)));

    int Nx, Ny, Nz;
    Nx = 32; Ny = 32; Nz = 8;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;

    Vector3d center;
    Vector3d l;

    d = 20;

    center << double(Nx - 1) / 2, double(Ny - 1) / 2, double(Nz - 1) / 2;
    l << Nx - 1, Ny - 1, Nz - 1;

    Structure s1(S.get_total_space(), l, center);



    S = S + s1;

    double lam = 700;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");




    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;


    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    Vector3i bind(4, 4, 4);

    SpacePara spacepara(&S, bind, "ONES");

    CoreStructure CStr(&spacepara, d);
    AProductCore Core(&CStr, lam, material, "LDR");
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    string save_position = ".\\32328-topogenerated-dataset-diffbind\\";

    //TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    //TestModel.update_E_in_structure();
    //TestModel.solve_E();
    //TestModel.output_to_file(save_position, 0, 0);
    //CStr.output_to_file(save_position, 0);
    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;

    int num_var = 500;
    int min_num = 1;       //Smallest number of focus
    int max_num = 1;       //Biggest number of focus
    bool sym = false;
    int start_num = 20;
    int max_evo = 10;
    Vector3d lower_bound, upper_bound;

    ofstream TotalTime;
    TotalTime.open(save_position + "TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();

    lower_bound << center(0) - 2, center(1) - 2, 4;
    upper_bound << center(0) - 2, center(1) - 2, 4;
    Evo_Focus(&spacepara, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);


    /*
    lower_bound << 1, 1, 1;
    upper_bound << Nx - 2, Ny - 2, Nz - 2;
    for (int i = 0; i <= num_var - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, "ONES", spacepara.get_geometryPara());
        start_num += max_evo;
        Evo_Focus(&spacepara_tmp, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);
    }

    for (int i = 0; i <= num_var - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, "RANDOM", spacepara.get_geometryPara());
        start_num += max_evo;
        Evo_Focus(&spacepara_tmp, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);
    }

    sym = true;
    max_num = 2;
    for (int i = 0; i <= num_var - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, "ONES", spacepara.get_geometryPara());
        start_num += max_evo;
        Evo_Focus(&spacepara_tmp, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);
    }

    for (int i = 0; i <= num_var - 1; i++) {
        SpacePara spacepara_tmp(&S, bind, "RANDOM", spacepara.get_geometryPara());
        start_num += max_evo;
        Evo_Focus(&spacepara_tmp, &CStr, &TestModel, save_position, start_num, max_evo, min_num, max_num, lower_bound, upper_bound, sym);
    }
    */

    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_end - t_start).count();
    TotalTime << duration / 1000 << endl;
    TotalTime.close();





    return 0;

}