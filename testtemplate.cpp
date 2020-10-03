
//General template for inverse design of 2D sym or sym
int main(){

    /*
    int Nx,Ny,Nz;
    Nx=81;Ny=81;Nz=31;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);


    Structure s1(S.get_total_space(), "Lens", 0);
    //Structure s0(S.get_total_space(), "0", 0);
    */

    
    int Nx, Ny, Nz;
    //Nx = 102; Ny = 102; Nz = 15;
    Nx = 83; Ny = 83; Nz = 13;
    //Nx = 103;
    //Ny = 103;
    //Nz = 16;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    //l << 49.0, 49.0, 0.0;
    //center << 24.5, 24.5, 14.0;
    l << 39.0, 39.0, 0.0;
    center << 19.5, 19.5, 11.0;
    //l << 79.0, 79.0, 11.0;
    //center << 39.5, 39.5, 5.5;
    direction<<0,0,-1;
    int times = 12;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);
    Structure s2(S.get_total_space(), &s1, 1);
    Structure s3(S.get_total_space(), &s1, 2);
    Structure s4(S.get_total_space(), &s1, 3);
    Structure s5(S.get_total_space(), &s1, direction, times, 2);
    Structure s6(S.get_total_space(), &s2, direction, times, 2);
    Structure s7(S.get_total_space(), &s3, direction, times, 2);
    Structure s8(S.get_total_space(), &s4, direction, times, 2);



    S = S + s1;
    S = S + s2;
    S = S + s3;
    S = S + s4;
    S = S + s5;
    S = S + s6;
    S = S + s7;
    S = S + s8;
    //S=S+s0;
    //S.show_something_about_Structures();

    
    double d = 25;

    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;

    double focus = 295;   //nm       

    //Vector3d r;
    //r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 50;

    list<string> ObjectFunctionNames{ "ExtSurfaceEExp_CPU" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter1{ focus, exponent, ratio };
    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter1 };
    string save_position = "";


    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);

    TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");



     
    
    return 0;
    
}


//Template for doing parameter sweep for a DDA str
int main() {



    double regionx, regiony, regionz;
    regionx = 1000;
    regiony = 1000;
    regionz = 300;

    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;
    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 2;

    list<string> ObjectFunctionNames{ "ExtSurfaceEExp_CPU" };

    double exponent = 2;
    double ratio = 4;


    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    string save_position = "";

    ofstream objectives;
    objectives.open("objectives.txt");
    for (double focus = 300; focus <= 380; focus += 5) {
        for (double d = 15; d <= 25; d += 2) {
            int Nx, Ny, Nz;
            Nx = int(round(regionx / d)) + 5;
            Ny = int(round(regiony / d)) + 5;
            Nz = int(round(regionz / d)) + 1;

            int N = 0;
            VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
            list<Structure> ln;

            Space S(&total_space, Nx, Ny, Nz, N, &ln);

            Vector3d l;
            Vector3d center;
            l << round(regionx / d), round(regiony / d), round(regionz / d);
            center << l(0) / 2, l(1) / 2, l(2) / 2;
            Structure s1(S.get_total_space(), "ONES", l, center, 1);
            S = S + s1;

            list<double> ObjectParameter1{ focus, exponent, ratio };
            list<list<double>*> ObjectParameters{ &ObjectParameter1 };

            EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);
            string tmp = to_string(int(focus)) + to_string(int(d));
            cout << "tmp: " << tmp << endl;
            int Name = stoi(tmp);
            cout << "Name: " << Name << endl;
            double obj = TestModel.CalTheObjForSingleStr(MAX_ITERATION_DDA, MAX_ERROR, Name);
            cout << "Objective: " << obj << endl;
            objectives << focus << " " << d << " " << obj << endl;

        }
    }
    objectives.close();










    return 0;

}

//Template for import str and calculate
int main() {

    int Nx, Ny, Nz;
    Nx = 83; Ny = 83; Nz = 13;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);


    Structure s1(S.get_total_space(), "55", 0);
    S = S + s1;

    double d = 25;
    double focus = 350;
    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;
    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 2;

    list<string> ObjectFunctionNames{ "ExtSurfaceEExp_CPU" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter1{ focus, exponent, ratio };
    list<list<double>*> ObjectParameters{ &ObjectParameter1 };



    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    string save_position = "";

    int Name = 666;

    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);
    TestModel.output_to_file();
    double obj = TestModel.CalTheObjForSingleStr(MAX_ITERATION_DDA, MAX_ERROR, Name);
    cout << "Objective: " << obj << endl;
    //objectives << focus << " " << d << " " << obj << endl;










    return 0;

}

//Model verification
int main() {

    /*
    int Nx,Ny,Nz;
    Nx=81;Ny=81;Nz=31;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Structure s1(S.get_total_space(), "Lens", 0);
    //Structure s0(S.get_total_space(), "0", 0);
    */


    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
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

    Structure s1(S.get_total_space(), "ONES", r, center, 1);



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
    Vector2cd material = Get_2_material("Air", "Au", lam, "nm");

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
    //string AMatrixMethod = "LDR";

    Model TestModel(&S, d, lam, n_K, E0, n_E0, material);
    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();
    TestModel.output_to_file();

    double x, y, z;
    y = center(1) * d;
    z = center(2) * d;


    N = TestModel.get_N();
    VectorXcd* P = TestModel.get_P();
    VectorXi* R = TestModel.get_R();
    double K = 2 * M_PI / lam;
    Vector3cd E_sum = Vector3cd::Zero();
    Vector3cd E_ext = Vector3cd::Zero();

    string name = "ModelEfieldscand=" + to_string(d) + ".txt";
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
            Matrix3cd A = TestModel.A_dic_generator(rx, ry, rz);
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

//DDAModel verification
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

    Structure s1(S.get_total_space(), "ONES", r, center, 1);



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
    Vector2cd material = Get_2_material("Air", "Au", lam, "nm");

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

    AProductCore Core(&S, d, lam, material);
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

//Template for EvoModel test
int main() {

    /*
    int Nx,Ny,Nz;
    Nx=81;Ny=81;Nz=31;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    Structure s1(S.get_total_space(), "Lens", 0);
    //Structure s0(S.get_total_space(), "0", 0);
    */


    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = 83; Ny = 83; Nz = 13;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    //l << 50.0, 50.0, 15.0;
    //center << 25.0, 25.0, 7.5;
    l << 80.0, 80.0, 12.0;
    center << 40.0, 40.0, 6.0;
    //l << 79.0, 79.0, 11.0;
    //center << 39.5, 39.5, 5.5;
    direction << 0, 0, -1;
    int times = 15;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);
    //Structure s2(S.get_total_space(), &s1, 1);
    //Structure s3(S.get_total_space(), &s1, 2);
    //Structure s4(S.get_total_space(), &s1, 3);
    //Structure s5(S.get_total_space(), &s1, direction, times, 2);
    //Structure s6(S.get_total_space(), &s2, direction, times, 2);
    //Structure s7(S.get_total_space(), &s3, direction, times, 2);
    //Structure s8(S.get_total_space(), &s4, direction, times, 2);



    S = S + s1;
    //S = S + s2;
    //S = S + s3;
    //S = S + s4;
    //S = S + s5;
    //S = S + s6;
    //S = S + s7;
    //S = S + s8;
    //S=S+s0;
    //S.show_something_about_Structures();


    double d = 25;

    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;

    double focus = 350;   //nm       

    //Vector3d r;
    //r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 200;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    //list<double> ObjectParameter1{ focus, exponent, ratio };
    list<double> ObjectParameter2{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter2 };
    string save_position = "";


    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);

    TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//Template for EvoDDAModel test
int main() {

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = 83; Ny = 83; Nz = 13;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    //l << 50.0, 50.0, 15.0;
    //center << 25.0, 25.0, 7.5;
    l << 80.0, 80.0, 12.0;
    center << 40.0, 40.0, 6.0;
    //l << 79.0, 79.0, 11.0;
    //center << 39.5, 39.5, 5.5;
    direction << 0, 0, -1;
    int times = 15;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);
    //Structure s2(S.get_total_space(), &s1, 1);
    //Structure s3(S.get_total_space(), &s1, 2);
    //Structure s4(S.get_total_space(), &s1, 3);
    //Structure s5(S.get_total_space(), &s1, direction, times, 2);
    //Structure s6(S.get_total_space(), &s2, direction, times, 2);
    //Structure s7(S.get_total_space(), &s3, direction, times, 2);
    //Structure s8(S.get_total_space(), &s4, direction, times, 2);



    S = S + s1;
    //S = S + s2;
    //S = S + s3;
    //S = S + s4;
    //S = S + s5;
    //S = S + s6;
    //S = S + s7;
    //S = S + s8;
    //S=S+s0;
    //S.show_something_about_Structures();


    double d = 25;

    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;

    double focus = 350;   //nm       

    //Vector3d r;
    //r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 10;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    //list<double> ObjectParameter1{ focus, exponent, ratio };
    list<double> ObjectParameter2{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter2 };
    string save_position = "";

    AProductCore Core(&S, d, lam, material);
    DDAModel Model1(&Core, n_K, E0, n_E0);

    list<DDAModel*> ModelList;
    ModelList.push_back(&Model1);

    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &Core, ModelList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//Template for angle sweep
int main() {

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = 83; Ny = 83; Nz = 13;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    l << 80.0, 80.0, 12.0;
    center << 40.0, 40.0, 6.0;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);



    S = S + s1;


    double d = 25;

    double lam = 500;

    double E0 = 1.0;

    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    double epsilon = 100;

    double focus = 350;   //nm       


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 100;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = "";
    double angle = 10 * M_PI / 180;
    Vector3d n_K;
    Vector3d n_E0;

    AProductCore Core(&S, d, lam, material);

    n_K << 0.0, 0.0, 1.0;
    n_E0 << 1.0, 0.0, 0.0;
    DDAModel Model0(&Core, n_K, E0, n_E0);
    n_K << sin(angle), 0.0, cos(angle);
    n_E0 << cos(angle), 0.0, -sin(angle);
    DDAModel Model1(&Core, n_K, E0, n_E0);
    n_K << -sin(angle), 0.0, cos(angle);
    n_E0 << cos(angle), 0.0, sin(angle);
    DDAModel Model2(&Core, n_K, E0, n_E0);
    //n_K << 0.0, sin(angle), cos(angle);
    //n_E0 << 1.0, 0.0, 0.0;
    //DDAModel Model3(&Core, n_K, E0, n_E0);
    //n_K << 0.0, sin(angle), -cos(angle);
    //n_E0 << 1.0, 0.0, 0.0;
    //DDAModel Model4(&Core, n_K, E0, n_E0);

    list<DDAModel*> ModelList;
    ModelList.push_back(&Model0);
    ModelList.push_back(&Model1);
    ModelList.push_back(&Model2);
    //ModelList.push_back(&Model3);
    //ModelList.push_back(&Model4);

    EvoDDAModel EModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &Core, ModelList);


    EModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//Template for traditional lens

int main() {

    int Nx, Ny, Nz;
    Nx = 83; Ny = 83; Nz = 34;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);


    Structure s1(S.get_total_space(), "Lens", 0);
    S = S + s1;

    double d = 25;
    double focus = 1600;
    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;
    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 2;

    list<string> ObjectFunctionNames{ "PointE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 40 * d,40 * d,focus };
    list<list<double>*> ObjectParameters{ &ObjectParameter };



    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    string save_position = "";

    int Name = 0;

    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);
    TestModel.output_to_file();
    double obj = TestModel.CalTheObjForSingleStr(MAX_ITERATION_DDA, MAX_ERROR, Name);
    cout << "Objective: " << obj << endl;
    //objectives << focus << " " << d << " " << obj << endl;










    return 0;

}