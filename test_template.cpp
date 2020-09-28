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

//Template for simple inverse design(3D in it)
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

    list<double> ObjectParameter1{ focus, exponent, ratio };
    //list<double> ObjectParameter2{center(0)*d,center(1)*d,focus};

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter1 };
    string save_position = "";


    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);

    TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}

//Template for simple inverse design(2D in it with same scope as previous 3D)
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
    l << 80.0, 80.0, 1.0;
    center << 40.0, 40.0, 12.0;
    //l << 79.0, 79.0, 11.0;
    //center << 39.5, 39.5, 5.5;
    direction << 0, 0, -1;
    int times = 13;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);
    //Structure s2(S.get_total_space(), &s1, 1);
    //Structure s3(S.get_total_space(), &s1, 2);
    //Structure s4(S.get_total_space(), &s1, 3);
    Structure s5(S.get_total_space(), &s1, direction, times, 2);
    //Structure s6(S.get_total_space(), &s2, direction, times, 2);
    //Structure s7(S.get_total_space(), &s3, direction, times, 2);
    //Structure s8(S.get_total_space(), &s4, direction, times, 2);



    S = S + s1;
    //S = S + s2;
    //S = S + s3;
    //S = S + s4;
    S = S + s5;
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
    list<double> ObjectParameter2{center(0)*d,center(1)*d,focus};

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter2 };
    string save_position = "./test-2D/";


    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);

    TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}


