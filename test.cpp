#include "definition.h"



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
    Nx = 102; Ny = 102; Nz = 15;
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
    l << 49.0, 49.0, 0.0;
    center << 24.5, 24.5, 14.0;
    //l << 39.0, 39.0, 0.0;
    //center << 19.5, 19.5, 11.0;
    //l << 79.0, 79.0, 11.0;
    //center << 39.5, 39.5, 5.5;
    direction<<0,0,-1;
    int times = 15;
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

    
    double d = 20;

    double lam = 500;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 100;

    double focus = 325;   //nm       

    //Vector3d r;
    //r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 100;

    list<string> ObjectFunctionNames{ "ExtSurfaceEExp_CPU" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter1{ focus, exponent, ratio };
    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter1 };
    string save_position = "./275-290-2D-sym-20/";


    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);

    TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");



     
    
    return 0;
    
}



