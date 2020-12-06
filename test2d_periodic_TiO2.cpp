#include "definition.h"

int main(){


    
    int Nx, Ny, Nz;
    Nx = 22; Ny =22; Nz = 10;
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
    l << (Nx-1.0), (Ny-1.0), 1.0;
    center << (Nx-1.0)/2.0, (Nx-1.0)/2.0, 0.0;
    direction<<0,0,1;
    int times = Nz;
    Structure s1(S.get_total_space(), "CIRCLE", l, center, 1);
    Structure s2(S.get_total_space(), &s1, direction, times, 2);



    S = S + s1;
    S = S + s2;
    //S=S+s0;
    //S.show_something_about_Structures();

    
    double d = 15;

    double lam = 542;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "TiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon = 0.2;

    double focus = 295;   //nm       

    //Vector3d r;
    //r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA = 50000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 500;

    list<string> ObjectFunctionNames{ "IntegratedE" };
    list<double> ObjectParameter1{ 5*d, 20*d, 5*d };
    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};
    bool HavePathRecord = false;
    bool HavePenalty = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter1 };
    string save_position = "./Test2d_integratedE_dielTiO2_d15_p330_period50_eps0.2_Adam_start_at_circle_beta_0.7_gamma_0.01_fullArea/";

    double Lm = 22*d;
    double Ln = 22*d;
    int Mmax = 50;
    int Nmax = 50;

    AProductCore Core(&S, d, lam, material,Mmax, Nmax, Lm, Ln, "FCD");
    DDAModel Model0(&Core, n_K, E0, n_E0);
    list<DDAModel*> ModelList;
    ModelList.push_back(&Model0);
    EvoDDAModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &Core, ModelList);

    TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");



     
    
    return 0;
    
}