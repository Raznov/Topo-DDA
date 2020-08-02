#include "definition.h"



int main(){

    /*
    int Nx,Ny,Nz;
    Nx=81;Ny=81;Nz=32;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);
    
    
    Structure s1(S.get_total_space(), "Lens", 0);
    //Structure s0(S.get_total_space(), "0", 0);
    */

    
    int Nx,Ny,Nz;
    Nx=81;Ny=81;Nz=23;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    l<<79.0,79.0,22.0;
    center<<39.5,39.5,11.0;
    //direction<<0,0,1;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);
    
    
    
    
    S=S+s1;
    //S=S+s0;
    

    double d=25;

    double lam=500;
    Vector3d n_K;
    n_K<<0.0,0.0,1.0;
    double E0=1.0;
    Vector3d n_E0;
    n_E0<<1.0,0.0,0.0;
    Vector2cd material=Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon=1;
    
    double focus = 1600;   //nm       

    //Vector3d r;
    //r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA=10000;
    double MAX_ERROR=0.00001;
    int MAX_ITERATION_EVO=100;
        
    list<string> ObjectFunctionNames{ "ExtSurfaceEExp_CPU"};
    
    double exponent = 2;
    int ratio = 4;

    list<double> ObjectParameter1{focus, exponent, ratio};
    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};
    
    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{&ObjectParameter1};
    string save_position="./550thick-ExtSurfaceEExp-focus1600-eps1-ratio4-WithPathRecord/";

    
    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material);

    TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");
    

    
    
    return 0;

}



//for a traditional lens
/*
int main(){
    int Nx,Ny,Nz,N;
    Nx=30;
    Ny=30;
    Nz=30;
    N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    //Vector3d center;
    //center<<20.0,20.0,20.0;
    //Structure s1(S.get_total_space(), "ONES", 10.0, center, 0);
    
    Structure s1(S.get_total_space(), "Lens", 0);
    S=S+s1;
    S.show_something_about_Structures();

    double d=25;

    double lam=500;
    Vector3d n_K;
    n_K<<0.0,0.0,1.0;
    double E0=1.0;
    Vector3d n_E0;
    n_E0<<1.0,0.0,0.0;
    Vector2cd material=Get_2_material("Air", "SiO2", lam, "nm");
    Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    cout<<test_model.get_N()<<endl;
    test_model.bicgstab(1000, 0.00001);
    //cout<<"fuck"<<endl;
    test_model.update_E_in_structure();
    test_model.solve_E();
    test_model.output_to_file();
}
*/

