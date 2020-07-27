#include "definition.h"

int main(){

    /*   
    int Nx,Ny,Nz;
    Nx=20;Ny=20;Nz=20;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3d center;
    center<<10.0,10.0,10.0;
    Structure s1(S.get_total_space(), "ZEROS", 5.0, center, 0);

    Vector3i direction;
    Vector3d l;
    l<<10.0,10.0,0.0;
    center<<10.5,10.5,0.0;
    direction<<0,0,1;
    Structure s2(S.get_total_space(), "RANDOM", l, center, 1);
    cout<<"fuck"<<endl;
    Structure s3(S.get_total_space(), &s2, direction, 2, 2);

    center<<5.0,5.0,5.0;
    Structure s4(S.get_total_space(), "ONES", 2.5, center, 0);

    //int times=2;
    //Structure s2(S.get_total_space(), &s1, direction, times);
    S=S+s1+s2+s3+s4;
    //S=S+s2;
    //S=S+s3;
    //S=S+s4;
    cout<<S.get_total_space()<<endl;
    cout<<S.get_ln_size()<<endl;
    S.show_something_about_Structures();
    */

    



    int Nx,Ny,Nz,N;
    Nx=200;
    Ny=200;
    Nz=4;
    N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    //Vector3d center;
    //center<<20.0,20.0,20.0;
    //Structure s1(S.get_total_space(), "ONES", 10.0, center, 0);
    
    Structure s1(S.get_total_space(), "Lens", 0);
    S=S+s1;
    double d=25;


    
    //For X-Y plane
    
    double plotz = 50;
    double plotx = 50;
    double ploty = 50;
    int NPx = int(round(plotx/d));
    int NPy = int(round(ploty/d));
    VectorXi PlotGeometry = VectorXi::Zero(3 * NPx * NPy);
    int position = 0;
    for(int i = 0; i<=NPx - 1; i++){
        for(int j = 0; j<=NPy - 1; j++){
            PlotGeometry(3*position) = i + int(round((Nx - NPx)/2));
            PlotGeometry(3*position+1) =  j + int(round((Ny - NPy)/2));
            PlotGeometry(3*position+2) = int(round(plotz/d));
            position = position + 1;
            //cout<<PlotGeometry(3*(i+j))<<" "<<PlotGeometry(3*(i+j)+1)<<" "<<PlotGeometry(3*(i+j)+2)<<" "<<endl;
        }
    }
    

    /*
    //For X-Z plane
    
    double plotz = 1000;
    double plotx = 1000;
    double ploty = 2460;
    int NPz = int(round(plotz/d));
    int NPx = int(round(plotx/d));
    VectorXi PlotGeometry = VectorXi::Zero(3 * NPz * NPx);
    int position = 0;
    for(int i = 0; i<=NPz - 1; i++){
        for(int j = 0; j<=NPx - 1; j++){
            PlotGeometry(3*position) = j + int(round((Nx - NPx)/2));
            PlotGeometry(3*position+1) = int(round(ploty/d)); 
            PlotGeometry(3*position+2) = i;
            position = position + 1;
            //cout<<PlotGeometry(3*(i+j))<<" "<<PlotGeometry(3*(i+j)+1)<<" "<<PlotGeometry(3*(i+j)+2)<<" "<<endl;
        }
    } 
    */
    
    

    //For Y-Z plane
    /*
    double plotz = 1000;
    double plotx = 2460;
    double ploty = 1000;
    int NPz = int(round(plotz/d));
    int NPy = int(round(ploty/d));
    VectorXi PlotGeometry = VectorXi::Zero(3 * NPz * NPy);
    int position = 0;
    for(int i = 0; i<=NPz - 1; i++){
        for(int j = 0; j<=NPy - 1; j++){
            PlotGeometry(3*position) = int(round(plotx/d)); 
            PlotGeometry(3*position+1) =  j + int(round((Ny - NPy)/2));
            PlotGeometry(3*position+2) = i;
            position = position + 1;
            //cout<<PlotGeometry(3*(i+j))<<" "<<PlotGeometry(3*(i+j)+1)<<" "<<PlotGeometry(3*(i+j)+2)<<" "<<endl;
        }
    }
    */
    
    
    double lam=550;
    Vector3d n_K;
    n_K<<0.0,0.0,-1.0;
    double E0=1.0;
    Vector3d n_E0;
    n_E0<<1.0,0.0,0.0;
    Vector2cd material=Get_2_material("Air", "SiO2", lam, "nm");
    Model test_model(&S, d, lam, n_K, E0, n_E0, material, &PlotGeometry);
    cout<<test_model.get_N()<<endl;
    test_model.bicgstab(1000, 0.00001);
    //cout<<"fuck"<<endl;
    test_model.solve_E();
    test_model.output_to_file();
    
    




    /*
    int Nx,Ny,Nz;
    Nx=41;Ny=41;Nz=4;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    l<<40.0,40.0,3.0;
    center<<20.0,20.0,1.5;
    //direction<<0,0,1;
    Structure s1(S.get_total_space(), "ZEROS", l, center, 1);
    //Structure s2(S.get_total_space(), &s2, direction, 3, 2);
    S=S+s1;
    S.show_something_about_Structures();
    double d=25;

    //For X-Y plane
    
    double plotz = 500;
    double plotx = 1000;
    double ploty = 1000;
    int NPx = int(round(plotx/d));
    int NPy = int(round(ploty/d));
    VectorXi PlotGeometry = VectorXi::Zero(3 * NPx * NPy);
    int position = 0;
    for(int i = 0; i<=NPx - 1; i++){
        for(int j = 0; j<=NPy - 1; j++){
            PlotGeometry(3*position) = i + int(round((Nx - NPx)/2));
            PlotGeometry(3*position+1) =  j + int(round((Ny - NPy)/2));
            PlotGeometry(3*position+2) = int(round(plotz/d));
            position = position + 1;
            //cout<<PlotGeometry(3*(i+j))<<" "<<PlotGeometry(3*(i+j)+1)<<" "<<PlotGeometry(3*(i+j)+2)<<" "<<endl;
        }
    }




    double lam=550;
    Vector3d n_K;
    n_K<<0.0,0.0,-1.0;
    double E0=1.0;
    Vector3d n_E0;
    n_E0<<1.0,0.0,0.0;
    Vector2cd material=Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon=1000;

    Vector3d r;
    r<<20.0*d, 20.0*d, 500.0;
    int MAX_ITERATION_DDA=10000;
    double MAX_ERROR=0.00001;
    int MAX_ITERATION_EVO=20;
        
    list<string> ObjectFunctionNames{ "PointE"};
    
    list<double> ObjectParameter1{r(0), r(1), r(2)};
    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};
    
    string PenaltyType = "No Penalty";

    list<list<double>*> ObjectParameters{&ObjectParameter1};
    string save_position="./404003-202015-25-785-01-PointE-202020/";
    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, PenaltyType, save_position, &S, d, lam, n_K, E0, n_E0, material, &PlotGeometry);
    TestModel.EvoOptimization(epsilon, MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO);
    
    
    //complex<double> diel=-2.54038+3.61194*1.0i;
    //double K=2*M_PI/lam;
    //cout<<endl<<endl;
    //complex<double> tmp=Get_Alpha(lam, K, d, diel);
    
    //cout<<(1.0/tmp)<<endl;
    
    */
    
    return 0;

}



//focus len design with GPU. Total time should between 1~2h
int main(){


    int Nx,Ny,Nz;
    Nx=81;Ny=81;Nz=12;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    Vector3d l;
    Vector3d center;
    l<<80.0,80.0,11.0;
    center<<40.5,40.5,6.0;
    //direction<<0,0,1;
    Structure s1(S.get_total_space(), "ONES", l, center, 1);

    S=S+s1;

    double d=25;

    double lam=500;
    Vector3d n_K;
    n_K<<0.0,0.0,1.0;
    double E0=1.0;
    Vector3d n_E0;
    n_E0<<1.0,0.0,0.0;
    Vector2cd material=Get_2_material("Air", "SiO2", lam, "nm");
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    double epsilon=100;
    
    double focus = 2500;   //nm       

    Vector3d r;
    r<<center(0)*d, center(1)*d, focus;


    int MAX_ITERATION_DDA=10000;
    double MAX_ERROR=0.00001;
    int MAX_ITERATION_EVO=100;
        
    list<string> ObjectFunctionNames{ "PointE"};
    
    list<double> ObjectParameter1{r(0), r(1), r(2)};
    //list<double> ObjectParameter2{r(0), r(1), r(2), l(0)*d, l(1)*d, d};
    
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;

    list<list<double>*> ObjectParameters{&ObjectParameter1};
    string save_position="./818112-4054056-25-500-focus2500/";
    //Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    //test_model.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    //test_model.solve_E();
    //test_model.output_to_file();



    double plotz = focus;
    double plotx = 2000;
    double ploty = 2000;
    int NPx = int(round(plotx/d));
    int NPy = int(round(ploty/d));
    VectorXi PlotGeometry = VectorXi::Zero(3 * NPx * NPy);
    int position = 0;
    for(int i = 0; i<=NPx - 1; i++){
        for(int j = 0; j<=NPy - 1; j++){
            PlotGeometry(3*position) = i + int(round((Nx - NPx)/2));
            PlotGeometry(3*position+1) =  j + int(round((Ny - NPy)/2));
            PlotGeometry(3*position+2) = int(round(plotz/d));
            position = position + 1;
            //cout<<PlotGeometry(3*(i+j))<<" "<<PlotGeometry(3*(i+j)+1)<<" "<<PlotGeometry(3*(i+j)+2)<<" "<<endl;
        }
    }

    
    EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material, &PlotGeometry);
    TestModel.EvoOptimization(epsilon, MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");
    
    
    //complex<double> diel=-2.54038+3.61194*1.0i;
    //double K=2*M_PI/lam;
    //cout<<endl<<endl;
    //complex<double> tmp=Get_Alpha(lam, K, d, diel);
    
    //cout<<(1.0/tmp)<<endl;
    
    
    return 0;

}


//A SiO2 sphere
int main(){


    int Nx,Ny,Nz;
    Nx=100;Ny=100;Nz=100;
    int N=0;
    VectorXi total_space=build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;
    double r;
    Vector3d center;
    r = 45.0;
    center<<45.0,45.0,45.0;
    
    Structure s1(S.get_total_space(), "ONES", r, center, 1);

    S=S+s1;

    double d=1.0;

    double lam=508;
    Vector3d n_K;
    n_K<<0.0,0.0,-1.0;
    double E0=1.0;
    Vector3d n_E0;
    n_E0<<1.0,0.0,0.0;
    Vector2cd material=Get_2_material("Air", "SiO2", lam, "nm");
    

    int MAX_ITERATION_DDA=10000;
    double MAX_ERROR=0.00001;

    
    Model test_model(&S, d, lam, n_K, E0, n_E0, material);
    test_model.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    test_model.solve_E();
    test_model.output_to_file();
    return 0;

}


