#ifndef PREPROCESSING_H_INCLUDED
#define PREPROCESSING_H_INCLUDED
#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <chrono>
#include <omp.h>
#include <string>
#include <cstdlib>
#include <map>
#include <list>
#include <tuple>

#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Core"
#include "eigen/Eigen/IterativeLinearSolvers"
#include "eigen/Eigen/Sparse"
#include "eigen/unsupported/Eigen/CXX11/Tensor"

//#include "cufftw.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cufft.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

//kernel wrappers:
void A2As(double *A, cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, int NxFFT, int NyFFT, int NzFFT);
void B2Bs(double *bDev, cufftDoubleComplex *bxDev, cufftDoubleComplex *byDev, cufftDoubleComplex *bzDev, int NxFFT, int NyFFT, int NzFFT);
void Conv(cufftDoubleComplex *Convx, cufftDoubleComplex *Convy, cufftDoubleComplex *Convz, cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, cufftDoubleComplex *bxDev, cufftDoubleComplex *byDev, cufftDoubleComplex *bzDev, int NxFFT, int NyFFT, int NzFFT);
void Conv2B(cufftDoubleComplex *Convx, cufftDoubleComplex *Convy, cufftDoubleComplex *Convz, double *bDev, int NxFFT, int NyFFT, int NzFFT);
void APtoESum(cufftDoubleComplex *A00, cufftDoubleComplex *A01, cufftDoubleComplex *A02, cufftDoubleComplex *A11, cufftDoubleComplex *A12, cufftDoubleComplex *A22, cufftDoubleComplex *PxDev, cufftDoubleComplex *PyDev, cufftDoubleComplex *PzDev, cufftDoubleComplex *ESumxDev, cufftDoubleComplex *ESumyDev, cufftDoubleComplex *ESumzDev, int NxFFT, int NyFFT, int NzFFT, int NxA, int NyA, int NzA, int index1, int index2, int index3, int deduction);

complex<double> Get_Alpha(double lam, double K, double d, complex<double> diel, Vector3d n_E0, Vector3d n_K);
complex<double> Get_Alpha_FCD(double lam,double K,double d,complex<double> diel);
complex<double> Get_material(string mat, double wl, string unit);                  //name of mat to get its diel function at certain wavlength              
Vector2cd Get_2_material(string sub, string mat, double wl, string unit);          //a wrapper for Get_material
double Average(VectorXcd* E, int N, double exponent);
double Get_Max(VectorXcd* E, int N);
double G(VectorXcd* E, int N, double exponent, double E0);
//ArrayXcd FFT(int nx,int ny,int nz,ArrayXcd in,int _direction);
VectorXi build_a_bulk(int Nx, int Ny, int Nz);

class Structure{
    private:
        VectorXd diel;                  //0~1
        VectorXi geometry;
        VectorXi geometry_dep;          //Same size with geometry. If para=2, stores the corresponding para position of each point. Has no meaning in para=0&1
        int para;                       //Parameter condition: 0-not parameter, 1-parameter, 2-duplicated from a 2D parameter, controlled by parameter.
    public:
        //--------------------------------------Dont support dependent para build-------------------------------------------------------
        Structure(VectorXi *total_space, VectorXd *diel_, VectorXi *geometry_, int para_);            //vanilla initialization

        //Sphere
        Structure(VectorXi *total_space, string initial_diel, double r, Vector3d center, int para_);  //r: actual radius/d. center: actual center/d. In charge of Sphere 
        
        //Circle
        Structure(VectorXi *total_space, string initial_diel, double r, Vector3i center, Vector3i direction, int para_); //build a circle, direction is its normalized direction in Cart. coord.

        //From file
        Structure(VectorXi* total_space, string FileName, int para_);                                                    //Read a structure from txt file
        
        //--------------------------------------Support dependent para build-------------------------------------------------------------
        //Bulk
        Structure(VectorXi *total_space, string initial_diel, Vector3d l, Vector3d center, int para_);    //Ractangular(both 2D and 3D). Para can only = 0&1 
        
        //Duplicate
        Structure(VectorXi *total_space, Structure *s, Vector3i direction, int times, int para_);                     //Initializa a Structure by duplicating a existing structure along a certain direction for several times. Direction is normalized and can only be alone x, y or z.
        //The original structure is not included. original structure + new structure = times * original structure. If set para=2 and original para=1, then depend on origin str as geometry_dep. If para=2 and original para=2, will copy origin geometry_dep.
        Structure(VectorXi* total_space, Structure* s, int dep_way);     //Special copy for setting up 2-fold and 4-fold symmetry dependence in xy plane. Para auto set to 2.
        //The original str at left down corner of the xy plane(left down corner is (0,0)). dep_way=1, 2, 3 corresponds to the other blocks in clock-wise. 


        //-------------------------------------Other member functions--------------------------------------------------------------------
        VectorXi *get_geometry();
        VectorXi* get_geometry_dep();
        VectorXd *get_diel();
        void cut(VectorXi *big, VectorXi *smalll, VectorXd *small_diel);
        int get_geometry_size();
        int get_para();
};

class Space{
    private:
        VectorXi *total_space;
        int Nx, Ny, Nz;
        int N;                        //total size of all the geometry inside the list(dipole size which is actual size/3)
        list<Structure> *ln;
    public:
        Space(VectorXi *total_space_, int Nx_, int Ny_, int Nz_, int N_, list<Structure> *ln_);
        VectorXi *get_total_space();
        int get_ln_size();
        tuple<int, int, int, int> get_Ns();
        list<Structure> *get_ln();
        void show_something_about_Structures() const;
        friend Space operator+(const Space &s1, Structure &s2);
        

};

//Abstract parent class for objective function.
class Objective;
class ObjectivePointE;
class ObjectiveSurfaceEExp;

class SiCi {
public:
    int numberSi;
    int numberCi;
    double disSi;
    double disCi;
    VectorXd Si;
    VectorXd Ci;
    SiCi();
    double get_Si(double x);
    double get_Ci(double y);

};


class Model{
    protected:
        int N;                        //Number of dipoles
        int Nx;                       //scope of space. Nx*Ny*Nz!=N
        int Ny;
        int Nz;
        int time;
        int ITERATION;
        double d;
        double E0;
        double K;
        double lam;
        double Error;
        VectorXi R;                      //Position of dipoles. Both R and RResult are unitless, so need to time d to get real number.
        VectorXi RDep;                   //Position of the dependent para points in space of the points in the same position as R
        list<list<int>> PositionDep;    //First D has the D of para(para=1). Second D is the positions of other points dependent on the para in the 1stD.
        VectorXi PositionPara;          //Position i(in R) of the parameters (3*i=x, 3*i+1=y, 3*i+2=z)
        bool RResultSwitch;               //0(false) for plot only E field on the structure points (fast), 1(true) for using RResult different from R to plot (slow but adjustable).
        VectorXi RResult;                //The position matrix for the EResult (where you want to plot the E field)
        list<int> para_nums;
        list<int> para_starts;
        list<int> para_dep_nums;
        list<int> para_dep_starts;
        Vector3d n_E0;
        Vector3d n_K;
        VectorXcd diel;                   //real diel after 0~1 corresponds to real numbers
        VectorXd diel_old;                //The 0~1 version of diel
        VectorXcd P;
        VectorXcd E;
        VectorXcd Einternal;              //E field on structure points
        VectorXcd EResult;                //E field on designated points
        Vector2cd material;
        VectorXcd al;                       // 1 over alpha instead of alpha.
        bool verbose;

        VectorXcd diel_max;                         //corresponds to the previous maximum obj
        VectorXd diel_old_max;
        VectorXcd P_max;
        VectorXcd al_max;

        int MAXm;                            //-MAXm<=m<=MAXm
        int MAXn;                             //-MAXn<=n<=MAXn
        double Lm;                         //desplacement vector for one period, Currently should only be in x and y direction; d included do not need to time d
        double Ln;
        
    public:
        Model(Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_);
        Model(Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_, VectorXi *RResult_);
        Model(Space* space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_, int MAXm_, int MAXn_, double Lm_, double Ln_);
        ~Model();
        Matrix3cd A_dic_generator(double x,double y,double z);
        Matrix3cd A_dic_generator(double x, double y, double z, int m, int n);
        VectorXcd Aproduct(VectorXcd &b);
        void bicgstab(int MAX_ITERATION,double MAX_ERROR);
        void change_E(VectorXcd E_);
        void reset_E();             //reset E to E0 
        //double get_step_length(VectorXd gradients, double epsilon);                                        
        void change_para_diel(VectorXd step);
        VectorXcd get_diel();
        Vector2cd get_material();
        VectorXcd get_alpha();
        tuple<list<int>, list<int>, list<int>, list<int>> get_para_info();
        int get_N();
        int get_Nx();
        int get_Ny();
        int get_Nz();
        VectorXcd* get_P();
        VectorXi* get_R();
        double get_d();
        Vector3d get_nE0();
        Vector3d get_nK();
        double get_E0();
        double get_wl();
        VectorXcd* get_Einternal();
      
        //Vector3d* get_nE0();
        void solve_E();                                                        //update the result E field on each dipole or on a designated space
        void update_E_in_structure();                                          //update the result E field on each dipole 
        void output_to_file();
        void output_to_file(string save_position, int iteration);              //especially used for EvoOptimization
        
        //FFT related variables;
        double *AHos;                              //A_dicDoubl
        double *ADev;                              // double, but actually double*2 course real and imag are both stored in this double
        cufftDoubleComplex *A00, *A01, *A02, *A11, *A12, *A22; //The components in A for FFT, only in device
        double *bHos;                              //b in Aproduct
        double *bDev;                              // double, but actually double*2 course real and imag are both stored in this double
        cufftDoubleComplex *bxDev, *byDev, *bzDev; //b components for FFT, only in device
        int NxFFT, NyFFT, NzFFT;                   //2*(Nx, Ny, Nz) - 1
        int NFFT;                                  //NxFFT*NyFFT*NzFFT
        cufftDoubleComplex *Convx, *Convy, *Convz; //convolution of A and b on device

        //FFT plan for all
        cufftHandle Plan;
};

class Objective;

class EvoModel : public Model{
    private:

        string save_position;
        
        list<list<double>*> *ObjectParameters;
        list<double> MajorObjectParameters;
        list<list<double>> MinorObjectParameters;
        
        list<string> *ObjectFunctionNames;
        string MajorObjectFunctionName;
        double MajorObjectFunctionResult;
        list<string> MinorObjectFunctionNames;
        list<double> MinorObjectFunctionResults;
        Objective* objective; 
        double origin;                               //Record the objective function for partial derivative (the value before change)   
        bool HavePenalty;
        double PenaltyFactor;

        double MaxObj;                                //Record the historical maximum obj func
        double epsilon_fix;
        double epsilon_tmp;                         //The epsilon used for calculation (can be different from the fixed input epsilon)
        bool HavePathRecord;
        int Stephold;
    public:
        EvoModel(list<string>* ObjectFunctionNames_, list<list<double>*>* ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, double PenaltyFactor_, string save_position_, Space* space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_);
        EvoModel(list<string>* ObjectFunctionNames_, list<list<double>*>* ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, double PenaltyFactor_, string save_position_, Space* space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_, VectorXi* RResult_);
        EvoModel(list<string>* ObjectFunctionNames_, list<list<double>*>* ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, double PenaltyFactor_, string save_position_, Space* space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_, int MAXm_, int MAXn_, double Lm_, double Ln_);
        
        //functions used to calculate partial derivatives                                 
        tuple<VectorXd, VectorXcd> devx_and_Adevxp(double epsilon);                       //partial derivative of obj to parameter and A to x times p
        VectorXcd devp(double epsilon);                       //partial derivative of obj to P. Size of P
        
        void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method);
        double CalTheObjForSingleStr(int MAX_ITERATION, double MAX_ERROR, int Name);                    //If you want to calculate the objective for single DDA structure.

        //The objective choosing function:
        double MajorObjective();
        list<double> MinorObjective();
        Objective* ObjectiveFactory(string ObjectName,  list<double> ObjectParameters);

        //functions that can be chosen as objective functions:
        double PointE(list<double> Parameter);                   //r is the real position in the Model coordinate, d should have been considered when defining r.
        double PointEWithPenalty(list<double> Parameter);
        double SurfaceEExp(list<double> Parameter);
        double SurfaceEExpWithPenalty(list<double> Parameter);
        double L1Norm();
        
        

};
//void EvoOptimization(Model *model, double epsilon, Vector3d r, int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO);

class AProductCore {
private:
    //---------------------------------Geometries-------------------------------
    Space* space;
    //---------------------------------Necessary values for A matrix-----------------------------------------------------
    int N;                        //Number of dipoles
    int Nx;                       //scope of space. Nx*Ny*Nz!=N
    int Ny;
    int Nz;
    double d;
    double K;
    double lam;
    VectorXi R;                      //Position of dipoles. Both R and RResult are unitless, so need to time d to get real number.

    //FFT related variables;
    double* AHos;                              //A_dicDoubl
    double* ADev;                              // double, but actually double*2 course real and imag are both stored in this double
    cufftDoubleComplex* A00, * A01, * A02, * A11, * A12, * A22; //The components in A for FFT, only in device
    double* bHos;                              //b in Aproduct
    double* bDev;                              // double, but actually double*2 course real and imag are both stored in this double
    cufftDoubleComplex* bxDev, * byDev, * bzDev; //b components for FFT, only in device
    int NxFFT, NyFFT, NzFFT;                   //2*(Nx, Ny, Nz) - 1
    int NFFT;                                  //NxFFT*NyFFT*NzFFT
    cufftDoubleComplex* Convx, * Convy, * Convz; //convolution of A and b on device

    //FFT plan for all
    cufftHandle Plan;

    //---------------------------------------Necessary for periodic A matrix----------------------------------------
    int MAXm;                            //-MAXm<=m<=MAXm
    int MAXn;                             //-MAXn<=n<=MAXn
    double Lm;                         //desplacement vector for one period, Currently should only be in x and y direction; d included do not need to time d
    double Ln;

    //--------------------------------Not necessary for A matrix but should be the same for diff DDAModel using the same A matrix------------------------------
    VectorXi RDep;                   //Position of the dependent para points in space of the points in the same position as R
    list<list<int>> PositionDep;    //First D has the D of para(para=1). Second D is the positions of other points dependent on the para in the 1stD.
    VectorXi PositionPara;          //Position i(in R) of the parameters (3*i=x, 3*i+1=y, 3*i+2=z)
    list<int> para_nums;
    list<int> para_starts;
    list<int> para_dep_nums;
    list<int> para_dep_starts;
    VectorXcd diel;                   //real diel after 0~1 corresponds to real numbers
    VectorXd diel_old;                //The 0~1 version of diel
    Vector2cd material;
    VectorXcd diel_max;                         //corresponds to the previous maximum obj
    VectorXd diel_old_max;

    string AMatrixMethod;
    SiCi* SiCiValue;


public:
    AProductCore(Space* space_, double d_, double lam_, Vector2cd material_, string AMatrixMethod_);
    AProductCore(Space* space_, double d_, double lam_, Vector2cd material_, int MAXm_, int MAXn_, double Lm_, double Ln_, string AMatrixMethod_);
    ~AProductCore();
    Matrix3cd A_dic_generator(double x, double y, double z);
    Matrix3cd A_dic_generator(double x, double y, double z, int m, int n);
    VectorXcd Aproduct(VectorXcd& b);                                          //without al*b because al is in DDAModel and can be diff for the same AMatrix
    void UpdateStr(VectorXd step);                                      //alpha not updated because in DDAModel, do not forget!
    void output_to_file();
    void output_to_file(string save_position, int iteration);
    
    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    tuple<list<int>, list<int>, list<int>, list<int>> get_para_info();
    VectorXi* get_R();
    double get_d();
    Space* get_space();
    double get_lam();
    VectorXcd* get_diel();
    list<list<int>>* get_PositionDep();
    VectorXi* get_PositionPara();          
    list<int>* get_para_nums();
    list<int>* get_para_starts();
    list<int>* get_para_dep_nums();
    list<int>* get_para_dep_starts();                 
    VectorXd* get_diel_old();               
    Vector2cd* get_material();
    VectorXcd* get_diel_max();                        
    VectorXd* get_diel_old_max();
    Matrix3cd FCD_inter(double x, double y, double z);
    Matrix3cd LDR_inter(double x,double y,double z);
};

class DDAModel {
private:
    //------------------------------------Get from AProductCore------------------------------ For the same AMatrix, these are all the same
    AProductCore* Core;

    //-----------------------------------Not from AProductCore------------------------------- For the same AMatrix, these can be diff for diff DDAModel
    bool RResultSwitch;               //0(false) for plot only E field on the structure points (fast), 1(true) for using RResult different from R to plot (slow but adjustable).
    VectorXi RResult;                //The position matrix for the EResult (where you want to plot the E field)
    double E0; 
    Vector3d n_E0;
    Vector3d n_K;
    VectorXcd P;
    VectorXcd P_store;                //Previous solution for P in the last iteration.
    VectorXcd E;
    VectorXcd Einternal;              //E field on structure points
    VectorXcd EResult;                //E field on designated points
    VectorXcd al;                       // 1 over alpha instead of alpha.
    bool verbose;
    VectorXcd P_max;
    VectorXcd al_max;

    //------------------different for different angles------------------
    int time;
    int ITERATION;
    double Error;

public:
    DDAModel(AProductCore* AProductCore_, Vector3d n_K_, double E0_, Vector3d n_E0_);
    DDAModel(AProductCore* AProductCore_, Vector3d n_K_, double E0_, Vector3d n_E0_, VectorXi* RResult_);
    int bicgstab(int MAX_ITERATION, double MAX_ERROR);
    void save_P();             //set P_store to P and set P to 0.
    void set_P();              //set P to P_store.
    void change_E(VectorXcd E_);
    void reset_E();             //reset E to E0                                
    void UpdateAlpha();                                //update alpha according to updated diel in AProductCore.
    void solve_E();                                                        //update the result E field on each dipole or on a designated space
    void update_E_in_structure();                                          //update the result E field on each dipole 
    VectorXcd Aproductwithalb(VectorXcd& b);                    //add the al*b term on base of AproductCore
    void output_to_file();
    void output_to_file(string save_position, int iteration, int ModelLabel);              //especially used for EvoOptimization
    VectorXcd* get_P();
    Vector3d get_nE0();
    Vector3d get_nK();
    double get_E0();
    VectorXcd* get_Einternal();
    AProductCore* get_Core();
    VectorXcd* get_al();
    VectorXcd* get_P_max();
    VectorXcd* get_al_max();
};

class ObjectiveDDAModel;

class EvoDDAModel {
private:
    AProductCore* Core;                           //The AProductCore
    list<DDAModel*> ModelList;                    //List of DDA models sharing the same AProductCore : "Core"
    int ModelNum;                                 //number of DDA model
    string save_position;

    list<list<double>*>* ObjectParameters;
    list<double> MajorObjectParameters;
    list<list<double>> MinorObjectParameters;

    list<string>* ObjectFunctionNames;
    string MajorObjectFunctionName;
    double MajorObjectFunctionResult;
    list<string> MinorObjectFunctionNames;
    list<double> MinorObjectFunctionResults;
    list<ObjectiveDDAModel*> ObjList;
    VectorXd Originarray;                               //Record the objective function for partial derivative (the value before change)   
    bool HavePenalty;
    double PenaltyFactor;

    double MaxObj;                                //Record the historical maximum obj func
    VectorXd MaxObjarray;                         //the individual objs for each model when the average obj is maximum(not necessaily the maximum individual objs)
    double epsilon_fix;
    double epsilon_tmp;                         //The epsilon used for calculation (can be different from the fixed input epsilon)
    bool HavePathRecord;
    int Stephold;

    int iteration;                              //current iteration.
    int MAX_ITERATION_EVO;                      //maximum number of iteration for optimizaiton.
public:
    EvoDDAModel(list<string>* ObjectFunctionNames_, list<list<double>*>* ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, double PenaltyFactor_, string save_position_, AProductCore* Core_, list<DDAModel*> ModelList_);
    
    //functions used to calculate partial derivatives                                 
    tuple<VectorXd, VectorXcd> devx_and_Adevxp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin);                       //partial derivative of obj to parameter and A to x times p
    VectorXcd devp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin);                       //partial derivative of obj to P. Size of P

    void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO_, string method);
    double CalTheObjForSingleStr(int MAX_ITERATION, double MAX_ERROR, int Name);                    //If you want to calculate the objective for single DDA structure.

    //The objective choosing function:
    ObjectiveDDAModel* ObjectiveFactory(string ObjectName, list<double> ObjectParameters, DDAModel* ObjDDAModel);

    double L1Norm();



};

//Abstract parent class for objective function.
class Objective {
public:
    bool Have_Devx;
    bool Have_Penalty;
    virtual void SingleResponse(int idx, bool deduction) = 0;
    virtual double GroupResponse() = 0;
    virtual void Reset() = 0;
    virtual double GetVal() = 0;
};


//Child classes for objective function


class ObjectivePointE : public Objective {
private:
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    EvoModel* model;
    Vector3cd E_sum;
    Vector3cd E_ext;
public:
    ObjectivePointE(list<double> parameters, EvoModel* model_, bool HavePenalty_);
    ObjectivePointE(list<double> parameters, DDAModel* model_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveSurfaceEExp : public Objective {
private:
    bool Have_Devx;
    bool Have_Penalty;
    int N;
    int Nobj;
    double z;
    double exponent;
    double d;
    int nz;
    VectorXcd* P;
    VectorXi* R;
    VectorXi Robj;
    VectorXcd* Einternal;
    EvoModel* model;
    VectorXcd E_sum;

public:
    ObjectiveSurfaceEExp(list<double> parameters, EvoModel* model_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveExtSurfaceEExp_CPU : public Objective {
private:
    bool Have_Devx;
    bool Have_Penalty;
    double d;
    int N;
    int Nobj;
    double exponent;                              //2, 4 or something else for E^?
    double ExtSurfaceEExpRz;
    int Nx;
    int Ny;
    int Nz;                                       //The entire(as big as focus in Nz, which is bigger then the str length)
    int ratio;                                 //Nx_obj = Nx/ratio; Ny_obj = Ny/ratio;
    vector<vector<vector<Matrix3cd>>> A_dic;

    VectorXcd* P;
    VectorXi* R;
    VectorXi Robj;
    EvoModel* model;
    VectorXcd E_sum;
    VectorXcd E_ext;
    double distance0;                                  //shortest distance between the plane and the str(corresponds to nz=0 in the A_dic (longest z corresponds to nz=max(nz)))

public:
    ObjectiveExtSurfaceEExp_CPU(list<double> parameters, EvoModel* model_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveExtSurfaceEMax : public Objective {
private:
    bool Have_Devx;
    bool Have_Penalty;
    double d;
    int N;
    int Nobj;
    double exponent;                              //2, 4 or something else for E^?
    double ExtSurfaceEExpRz;
    int Nx;
    int Ny;
    int Nz;                                       //The entire(as big as focus in Nz, which is bigger then the str length)
    int ratio;                                 //Nx_obj = Nx/ratio; Ny_obj = Ny/ratio;
    vector<vector<vector<Matrix3cd>>> A_dic;

    VectorXcd* P;
    VectorXi* R;
    VectorXi Robj;
    EvoModel* model;
    VectorXcd E_sum;
    VectorXcd E_ext;
    double distance0;                                  //shortest distance between the plane and the str(corresponds to nz=0 in the A_dic (longest z corresponds to nz=max(nz)))

public:
    ObjectiveExtSurfaceEMax(list<double> parameters, EvoModel* model_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveExtSurfaceEExp_CPU_Old : public Objective {
private:
    bool Have_Devx;
    bool Have_Penalty;
    double d;
    int N;
    int Nobj;
    double exponent;                              //2, 4 or something else for E^?
    double ExtSurfaceEExpRz;
    int Nx;
    int Ny;
    int Nz;                                       //The entire(as big as focus in Nz, which is bigger then the str length)
    int ratio;                                 //Nx_obj = Nx/ratio; Ny_obj = Ny/ratio;
    vector<vector<vector<Matrix3cd>>> A_dic;

    VectorXcd* P;
    VectorXi* R;
    VectorXi Robj;
    EvoModel* model;
    VectorXcd E_sum;
    VectorXcd E_ext;
    int distance0;                                  //shortest distance between the plane and the str(corresponds to nz=0 in the A_dic (longest z corresponds to nz=max(nz)))

public:
    ObjectiveExtSurfaceEExp_CPU_Old(list<double> parameters, EvoModel* model_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveExtSurfaceEExp : public Objective {
private:
    bool Have_Devx;
    bool Have_Penalty;
    double d;
    int N;
    int Nobj;
    int NxA;
    int NyA;
    int NzA;
    int NA;
    double exponent;                              //2, 4 or something else for E^?
    double ExtSurfaceEExpRz;
    int Nx;
    int Ny;
    int Nz;                                       //The entire(as big as focus in Nz, which is bigger then the str length)
    double* AHos;
    double* ADev;                              // double, but actually double*2 course real and imag are both stored in this double
    cufftDoubleComplex* A00, * A01, * A02, * A11, * A12, * A22; //The components in A for FFT, only in device
    double* PHos;                              //b in Aproduct
    double* PDev;                              // double, but actually double*2 course real and imag are both stored in this double
    cufftDoubleComplex* PxDev, * PyDev, * PzDev; //b components for FFT, only in device
    double* ESumHos;
    double* ESumDev;
    cufftDoubleComplex* ESumxDev, * ESumyDev, * ESumzDev;

    VectorXcd* P;
    VectorXi* R;
    VectorXi Robj;
    EvoModel* model;
    VectorXcd E_sum;
    VectorXcd E_ext;
    int distance0;                                  //shortest distance between the plane and the str(corresponds to nz=0 in the A_dic (longest z corresponds to nz=max(nz)))

public:
    ObjectiveExtSurfaceEExp(list<double> parameters, EvoModel* model_, bool HavePenalty_);
    ~ObjectiveExtSurfaceEExp();
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveG : public Objective {
private:
    bool Have_Devx;
    bool Have_Penalty;
    double d;
    double E0;                                    //input electric field
    int N;
    int Nobj;
    double exponent;                              //2, 4 or something else for E^?
    double ExtSurfaceEExpRz;
    int Nx;
    int Ny;
    int Nz;                                       //The entire(as big as focus in Nz, which is bigger then the str length)
    int ratio;                                 //Nx_obj = Nx/ratio; Ny_obj = Ny/ratio;
    vector<vector<vector<Matrix3cd>>> A_dic;

    VectorXcd* P;
    VectorXi* R;
    VectorXi Robj;
    EvoModel* model;
    VectorXcd E_sum;
    VectorXcd E_ext;
    double distance0;                                  //shortest distance between the plane and the str(corresponds to nz=0 in the A_dic (longest z corresponds to nz=max(nz)))

public:
    ObjectiveG(list<double> parameters, EvoModel* model_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};


class ObjectiveDDAModel {
public:
    bool Have_Devx;
    bool Have_Penalty;
    virtual void SingleResponse(int idx, bool deduction) = 0;
    virtual double GroupResponse() = 0;
    virtual void Reset() = 0;
    virtual double GetVal() = 0;
};


//Child classes for objective function


class ObjectivePointEDDAModel : public ObjectiveDDAModel {
private:
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vector3cd E_sum;
    Vector3cd E_ext;
public:
    ObjectivePointEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveIntegratedEDDAModel : public ObjectiveDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    VectorXcd* P;
    VectorXcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    VectorXcd E;
    double E_int;
public:
    ObjectiveIntegratedEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};


class ObjectiveMDipoleDDAModel : public ObjectiveDDAModel {
private:
    double d;
    int N;
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    VectorXcd* P;
    VectorXi* R;
    VectorXcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vector3cd Dipole_sum;
public:
    ObjectiveMDipoleDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveTDipoleDDAModel : public ObjectiveDDAModel {
private:
    double d;
    int N;
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    VectorXcd* P;
    VectorXi* R;
    VectorXcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vector3cd Dipole_sum;
public:
    ObjectiveTDipoleDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};


class ObjectiveEDipoleDDAModel : public ObjectiveDDAModel {
private:
    double d;
    int N;
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    VectorXcd* P;
    VectorXi* R;
    VectorXcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    Vector3cd Dipole_sum;
public:
    ObjectiveEDipoleDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveMDipole_TDipoleDDAModel : public ObjectiveDDAModel {
private:
    double d;
    int N;
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    VectorXcd* P;
    VectorXi* R;
    VectorXcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    double K = 0.0;
    Vector3cd MDipole_sum;
    Vector3cd TDipole_sum;
public:
    ObjectiveMDipole_TDipoleDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveEDipole_TDipoleDDAModel : public ObjectiveDDAModel {
private:
    double d;
    int N;
    double x;
    double y;
    double z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    VectorXcd* P;
    VectorXi* R;
    VectorXcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    double K = 0.0;
    Vector3cd EDipole_sum;
    Vector3cd TDipole_sum;
public:
    ObjectiveEDipole_TDipoleDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

















#endif // PREPROCESSING_H_INCLUDED