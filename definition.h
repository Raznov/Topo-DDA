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
#include <set>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>

//#include "cufftw.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "cufft.h"
//#include "direct.h"




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
complex<double> Get_Alpha_FCD(double lam, double K, double d, complex<double> diel);
complex<double> Get_material(string mat, double wl, string unit);                  //name of mat to get its diel function at certain wavlength              
Vector2cd Get_2_material(string sub, string mat, double wl, string unit);          //a wrapper for Get_material
VectorXcd Get_X_material(list<string> mat_l, double wl, string unit);
double Average(VectorXcd* E, int N, double exponent);
double Get_Max(VectorXcd* E, int N);
double G(VectorXcd* E, int N, double exponent, double E0);
//ArrayXcd FFT(int nx,int ny,int nz,ArrayXcd in,int _direction);
VectorXi build_a_bulk(int Nx, int Ny, int Nz);
bool CheckPerp(Vector3d v1, Vector3d v2);
Vector3d nEPerpinXZ(double theta, double phi);
MatrixXi find_scope_3_dim(VectorXi* x);
VectorXd initial_diel_func(string initial_diel, int N);
double initial_diel_func(string initial_diel);
list<double> makelist(double start, double end, double interval);
list<double> makelist(double start, double end, int number);

double exp_update(const double x, const double x_max, const double y_min, const double y_max);
double piecewise_update(const double x, const double x_max, const double y_min, const double y_max);
double linear_update(const double x, const double x_max, const double y_min, const double y_max);

int makedirect(string name);

class FilterOption {
    private:
        double beta;
        double ita;
        string beta_type;
        double beta_min;
        double beta_max;
        double rfilter;
        bool fixit;
        int MAX_ITERATION_FIXED;

    public:
        FilterOption(double beta_min_, double beta_max_, double ita_, string beta_type_, double rfilter_, bool fixit_=false, int MAX_ITERATION_FIXED_=100);
        double get_beta();
        double get_ita();
        double get_rfilter();
        void update_beta(const int iteration, const int Max_iteration);
        double SmoothDensity(double input);
};

class Structure{
    private:
        //VectorXd diel;                  //0~1
        VectorXi geometry;
        //VectorXi geometry_dep;          //Same size with geometry. If para=2, stores the corresponding para position of each point. Has no meaning in para=0&1
        //int para;                       //Parameter condition: 0-not parameter, 1-parameter, 2-duplicated from a 2D parameter, controlled by parameter.
    public:
        //--------------------------------------Dont support dependent para build-------------------------------------------------------
        Structure(VectorXi *total_space, VectorXi *geometry_);            //vanilla initialization

        //Sphere
        Structure(VectorXi *total_space, double r, Vector3d center);  //r: actual radius/d. center: actual center/d. In charge of Sphere 
        
        //Circle
        //Structure(VectorXi *total_space, string initial_diel, double r, Vector3i center, Vector3i direction, int para_); //build a circle, direction is its normalized direction in Cart. coord.

        //From file
        //Structure(VectorXi* total_space, string FileName, int para_);                                                    //Read a structure from txt file
        
        //--------------------------------------Support dependent para build-------------------------------------------------------------
        //Bulk
        Structure(VectorXi *total_space, Vector3d l, Vector3d center);    //Ractangular(both 2D and 3D). Para can only = 0&1 
        
        //Duplicate
        //Structure(VectorXi *total_space, Structure *s, Vector3i direction, int times, int para_);                     //Initializa a Structure by duplicating a existing structure along a certain direction for several times. Direction is normalized and can only be alone x, y or z.
        //The original structure is not included. original structure + new structure = times * original structure. If set para=2 and original para=1, then depend on origin str as geometry_dep. If para=2 and original para=2, will copy origin geometry_dep.
        //Structure(VectorXi* total_space, Structure* s, int dep_way);     //Special copy for setting up 2-fold and 4-fold symmetry dependence in xy plane. Para auto set to 2.
        //The original str at left down corner of the xy plane(left down corner is (0,0)). dep_way=1, 2, 3 corresponds to the other blocks in clock-wise. 


        //-------------------------------------Other member functions--------------------------------------------------------------------
        VectorXi *get_geometry();
        //VectorXi* get_geometry_dep();
        //VectorXd *get_diel();
        void cut(VectorXi *big, VectorXi *smalll);
        int get_geometry_size();
        //int get_para();
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

struct WeightPara {
    double weight;
    int position;
};

class SpacePara {
private:
    Space* space;
    VectorXi geometry;                //3N dimension
    VectorXi geometryPara;            //N dimension. N=number of dipoles. Each position stores the para index in VectorXi Para : 0->Para[0]...
    VectorXd Para;                    //P dimension. P=number of parameters. Same as Para_origin if Filter=False. Filtered and biased para if Filter=True.
    VectorXd Para_origin;             //Un-filtered, unbiased para. No use when Filter=False.
    VectorXd Para_filtered;           //Filtered but unbiased para. No use when Filter=False.
    MatrixXi scope;                   //[[xmin, xmax],[ymin, ymax],[zmin, zmax]]
    Vector3i bind;
    VectorXi FreeparatoPara;          //Position of free parameters inside Para. dimension<=P. FreeparatoPara[i] is the index of a free parameter inside Para.
    //vector<list<int>> Paratogeometry;  //P dimension. Each position stores a list of corresponding dipole index for parameter for this specific position.
    bool Filter;                      //True for with Filter. False for without Filter. Defualt as False for most initilizations.
    FilterOption* Filterstats;        //Only used when Filter=True
    vector<vector<WeightPara>> FreeWeight;

public:
    SpacePara(Space* space_, string initial_diel, VectorXi geometry_, VectorXd diel_); //Can freeze part of the parameter space
    SpacePara(Space* space_, Vector3i bind_, VectorXi* geometryPara_, VectorXd* Para_, VectorXi* FreeparatoPara_);
    SpacePara(Space* space_, Vector3i bind_, string initial_diel); //l, center similar to bulk build in Structure class. Every 'bind' nearby dipoles correspond 
                                                                    //to 1 parameter in this bulk. bind=(2,2,2): 2*2*2; bind=(1,1,3):1*1*3
    SpacePara(Space* space_, Vector3i bind_, string initial_diel1, string initial_diel2); //One constant layer at bottom. One design region on top. divide_pos is the divide in geometry array of the two parts.
    
    SpacePara(Space* space_, Vector3i bind_, string initial_diel, VectorXi* geometryPara_);
    SpacePara(Space* space_, Vector3i bind_, string initial_diel_center, string initial_diel_ring, double r, string type);   //ONly for 2d cylinder or spheres. r is raidus/d.

    SpacePara(Space* space_, Vector3i bind_, string initial_diel_background, list<string>* initial_diel_list, list<double>* r_list, list<Vector2d>* center_list);
    //Build 2d cylinders with diel in the list, rest of the diel is the backgroudn diel.

    SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2);  //random rect in a region with extruded 2D geometry
    SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2);  //random rect with 3D
    SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, VectorXi* geometryPara_);
    SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2, VectorXi* geometryPara_);

    SpacePara(Space* space_, Vector3i bind_, string initial_diel, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_);
    SpacePara(Space* space_, Vector3i bind_, string initial_diel, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_, bool Filter_, FilterOption* Filterstats_=NULL); //Should be exactly the same with the previous one except of Filter
    void ChangeBind(Vector3i bind_);                                  //Change bind number
    VectorXi cut(VectorXi* big, VectorXi* smalll);

    Space* get_space();
    VectorXi get_geometry();
    VectorXi* get_geometryPara();
    VectorXd* get_Para();
    VectorXd* get_Para_origin();
    VectorXd* get_Para_filtered();
    Vector3i* get_bind();
    VectorXi* get_Free();
    //vector<list<int>>* get_Paratogeometry();
    bool get_Filter();
    FilterOption* get_Filterstats();
    vector<vector<WeightPara>>* get_FreeWeight();
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

class CoreStructure {
private:
    //---------------------------------Geometries, not related to wavelength-------------------------------
    Space* space;
    SpacePara* spacepara;
    int N;                        //Number of dipoles
    int Nx;                       //scope of space. Nx*Ny*Nz!=N
    int Ny;
    int Nz;
    double d;
    VectorXi R;                      //Position of dipoles. Both R and RResult are unitless, so need to time d to get real number.

    //VectorXi RDep;                   //Position of the dependent para points in space of the points in the same position as R
    //list<list<int>> PositionDep;    //First D has the D of para(para=1). Second D is the positions of other points dependent on the para in the 1stD.
    //VectorXi PositionPara;          //Position i(in R) of the parameters (3*i=x, 3*i+1=y, 3*i+2=z)
    //list<int> para_nums;
    //list<int> para_starts;
    //list<int> para_dep_nums;
    //list<int> para_dep_starts;
    VectorXd diel_old;                //The 0~1 version of diel, 3*N
    VectorXd diel_old_max;
public:
    CoreStructure(SpacePara* spacepara_, double d_);
    void UpdateStr(VectorXd step, int current_it, int Max_it);
    void UpdateStr(SpacePara* spacepara_);
    void UpdateStrSingle(int idx, double value);
    void output_to_file();
    void output_to_file(string save_position, int iteration, string mode = "normal");

    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    VectorXi* get_R();
    double get_d();
    SpacePara* get_spacepara();
    //list<list<int>>* get_PositionDep();
    //VectorXi* get_PositionPara();
    //list<int>* get_para_nums();
    //list<int>* get_para_starts();
    //list<int>* get_para_dep_nums();
    //list<int>* get_para_dep_starts();
    VectorXd* get_diel_old();
    VectorXd* get_diel_old_max();

};

class AProductCore {
private:
    //---------------------------------Necessary values for A matrix-----------------------------------------------------
    CoreStructure* CStr;
    double K;
    double lam;
    
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
    
    //VectorXcd diel;                   //real diel after 0~1 corresponds to real numbers
    VectorXcd material;
    //VectorXcd diel_max;                         //corresponds to the previous maximum obj
    
    //------------------------------For FCD and LDR choice-----------------------
    string AMatrixMethod;
    SiCi* SiCiValue;

public:
    AProductCore(CoreStructure* CStr_, double lam_, VectorXcd material_, string AMatrixMethod_);
    AProductCore(CoreStructure* CStr_, double lam_, VectorXcd material_, int MAXm_, int MAXn_, double Lm_, double Ln_, string AMatrixMethod_);
    ~AProductCore();
    Matrix3cd A_dic_generator(double x, double y, double z);
    Matrix3cd A_dic_generator(double x, double y, double z, int m, int n);
    VectorXcd Aproduct(VectorXcd& b);                                          //without al*b because al is in DDAModel and can be diff for the same AMatrix
    //void UpdateStr(VectorXd step);                                      //alpha not updated because in DDAModel, do not forget!
    //void output_to_file();
    //void output_to_file(string save_position, int iteration);
    
    CoreStructure* get_CStr();
    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    double get_Lm();
    double get_Ln();
    tuple<list<int>, list<int>, list<int>, list<int>> get_para_info();
    VectorXi* get_R();
    double get_d();
    double get_lam();
    //VectorXcd* get_diel();        
    VectorXd* get_diel_old();               
    VectorXcd* get_material();
    //VectorXcd* get_diel_max();                        
    VectorXd* get_diel_old_max();
    Matrix3cd FCD_inter(double x, double y, double z);
    Matrix3cd LDR_inter(double x, double y, double z);
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
    VectorXcd E;
    VectorXcd Einternal;              //E field on structure points
    VectorXcd EResult;                //E field on designated points
    VectorXcd al;                       // 1 over alpha instead of alpha.
    VectorXcd diel;                     //Real dielectric from diel_old. Needed to calculate the Lorentz factor.
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
    void bicgstab(int MAX_ITERATION, double MAX_ERROR);
    void bicgstab(int MAX_ITERATION, double MAX_ERROR, int EVOITERATION);  //FOR DEBUG ONLY. OUTPUT SOME VALUE AT CERTAIN EVO ITERATION.
    void change_E(VectorXcd E_);
    void reset_E();             //reset E to E0                                
    void UpdateAlpha();                                //update alpha according to updated diel in AProductCore.
    void UpdateAlphaSingle(int idx);
    void solve_E();                                                        //update the result E field on each dipole or on a designated space
    void update_E_in_structure();                                          //update the result E field on each dipole 
    VectorXcd Aproductwithalb(VectorXcd& b);                    //add the al*b term on base of AproductCore
    void output_to_file();
    void output_to_file(string save_position, int iteration, int ModelLabel);              //especially used for EvoOptimization
    void output_to_file(string save_position, int iteration);             //For simplify output
    void output_to_file(string save_position, double wavelength, int iteration);
    void InitializeP(VectorXcd& Initializer);
    VectorXcd* get_P();
    Vector3d get_nE0();
    Vector3d get_nK();
    double get_E0();
    VectorXcd* get_Einternal();
    AProductCore* get_Core();
    VectorXcd* get_al();
    VectorXcd* get_P_max();
    VectorXcd* get_al_max();
    int get_ITERATION();

    //-----------------From AProductCore-----------------------

    int get_N();
    int get_Nx();
    int get_Ny();
    int get_Nz();
    VectorXi* get_R();
    double get_d();
    SpacePara* get_spacepara();
    double get_lam();
    VectorXd* get_diel_old();
    VectorXcd* get_material();
    VectorXd* get_diel_old_max();
};

class ObjectiveDDAModel;


class EvoDDAModel {
private:
    double output_time;
    CoreStructure* CStr;
    list<DDAModel*> ModelList;                    //List of DDA models sharing the same AProductCore : "Core"
    int ModelNum;                                 //number of DDA model
    string save_position;

    list<VectorXcd> PforOrigin;
    list<VectorXcd> PforAdjoint;
    list<VectorXcd> PforOriginMax;
    list<VectorXcd> PforAdjointMax;

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
    double PreviousObj;                            //The previous obj
    int CutoffHold;
    VectorXd MaxObjarray;                         //the individual objs for each model when the average obj is maximum(not necessaily the maximum individual objs)
    double epsilon_fix;
    double epsilon_tmp;                         //The epsilon used for calculation (can be different from the fixed input epsilon)
    bool HavePathRecord;
    bool HaveOriginHeritage;
    bool HaveAdjointHeritage;
    int Stephold;

    VectorXd gradientsquare;                    //cumulative summation of gradients square. Used in Adagrad.
public:
    EvoDDAModel(list<string>* ObjectFunctionNames_, list<list<double>*>* ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, bool HaveOriginHeritage_, bool HaveAdjointHeritage_, double PenaltyFactor_, string save_position_, CoreStructure* CStr_, list<DDAModel*> ModelList_);
    
    //functions used to calculate partial derivatives                                 
    tuple<VectorXd, VectorXcd> devx_and_Adevxp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin);                       //partial derivative of obj to parameter and A to x times p
    tuple<VectorXd, VectorXcd> devx_and_Adevxp_tmp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin);
    VectorXcd devp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin);                       //partial derivative of obj to P. Size of P

    void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, double start_num=0);
    void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, VectorXd* V_, VectorXd* S_);
    double CalTheObjForSingleStr(int MAX_ITERATION, double MAX_ERROR, int Name);                    //If you want to calculate the objective for single DDA structure.

    //The objective choosing function:
    ObjectiveDDAModel* ObjectiveFactory(string ObjectName, list<double> ObjectParameters, DDAModel* ObjDDAModel);

    double get_output_time();
    double L1Norm();
    VectorXd gradients_filtered(VectorXd gradients, int current_it, int Max_it);



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

class ObjectivePointListEDDAModel : public ObjectiveDDAModel {
private:
    VectorXd x;
    VectorXd y;
    VectorXd z;      // Here x, y, z are absolute coordinates. (No need to multiply d).
    int PNum;
    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    MatrixXcd E_sum;
    MatrixXcd E_ext;
public:
    ObjectivePointListEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectivePointIDDAModel : public ObjectiveDDAModel {
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
    ObjectivePointIDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
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
    VectorXi* R;
public:
    ObjectiveIntegratedEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class ObjectiveMidAvgEDDAModel : public ObjectiveDDAModel {
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
    double E_avg;
    double r;                 //radius of the middle regions (for 2D)
    VectorXi* R;
    double centerx;
    double centery;
public:
    ObjectiveMidAvgEDDAModel(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class Objectivescattering2D : public ObjectiveDDAModel {
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

    Vector3d n_K_s;
    Matrix3d FconstM;

public:
    Objectivescattering2D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
    double FTUCnsquare();
};

class Objectivescattering0D : public ObjectiveDDAModel {
private:

    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    EvoDDAModel* evomodel;
    double K;
    double E0;
    int Paralength;
    list<Vector3d> n_K_s_l;
    list<Vector3cd> PSum_l;
    list<Matrix3d> FconstM_l;
    //double ATUC;
    //Matrix3d FconstM;


public:
    Objectivescattering0D(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
    double FTUCnsquare();
};

class ObjectiveAbs : public ObjectiveDDAModel {
private:
    double d;
    int N;
    int Nx;
    int Ny;
    int Nz;
    double K;
    double K3;
    double E0;
    VectorXcd* P;
    VectorXcd* al;
    DDAModel* model;
    EvoDDAModel* evomodel;
    VectorXcd E;
    double Cabs;
    VectorXi* R;
public:
    ObjectiveAbs(list<double> parameters, DDAModel* model_, EvoDDAModel* evomodel_, bool HavePenalty_);
    void SingleResponse(int idx, bool deduction);
    double GroupResponse();
    double GetVal();
    void Reset();
};

class FOMscattering2D{
private:
    VectorXd FOMParameters;
    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    double K;
    int Paralength;
    list<Vector3d> n_K_s_l;
    double ATUC;
    double E0;
    Vector3d n_E0;
    //Matrix3d FconstM;
    

public:
    FOMscattering2D(list<double> parameters, DDAModel* model_);
    list<double> GetVal();                                                        //Return list of far field abs(Es) at specified directions
    Vector3cd FTUC(Vector3d n_K_s);
};

class FOMreflect2D {
private:
    VectorXd FOMParameters;
    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    double K;
    int Paralength;
    list<Vector3d> n_K_s_l;
    double ATUC;
    double E0;
    Vector3d n_E0;
    //Matrix3d FconstM;


public:
    FOMreflect2D(list<double> parameters, DDAModel* model_);
    list<double> GetVal();                                                        //Return list of far field abs(Es) at specified directions
    Vector3cd FTUC(Vector3d n_K_s);
};

class FOMscattering0D {
private:

    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    double K;
    double E0;
    int Paralength;
    list<Vector3d> n_K_s_l;
    //double ATUC;
    //Matrix3d FconstM;


public:
    FOMscattering0D(list<double> parameters, DDAModel* model_);
    list<double> GetVal();                                                      //Return list of dCsca/dOmega at specified directions
    double FTUCnsquare(Vector3d n_K_s);
};

/*
class FOMAbs {
private:

    double d;
    int N;
    VectorXcd* P;
    VectorXi* R;
    DDAModel* model;
    double K;
    double E0;
    int Paralength;
    list<Vector3d> n_K_s_l;
    //double ATUC;
    //Matrix3d FconstM;


public:
    FOMAbs(list<double> parameters, DDAModel* model_);
    list<double> GetVal();                                                      //Return list of dCsca/dOmega at specified directions
    double FTUCnsquare(Vector3d n_K_s);
};
*/










void Evo_Focus(SpacePara* spacepara_tmp, CoreStructure* CStr, DDAModel* TestModel, string save_position, int start_num, int max_evo,
    int min_num, int max_num, Vector3d lower_bound, Vector3d upper_bound, bool sym  //Parameters for focus generation
);
void Evo_single(string save_position, Vector3i bind, Vector3d l, int MAX_ITERATION_EVO, Vector3d move_focus);




void eval_FOM(string name, DDAModel* TestModel, list<double> theta, list<double> phi);
void eval_FOM_2Dperiod(string name, DDAModel* TestModel, list<double> theta, list<double> phi);

#endif // PREPROCESSING_H_INCLUDED