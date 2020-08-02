#ifndef PREPROCESSING_H_INCLUDED
#define PREPROCESSING_H_INCLUDED
#define quote(x) #x
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

complex<double> Get_Alpha(double lam,double K,double d,complex<double> diel);
complex<double> Get_Alpha_FLTRCD(double lam,double K,double d,complex<double> diel);
complex<double> Get_material(string mat, double wl, string unit);                  //name of mat to get its diel function at certain wavlength              
Vector2cd Get_2_material(string sub, string mat, double wl, string unit);          //a wrapper for Get_material
double Average(VectorXcd* E, int N, double exponent);
//ArrayXcd FFT(int nx,int ny,int nz,ArrayXcd in,int _direction);
VectorXi build_a_bulk(int Nx, int Ny, int Nz);

class Structure{
    private:
        VectorXd diel;                  //0~1
        VectorXi geometry;
        int para;                       //Parameter condition: 0-not parameter, 1-parameter, 2-duplicated from a 2D parameter, controlled by parameter.
    public:
        Structure(VectorXi *total_space, VectorXd *diel_, VectorXi *geometry_, int para_);
        Structure(VectorXi *total_space, string initial_diel, double r, Vector3d center, int para_);  //r: actual radius/d. center: actual center/d. In charge of Sphere 
        Structure(VectorXi *total_space, string initial_diel, double r, Vector3i center, Vector3i direction, int para_); //build a circle, direction is its normalized direction in Cart. coord.
        Structure(VectorXi *total_space, string initial_diel, Vector3d l, Vector3d center, int para_);    //Ractangular(both 2D and 3D). 
        Structure(VectorXi *total_space, Structure *s, Vector3i direction, int times, int para_);                     //Initializa a Structure by duplicating a existing structure along a certain direction for several times. Direction is normalized and can only be alone x, y or z.
        Structure(VectorXi *total_space, string FileName, int para_);                                                    //Read a structure from txt file
        VectorXi *get_geometry();
        VectorXd *get_diel();
        void cut(VectorXi *big, VectorXi *small, VectorXd *small_diel);
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
        double ERROR;
        VectorXi R;                      //Position of dipoles. Both R and RResult are unitless, so need to time d to get real number.
        bool RResultSwitch;               //0(false) for plot only E field on the structure points (fast), 1(true) for using RResult different from R to plot (slow but adjustable).
        VectorXi RResult;                //The position matrix for the EResult (where you want to plot the E field)
        list<int> para_nums;
        list<int> para_starts;
        list<int> para_dep_nums;
        list<int> para_dep_starts;
        Vector3d n_E0;
        Vector3d n_K;
        VectorXcd diel;                   //real diel after 0~1 corresponds to real numbers
        VectorXcd diel_max;                         //corresponds to the previous maximum obj
        VectorXd diel_old;                //The 0~1 version of diel
        VectorXcd P;
        VectorXcd E;
        VectorXcd Einternal;              //E field on structure points
        VectorXcd EResult;                //E field on designated points
        Vector2cd material;
        VectorXcd al;                       // 1 over alpha instead of alpha.
        bool verbose;
    public:
        Model(Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_);
        Model(Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_, VectorXi *RResult_);
        ~Model();
        Matrix3cd A_dic_generator(double x,double y,double z);
        VectorXcd Aproduct(VectorXcd &b);
        void bicgstab(int MAX_ITERATION,double MAX_ERROR);
        void change_E(VectorXcd E_);
        void reset_E();             //reset E to E0 
        double get_step_length(VectorXd gradients, double epsilon);                                        
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
    public:
        EvoModel(list<string> *ObjectFunctionNames_, list<list<double>*> *ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, double PenaltyFactor_, string save_position_, Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_);
        EvoModel(list<string> *ObjectFunctionNames_, list<list<double>*> *ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, double PenaltyFactor_, string save_position_, Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_, VectorXi *RResult_);
        
        
        //functions used to calculate partial derivatives                                 
        tuple<VectorXd, VectorXcd> devx_and_Adevxp(double epsilon);                       //partial derivative of obj to parameter and A to x times p
        VectorXcd devp(double epsilon);                       //partial derivative of obj to P. Size of P
        
        void EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method);

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









//Abstract parent class for objective function.
class Objective{
    public:
        bool Have_Devx;
        bool Have_Penalty;
        virtual void SingleResponse(int idx, bool deduction) = 0;
        virtual double GroupResponse() = 0;
        virtual void Reset() = 0;
        virtual double GetVal() = 0;
};


//Child classes for objective function


class ObjectivePointE : public Objective{
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
      void SingleResponse(int idx, bool deduction);
      double GroupResponse();
      double GetVal();
      void Reset();
};

class ObjectiveSurfaceEExp : public Objective{
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

class ObjectiveExtSurfaceEExp_CPU : public Objective{
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
      ObjectiveExtSurfaceEExp_CPU(list<double> parameters, EvoModel* model_, bool HavePenalty_);
      void SingleResponse(int idx, bool deduction);
      double GroupResponse();
      double GetVal();
      void Reset();
};



class ObjectiveExtSurfaceEExp : public Objective{
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
        double *ADev;                              // double, but actually double*2 course real and imag are both stored in this double
        cufftDoubleComplex *A00, *A01, *A02, *A11, *A12, *A22; //The components in A for FFT, only in device
        double *PHos;                              //b in Aproduct
        double *PDev;                              // double, but actually double*2 course real and imag are both stored in this double
        cufftDoubleComplex *PxDev, *PyDev, *PzDev; //b components for FFT, only in device
        double *ESumHos;                              
        double *ESumDev;                              
        cufftDoubleComplex *ESumxDev, *ESumyDev, *ESumzDev; 

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
























#endif // PREPROCESSING_H_INCLUDED