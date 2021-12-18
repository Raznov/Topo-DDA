#include "definition.h"

#include "definition.h"
#define PI 3.14159265


int main() {

    string save_position = "./1um1um1um-SiO2-phi0theta0-lam500-theta90/";       //output file
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



