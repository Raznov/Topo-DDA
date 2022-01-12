#include "definition.h"
#define PI 3.14159265



int main() {

    ofstream TotalTime;
    TotalTime.open("TotalTime.txt");
    high_resolution_clock::time_point t_start = high_resolution_clock::now();



    Vector3d l;
    Vector3d center;
    l << 21.0, 21.0, 9.0;
    //l << 40.0, 40.0, 8.0;
    center << l(0) / 2, l(1) / 2, l(2) / 2;

    int Nx, Ny, Nz;
    //Nx = 103; Ny = 103; Nz = 16;
    Nx = round(l(0) + 1); Ny = round(l(1) + 1); Nz = round(l(2) + 1);
    cout << center << endl;
    //Nx = 23; Ny = 23; Nz = 10;
    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    Vector3i direction;

    //l << 20.0, 20.0, 9.0;
    //center << 10.0, 10.0, 4.5;
    Structure s1(S.get_total_space(), l, center);



    S = S + s1;
    double d = 15;
    double r = 150 / d;




    Vector3i bind(1, 1, 10);
    VectorXi* s1geometry = s1.get_geometry();
    list<VectorXi*> FPGeometryl{ s1geometry };
    list<VectorXi*> BPGeometryl{};
    list<double> BParal{};
    bool Filter = true;
    double r_f = 2.0;
    FilterOption filteropt(0.0, 80.0, 0.5, "piecewise", r_f);
    SpacePara spacepara(&S, bind, "RANDOM", FPGeometryl, BPGeometryl, BParal, Filter, &filteropt);

    //SpacePara spacepara(&S, bind, "RANDOM");


    double E0 = 1.0;


    double epsilon = 0.2;
    //double epsilon = 1;

    //double focus = (l(2) + 2) * d;   //nm       
    double focus = (l(2) + 2) * d;
    cout << focus << endl;

    int MAX_ITERATION_DDA = 100000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 100;

    list<string> ObjectFunctionNames{ "IntegratedE" };

    double exponent = 2;
    double ratio = 4;

    list<double> ObjectParameter{ 70.0 };

    bool HavePathRecord = false;
    bool HavePenalty = false;
    bool HaveOriginHeritage = false;
    bool HaveAdjointHeritage = false;
    double PenaltyFactor = 1;
    list<list<double>*> ObjectParameters{ &ObjectParameter };
    string save_position = "./p330-lam542-TiO2-InE-2layer-rand-filter-beta80/";

    Vector3d n_K;
    Vector3d n_E0;



    list<DDAModel> ModelList;
    list<DDAModel*> ModelpointerList;

    ofstream AngleInfo(save_position + "AngleInfo.txt");
    ofstream nEInfo(save_position + "nEInfo.txt");

    int theta_num = 1;
    VectorXd theta(theta_num);
    theta << 0;
    int phi_num = 1;
    VectorXd phi(phi_num);
    phi << 0;
    int lam_num = 1;
    VectorXd lam(lam_num);
    lam << 542;

    CoreStructure CStr(&spacepara, d);
    list<AProductCore> CoreList;
    list<AProductCore*> CorePointList;
    Vector2cd material;
    material = Get_2_material("Air", "TiO2", lam(0), "nm");
    //AProductCore Core1(&CStr, lam(0), material, "LDR");

    int m, n;
    double Lm, Ln;
    m = 50;
    n = 50;
    Lm = 22 * d;
    Ln = 22 * d;
    AProductCore Core1(&CStr, lam(0), material, m, n, Lm, Ln, "FCD");

    CorePointList.push_back(&Core1);

    /*
    for (int k = 0; k <= lam_num - 1; k++) {
        Vector2cd material = Get_2_material("Air", "SiO2", lam(k), "nm");
        AProductCore Core_tmp(&CStr, lam(k), material);
        CoreList.push_back(Core_tmp);
    }
    */

    list<AProductCore*>::iterator it = CorePointList.begin();
    for (int k = 0; k <= lam_num - 1; k++) {
        AProductCore* Core = (*it);
        for (int i = 0; i <= theta_num - 1; i++) {
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


    cout << "nK" << n_K << endl;
    cout << "nE0" << n_E0 << endl;

    AngleInfo.close();
    nEInfo.close();
    cout << "Number of DDA Model : " << ModelList.size() << endl;


    ofstream Common;
    Common.open(save_position + "Commondata.txt");
    Common << CStr.get_Nx() << endl << CStr.get_Ny() << endl << CStr.get_Nz() << endl << CStr.get_N() << endl;
    Common << (spacepara.get_geometry()) << endl;
    Common << d << endl;
    Common << n_E0 << endl;
    Common << n_K << endl;


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