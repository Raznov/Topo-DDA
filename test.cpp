#include "definition.h"

int main() {

    int Nx, Ny, Nz;
    Nx = 56; Ny = 56; Nz = 56;

    int N = 0;
    VectorXi total_space = build_a_bulk(Nx, Ny, Nz);
    list<Structure> ln;
    Space S(&total_space, Nx, Ny, Nz, N, &ln);

    double d;
    double r;
    Vector3d center;

    d = 2;
    r = 55 / d;
    center << Nx / 2, Ny / 2, Nz / 2;

    Structure s1(S.get_total_space(), "ONES", r, center, 1);



    S = S + s1;

    double ratioESItoG = 1 / (2.998 * pow(10, 4));
    double rationmtocm = 1 / (pow(10, 7));
    double ratioESItoGnm = (1 / (2.998 * pow(10, 4))) / (pow(10, 7));



    double lam = 525;
    Vector3d n_K;
    n_K << 0.0, 0.0, 1.0;
    double E0 = 1.0;
    Vector3d n_E0;
    n_E0 << 1.0, 0.0, 0.0;
    Vector2cd material = Get_2_material("Air", "Au", lam, "nm");

    double epsilon = 100;

    double focus = 350;   //nm       


    int MAX_ITERATION_DDA = 10000;
    double MAX_ERROR = 0.00001;
    int MAX_ITERATION_EVO = 100;

    //d = rationmtocm * d;
    //lam = rationmtocm * lam;
    //E0 = ratioESItoGnm * E0;
    //focus = rationmtocm * focus;


    list<string> ObjectFunctionNames{ "PointEratio" };

    double exponent = 2;
    double ratio = 4;

    //list<double> ObjectParameter1{ focus, exponent, ratio };
    list<double> ObjectParameter2{ center(0) * d,center(1) * d,focus };

    bool HavePathRecord = true;
    bool HavePenalty = false;
    double PenaltyFactor = 0.0001;
    list<list<double>*> ObjectParameters{ &ObjectParameter2 };
    string save_position = "";

    AProductCore Core(&S, d, lam, material);
    DDAModel TestModel(&Core, n_K, E0, n_E0);

    TestModel.bicgstab(MAX_ITERATION_DDA, MAX_ERROR);
    TestModel.update_E_in_structure();
    TestModel.solve_E();
    TestModel.output_to_file();

    double x, y, z;
    y = center(1) * d;
    z = center(2) * d;


    N = Core.get_N();
    VectorXcd* P = TestModel.get_P();
    VectorXi* R = Core.get_R();
    double K = 2 * M_PI / lam;
    Vector3cd E_sum = Vector3cd::Zero();
    Vector3cd E_ext = Vector3cd::Zero();

    string name = "DDAModelEfieldscand=" + to_string(d) + ".txt";
    ofstream fout(name);
    for (x = center(0) * d; x <= center(0) * d + 350; x += 2) {
        cout << "x" << x - center(0) * d << endl;
        E_ext(0) = E0 * n_E0(0) * (cos(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_ext(1) = E0 * n_E0(1) * (cos(K * (n_K(0) * y + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_ext(2) = E0 * n_E0(2) * (cos(K * (n_K(0) * z + n_K(1) * y + n_K(2) * z)) + sin(K * (n_K(0) * x + n_K(1) * y + n_K(2) * z)) * 1i);
        E_sum(0) = E_ext(0);
        E_sum(1) = E_ext(1);
        E_sum(2) = E_ext(2);
        for (int i = 0; i <= N - 1; i++) {
            double rx = x - d * (*R)(3 * i);                  //R has no d in it, so needs to time d
            double ry = y - d * (*R)(3 * i + 1);
            double rz = z - d * (*R)(3 * i + 2);
            Matrix3cd A = Core.A_dic_generator(rx, ry, rz);
            E_sum(0) -= (A(0, 0) * (*P)(3 * i) + A(0, 1) * (*P)(3 * i + 1) + A(0, 2) * (*P)(3 * i + 2));
            E_sum(1) -= (A(1, 0) * (*P)(3 * i) + A(1, 1) * (*P)(3 * i + 1) + A(1, 2) * (*P)(3 * i + 2));
            E_sum(2) -= (A(2, 0) * (*P)(3 * i) + A(2, 1) * (*P)(3 * i + 1) + A(2, 2) * (*P)(3 * i + 2));
        }
        double normE = E_sum.norm();
        cout << "E" << normE << endl;
        fout << x - center(0) * d << " " << normE << endl;
    }
    fout.close();

    //EvoModel TestModel(&ObjectFunctionNames, &ObjectParameters, epsilon, HavePathRecord, HavePenalty, PenaltyFactor, save_position, &S, d, lam, n_K, E0, n_E0, material, AMatrixMethod);

    //TestModel.EvoOptimization(MAX_ITERATION_DDA, MAX_ERROR, MAX_ITERATION_EVO, "Adam");





    return 0;

}





