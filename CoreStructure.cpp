#include "definition.h"

CoreStructure::CoreStructure(Space* space_, double d_) {

    space = space_;
    d = d_;

    cout << "(d=" << d << ") " << endl;

    tie(Nx, Ny, Nz, N) = (*space_).get_Ns();
    list<Structure>* ln = (*space_).get_ln();
    R = VectorXi::Zero(3 * N);
    RDep = VectorXi::Zero(3 * N);
    VectorXd diel_tmp = VectorXd::Zero(3 * N);

    //-----------------------------------------------------------------Input strs-------------------------------------------------------------
    list<Structure>::iterator it = (*ln).begin();

    int PositionParaN = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = ((*((*it).get_geometry())).size());
        if ((*it).get_para() == 1) {
            PositionParaN += int(round(n2 / 3));
        }
        it++;
    }
    PositionPara = VectorXi::Zero(PositionParaN);

    //Build R, RDep
    int n1 = 0;
    it = (*ln).begin();
    int PositionParaPos = 0;
    VectorXi RPara = VectorXi::Zero(3 * N);
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = ((*((*it).get_geometry())).size());
        if ((*it).get_para() == 1) {
            para_nums.push_back(n2);
            para_starts.push_back(n1);
            for (int i = int(round(n1 / 3)); i <= int(round(n1 / 3)) + int(round(n2 / 3)) - 1; i++) {
                PositionPara(PositionParaPos) = i;
                PositionParaPos += 1;
            }
        }
        if ((*it).get_para() == 2) {
            para_dep_nums.push_back(n2);
            para_dep_starts.push_back(n1);
        }
        for (int j = 0; j <= n2 - 1; j++) {
            R(n1 + j) = (*((*it).get_geometry()))(j);
            diel_tmp(n1 + j) = (*((*it).get_diel()))(j);
            RDep(n1 + j) = (*((*it).get_geometry_dep()))(j);
            RPara(n1 + j) = (*it).get_para();
        }
        it++;
        n1 = n1 + n2;
    }


    for (int i = 0; i <= PositionPara.size() - 1; i++) {
        list<int> tmpPos;
        int Parax = R(3 * PositionPara(i));
        int Paray = R(3 * PositionPara(i) + 1);
        int Paraz = R(3 * PositionPara(i) + 2);
        for (int j = 0; j <= N - 1; j++) {
            int jx = R(3 * j);
            int jy = R(3 * j + 1);
            int jz = R(3 * j + 2);
            int Depx = RDep(3 * j);
            int Depy = RDep(3 * j + 1);
            int Depz = RDep(3 * j + 2);
            if (RPara(3 * j) == 2) {
                if (Depx == Parax && Depy == Paray && Depz == Paraz) {
                    tmpPos.push_back(j);
                }
            }
        }
        PositionDep.push_back(tmpPos);
    }

    int PositionDepN = 0;
    list<list<int>>::iterator it1PositionDep = PositionDep.begin();
    for (int i = 0; i <= int(PositionDep.size()) - 1; i++) {
        PositionDepN += (*it1PositionDep).size();
        it1PositionDep++;
    }
    if (PositionDepN + PositionParaN != N) {
        cout << "PositionDepN = " << PositionDepN << endl;
        cout << "PositionParaN = " << PositionParaN << endl;
        cout << "In CoreStructure::CoreStructure(Space* space_, double d_) : PositionDepN + PositionParaN! = N" << endl;
    }

    //---------------------------------------------------initial diel------------------------------------
    diel_old = VectorXd::Zero(3 * N);
    diel_old_max = diel_old;
    for (int i = 0; i <= 3 * N - 1; i++) {
        diel_old(i) = diel_tmp(i);
    }
}

void CoreStructure::UpdateStr(VectorXd step) {

    cout << "step in UpdateStr" << step.mean() << endl;

    int para_size = para_nums.size();
    int para_dep_size = para_dep_nums.size();
    //cout << para_size << ' ' << para_dep_size << endl;
    if (para_dep_size != 0) {
        if (PositionPara.size() != PositionDep.size()) {
            cout << "In Model::change_para_diel(VectorXd step) : PositionPara.size() != PositionDep.size()" << endl;
        }

        list<list<int>>::iterator it1 = PositionDep.begin();
        for (int i = 0; i <= PositionPara.size() - 1; i++) {
            int position1 = PositionPara(i);
            diel_old(3 * position1) += step(i);
            if (diel_old(3 * position1) >= 1) {
                diel_old(3 * position1) = 1;
            }
            if (diel_old(3 * position1) <= 0) {
                diel_old(3 * position1) = 0;
            }

            diel_old(3 * position1 + 1) = diel_old(3 * position1);
            diel_old(3 * position1 + 2) = diel_old(3 * position1);

            list<int>::iterator it2 = (*it1).begin();
            for (int j = 0; j <= (*it1).size() - 1; j++) {
                int position2 = (*it2);
                diel_old(3 * position2) = diel_old(3 * position1);
                diel_old(3 * position2 + 1) = diel_old(3 * position1);
                diel_old(3 * position2 + 2) = diel_old(3 * position1);

                it2++;
            }
            it1++;
        }
    }
    else {

        list<int>::iterator it1 = para_nums.begin();
        list<int>::iterator it2 = para_starts.begin();
        int position = 0;
        for (int i = 0; i <= para_size - 1; i++) {
            int para_begin = round((*it2) / 3);
            int para_number = round((*it1) / 3);

            for (int j = 0; j <= para_number - 1; j++) {
                int position1 = (j + para_begin);
                diel_old(3 * position1) += step(position);
                if (diel_old(3 * position1) >= 1) {
                    diel_old(3 * position1) = 1;
                }
                if (diel_old(3 * position1) <= 0) {
                    diel_old(3 * position1) = 0;
                }

                diel_old(3 * position1 + 1) = diel_old(3 * position1);
                diel_old(3 * position1 + 2) = diel_old(3 * position1);

                position = position + 1;

            }
            it1++;
            it2++;
        }

    }

}

void CoreStructure::output_to_file() {

    ofstream fout("CoreStructure.txt");
    fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
    fout << R << endl;
    fout << diel_old << endl;
    fout << d << endl;
    fout.close();
}

void CoreStructure::output_to_file(string save_position, int iteration) {

    string name;
    //name=save_position+"AProductCoreit" + to_string(iteration) + ".txt";
    name = save_position + "CoreStructure" + to_string(iteration) + ".txt";
    ofstream fout(name);
    fout << Nx << endl << Ny << endl << Nz << endl << N << endl;
    fout << R << endl;
    fout << diel_old << endl;
    fout << d << endl;
    fout.close();
}

int CoreStructure::get_N() {
    return N;
}
int CoreStructure::get_Nx() {
    return Nx;
}
int CoreStructure::get_Ny() {
    return Ny;
}
int CoreStructure::get_Nz() {
    return Nz;
}
tuple<list<int>, list<int>, list<int>, list<int>> CoreStructure::get_para_info() {
    return make_tuple(para_nums, para_starts, para_dep_nums, para_dep_starts);
}
VectorXi* CoreStructure::get_R() {
    return &R;
}
double CoreStructure::get_d() {
    return d;
}
Space* CoreStructure::get_space() {
    return space;
}
list<list<int>>* CoreStructure::get_PositionDep() {
    return &PositionDep;
}
VectorXi* CoreStructure::get_PositionPara() {
    return &PositionPara;
}
list<int>* CoreStructure::get_para_nums() {
    return &para_nums;
}
list<int>* CoreStructure::get_para_starts() {
    return &para_starts;
}
list<int>* CoreStructure::get_para_dep_nums() {
    return &para_dep_nums;
}
list<int>* CoreStructure::get_para_dep_starts() {
    return &para_dep_starts;
}
VectorXd* CoreStructure::get_diel_old() {
    return &diel_old;
}
VectorXd* CoreStructure::get_diel_old_max() {
    return &diel_old_max;
}