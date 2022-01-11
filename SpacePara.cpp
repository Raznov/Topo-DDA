#include "definition.h"

double calweight(int xo, int yo, int x, int y, double r) {
    return r - sqrt(double(x - xo) * double(x - xo) + double(y - yo) * double(y - yo));
}

bool circlerange(int xo, int yo, int x, int y, double r) {
    if (double(x - xo) * double(x - xo) + double(y - yo) * double(y - yo) <= r * r - 0.01) {
        return true;
    }
    else {
        return false;
    }
}

int Get3divSize(VectorXi* geometry) {
    int N = (*geometry).size();
    if (N % 3 != 0) {
        cout << "int Get3divSize(VectorXi geometry):Geometry size must be time of 3" << endl;
        throw N;
    }
    return int(round(N / 3));
}

set<vector<int>> Get3divSet(VectorXi* geometry) {
    set<vector<int>> result;
    int N = Get3divSize(geometry);
    for (int i = 0; i <= N - 1; i++) {
        vector<int> tmp{ (*geometry)(3 * i), (*geometry)(3 * i + 1), (*geometry)(3 * i + 2) };
        if (result.count(tmp)) {
            cout << "set<vector<int>> Get3divSet(VectorXi* geometry): Element is present in the set" << endl;
            throw i;
        }
        result.insert(tmp);
    }
    return result;
}

list<set<vector<int>>> Get3divSetList(list<VectorXi*> geometry) {
    list<set<vector<int>>> result;
    list<VectorXi*>::iterator it = geometry.begin();
    while (it != geometry.end()) {
        set<vector<int>> tmp = Get3divSet(*it);
        result.push_back(tmp);
        it++;
    }
    return result;
}

bool CheckOverlap(VectorXi* geometry1, VectorXi* geometry2) {
    int n1 = Get3divSize(geometry1);
    int n2 = Get3divSize(geometry2);
    for (int i = 0; i <= n1 - 1; i++) {
        int x1 = (*geometry1)(3 * i);
        int y1 = (*geometry1)(3 * i + 1);
        int z1 = (*geometry1)(3 * i + 2);
        for (int j = 0; j <= n2 - 1; j++) {
            int x2 = (*geometry2)(3 * j);
            int y2 = (*geometry2)(3 * j + 1);
            int z2 = (*geometry2)(3 * j + 2);
            if (x1 == x2 && y1 == y2 && z1 == z2) {
                return false;
            }
        }
    }

    return true;
}

bool CheckOverlapList(list<VectorXi*> geometry) {
    list<VectorXi*>::iterator it = geometry.begin();
    int i = 0;
    while (next(it) != geometry.end()) {
        list<VectorXi*>::iterator it1 = next(it);
        int j = i + 1;
        while (it1 != geometry.end()) {
            if (!CheckOverlap(*it, *it1)) {
                cout << "geometry " << i << " and geometry " << j << " overlap" << endl;
                return false;
            }
            it1++;
            j++;
        }
        it++;
        i++;
    }
    return true;
}

VectorXi ConnectGeometry(list<VectorXi*> geometry) {
    if (!CheckOverlapList(geometry)) {
        cout << "CheckOverlapList Fail" << endl;
        throw false;
    }

    list<VectorXi*>::iterator it = geometry.begin();
    int N = 0;
    while (it != geometry.end()) {
        N += Get3divSize(*it);
        it++;
    }
    
    VectorXi result = VectorXi::Zero(3 * N);
    it = geometry.begin();
    int i = 0;
    while (it != geometry.end()) {
        for (int j = 0; j <= Get3divSize(*it) - 1; j++) {
            result(3 * i) = (*(*it))(3 * j);
            result(3 * i + 1) = (*(*it))(3 * j + 1);
            result(3 * i + 2) = (*(*it))(3 * j + 2);
            i++;
        }
        it++;
    }
    return result;

}

VectorXi SpacePara::cut(VectorXi* big, VectorXi* smalll) {

    int number_origin = round((*smalll).size() / 3);
    MatrixXi big_scope = find_scope_3_dim(big);
    //cout<<"big_scope "<<big_scope<<endl;
    list<int> positions_in;
    int number_out = 0;
    //cout<<"small_scope "<<find_scope_3_dim(small)<<endl;
    for (int i = 0; i <= number_origin - 1; i++) {
        if (((*smalll)(3 * i) < big_scope(0, 0)) || ((*smalll)(3 * i) > big_scope(0, 1)) ||
            ((*smalll)(3 * i + 1) < big_scope(1, 0)) || ((*smalll)(3 * i + 1) > big_scope(1, 1)) ||
            ((*smalll)(3 * i + 2) < big_scope(2, 0)) || ((*smalll)(3 * i + 2) > big_scope(2, 1))) {
            number_out += 1;
        }
        else {
            positions_in.push_back(i);
        }
    }
    int number_in = positions_in.size();
    VectorXi geometry = VectorXi::Zero(3 * number_in);
    for (int i = 0; i <= number_in - 1; i++) {
        int j = positions_in.front();
        positions_in.pop_front();
        geometry(3 * i) = (*smalll)(3 * j);
        geometry(3 * i + 1) = (*smalll)(3 * j + 1);
        geometry(3 * i + 2) = (*smalll)(3 * j + 2);
    }
    if (number_out == 0) {
        cout << "The geometry you built is entirely in the space." << endl;
        cout << "number_origin " << number_origin << endl;
        cout << "number_real " << number_in << endl;
    }
    else {
        cout << "The geometry you built is at least partially outside of space." << endl;
        cout << "number_origin " << number_origin << endl;
        cout << "number_out " << number_out << endl;
        cout << "number_real " << number_in << endl;
    }

    return geometry;
}

SpacePara::SpacePara(Space* space_, string initial_diel, VectorXi geometry_, VectorXd diel_) {
    Filter = false;
    space = space_;
    bind << 1, 1, 1;                                    //---------Only work for 1,1,1
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);
    cout << N << endl;
    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    MatrixXi scopefrozen = find_scope_3_dim(&geometry_);            //Only work for rect input with same height
    int Nxfrozen, Nyfrozen, Nzfrozen, Nfrozen;
    Nxfrozen = ceil(double(scopefrozen(0, 1) - scopefrozen(0, 0) + 1) / bind(0));
    Nyfrozen = ceil(double(scopefrozen(1, 1) - scopefrozen(1, 0) + 1) / bind(1));
    Nzfrozen = ceil(double(scopefrozen(2, 1) - scopefrozen(2, 0) + 1) / bind(2));
    Nfrozen = Nxfrozen * Nyfrozen * Nzfrozen;

    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    FreeparatoPara = VectorXi::Zero(Npara - Nfrozen);

    Vector3i relativepos;
    relativepos << int((Nparax - Nxfrozen) / 2), int((Nparay - Nyfrozen) / 2), int((Nparaz - Nzfrozen) / 2);  //Does not neceesarily /2=int, can have a bit deviation

    cout << "Nparax" << Nparax << endl;
    cout << "Nxfrozen" << Nxfrozen << endl;
    cout << "relativepos" << relativepos << endl;

    Para = initial_diel_func(initial_diel, Npara);
    
    int nfree = 0;
    for (int i = 0; i <= Nparax - 1; i++) {
        for (int j = 0; j <= Nparay - 1; j++) {
            for (int k = 0; k <= Nparaz - 1; k++) {
                int pos = k + Nparaz * (j + Nparay * i);
                
                if (relativepos(0) <= i && i <= Nparax - 1 - relativepos(0) 
                    && relativepos(1) <= j && j <= Nparay - 1 - relativepos(1) 
                    && relativepos(2) <= k && k <= Nparaz - 1 - relativepos(2)) 
                {
                    int pos_tmp = (k - relativepos(2)) + Nzfrozen * ((j - relativepos(1)) + Nyfrozen * (i - relativepos(0)));
                    //cout << i << endl;
                    //cout << j << endl;
                    //cout << k << endl;
                    //cout << pos_tmp << endl;
                    Para(pos) = (diel_)(3 * pos_tmp);
                }
                else {
                    FreeparatoPara(nfree)=pos;
                    nfree += 1;
                }
            }
        }
    }

    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }

    /*
    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
    */
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, VectorXi* geometryPara_, VectorXd* Para_, VectorXi* FreeparatoPara_) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = *geometryPara_;
    Para = *Para_;
    
    FreeparatoPara = *FreeparatoPara_;

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1+j)= (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    //cout << "scope" << endl << scope << endl;
    //cout << "(Npara, Nparay, Nparaz) " << "(" << Nparax << ", " << Nparay << ", " << Nparaz << ")" << endl;
    //cout << "Npara: " << Npara << endl;
    Para = initial_diel_func(initial_diel, Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }
    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos= paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel1, string initial_diel2) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);
    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= 0; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        cout << n2 << endl;
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }
    //cout << n1 << endl;
    it = (*ln).begin();
    geometryPara = VectorXi::Zero(N);
    MatrixXi scopefree;
    scopefree = find_scope_3_dim((*it).get_geometry());
    int Nparafreex, Nparafreey, Nparafreez, Nparafree;
    Nparafreex = ceil(double(scopefree(0, 1) - scopefree(0, 0) + 1) / bind(0));
    Nparafreey = ceil(double(scopefree(1, 1) - scopefree(1, 0) + 1) / bind(1));
    Nparafreez = ceil(double(scopefree(2, 1) - scopefree(2, 0) + 1) / bind(2));
    Nparafree = Nparafreex * Nparafreey * Nparafreez;
    VectorXd para_free = initial_diel_func(initial_diel1, Nparafree);

    FreeparatoPara = VectorXi::Zero(Nparafree);
    for (int i = 0; i <= Nparafree - 1; i++) {
        FreeparatoPara(i) = i;
    }
    //cout << scopefree << endl;
    for (int i = 0; i <= int(round(n1/3)) - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scopefree(0, 0)) / bind(0));
        int paray = floor((y - scopefree(1, 0)) / bind(1));
        int paraz = floor((z - scopefree(2, 0)) / bind(2));
        int pos = paraz + Nparafreez * (paray + Nparafreey * parax);
        geometryPara(i) = pos;
    }

    int n1free = n1;                                //n1free does not equal to Nparafree
    //n1 = 0;
    it = (*ln).begin();
    it++;
    for (int i = 1; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        //cout << n2 << endl;
        //cout << ((*it).get_geometry_size()) << endl;
        //cout << n1 << endl;
        //cout << 3 * N << endl;
        for (int j = 0; j <= n2 - 1; j++) {
            //cout << j << endl;
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }
    int Npararest = int(round((n1 - n1free) / 3));
    VectorXd para_rest = initial_diel_func(initial_diel2, Npararest);
    Para = VectorXd::Zero(Nparafree + Npararest);
    for (int i = 0; i <= Nparafree - 1; i++) {
        Para(i) = para_free(i);
    }
    for (int i = Nparafree; i <= Nparafree + Npararest - 1; i++) {
        Para(i) = para_rest(i - Nparafree);
    }

    for (int i = int(round(n1free / 3)); i <= int(round(n1 / 3)) - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int pos = (i - int(round(n1free / 3))) + Nparafree;
        geometryPara(i) = pos;
    }




}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel, VectorXi* geometryPara_) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = *geometryPara_;
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    //cout << "scope" << endl << scope << endl;
    //cout << "(Npara, Nparay, Nparaz) " << "(" << Nparax << ", " << Nparay << ", " << Nparaz << ")" << endl;
    //cout << "Npara: " << Npara << endl;
    Para = initial_diel_func(initial_diel, Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel_center, string initial_diel_ring, double r, string type) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;

    double xcenter, ycenter, zcenter;
    xcenter = double(scope(0, 1) + scope(0, 0)) / 2;
    ycenter = double(scope(1, 1) + scope(1, 0)) / 2;
    zcenter = double(scope(2, 1) + scope(2, 0)) / 2;
    
    cout << "xcenter" << xcenter << endl;
    cout << "ycenter" << ycenter << endl;
    cout << "zcenter" << zcenter << endl;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Nparax - 1; i++) {
        for (int j = 0; j <= Nparay - 1; j++) {
            for (int k = 0; k <= Nparaz - 1; k++) {
                int pos = k + Nparaz * (j + Nparay * i);
                if (type == "CYLINDER") {
                    double x, y;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;
                    double rx, ry;
                    rx = x - xcenter;
                    ry = y - ycenter;
                    if (rx * rx + ry * ry < (r+0.1) * (r+0.1)) {
                        Para(pos) = initial_diel_func(initial_diel_center);
                    }
                    else {
                        Para(pos) = initial_diel_func(initial_diel_ring);
                    }
                }

                if (type == "SPHERE") {
                    double x, y, z;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;
                    z = bind(2) * (2 * k + 1) / 2;
                    double rx, ry, rz;
                    rx = x - xcenter;
                    ry = y - ycenter;
                    rz = z - zcenter;
                    if (rx * rx + ry * ry + rz * rz < (r + 0.1) * (r + 0.1)) {
                        Para(pos) = initial_diel_func(initial_diel_center);
                    }
                    else {
                        Para(pos) = initial_diel_func(initial_diel_ring);
                    }
                }

                

            }
        }
    }

    

    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel_background, list<string>* initial_diel_list, list<double>* r_list, list<Vector2d>* center_list){
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func(initial_diel_background);
    }

    list<string>::iterator it_diel = (*initial_diel_list).begin();
    list<double>::iterator it_r = (*r_list).begin();
    list<Vector2d>::iterator it_center = (*center_list).begin();

    for (int itnum = 0; itnum <= (*initial_diel_list).size() - 1; itnum++) {
        string initial_diel = (*it_diel);
        double r = (*it_r);
        Vector2d center = (*it_center);
        double xcenter, ycenter;
        xcenter = center(0);
        ycenter = center(1);
        cout << "xcenter" << xcenter << endl;
        cout << "ycenter" << ycenter << endl;
 

        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;
                    double rx, ry;
                    rx = x - xcenter;
                    ry = y - ycenter;
                    if (rx * rx + ry * ry < (r+0.1) * (r+0.1)) {
                        Para(pos) = initial_diel_func(initial_diel);
                    }
                }
            }
        }

        it_diel++;
        it_r++;
        it_center++;


    }
    

    
    
    


    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func("ZEROS");
    }

    int xcenter, ycenter;
    int lx, ly;


    for (int m = 0; m < number; m++) {
        xcenter = round(((double)rand() / RAND_MAX) * (double(scope(0, 1) - scope(0, 0) + 1))) + scope(0, 0);
        ycenter = round(((double)rand() / RAND_MAX) * (double(scope(1, 1) - scope(1, 0) + 1))) + scope(1, 0);
        lx = round(((double)rand() / RAND_MAX) * (limitx2 - limitx1) + limitx1);
        ly = round(((double)rand() / RAND_MAX) * (limity2 - limity1) + limity1);

        cout << "xcenter=" << xcenter << " ycenter=" << ycenter << endl;
        cout << "lx=" << lx << " " << "ly=" << ly << endl;
        


        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;

                    if ((abs(x - xcenter) <= lx) && (abs(y - ycenter) <= ly)) {
                        Para(pos) = initial_diel_func("ONES");
                    }
                }
            }
        }        
    }

    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }
  
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = VectorXi::Zero(N);
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func("ZEROS");
    }

    int xcenter, ycenter, zcenter;
    int lx, ly, lz;


    for (int m = 0; m < number; m++) {
        xcenter = round(((double)rand() / RAND_MAX) * (double(scope(0, 1) - scope(0, 0) + 1))) + scope(0, 0);
        ycenter = round(((double)rand() / RAND_MAX) * (double(scope(1, 1) - scope(1, 0) + 1))) + scope(1, 0);
        zcenter = round(((double)rand() / RAND_MAX) * (double(scope(2, 1) - scope(2, 0) + 1))) + scope(2, 0);
        lx = round(((double)rand() / RAND_MAX) * (limitx2 - limitx1) + limitx1);
        ly = round(((double)rand() / RAND_MAX) * (limity2 - limity1) + limity1);
        lz = round(((double)rand() / RAND_MAX) * (limitz2 - limitz1) + limitz1);


        cout << "xcenter=" << xcenter << " ycenter=" << ycenter << " zcenter=" << zcenter << endl;
        cout << "lx=" << lx << " " << "ly=" << ly << " " << "lz=" << lz << endl;



        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y, z;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;
                    z = bind(2) * (2 * k + 1) / 2;

                    if ((abs(x - xcenter) <= lx) && (abs(y - ycenter) <= ly) && (abs(z - zcenter) <= lz)) {
                        Para(pos) = initial_diel_func("ONES");
                    }
                }
            }
        }
    }

    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, VectorXi* geometryPara_) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = *geometryPara_;                //CAN IMPROVE
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func("ZEROS");
    }

    int xcenter, ycenter;
    int lx, ly;

    for (int m = 0; m < number; m++) {
        xcenter = round(((double)rand() / RAND_MAX) * (double(scope(0, 1) - scope(0, 0) + 1))) + scope(0, 0);
        ycenter = round(((double)rand() / RAND_MAX) * (double(scope(1, 1) - scope(1, 0) + 1))) + scope(1, 0);
        lx = round(((double)rand() / RAND_MAX) * (limitx2 - limitx1) + limitx1);
        ly = round(((double)rand() / RAND_MAX) * (limity2 - limity1) + limity1);

        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;

                    if ((abs(x - xcenter) <= lx) && (abs(y - ycenter) <= ly)) {
                        Para(pos) = initial_diel_func("ONES");
                    }
                }
            }
        }
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2, VectorXi* geometryPara_) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= (*ln).size() - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }

    scope = find_scope_3_dim(&geometry);
    geometryPara = *geometryPara_;
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil(double(scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil(double(scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil(double(scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;
    Para = VectorXd::Zero(Npara);
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    for (int i = 0; i <= Npara - 1; i++) {                          //First inital all to background
        Para(i) = initial_diel_func("ZEROS");
    }

    int xcenter, ycenter, zcenter;
    int lx, ly, lz;


    for (int m = 0; m < number; m++) {
        xcenter = round(((double)rand() / RAND_MAX) * (double(scope(0, 1) - scope(0, 0) + 1))) + scope(0, 0);
        ycenter = round(((double)rand() / RAND_MAX) * (double(scope(1, 1) - scope(1, 0) + 1))) + scope(1, 0);
        zcenter = round(((double)rand() / RAND_MAX) * (double(scope(2, 1) - scope(2, 0) + 1))) + scope(2, 0);
        lx = round(((double)rand() / RAND_MAX) * (limitx2 - limitx1) + limitx1);
        ly = round(((double)rand() / RAND_MAX) * (limity2 - limity1) + limity1);
        lz = round(((double)rand() / RAND_MAX) * (limitz2 - limitz1) + limitz1);


        //cout << "xcenter=" << xcenter << " ycenter=" << ycenter << " zcenter=" << zcenter << endl;
        //cout << "lx=" << lx << " " << "ly=" << ly << " " << "lz=" << lz << endl;



        for (int i = 0; i <= Nparax - 1; i++) {
            for (int j = 0; j <= Nparay - 1; j++) {
                for (int k = 0; k <= Nparaz - 1; k++) {
                    int pos = k + Nparaz * (j + Nparay * i);
                    double x, y, z;      //center of one 2D para region
                    x = bind(0) * (2 * i + 1) / 2;
                    y = bind(1) * (2 * j + 1) / 2;
                    z = bind(2) * (2 * k + 1) / 2;

                    if ((abs(x - xcenter) <= lx) && (abs(y - ycenter) <= ly) && (abs(z - zcenter) <= lz)) {
                        Para(pos) = initial_diel_func("ONES");
                    }
                }
            }
        }
    }



}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_) {
    Filter = false;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);
    
    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }
    geometryPara = VectorXi::Zero(N);
    scope = find_scope_3_dim(&geometry);

    /*
    list<VectorXi*> ParaCheck;
    for (list<VectorXi*>::iterator ittmp = FParaGeometry_.begin(); ittmp != FParaGeometry_.end(); ittmp++) {
        ParaCheck.push_back(*ittmp);
    }
    for (list<VectorXi*>::iterator ittmp = BParaGeometry_.begin(); ittmp != BParaGeometry_.end(); ittmp++) {
        ParaCheck.push_back(*ittmp);
    }
    cout << ParaCheck.size() << endl;
    if (!CheckOverlapList(ParaCheck)) {
        cout << "SpacePara::SpacePara(Space* space_, string initial_diel, Vector3i bind_, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_)" << endl;
        cout << "ParaCheck Fail" << endl;
        throw 1;
    }
    */

    VectorXi FParaGeometry = ConnectGeometry(FParaGeometry_);
    set<vector<int>> FParaGeometrySet = Get3divSet(&FParaGeometry);
    MatrixXi FParascope = find_scope_3_dim(&FParaGeometry);
    int NFparax, NFparay, NFparaz, NFpara;
    NFparax = ceil(double(FParascope(0, 1) - FParascope(0, 0) + 1) / bind(0));
    NFparay = ceil(double(FParascope(1, 1) - FParascope(1, 0) + 1) / bind(1));
    NFparaz = ceil(double(FParascope(2, 1) - FParascope(2, 0) + 1) / bind(2));
    NFpara = NFparax * NFparay * NFparaz;
    
    VectorXd Para1 = initial_diel_func(initial_diel, NFpara);
    FreeparatoPara = VectorXi::Zero(NFpara);                              //Points to the position of free para in Para. In this function, actually it is the first NFpara elements in Para.
    for (int i = 0; i <= NFpara - 1; i++) {
        FreeparatoPara(i) = i;
    }
    
    if (int(BParaGeometry_.size()) != int(BPara_.size())) {
        cout << "SpacePara::SpacePara(Space* space_, string initial_diel, Vector3i bind_, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_): BParaGeometry_.size() != BPara_.size()" << endl;
        throw 1;
    }
    list<VectorXi*>::iterator it1 = BParaGeometry_.begin();
    list<int> BParaDividePos;
    int Npara = NFpara;
    while (it1 != BParaGeometry_.end()) {
        BParaDividePos.push_back(Npara);       //The first pos for BPara is NFpara. So this line is in front of the update.
        Npara += Get3divSize(*it1);
        it1++;
    }
    Para = VectorXd::Zero(Npara);
    for (int i = 0; i <= NFpara - 1; i++) {
        Para(i) = Para1(i);
    }

    it1 = BParaGeometry_.begin();
    list<double>::iterator it2 = BPara_.begin();
    int Parapos = NFpara;
    while (it2 != BPara_.end()) {
        int tmpN = Get3divSize(*it1);
        for (int i = 0; i <= tmpN - 1; i++) {
            Para(Parapos) = (*it2);
            Parapos++;
        }
        it1++;
        it2++;
    }

    list<set<vector<int>>> BParaGeometrySetList=Get3divSetList(BParaGeometry_);
    
    VectorXi BParaCurrentPos = VectorXi::Zero(BParaGeometry_.size());
    list<int>::iterator it_tmp2 = BParaDividePos.begin();

    //cout << BParaGeometry_.size() - 1 << endl;
    //cout << int(BParaGeometry_.size()) - 1 << endl;                       These two are actually different when size=0

    for (int i = 0; i <= int(BParaGeometry_.size()) - 1; i++) {
        BParaCurrentPos(i) = (*it_tmp2);
        it_tmp2++;
    }

    for (int i = 0; i <= N - 1; i++) {
        int x = geometry(3 * i);
        int y = geometry(3 * i + 1);
        int z = geometry(3 * i + 2);
        vector<int> tmp{ x,y,z };
        if (FParaGeometrySet.count(tmp)) {
            int parax = floor((double(x) - FParascope(0, 0)) / bind(0));
            int paray = floor((double(y) - FParascope(1, 0)) / bind(1));
            int paraz = floor((double(z) - FParascope(2, 0)) / bind(2));
            int pos = paraz + NFparaz * (paray + NFparay * parax);
            geometryPara(i) = pos;                                     //The first NFpara elements in Para is the elements in FreeParatoPara 
        }
        else {
            int j = 0;
            list<set<vector<int>>>::iterator it_tmp1 = BParaGeometrySetList.begin();
            
            while (it_tmp1 != BParaGeometrySetList.end()) {
                if ((*it_tmp1).count(tmp)) {
                    geometryPara(i) = BParaCurrentPos(j);
                    BParaCurrentPos(j) = BParaCurrentPos(j) + 1;
                }
                it_tmp1++;
                j++;
            }
            
        }

        
        
       
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_, bool Filter_, FilterOption* Filterstats_) {
    Filter = Filter_;
    space = space_;
    bind = bind_;
    VectorXi* total_space = (*space).get_total_space();
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    geometry = VectorXi::Zero(3 * N);

    list<Structure>* ln = (*space).get_ln();
    list<Structure>::iterator it = (*ln).begin();
    int n1 = 0;
    for (int i = 0; i <= int((*ln).size()) - 1; i++) {
        int n2 = 3 * ((*it).get_geometry_size());
        for (int j = 0; j <= n2 - 1; j++) {
            geometry(n1 + j) = (*((*it).get_geometry()))(j);
        }
        n1 = n1 + n2;
        it++;
    }
    geometryPara = VectorXi::Zero(N);
    scope = find_scope_3_dim(&geometry);

    VectorXi FParaGeometry = ConnectGeometry(FParaGeometry_);
    set<vector<int>> FParaGeometrySet = Get3divSet(&FParaGeometry);
    MatrixXi FParascope = find_scope_3_dim(&FParaGeometry);
    int NFparax, NFparay, NFparaz, NFpara;
    NFparax = ceil(double(FParascope(0, 1) - FParascope(0, 0) + 1) / bind(0));
    NFparay = ceil(double(FParascope(1, 1) - FParascope(1, 0) + 1) / bind(1));
    NFparaz = ceil(double(FParascope(2, 1) - FParascope(2, 0) + 1) / bind(2));
    NFpara = NFparax * NFparay * NFparaz;

    VectorXd Para1 = initial_diel_func(initial_diel, NFpara);
    FreeparatoPara = VectorXi::Zero(NFpara);                              //Points to the position of free para in Para. In this function, actually it is the first NFpara elements in Para.
    for (int i = 0; i <= NFpara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    if (int(BParaGeometry_.size()) != int(BPara_.size())) {
        cout << "SpacePara::SpacePara(Space* space_, string initial_diel, Vector3i bind_, list<VectorXi*> FParaGeometry_, list<VectorXi*> BParaGeometry_, list<double> BPara_): BParaGeometry_.size() != BPara_.size()" << endl;
        throw 1;
    }
    list<VectorXi*>::iterator it1 = BParaGeometry_.begin();
    list<int> BParaDividePos;
    int Npara = NFpara;
    while (it1 != BParaGeometry_.end()) {
        BParaDividePos.push_back(Npara);       //The first pos for BPara is NFpara. So this line is in front of the update.
        Npara += Get3divSize(*it1);
        it1++;
    }
    Para = VectorXd::Zero(Npara);
    for (int i = 0; i <= NFpara - 1; i++) {
        Para(i) = Para1(i);
    }

    it1 = BParaGeometry_.begin();
    list<double>::iterator it2 = BPara_.begin();
    int Parapos = NFpara;
    while (it2 != BPara_.end()) {
        int tmpN = Get3divSize(*it1);
        for (int i = 0; i <= tmpN - 1; i++) {
            Para(Parapos) = (*it2);
            Parapos++;
        }
        it1++;
        it2++;
    }

    list<set<vector<int>>> BParaGeometrySetList = Get3divSetList(BParaGeometry_);

    VectorXi BParaCurrentPos = VectorXi::Zero(BParaGeometry_.size());
    list<int>::iterator it_tmp2 = BParaDividePos.begin();

    //cout << BParaGeometry_.size() - 1 << endl;
    //cout << int(BParaGeometry_.size()) - 1 << endl;                       These two are actually different when size=0

    for (int i = 0; i <= int(BParaGeometry_.size()) - 1; i++) {
        BParaCurrentPos(i) = (*it_tmp2);
        it_tmp2++;
    }

    for (int i = 0; i <= N - 1; i++) {
        int x = geometry(3 * i);
        int y = geometry(3 * i + 1);
        int z = geometry(3 * i + 2);
        vector<int> tmp{ x,y,z };
        if (FParaGeometrySet.count(tmp)) {
            int parax = floor((double(x) - FParascope(0, 0)) / bind(0));
            int paray = floor((double(y) - FParascope(1, 0)) / bind(1));
            int paraz = floor((double(z) - FParascope(2, 0)) / bind(2));
            int pos = paraz + NFparaz * (paray + NFparay * parax);
            geometryPara(i) = pos;                                     //The first NFpara elements in Para is the elements in FreeParatoPara 
        }
        else {
            int j = 0;
            list<set<vector<int>>>::iterator it_tmp1 = BParaGeometrySetList.begin();

            while (it_tmp1 != BParaGeometrySetList.end()) {
                if ((*it_tmp1).count(tmp)) {
                    geometryPara(i) = BParaCurrentPos(j);
                    BParaCurrentPos(j) = BParaCurrentPos(j) + 1;
                }
                it_tmp1++;
                j++;
            }
        }
    }
    if (Filter == true) {
        Para_origin = Para;
        Para_filtered = Para;
        Filterstats = Filterstats_;
        if (Filterstats == NULL) {
            cout << "ERROR: SpacePara::SpacePara: Filter==true then Filterstats must be passed in." << endl;
            throw 1;
        }
        vector<vector<int>> Paratogeometry(Npara);
        for (int i = 0; i <= N - 1; i++) {
            (Paratogeometry[geometryPara(i)]).push_back(i);
        }

        //Only works when freepara is 2D binding (Only do the filter in 2D)
        FreeWeight = vector<vector<WeightPara>>(NFpara);
        for (int i = 0; i <= NFpara - 1; i++) {
            int poso = Paratogeometry[FreeparatoPara(i)][0];                 //As 2D extrusion is assumed, different z does not matter
            int xo = geometry(3 * poso);
            int yo = geometry(3 * poso + 1);
            int zo = geometry(3 * poso + 2);

            double rfilter = (*Filterstats).get_rfilter();

            for (int j = 0; j <= NFpara - 1; j++) {
                //bool inornot = false;
                for (int k = 0; k <= Paratogeometry[FreeparatoPara(j)].size() - 1; k++) {
                    int posr = Paratogeometry[FreeparatoPara(j)][k];
                    
                    int xr = geometry(3 * posr);
                    int yr = geometry(3 * posr + 1);
                    int zr = geometry(3 * posr + 2);
                    if ((zo==zr) && (circlerange(xo, yo, xr, yr, rfilter))) {
                        //1. Same xy plane 2. inside the circle in xy plane 
                        //para>=2 wont be in NFpara
                        int posweight = FreeparatoPara(j);
                        double weight = calweight(xo, yo, xr, yr, rfilter);
                        FreeWeight[i].push_back(WeightPara{ weight,posweight });
                        break;
                        //As soon as one in the entire z direction is verified, no need for looking at others for 1 j.
                    }

                }
                
            }
        }
    }


}

void SpacePara::ChangeBind(Vector3i bind_) {
    bind = bind_;
    VectorXi geometryPara_before = geometryPara;
    VectorXd Para_before = Para;
    int Nx, Ny, Nz, N;
    tie(Nx, Ny, Nz, N) = (*space).get_Ns();
    int Nparax, Nparay, Nparaz, Npara;
    Nparax = ceil((scope(0, 1) - scope(0, 0) + 1) / bind(0));
    Nparay = ceil((scope(1, 1) - scope(1, 0) + 1) / bind(1));
    Nparaz = ceil((scope(2, 1) - scope(2, 0) + 1) / bind(2));
    Npara = Nparax * Nparay * Nparaz;

    vector<vector<int>> Paratogeometry(Npara, vector<int>(2, 0));
    for (int i = 0; i <= N - 1; i++) {
        double x = geometry(3 * i);
        double y = geometry(3 * i + 1);
        double z = geometry(3 * i + 2);
        int parax = floor((x - scope(0, 0)) / bind(0));
        int paray = floor((y - scope(1, 0)) / bind(1));
        int paraz = floor((z - scope(2, 0)) / bind(2));
        int pos = paraz + Nparaz * (paray + Nparay * parax);
        geometryPara(i) = pos;
        Paratogeometry[pos][0] += 1;
        Paratogeometry[pos][1] += Para_before(geometryPara_before(i));
    }

    for (int i = 0; i <= Npara - 1; i++) {
        Para(i) = Paratogeometry[i][1] / Paratogeometry[i][0];
    }
    FreeparatoPara = VectorXi::Zero(Npara);
    for (int i = 0; i <= Npara - 1; i++) {
        FreeparatoPara(i) = i;
    }

    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }

    


}

Space* SpacePara::get_space() {
    return space;
}

VectorXi SpacePara::get_geometry() {
    return geometry;
}

VectorXi* SpacePara::get_geometryPara() {
    return &geometryPara;
}

VectorXd* SpacePara::get_Para() {
    return &Para;
}

VectorXd* SpacePara::get_Para_origin() {
    if (!Filter) {
        cout << "ERROR: SpacePara::get_Para_origin()--Filter can not be false" << endl;
        throw 1;
    }
    return &Para_origin;
}

VectorXd* SpacePara::get_Para_filtered() {
    if (!Filter) {
        cout << "ERROR: SpacePara::get_Para_filtered()--Filter can not be false" << endl;
        throw 1;
    }
    return &Para_filtered;
}

Vector3i* SpacePara::get_bind() {
    return &bind;
}

VectorXi* SpacePara::get_Free() {
    return &FreeparatoPara;
}

bool SpacePara::get_Filter() {
    return Filter;
}

FilterOption* SpacePara::get_Filterstats() {
    if (!Filter) {
        cout << "ERROR: SpacePara::get_Filterstats()--Filter can not be false" << endl;
        throw 1;
    }
    return Filterstats;
}

vector<vector<WeightPara>>* SpacePara::get_FreeWeight() {
    if (!Filter) {
        cout << "ERROR: SpacePara::get_FreeWeight()---Filter can not be false" << endl;
        throw 1;
    }
    return &FreeWeight;
}
