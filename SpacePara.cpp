#include "definition.h"

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

    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, VectorXi* geometryPara_, VectorXd* Para_, VectorXi* FreeparatoPara_) {
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
    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel) {
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
    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel, VectorXi* geometryPara_) {
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
    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel_center, string initial_diel_ring, double r, string type) {
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
    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, string initial_diel_background, list<string>* initial_diel_list, list<double>* r_list, list<Vector2d>* center_list){
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
    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2) {
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
    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
  
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2) {
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
    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }

}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, VectorXi* geometryPara_) {
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

    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
    }
}

SpacePara::SpacePara(Space* space_, Vector3i bind_, int number, double limitx1, double limitx2, double limity1, double limity2, double limitz1, double limitz2, VectorXi* geometryPara_) {
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

    Paratogeometry.resize(Para.size());
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[geometryPara(i)]).push_back(i);
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

Vector3i* SpacePara::get_bind() {
    return &bind;
}

VectorXi* SpacePara::get_Free() {
    return &FreeparatoPara;
}

vector<list<int>>* SpacePara::get_Paratogeometry() {
    return &Paratogeometry;
}