#include "definition.h"


Structure::Structure(VectorXi *total_space, VectorXi *geometry_, bool para_cond_){
    geometry = *geometry_;
    para_cond = para_cond_;
}
Structure::Structure(VectorXi* total_space, double r, Vector3d center, bool para_cond_) {
    para_cond = para_cond_;
    int N = round((*total_space).size() / 3);
    list<int> positions;
    for (int i = 0; i <= N - 1; i++) {
        double x = (*total_space)(3 * i) - center(0);
        double y = (*total_space)(3 * i + 1) - center(1);
        double z = (*total_space)(3 * i + 2) - center(2);
        if (x * x + y * y + z * z < (r + 0.1) * (r + 0.1)) {
            positions.push_back(i);
        }
    }
    int N_want = positions.size();
    VectorXi geometry_tmp = VectorXi::Zero(3 * N_want);
    for (int i = 0; i <= N_want - 1; i++) {
        int j = positions.front();
        positions.pop_front();
        geometry_tmp(3 * i) = (*total_space)(3 * j);
        geometry_tmp(3 * i + 1) = (*total_space)(3 * j + 1);
        geometry_tmp(3 * i + 2) = (*total_space)(3 * j + 2);
    }
    cut(total_space, &geometry_tmp);
}
Structure::Structure(VectorXi* total_space, double r, double h, Vector3d center, bool para_cond_) {
    para_cond = para_cond_;
    int N = round((*total_space).size() / 3);
    list<int> positions;
    for (int i = 0; i <= N - 1; i++) {
        double x = (*total_space)(3 * i) - center(0);
        double y = (*total_space)(3 * i + 1) - center(1);
        double z = (*total_space)(3 * i + 2) - center(2);
        if ((x * x + y * y < (r + 0.1) * (r + 0.1)) && (abs(z)<(h/2)+0.1)) {
            positions.push_back(i);
        }
    }
    int N_want = positions.size();
    VectorXi geometry_tmp = VectorXi::Zero(3 * N_want);
    for (int i = 0; i <= N_want - 1; i++) {
        int j = positions.front();
        positions.pop_front();
        geometry_tmp(3 * i) = (*total_space)(3 * j);
        geometry_tmp(3 * i + 1) = (*total_space)(3 * j + 1);
        geometry_tmp(3 * i + 2) = (*total_space)(3 * j + 2);
    }
    cut(total_space, &geometry_tmp);
}


Structure::Structure(VectorXi *total_space, Vector3d l, Vector3d center, bool para_cond_){
    para_cond = para_cond_;
    int N=round((*total_space).size()/3);
    list<int> positions; 
    for(int i=0;i<=N-1;i++){
        double x=double((*total_space)(3*i))-center(0);
        double y=double((*total_space)(3*i+1))-center(1);
        double z=double((*total_space)(3*i+2))-center(2);
        if((abs(x)<=l(0)/2+0.001)&&(abs(y)<=l(1)/2 + 0.001)&&(abs(z)<=l(2)/2 + 0.001)){
            positions.push_back(i);
        }   
    }
    int N_want=positions.size();
    VectorXi geometry_tmp=VectorXi::Zero(3*N_want);
    for(int i=0;i<=N_want-1;i++){
        int j=positions.front();
        positions.pop_front();
        geometry_tmp(3*i)=(*total_space)(3*j);
        geometry_tmp(3*i+1)=(*total_space)(3*j+1);
        geometry_tmp(3*i+2)=(*total_space)(3*j+2);
    }
    cut(total_space, &geometry_tmp);

}
Structure::Structure(VectorXi* total_space, Vector3d l, Vector3d center, Structure* Str, bool para_cond_) {
    
    para_cond = para_cond_;
    int N = round((*total_space).size() / 3);
    vector<int> positions;
    set<vector<int>> Strpoints;
    VectorXi* Strgeo = (*Str).get_geometry();
    for (int i = 0; i < (*Str).get_geometry_size(); i++) {
        Strpoints.insert(vector<int>{(*Strgeo)(3 * i), (*Strgeo)(3 * i + 1), (*Strgeo)(3 * i + 2)});
    }

    for (int i = 0; i <= N - 1; i++) {
        int xi = (*total_space)(3 * i);
        int yi = (*total_space)(3 * i + 1);
        int zi = (*total_space)(3 * i + 2);
        if (!Strpoints.count(vector<int>{xi, yi, zi})) {
            double x = double((*total_space)(3 * i)) - center(0);
            double y = double((*total_space)(3 * i + 1)) - center(1);
            double z = double((*total_space)(3 * i + 2)) - center(2);
            if ((abs(x) <= l(0) / 2 + 0.001) && (abs(y) <= l(1) / 2 + 0.001) && (abs(z) <= l(2) / 2 + 0.001)) {
                positions.push_back(i);
            }
        }
    }
    int N_want = positions.size();
    VectorXi geometry_tmp = VectorXi::Zero(3 * N_want);
    for (int i = 0; i <= N_want - 1; i++) {
        int j = positions[i];
        geometry_tmp(3 * i) = (*total_space)(3 * j);
        geometry_tmp(3 * i + 1) = (*total_space)(3 * j + 1);
        geometry_tmp(3 * i + 2) = (*total_space)(3 * j + 2);
    }
    cut(total_space, &geometry_tmp);

}
VectorXi* Structure::get_geometry() {
    return &geometry;
}
void Structure::cut(VectorXi* big, VectorXi* smalll) {

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
    geometry = VectorXi::Zero(3 * number_in);
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
}
int Structure::get_geometry_size() {
    return round(geometry.size() / 3);
}
bool Structure::para_or_not() {
    return para_cond;
}
