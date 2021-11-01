#include "definition.h"
VectorXi *Space::get_total_space(){
    return total_space;
}
int Space::get_ln_size(){
    return int((*ln).size());
}
list<Structure> *Space::get_ln(){
    return ln;
}
tuple<int, int, int, int> Space::get_Ns(){
    return make_tuple(Nx, Ny, Nz, N);
}
Space::Space(VectorXi *total_space_, int Nx_, int Ny_, int Nz_, int N_, list<Structure> *ln_){
    total_space=total_space_;
    Nx=Nx_;
    Ny=Ny_;
    Nz=Nz_;
    N=N_;
    ln=ln_;
}
void Space::show_something_about_Structures() const{ 
    ofstream fout("Space.txt");
    fout<<Nx<<endl<<Ny<<endl<<Nz<<endl<<N<<endl;
    list<Structure>::iterator it=(*ln).begin();
    for(int i=0;i<=int((*ln).size())-1;i++){
        fout<<*((*it).get_geometry())<<endl;
        it++;
    }

    /*
    it=(*ln).begin();
    for(int i=0;i<=int((*ln).size())-1;i++){
        fout<<*((*it).get_diel())<<endl;
        it++;
    }

    it=(*ln).begin();
    for(int i=0;i<=int((*ln).size())-1;i++){
        int n=(*((*it).get_diel())).size();
        for(int j=0;j<=n-1;j++){
            fout<<((*it).get_para())<<endl;
        }
        it++;
    }
    */



    fout.close();
}


Space operator+(const Space &s1, Structure &s2){
    list<Structure> *tmp=s1.ln;
    (*tmp).push_back(s2);
    return Space(s1.total_space, s1.Nx, s1.Ny, s1.Nz, s1.N+s2.get_geometry_size(), tmp);
}

