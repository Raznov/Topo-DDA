#include "definition.h"

VectorXi *Structure::get_geometry(){
    return &geometry;
}

VectorXd *Structure::get_diel(){
    return &diel;
}

int Structure::get_para(){
    return para;
}
int Structure::get_geometry_size(){
    if(geometry.size()==diel.size()){
        return round(geometry.size()/3);
    }
    else{
        cout<<"geometry size not equal to diel size"<<endl;
        return round(geometry.size()/3);
    }
}

//Find the Max and Min of the input geometry in each direction
MatrixXi find_scope_3_dim(VectorXi *x){                          
    int N=round((*x).size()/3);
    MatrixXi result(3,2);
    result(0,0)=(*x)(0);
    result(0,1)=(*x)(0);
    result(1,0)=(*x)(1);
    result(1,1)=(*x)(1);
    result(2,0)=(*x)(2);
    result(2,1)=(*x)(2);

    for(int i=0;i<=N-1;i++){
        if(result(0,0)>=(*x)(3*i)){
            result(0,0)=(*x)(3*i);
        }
        if(result(0,1)<=(*x)(3*i)){
            result(0,1)=(*x)(3*i);
        }
        if(result(1,0)>=(*x)(3*i+1)){
            result(1,0)=(*x)(3*i+1);
        }
        if(result(1,1)<=(*x)(3*i+1)){
            result(1,1)=(*x)(3*i+1);
        }
        if(result(2,0)>=(*x)(3*i+2)){
            result(2,0)=(*x)(3*i+2);
        }
        if(result(2,1)<=(*x)(3*i+2)){
            result(2,1)=(*x)(3*i+2);
        }
    }
    return result;
}

void Structure::cut(VectorXi *big, VectorXi *smalll, VectorXd *small_diel){
    
    int number_origin=round((*smalll).size()/3);
    MatrixXi big_scope=find_scope_3_dim(big);
    //cout<<"big_scope "<<big_scope<<endl;
    list<int> positions_in;
    int number_out=0;
    //cout<<"small_scope "<<find_scope_3_dim(small)<<endl;
    for(int i=0;i<=number_origin-1;i++){
        if(((*smalll)(3*i)<big_scope(0,0))||((*smalll)(3*i)>big_scope(0,1))||
           ((*smalll)(3*i+1)<big_scope(1,0))||((*smalll)(3*i+1)>big_scope(1,1))||
           ((*smalll)(3*i+2)<big_scope(2,0))||((*smalll)(3*i+2)>big_scope(2,1))){
               number_out+=1;
           }
        else{
            positions_in.push_back(i);
        }
    }
    int number_in=positions_in.size();
    geometry=VectorXi::Zero(3*number_in);
    diel=VectorXd::Zero(3*number_in);
    for(int i=0;i<=number_in-1;i++){
        int j=positions_in.front();
        positions_in.pop_front();
        geometry(3*i)=(*smalll)(3*j);
        geometry(3*i+1)=(*smalll)(3*j+1);
        geometry(3*i+2)=(*smalll)(3*j+2);
        diel(3*i)=(*small_diel)(3*j);
        diel(3*i+1)=(*small_diel)(3*j+1);
        diel(3*i+2)=(*small_diel)(3*j+2);
    }
    if(number_out==0){
        cout<<"The geometry you built is entirely in the space."<<endl;
        cout<<"number_origin "<<number_origin<<endl;
        cout<<"number_real "<<number_in<<endl;
    }
    else{
        cout<<"The geometry you built is at least partially outside of space."<<endl;
        cout<<"number_origin "<<number_origin<<endl;
        cout<<"number_out "<<number_out<<endl;
        cout<<"number_real "<<number_in<<endl;
    }
}

VectorXd initial_diel_func(string initial_diel, int N){
    VectorXd diel;
    if(initial_diel.compare("ZEROS")==0){
        diel=VectorXd::Zero(N);
    }
    else if(initial_diel.compare("ONES")==0){
        diel=VectorXd::Ones(N);
    }
    else if(initial_diel.compare("RANDOM")==0){
        int n=round(N/3);
        diel=VectorXd::Zero(N);
        for(int i=0;i<=n-1;i++){
            double r=((double)rand()/(RAND_MAX));
            diel(3*i)=r;
            diel(3*i+1)=r;
            diel(3*i+2)=r;
        }
        
    }
    else{
        diel=VectorXd::Zero(N);
        cout<<"The initial type given does not match any of the built in method"<<endl;
    }
    return diel;
}

Structure::Structure(VectorXi *total_space, VectorXd *diel_, VectorXi *geometry_, int para_){
    if((*diel_).size()!=(*geometry_).size()){
        cout<<"diel and geometry should have the same size which is not true"<<endl;
        geometry=*geometry_;
        diel=*diel_;
        para=para_;
    }
    else{
        cut(total_space, geometry_, diel_);             
        para=para_;
    if(diel.size()!=geometry.size()){
        cout<<"final size not equal"<<endl;
    }
    }
}

Structure::Structure(VectorXi *total_space, string initial_diel, double r, Vector3d center, int para_){
    int N=round((*total_space).size()/3);
    list<int> positions; 
    for(int i=0;i<=N-1;i++){
        double x=(*total_space)(3*i)-center(0);
        double y=(*total_space)(3*i+1)-center(1);
        double z=(*total_space)(3*i+2)-center(2);
        if(x*x+y*y+z*z<(r+0.1)*(r+0.1)){
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
    VectorXd diel_tmp=initial_diel_func(initial_diel, 3*N_want);
    cut(total_space, &geometry_tmp, &diel_tmp);
    para=para_;
    
    if(diel.size()!=geometry.size()){
        cout<<"diel and geometry should have the same size which is not true"<<endl;
    }

    

}

Structure::Structure(VectorXi *total_space, string initial_diel, double r, Vector3i center, Vector3i direction, int para_){
    int N=round((*total_space).size()/3);
    list<int> positions; 
    for(int i=0;i<=N-1;i++){
        int x=(*total_space)(3*i)-center(0);
        int y=(*total_space)(3*i+1)-center(1);
        int z=(*total_space)(3*i+2)-center(2);
        if((x*x+y*y+z*z<(r+0.1)*(r+0.1))&&(x*direction(0)+y*direction(1)+z*direction(2)==0)){
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
    VectorXd diel_tmp=initial_diel_func(initial_diel, 3*N_want);
    cut(total_space, &geometry_tmp, &diel_tmp);
    para=para_;
    
    if(diel.size()!=geometry.size()){
        cout<<"diel and geometry should have the same size which is not true"<<endl;
    }
    
}

Structure::Structure(VectorXi *total_space, string initial_diel, Vector3d l, Vector3d center, int para_){
    int N=round((*total_space).size()/3);
    list<int> positions; 
    for(int i=0;i<=N-1;i++){
        double x=(*total_space)(3*i)-center(0);
        double y=(*total_space)(3*i+1)-center(1);
        double z=(*total_space)(3*i+2)-center(2);
        if((abs(x)<=l(0)/2)&&(abs(y)<=l(1)/2)&&(abs(z)<=l(2)/2)){
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
    VectorXd diel_tmp=initial_diel_func(initial_diel, 3*N_want);
    cut(total_space, &geometry_tmp, &diel_tmp);
    para=para_;      
    
    if(diel.size()!=geometry.size()){
        cout<<"diel and geometry should have the same size which is not true"<<endl;
    }

}

Structure::Structure(VectorXi *total_space, Structure *s, Vector3i direction, int times, int para_){
    VectorXi *geometry_=(*s).get_geometry();
    VectorXd *diel_=(*s).get_diel();
    int N=round((*geometry_).size()/3);
    if(N!=round((*diel_).size()/3)){
        cout<<"Original Structure geometry diel not same length when duplicating"<<endl;
    }
    MatrixXi period=find_scope_3_dim(geometry_);
    int periodx=period(0,1)-period(0,0)+1;
    int periody=period(1,1)-period(1,0)+1;
    int periodz=period(2,1)-period(2,0)+1;
    int N_new=N*(times-1);                                                       //This build the duplicated geometry without original one
    VectorXi geometry_tmp=VectorXi::Zero(3*N_new);
    VectorXd diel_tmp=VectorXd::Zero(3*N_new);                                         //min(times)=2
    for(int i=0;i<=times-2;i++){                                                        //Some cautious about how this is stored.
        for(int j=0;j<=N-1;j++){
            geometry_tmp(3*(j+i*N))=(*geometry_)(3*j)+periodx*(i+1)*direction(0);
            geometry_tmp(3*(j+i*N)+1)=(*geometry_)(3*j+1)+periody*(i+1)*direction(1);
            geometry_tmp(3*(j+i*N)+2)=(*geometry_)(3*j+2)+periodz*(i+1)*direction(2);

            diel_tmp(3*(j+i*N))=(*diel_)(3*j);
            diel_tmp(3*(j+i*N)+1)=(*diel_)(3*j+1);
            diel_tmp(3*(j+i*N)+2)=(*diel_)(3*j+2);
        }
    }
    //cout<<geometry_tmp<<endl;
    cut(total_space, &geometry_tmp, &diel_tmp);
    para=para_;
    
}

Structure::Structure(VectorXi *total_space, string FileName, int para_){
    string GeometryName = FileName + "Geometry.txt";
    string DielName =FileName + "Diel.txt";
    ifstream GeometryFile;
    ifstream DielFile;
    GeometryFile.open(GeometryName);
    DielFile.open(DielName);


    list<int> GeometryList;
    int a;
    list<double> DielList;
    double b;
    while(GeometryFile>>a){
        GeometryList.push_back(a);
    }
    while(DielFile>>b){
        DielList.push_back(b);
    }
    GeometryFile.close();
    DielFile.close();
    if(GeometryList.size() != DielList.size()){
        cout<<"File input geometry does not the same Geometry and Diel size."<<endl;
    }


    int GDSize = GeometryList.size();
    VectorXi TmpGeometry = VectorXi::Zero(GDSize);
    VectorXd TmpDiel = VectorXd::Zero(GDSize);
    list<int>::iterator ItGeometryList = GeometryList.begin();
    list<double>::iterator ItDielList = DielList.begin();
    for(int i=0; i<=GDSize-1; i++){
        TmpGeometry(i) = (*ItGeometryList);
        TmpDiel(i) = (*ItDielList);
        ItGeometryList++;
        ItDielList++;
    }
    
    cut(total_space, &TmpGeometry, &TmpDiel);
    if(diel.size()!=geometry.size()){
        cout<<"diel and geometry should have the same size which is not true"<<endl;
    }
    para = para_;
}   