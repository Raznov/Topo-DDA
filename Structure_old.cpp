#include "definition.h"
/*
class Structure{
    private:
        int N;                           //N=Nx*Ny*Nz
        int Nx;
        int Ny;
        int Nz;
        int time;
        int ITERATION;
        double ERROR;
        double d;                        //distance between neighboring diples
        VectorXd R;                      //coordinates of dipoles
        double lam;                       //wavelength
        double K;                        //K value
        Vector3d n_K;                      //wave vector direction
        
        double E0;                       //Electric field intensity
        Vector3d n_E0;                   //Electric field direction
        VectorXcd E;                      //Electric field at each dipole
        VectorXcd P;                     //Polarization
        VectorXcd diel;                  //diel function of each dipole(the actual complex dielectric function)
        VectorXcd al;                    //alpha
        //VectorXd parameter;              //the parameter epsilon which diel=diel_sub+epsilon*(diel_mat-diel_sub)
        string ptype;                    //name of parameter type, 2D or 3D. Determines the way of building the geometry.
        Vector2cd material;              //[0]for substrate, [1]for material interested

        bool verbose;
    public:
        Structure(int Nx_, int Ny_, int Nz_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, 
        VectorXd parameter_, string ptype_, Vector2cd material_);//material is got from name to number outside
        Matrix3cd A_dic_generator(double x,double y,double z);
        VectorXcd Aproduct(vector<vector<vector<Matrix3cd>>> &A_dic,VectorXcd &b);
        void bicgstab(int MAX_ITERATION,double MAX_ERROR);
        void change_E(VectorXcd E_);
        VectorXcd get_diel();
        Vector2cd get_material();
        VectorXcd get_alpha();
};
*/
Structure::Structure(int Nx_, int Ny_, int Nz_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, VectorXd parameter_, string ptype_, Vector2cd material_){
    Nx=Nx_;
    Ny=Ny_;
    Nz=Nz_;
    N=Nx*Ny*Nz;
    time=0;
    ITERATION=0;
    ERROR=0.0;
    d=d_;
    lam=lam_;
    K=2*M_PI/lam;
    E0=E0_;
    n_K=n_K_;
    n_E0=n_E0_;
    R=VectorXd::Zero(3*N);
    diel=VectorXcd::Zero(3*N);
    ptype=ptype_;
    material=material_;
    if(ptype.compare("2D")==0){
        for(int x=0;x<=Nx-1;x++){
            for(int y=0;y<=Ny-1;y++){
                for(int z=0;z<=Nz-1;z++){
                    int position=z+Nz*(y+Ny*x);
                    R(3*position)=x;
                    R(3*position+1)=y;
                    R(3*position+2)=z;
                    diel(3*position)=material(0)+parameter_(y+Ny*x)*(material(1)-material(0));
                    diel(3*position+1)=material(0)+parameter_(y+Ny*x)*(material(1)-material(0));
                    diel(3*position+2)=material(0)+parameter_(y+Ny*x)*(material(1)-material(0));
                }
            }
        }
    }
    else if(ptype.compare("3D")==0){
        for(int x=0;x<=Nx-1;x++){
            for(int y=0;y<=Ny-1;y++){
                for(int z=0;z<=Nz-1;z++){
                    int position=z+Nz*(y+Ny*x);
                    R(3*position)=x;
                    R(3*position+1)=y;
                    R(3*position+2)=z;
                    diel(3*position)=material(0)+parameter_(position)*(material(1)-material(0));
                    diel(3*position+1)=material(0)+parameter_(position)*(material(1)-material(0));
                    diel(3*position+2)=material(0)+parameter_(position)*(material(1)-material(0));
                }
            }
        }
    }
    else{
        cout<<"ERROR: neither 2D nor 3D when initializing the geometry."<<endl;
        R=R;
        diel=diel;
    }

    P = VectorXcd::Zero(N*3);
    E = VectorXcd::Zero(N*3);
    for (int i=0;i<N;i++) {
        E(3*i) = E0*n_E0(0)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);
        E(3*i+1) = E0*n_E0(1)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);
        E(3*i+2) = E0*n_E0(2)*(cos(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))+sin(K*d*(n_K(0)*R(3*i)+n_K(1)*R(3*i+1)+n_K(2)*R(3*i+2)))*1i);                                            
    }
    al = VectorXcd::Zero(N*3);
    for (int i=0;i<N*3;i++) {
        al(i)=1.0/Get_Alpha(lam,K,d,diel(i));
    }  
    verbose = true;  
}

Matrix3cd Structure::A_dic_generator(double x,double y,double z){
    Matrix3cd result(3,3);
    double xsquare = x*x; double ysquare = y*y; double zsquare = z*z;
    double rnorm = sqrt(xsquare+ysquare+zsquare);

    if (rnorm==0.0) {
        result(0,0) = 0.0;
        result(1,1) = 0.0;
        result(2,2) = 0.0;
        return result;
    }
    else {
        double xy = x*y; double yz = y*z; double zx = z*x;
        double rsquare = rnorm*rnorm;
        double rcubic = rnorm*rnorm*rnorm;
        std::complex<double> const1 = 1i*K*rnorm;
        const1 = exp(const1)/rcubic;
        std::complex<double> const2 = (1.0-K*rnorm*1i);
        const2 = const2/rsquare;
        std::complex<double> const3 = K*K-3.0*const2;

        result(0,0) = const2*(rsquare-3*xsquare)-K*K*(ysquare+zsquare);
        result(0,1) = xy*const3;
        result(0,2) = zx*const3;
        result(1,1) = const2*(rsquare-3*ysquare)-K*K*(zsquare+xsquare);
        result(1,2) = yz*const3;
        result(2,2) = const2*(rsquare-3*zsquare)-K*K*(xsquare+ysquare);

        result(1,0) = result(0,1);
        result(2,0) = result(0,2);
        result(2,1) = result(1,2);

        result = const1*result;
        return result;
    }
}

VectorXcd Structure::Aproduct(vector<vector<vector<Matrix3cd>>> &A_dic, VectorXcd &b){
    //! in geometry use np.meshgrid with 'ij'  
    int nx=2*Nx-1;
    int ny=2*Ny-1;
    int nz=2*Nz-1;
    int n=nx*ny*nz;
    ArrayXcd b_x_in=ArrayXcd::Zero(n);
    ArrayXcd b_y_in=ArrayXcd::Zero(n);
    ArrayXcd b_z_in=ArrayXcd::Zero(n);
    ArrayXcd A_00_in=ArrayXcd::Zero(n);
    ArrayXcd A_01_in=ArrayXcd::Zero(n);
    ArrayXcd A_02_in=ArrayXcd::Zero(n);
    ArrayXcd A_11_in=ArrayXcd::Zero(n);
    ArrayXcd A_12_in=ArrayXcd::Zero(n);
    ArrayXcd A_22_in=ArrayXcd::Zero(n);
    for(int i=0;i<=N-1;i++){
        int x=round(R(3*i));
        int y=round(R(3*i+1));
        int z=round(R(3*i+2));
        int position=z+nz*(y+ny*x);
        b_x_in(position)=b(3*i);
        b_y_in(position)=b(3*i+1);
        b_z_in(position)=b(3*i+2);
    }

    for(int x=0;x<=nx-1;x++){
        for(int y=0;y<=ny-1;y++){
            for(int z=0;z<=nz-1;z++){
                int position=z+nz*(y+ny*x);
                A_00_in(position)=A_dic[x][y][z](0,0);
                A_01_in(position)=A_dic[x][y][z](0,1);
                A_02_in(position)=A_dic[x][y][z](0,2);
                A_11_in(position)=A_dic[x][y][z](1,1);
                A_12_in(position)=A_dic[x][y][z](1,2);
                A_22_in(position)=A_dic[x][y][z](2,2);
            }
        }
    }
    
    ArrayXcd b_x_out=FFT(nx,ny,nz,b_x_in,0);
    ArrayXcd b_y_out=FFT(nx,ny,nz,b_y_in,0);
    ArrayXcd b_z_out=FFT(nx,ny,nz,b_z_in,0);
    ArrayXcd A_00_out=FFT(nx,ny,nz,A_00_in,0);
    ArrayXcd A_01_out=FFT(nx,ny,nz,A_01_in,0);
    ArrayXcd A_02_out=FFT(nx,ny,nz,A_02_in,0);
    ArrayXcd A_11_out=FFT(nx,ny,nz,A_11_in,0);
    ArrayXcd A_12_out=FFT(nx,ny,nz,A_12_in,0);
    ArrayXcd A_22_out=FFT(nx,ny,nz,A_22_in,0);
    
    ArrayXcd a_x_in=ArrayXcd::Zero(n);
    ArrayXcd a_y_in=ArrayXcd::Zero(n);
    ArrayXcd a_z_in=ArrayXcd::Zero(n);

    for(int i=0;i<=n-1;i++){
        a_x_in(i)=A_00_out(i)*b_x_out(i)+A_01_out(i)*b_y_out(i)+A_02_out(i)*b_z_out(i);
        a_y_in(i)=A_01_out(i)*b_x_out(i)+A_11_out(i)*b_y_out(i)+A_12_out(i)*b_z_out(i);
        a_z_in(i)=A_02_out(i)*b_x_out(i)+A_12_out(i)*b_y_out(i)+A_22_out(i)*b_z_out(i);
    }

    ArrayXcd a_x_out=FFT(nx,ny,nz,a_x_in,1);
    ArrayXcd a_y_out=FFT(nx,ny,nz,a_y_in,1);
    ArrayXcd a_z_out=FFT(nx,ny,nz,a_z_in,1);

    VectorXcd result(3*N);
    for(int i=0;i<=N-1;i++){
        int x=round(R(3*i));
        int y=round(R(3*i+1));
        int z=round(R(3*i+2));
        int position=z+nz*(y+ny*x);
        result(3*i)=a_x_out(position)+b(3*i)*al(3*i);
        result(3*i+1)=a_y_out(position)+b(3*i+1)*al(3*i+1);
        result(3*i+2)=a_z_out(position)+b(3*i+2)*al(3*i+2);   ///FFT can only be used for periodic function so the alpha term has to be calculated individuelly
    
    }
    return result;

}


void Structure::bicgstab(int MAX_ITERATION,double MAX_ERROR){
    if (verbose) {
        cout << "--------------Calculation start. Iterative method used: BICGSTAB---------------" << endl;
        cout << endl;
    }
    
    high_resolution_clock::time_point t_start = high_resolution_clock::now();
    vector<vector<vector<Matrix3cd>>> A_dic(2*Nx-1,vector<vector<Matrix3cd>>(2*Ny-1,vector<Matrix3cd>(2*Nz-1,Matrix3cd::Zero())));    
      
    for (int i=0;i<=Nx-1;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*i; double y=d*j; double z=d*k;
                A_dic[i][j][k] = this->A_dic_generator(x,y,z);
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*j; double z=d*k;
                A_dic[i][j][k] = this->A_dic_generator(x,y,z);             
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*i; double y=d*(j-(2*Ny-1)); double z=d*k;
                A_dic[i][j][k] = this->A_dic_generator(x,y,z);      
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*i; double y=d*j; double z=d*(k-(2*Nz-1));
                A_dic[i][j][k] = this->A_dic_generator(x,y,z);
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=0;k<=Nz-1;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*(j-(2*Ny-1)); double z=d*k;
                A_dic[i][j][k] = this->A_dic_generator(x,y,z);
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=0;j<=Ny-1;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*j; double z=d*(k-(2*Nz-1));
                A_dic[i][j][k] = this->A_dic_generator(x,y,z);
            }
        }
    }
    for (int i=0;i<=Nx-1;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*i; double y=d*(j-(2*Ny-1)); double z=d*(k-(2*Nz-1));
                A_dic[i][j][k] = this->A_dic_generator(x,y,z);
            }
        }
    }
    for (int i=Nx;i<=2*Nx-2;i++) {
        for (int j=Ny;j<=2*Ny-2;j++) {
            for (int k=Nz;k<=2*Nz-2;k++) {
                double x=d*(i-(2*Nx-1)); double y=d*(j-(2*Ny-1)); double z=d*(k-(2*Nz-1));
                A_dic[i][j][k] = this->A_dic_generator(x,y,z);
            }
        }
    }

    high_resolution_clock::time_point t_init = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t_init-t_start).count();

    if (verbose) {
        cout << "--------------Initialization finished. Duration: " << duration/1000.0 << "s.-------------------" << endl;
    }
    fftw_init_threads();                                                       //////////Initialize the multi-thread
    fftw_plan_with_nthreads(omp_get_max_threads());
    VectorXcd p = VectorXcd::Zero(N*3);
    VectorXcd t = VectorXcd::Zero(N*3);
    VectorXcd w = VectorXcd::Zero(N*3);
    VectorXcd r = VectorXcd::Zero(N*3);
    VectorXcd r0 = VectorXcd::Zero(N*3);
    VectorXcd rl = VectorXcd::Zero(N*3);
    VectorXcd y = VectorXcd::Zero(N*3);
    VectorXcd u = VectorXcd::Zero(N*3);
    VectorXcd z = VectorXcd::Zero(N*3);
    VectorXcd x = VectorXcd::Zero(N*3);
    std::complex<double> alpha;
    std::complex<double> beta;
    std::complex<double> eta;
    std::complex<double> zeta;

    VectorXcd Ax0 = this->Aproduct(A_dic,P);

    r = E-Ax0;
    r0 = r;
    p = r;
    VectorXcd Ap0 = this->Aproduct(A_dic,p);


    alpha = r0.dot(r)/r0.dot(Ap0);
    t = r-alpha*Ap0;
    VectorXcd At0 = this->Aproduct(A_dic,t);

    zeta = At0.dot(t)/At0.dot(At0);
    u = zeta*Ap0;
    z = zeta*r-alpha*u;
    P = P+alpha*p+z;                                    //this will directly change P in this.
    rl = r;
    r = t-zeta*At0;

    VectorXcd Ap = VectorXcd::Zero(N*3);
    VectorXcd At = VectorXcd::Zero(N*3);
    for (int it=0;it<=MAX_ITERATION-1;it++) {
        if (verbose && (it+1)%20==0) {
            cout << "                Iter " << it+1 << ", error=" << r.norm()/E.norm() << endl;
        }
        beta = (alpha/zeta)*r0.dot(r)/r0.dot(rl);
        p = r+beta*(p-u);
        Ap = this->Aproduct(A_dic,p);
        alpha = r0.dot(r)/r0.dot(Ap);
        t = r-alpha*Ap;
        At = this->Aproduct(A_dic,t);

        zeta = At.dot(t)/At.dot(At);
        u = zeta*Ap;
        z = zeta*r-alpha*u;
        P = P+alpha*p+z;
        rl = r;
        r = t-zeta*At;
        if (r.norm()/E.norm()<=MAX_ERROR) {
            ERROR = r.norm()/E.norm();
            ITERATION = it+1;
            if (verbose) {
                high_resolution_clock::time_point t_end = high_resolution_clock::now();
                duration = duration_cast<milliseconds>(t_end-t_init).count();
                time = duration_cast<milliseconds>(t_end-t_start).count();
                cout << "--------------Calculation finished. Duration: " << duration/1000.0 << "s.-------------------" << endl;
                cout << "              Error: "<<ERROR<<endl;
                cout << "              Iteration: "<<ITERATION<<endl;
                cout << endl;
            }
            return;
        }
    }
    high_resolution_clock::time_point t_end = high_resolution_clock::now();
    time = duration_cast<milliseconds>(t_end-t_start).count();
    cout<<"                ERROR:does not converge in "<<MAX_ITERATION<<" iterations"<<endl;
    return;
}































VectorXcd Structure::get_diel(){
    return diel;
}

Vector2cd Structure::get_material(){
    return material;
}

VectorXcd Structure::get_alpha(){
    return al;
}