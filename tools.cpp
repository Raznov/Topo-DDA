#include "definition.h"


VectorXi build_a_bulk(int Nx, int Ny, int Nz){
    VectorXi result=VectorXi::Zero(3*Nx*Ny*Nz);
    for(int x=0;x<=Nx-1;x++){
        for(int y=0;y<=Ny-1;y++){
            for(int z=0;z<=Nz-1;z++){
                int position=z+Nz*(y+Ny*x);
                result(3*position)=x;
                result(3*position+1)=y;
                result(3*position+2)=z;
            }
        }
    }   
    return result; 
}


complex<double> Get_material(string mat, double wl, string unit){
    map<string,double> unit_dic;
    unit_dic.insert(pair<string,double>("nm",1.0e9));
    unit_dic.insert(pair<string,double>("um",1.0e6));
    unit_dic.insert(pair<string,double>("m",1.0));
    map<string,string> diel_dic;
    diel_dic.insert(pair<string,string>("Ag","diel/Ag (Silver) - CRC (raw)"));
    diel_dic.insert(pair<string,string>("Al","diel/Al (Aluminium) - Palik (raw)"));
    diel_dic.insert(pair<string,string>("Au","diel/Au (Gold) - Johnson and Christy (raw)"));
    diel_dic.insert(pair<string,string>("Si","diel/Si (Silicon) - Palik (raw)"));
    diel_dic.insert(pair<string,string>("SiO2","diel/SiO2 (Glass) - Palik (raw)"));
    diel_dic.insert(pair<string,string>("Air","diel/Air"));
    diel_dic.insert(pair<string,string>("1.5","diel/Diel1.5"));
    diel_dic.insert(pair<string,string>("2.0","diel/Diel2.0"));
    diel_dic.insert(pair<string,string>("2.5","diel/Diel2.5"));
    diel_dic.insert(pair<string,string>("3.0","diel/Diel3.0"));
    diel_dic.insert(pair<string,string>("3.5","diel/Diel3.5"));
    diel_dic.insert(pair<string,string>("4.0","diel/Diel4.0"));
    diel_dic.insert(pair<string,string>("4.5","diel/Diel4.5"));
    diel_dic.insert(pair<string,string>("5.0","diel/Diel5.0"));
    wl=wl/unit_dic[unit];
    string real="Re_eps.txt";
    string imag="Im_eps.txt";
    string mat_real_name=diel_dic[mat]+real;
    string mat_imag_name=diel_dic[mat]+imag; 
    double mat_real,mat_imag,a,b,up,down,up_value,down_value;
    
    up=1.0;
    down=0.0;

    ifstream mat_real_file;
    mat_real_file.open(mat_real_name);
    while(mat_real_file>>a>>b){
        if(a<=wl&&a>down){
            down=a;
            down_value=b;
        }
        else if(a>wl&&a<up){
            up=a;
            up_value=b;
        }

    }
    mat_real_file.close();
    if(up==1.0||down==0.0){
        cout<<"ERROR, wavelength out of range of the diel file provided."<<endl;
    }
    else{
        mat_real=down_value+(wl-down)*(up_value-down_value)/(up-down); 
        cout << "real:" << mat_real << endl;
    }
    up=1.0;
    down=0.0;

    ifstream mat_imag_file;
    mat_imag_file.open(mat_imag_name);
    while(mat_imag_file>>a>>b){
        if(a<=wl&&a>down){
            down=a;
            down_value=b;
        }
        else if(a>wl&&a<up){
            up=a;
            up_value=b;
        }

    }
    mat_imag_file.close();
    if(up==1.0||down==0.0){
        cout<<"ERROR, wavelength out of range of the diel file provided."<<endl;
    }
    else{
        mat_imag=down_value+(wl-down)*(up_value-down_value)/(up-down); 
        cout << "imag:" << mat_imag << endl;
    }
    
    complex<double> result=mat_real+1.0i*mat_imag;
    return result;
}

Vector2cd Get_2_material(string sub, string mat, double wl, string unit){
    Vector2cd result;
    result(0)=Get_material(sub,wl,unit);
    result(1)=Get_material(mat,wl,unit);
    return result;
}

complex<double> Get_Alpha(double lam, double K, double d, complex<double> diel, Vector3d n_E0, Vector3d n_K){
    double b1 = -1.891531;
    b1 = -1.891531;
    double b2 = 0.1648469;
    double b3 = -1.7700004;
    //cout<<"lam"<<lam<<"K"<<K<<"d"<<d<<endl;
    std::complex<double> a1=(3*pow(d,3)/(4*M_PI))*(diel-1.0)/(diel+2.0);
    //cout<<a1<<endl;
    std::complex<double> result=0.0+(2.0/3.0)*pow(K*d,3)*1.0i;
    double S = pow(n_E0(0) * n_K(0), 2) + pow(n_E0(1) * n_K(1), 2) + pow(n_E0(2) * n_K(2), 2);
    result = 1.0 + (a1 / pow(d, 3)) * ((b1 + diel * b2 + diel * b3 * S) * pow(K * d, 2) - result);
    //cout<<result<<endl;
    result=a1/result;
    return result;
}

/*
ArrayXcd FFT(int nx,int ny,int nz,ArrayXcd in,int _direction){
    //0 for forward, 1 for backward
    //FFTW_FORWARD and FFTW_BACKWARD are FFTW defined integer
    //normalize is needed because backward needs to be normalized to be the same as original

    int direction[2]={FFTW_FORWARD,FFTW_BACKWARD};
    int N=nx*ny*nz;
    int normalize[2]={1,N};
    fftw_complex* in1;
    fftw_complex* out;
    in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for(int i=0;i<=N-1;i++){
        in1[i][0]=in(i).real();
        in1[i][1]=in(i).imag();
    }
    fftw_plan p1;
    p1=fftw_plan_dft_3d(nx,ny,nz,in1,out,direction[_direction],FFTW_ESTIMATE);
    fftw_execute(p1);
    fftw_destroy_plan(p1);
    ArrayXcd result(N);
    for(int i=0;i<=N-1;i++){
        result(i)=out[i][0]+out[i][1]*1.0i;
    }
    result=result/normalize[_direction];
    fftw_free(in1);
    fftw_free(out);

    return result;

}
*/