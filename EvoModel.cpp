#include "definition.h"

EvoModel::EvoModel(list<string> *ObjectFunctionNames_, list<list<double>*> *ObjectParameters_, bool HavePenalty_, double PenaltyFactor_, string save_position_, Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_) : Model(space_, d_, lam_,  n_K_,  E0_,  n_E0_,  material_){
    ObjectFunctionNames = ObjectFunctionNames_;
    save_position = save_position_;
    ObjectParameters = ObjectParameters_;
    HavePenalty = HavePenalty_;
    PenaltyFactor = PenaltyFactor_;
    origin = 0;
    
    

    list<string>::iterator it0 = (*ObjectFunctionNames).begin();
    MajorObjectFunctionName = (*it0);
    it0++;
    
    for(int i = 1; i <= (*ObjectFunctionNames).size()-1; i++){
        MinorObjectFunctionNames.push_back(*it0);
        it0++;
        
    }
    
    list<list<double>*>::iterator it1 = (*ObjectParameters).begin();    
    for(int i = 0; i <= (*ObjectParameters).size()-1; i ++){
        if(i == 0){
            MajorObjectParameters = *(*it1);
        }
        else{
            MinorObjectParameters.push_back(*(*it1));
        }
        it1++;
    }

    list<double>::iterator it2 = MajorObjectParameters.begin();
    for(int i=0; i<= MajorObjectParameters.size()-1; i++){
        cout<<" "<<*it2<<" ";
        it2++;
    }
    cout<<endl;
    list<list<double>>::iterator it3 = MinorObjectParameters.begin();
    
    for(int i = 1; i<= (*ObjectParameters).size()-1; i++){
        
        list<double>::iterator it4 = (*it3).begin();
        for(int j = 0; j<= (*it3).size()-1; j++){
            cout<<" "<<*it4<<" ";
            it4++;
        }
        it3++;
        cout<<endl;
    }
    objective = ObjectiveFactory(MajorObjectFunctionName, MajorObjectParameters);
    
}
EvoModel::EvoModel(list<string> *ObjectFunctionNames_, list<list<double>*> *ObjectParameters_, bool HavePenalty_, double PenaltyFactor_, string save_position_, Space *space_, double d_, double lam_, Vector3d n_K_, double E0_, Vector3d n_E0_, Vector2cd material_, VectorXi *RResult_) : Model(space_, d_, lam_,  n_K_,  E0_,  n_E0_,  material_, RResult_){
    ObjectFunctionNames = ObjectFunctionNames_;
    save_position = save_position_;
    ObjectParameters = ObjectParameters_;
    HavePenalty = HavePenalty_;
    PenaltyFactor = PenaltyFactor_;
    origin = 0;

    list<string>::iterator it0 = (*ObjectFunctionNames).begin();
    MajorObjectFunctionName = (*it0);
    it0++;
    
    for(int i = 1; i <= (*ObjectFunctionNames).size()-1; i++){
        MinorObjectFunctionNames.push_back(*it0);
        it0++;
        
    }
    
    list<list<double>*>::iterator it1 = (*ObjectParameters).begin();    
    for(int i = 0; i <= (*ObjectParameters).size()-1; i ++){
        if(i == 0){
            MajorObjectParameters = *(*it1);
        }
        else{
            MinorObjectParameters.push_back(*(*it1));
        }
        it1++;
    }

    list<double>::iterator it2 = MajorObjectParameters.begin();
    for(int i=0; i<= MajorObjectParameters.size()-1; i++){
        cout<<" "<<*it2<<" ";
        it2++;
    }
    cout<<endl;
    list<list<double>>::iterator it3 = MinorObjectParameters.begin();
    
    for(int i = 1; i<= (*ObjectParameters).size()-1; i++){
        cout<<"FUCK"<<endl;
        list<double>::iterator it4 = (*it3).begin();
        for(int j = 0; j<= (*it3).size()-1; j++){
            cout<<" "<<*it4<<" ";
            it4++;
        }
        it3++;
        cout<<endl;
    }
    objective = ObjectiveFactory(MajorObjectFunctionName, MajorObjectParameters);
    if(HavePenalty) {
        cout << "Currently using L1 penalty. Penalty factor C is " << PenaltyFactor << endl;
    }
    
}



tuple<VectorXd, VectorXcd> EvoModel::devx_and_Adevxp(double epsilon){
    double origin=objective->GetVal();
    VectorXcd Adevxp=VectorXcd::Zero(3*N);

    int para_size=para_nums.size();
    int para_dep_size=para_dep_nums.size();
    
    int n_para=0;
    list<int>::iterator it=para_nums.begin();
    for(int i=0;i<=para_size-1;i++){
        n_para=n_para+(*it);
        it++;
    }
    n_para=round(n_para/3);
    VectorXd devx=VectorXd::Zero(n_para);

    if(para_dep_size!=0){
        if(para_dep_size!=para_size){
            cout<<"ERROR: para_dep_size not equal para_size"<<endl;
        }
        list<int>::iterator it1=para_nums.begin();
        list<int>::iterator it2=para_starts.begin();
        list<int>::iterator it3=para_dep_nums.begin();
        list<int>::iterator it4=para_dep_starts.begin();
        int position=0;
        for(int i=0;i<=para_size-1;i++){
            int times=round((*it3)/(*it1));
            int para_begin=round((*it2)/3);
            int para_number=round((*it1)/3);
            int para_dep_begin=round((*it4)/3);
            for(int j=0;j<=para_number-1;j++){
                int sign=0;
                if(diel_old(j)>=epsilon){
                    sign=-1;
                }
                else{
                    sign=1;
                }
                int position1=(j+para_begin);
                
                diel_old(3*position1)+=sign*epsilon;
                diel_old(3*position1+1)+=sign*epsilon;
                diel_old(3*position1+2)+=sign*epsilon;
                
                
                if(objective->Have_Devx) objective->SingleResponse(position1, true);
                
                diel(3*position1)=material(0)+diel_old(3*position1)*(material(1)-material(0));
                diel(3*position1+1)=diel(3*position1);
                diel(3*position1+2)=diel(3*position1);
                
                if(objective->Have_Devx) objective->SingleResponse(position1, false);
                
                Adevxp(3*position1)=((1.0/Get_Alpha(lam,K,d,diel(3*position1)))-al(3*position1))/(sign*epsilon);
                Adevxp(3*position1+1)=((1.0/Get_Alpha(lam,K,d,diel(3*position1+1)))-al(3*position1+1))/(sign*epsilon);
                Adevxp(3*position1+2)=((1.0/Get_Alpha(lam,K,d,diel(3*position1+2)))-al(3*position1+2))/(sign*epsilon);

                for(int k=0;k<=times-1;k++){
                    int position2=j+para_dep_begin+k*para_number;
                    
                    diel_old(3*position2)+=sign*epsilon;
                    diel_old(3*position2+1)+=sign*epsilon;
                    diel_old(3*position2+2)+=sign*epsilon;
                    
                    if(objective->Have_Devx) objective->SingleResponse(position2, true);
                    
                    diel(3*position2)=material(0)+diel_old(3*position2)*(material(1)-material(0));
                    diel(3*position2+1)=diel(3*position2);
                    diel(3*position2+2)=diel(3*position2);
                    
                    if(objective->Have_Devx) objective->SingleResponse(position2, false);
                    
                    Adevxp(3*position2)=((1.0/Get_Alpha(lam,K,d,diel(3*position2)))-al(3*position2))/(sign*epsilon);
                    Adevxp(3*position2+1)=((1.0/Get_Alpha(lam,K,d,diel(3*position2+1)))-al(3*position2+1))/(sign*epsilon);
                    Adevxp(3*position2+2)=((1.0/Get_Alpha(lam,K,d,diel(3*position2+2)))-al(3*position2+2))/(sign*epsilon);
                }
                devx(position)=(objective->GroupResponse()-origin)/(sign*epsilon);
                
                position=position+1;

                diel_old(3*position1)-=sign*epsilon;
                diel_old(3*position1+1)-=sign*epsilon;
                diel_old(3*position1+2)-=sign*epsilon;
                
                if(objective->Have_Devx) objective->SingleResponse(position1, true);
                
                diel(3*position1)=material(0)+diel_old(3*position1)*(material(1)-material(0));
                diel(3*position1+1)=diel(3*position1);
                diel(3*position1+2)=diel(3*position1);
                
                if(objective->Have_Devx) objective->SingleResponse(position1, false);

                for(int k=0;k<=times-1;k++){
                    int position2=j+para_dep_begin+k*para_number;
                    diel_old(3*position2)-=sign*epsilon;
                    diel_old(3*position2+1)-=sign*epsilon;
                    diel_old(3*position2+2)-=sign*epsilon;
                    
                    if(objective->Have_Devx) objective->SingleResponse(position2, true);
                    
                    diel(3*position2)=material(0)+diel_old(3*position2)*(material(1)-material(0));
                    diel(3*position2+1)=diel(3*position2);
                    diel(3*position2+2)=diel(3*position2);
                    
                    if(objective->Have_Devx) objective->SingleResponse(position2, false);
                    
                }
            }
            it1++;
            it2++;
            it3++;
            it4++;
        }
        for(int i=0;i<=3*N-1;i++){
            Adevxp(i)=Adevxp(i)*P(i);
        }

        return make_tuple(devx, Adevxp);
    }
    else{
        list<int>::iterator it1=para_nums.begin();
        list<int>::iterator it2=para_starts.begin();
        int position=0;
        for(int i=0;i<=para_size-1;i++){
            int para_begin=round((*it2)/3);
            int para_number=round((*it1)/3);
            for(int j=0;j<=para_number-1;j++){
                int sign=0;
                if(diel_old(j)>=epsilon){
                    sign=-1;
                }
                else{
                    sign=1;
                }
                int position1=(j+para_begin);
                diel_old(3*position1)+=sign*epsilon;
                diel_old(3*position1+1)+=sign*epsilon;
                diel_old(3*position1+2)+=sign*epsilon;
                
                if(objective->Have_Devx) objective->SingleResponse(position1, true);
                
                diel(3*position1)=material(0)+diel_old(3*position1)*(material(1)-material(0));
                //cout<<"material0"<<material(0)<<"material1"<<material(1)<<endl;
                diel(3*position1+1)=diel(3*position1);
                diel(3*position1+2)=diel(3*position1);
                
                if(objective->Have_Devx) objective->SingleResponse(position1, false);
                
                Adevxp(3*position1)=((1.0/Get_Alpha(lam,K,d,diel(3*position1)))-al(3*position1))/(sign*epsilon);
                Adevxp(3*position1+1)=((1.0/Get_Alpha(lam,K,d,diel(3*position1+1)))-al(3*position1+1))/(sign*epsilon);
                Adevxp(3*position1+2)=((1.0/Get_Alpha(lam,K,d,diel(3*position1+2)))-al(3*position1+2))/(sign*epsilon);
                //cout<<"diel"<<diel(3*position1)<<endl;
                //cout<<"lam"<<lam<<"K"<<K<<"d"<<d<<endl;
                //cout<<"1/alpha"<<(1.0/Get_Alpha(lam,K,d,diel(3*position1)))<<endl;
                devx(position)=(objective->GroupResponse()-origin)/(sign*epsilon);
                position=position+1;

                diel_old(3*position1)-=sign*epsilon;
                diel_old(3*position1+1)-=sign*epsilon;
                diel_old(3*position1+2)-=sign*epsilon;
                
                if(objective->Have_Devx) objective->SingleResponse(position1, true);
                
                diel(3*position1)=material(0)+diel_old(3*position1)*(material(1)-material(0));
                diel(3*position1+1)=diel(3*position1);
                diel(3*position1+2)=diel(3*position1);
                
                if(objective->Have_Devx) objective->SingleResponse(position1, false);
            }
            it1++;
            it2++;
        }
        
    }
    for(int i=0;i<=3*N-1;i++){
            //cout<<" i: "<<i<<" Adevx: "<<Adevxp(i)<<" P: "<<P(i)<<endl;
            Adevxp(i)=Adevxp(i)*P(i);
        }
    return make_tuple(devx, Adevxp);
}

VectorXcd EvoModel::devp(double epsilon){
    //move origin=objective0->GetVal() outside because it is the same for one partial derivative of the entire structure
    VectorXcd result=VectorXcd::Zero(P.size());
    for(int i=0;i<=P.size()-1;i++){
        int position = i/3;
        
        objective->SingleResponse(position, true);
        
        P(i)=P(i)+epsilon;
        
        objective->SingleResponse(position, false);
        
        result(i)+=(objective->GroupResponse()-origin)/epsilon;
        
        objective->SingleResponse(position, true);
        
        P(i)=P(i)-epsilon;
        
        complex<double> tmp=epsilon*1.0i;
        
        P(i)=P(i)+tmp;
        
        objective->SingleResponse(position, false);
        
        complex<double> tmpRes = (objective->GroupResponse()-origin)/tmp;
        result(i)+=tmpRes;
        
        objective->SingleResponse(position, true);
        
        P(i)=P(i)-tmp;
        
        objective->SingleResponse(position, false);
    }
    cout << "Devp: " << result.sum() << endl;
    return result;
}



void EvoModel::EvoOptimization(double epsilon, int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method){
    ofstream convergence;
    convergence.open(save_position+"convergence.txt");
    
    
    //Parameters for Adam Optimizer.
    double beta1 = 0.9;
    double beta2 = 0.99;
    VectorXd V;
    VectorXd S;
    
    
    double epsilon_partial=0.001;
    for(int iteration=0;iteration<=MAX_ITERATION_EVO-1;iteration++){
        //solve DDA
        cout<<"###EVO ITERATION "<<iteration<<endl;
        //get object function value
        cout<<"###START ORIGINAL PROBLEM"<<endl;
        
        this->bicgstab(MAX_ITERATION, MAX_ERROR);  
        //this->solve_E(); 

        this->update_E_in_structure();
        if(iteration==MAX_ITERATION_EVO-1){                                    //useless fix, not gonna to use RResultswithc = true feature in the future
            this->solve_E(); 
        }
        
        this->output_to_file(save_position+"Model_output/", iteration);


        high_resolution_clock::time_point t0 = high_resolution_clock::now();
        
        
        double obj = objective->GetVal();
        

        convergence << obj << " ";
        cout<<"objective function at iteration "<<iteration<<" is "<<obj<<endl;
              
        if((*ObjectFunctionNames).size()>1){
            list<double> obj_minor =  this->MinorObjective();
            list<double>::iterator it_obj_minor = obj_minor.begin();
            for(int i=0; i<=obj_minor.size()-1; i++){
                convergence << *it_obj_minor << " ";
                it_obj_minor++;
            }
        }
        
        convergence << "\n";

        
        
        //get partial derivative of current model
        high_resolution_clock::time_point t1 = high_resolution_clock::now();    

        cout<<"-------------------------------Time consumption of one obj with no pre-stored A"<<duration_cast<milliseconds>(t1-t0).count();

        cout<<"###START PARTIAL DERIVATIVE"<<endl;
        VectorXd devx;
        VectorXcd Adevxp;
        VectorXcd devp;

        tie(devx, Adevxp)=this->devx_and_Adevxp(epsilon_partial);
        
        origin=obj;                              //Origin equals to the objective function of current structure
        devp=this->devp(epsilon_partial);
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout<<"------------------------finished in "<<duration/1000<<" s"<<endl;

        //Solving adjoint problem
        cout<<"###START ADJOINT PROBLEM"<<endl;
        this->change_E(devp);
        this->bicgstab(MAX_ITERATION, MAX_ERROR);
        VectorXcd lambdaT=P;
        this->reset_E();                                  //reset E to initial value



        //times lambdaT and Adevxp together
        VectorXcd mult_result;
        list<int> para_nums, para_starts, para_dep_nums, para_dep_starts;
        tie(para_nums, para_starts, para_dep_nums, para_dep_starts)=this->get_para_info();
        int para_size=para_nums.size();
        int para_dep_size=para_dep_nums.size();
        int n_para=0;
        list<int>::iterator it=para_nums.begin();
        for(int i=0;i<=para_size-1;i++){
            n_para=n_para+(*it);
            it++;
        }
        n_para=round(n_para/3);                            //Total number of parameters
        mult_result=VectorXcd::Zero(n_para);               //multiplication result has the length of parameter
        if(para_dep_size!=0){
            if(para_dep_size!=para_size){
                cout<<"ERROR: para_dep_size not equal para_size"<<endl;
            }
            list<int>::iterator it1=para_nums.begin();
            list<int>::iterator it2=para_starts.begin();
            list<int>::iterator it3=para_dep_nums.begin();
            list<int>::iterator it4=para_dep_starts.begin();
            int position=0;
            for(int i=0;i<=para_size-1;i++){
                int times=round((*it3)/(*it1));
                int para_begin=round((*it2)/3);
                int para_number=round((*it1)/3);
                int para_dep_begin=round((*it4)/3);
                for(int j=0;j<=para_number-1;j++){
                    int position1=(j+para_begin);
                    mult_result(position)+=lambdaT(3*position1)*Adevxp(3*position1);
                    mult_result(position)+=lambdaT(3*position1+1)*Adevxp(3*position1+1);
                    mult_result(position)+=lambdaT(3*position1+2)*Adevxp(3*position1+2);
                    for(int k=0;k<=times-1;k++){
                        int position2=j+para_dep_begin+k*para_number;
                        mult_result(position)+=lambdaT(3*position2)*Adevxp(3*position2);
                        mult_result(position)+=lambdaT(3*position2+1)*Adevxp(3*position2+1);
                        mult_result(position)+=lambdaT(3*position2+2)*Adevxp(3*position2+2);
                    }
                    position=position+1;
                }
                it1++;
                it2++;
                it3++;
                it4++;
            }
        }
        else{
            list<int>::iterator it1=para_nums.begin();
            list<int>::iterator it2=para_starts.begin();
            int position=0;
            for(int i=0;i<=para_size-1;i++){
                int para_begin=round((*it2)/3);
                int para_number=round((*it1)/3);
                for(int j=0;j<=para_number-1;j++){
                    int position1=(j+para_begin);
                    mult_result(position)+=lambdaT(3*position1)*Adevxp(3*position1);
                    mult_result(position)+=lambdaT(3*position1+1)*Adevxp(3*position1+1);
                    mult_result(position)+=lambdaT(3*position1+2)*Adevxp(3*position1+2);
                    position=position+1;
                }
                it1++;
                it2++;
            }
        }
    

        //The final gradients for this iteration

        VectorXd gradients=VectorXd::Zero(n_para);
        VectorXd mult_result_real=VectorXd::Zero(n_para);
        for(int i=0;i<=n_para-1;i++){
            complex<double> tmp=mult_result(i);
            mult_result_real(i)=tmp.real();
        }

      
            
        gradients=devx-mult_result_real;              //What's the legitimacy in here to ignore the imag part?
        
        if(method == "Adam"){
            cout << "#####################Using Adam Optimizer.###################" << endl;
            if(iteration == 0){
                V = (1-beta1)*gradients/(1-pow(beta1,iteration+1));
                S = (1-beta2)*(gradients.array().pow(2).matrix())/(1-pow(beta2,iteration+1));
            }
            else{
                V = beta1*V + (1-beta1)*gradients/(1-pow(beta1,iteration+1));
                S = beta2*S + (1-beta2)*(gradients.array().pow(2).matrix())/(1-pow(beta2,iteration+1));          
            }
            for(int i=0;i<=n_para-1;i++){
                gradients(i) = V(i)/(sqrt(S(i))+0.00000001);
            } 
        } 
        cout << "gradients: " << gradients.mean() << endl;
        
        //cout<<"gradients1-3"<<endl<<gradients(0)<<endl<<gradients(1)<<endl<<gradients(2)<<endl;
        //double step_len = this->get_step_length(gradients,epsilon);
        VectorXd step=epsilon*gradients;               //Find the maximum. If -1 find minimum
        cout << "epsilon = " << epsilon << endl;
        //cout<<"devp1"<<devp<<endl;
        //cout<<"lambdaT1"<<endl<<lambdaT<<endl;
        //cout<<"Adevxp1"<<endl<<Adevxp<<endl;
        //cout<<"devx"<<endl<<devx<<endl;
        //cout<<"mult_result_real"<<endl<<mult_result_real<<endl;
        //cout<<"gradients"<<endl<<gradients<<endl;
        //cout<<"step"<<endl<<step<<endl;

           
        this->change_para_diel(step);
        
    }
    
    convergence.close();

}

double EvoModel::MajorObjective(){
    if(MajorObjectFunctionName == "PointE"){
        return this -> PointEWithPenalty(MajorObjectParameters);
    }
    else if(MajorObjectFunctionName == "SurfaceEExp"){
        return this -> SurfaceEExpWithPenalty(MajorObjectParameters);
    }
    else{
        cout<<"ERROR: NO suitable major objective function"<<endl;
        return 0.0;
    }
}

Objective* EvoModel::ObjectiveFactory(string ObjectName, list<double> ObjectParameters){
    if (HavePenalty) {
        cout << "Using L1 Penalty with Penalty Factor " << PenaltyFactor << endl;
    }
    if (MajorObjectFunctionName == "PointE"){
        return new ObjectivePointE(ObjectParameters, this, HavePenalty);
    }
    if (MajorObjectFunctionName == "SurfaceEExp"){
        return new ObjectiveSurfaceEExp(ObjectParameters, this, HavePenalty);
    }
    if (MajorObjectFunctionName == "ExtSurfaceEExp"){
        return new ObjectiveExtSurfaceEExp(ObjectParameters, this, HavePenalty);
    }
    else{
        // NOT FINALIZED. SHOULD RAISE AN EXCEPTION HERE.
        cout << "NOT A LEGIT OBJECTIVE NAME!" << endl;
        return new ObjectivePointE(ObjectParameters, this, HavePenalty);
    }

}

list<double> EvoModel::MinorObjective(){
    list<double> TemporalMinorObjectResults;
    list<string>::iterator it = MinorObjectFunctionNames.begin();
    list<list<double>>::iterator it1 = MinorObjectParameters.begin();
    for(int i = 0; i <= MinorObjectFunctionNames.size()-1; i++){
        if((*it) == "PointE"){
            TemporalMinorObjectResults.push_back(this -> PointE(*it1));
        }
        else if((*it) == "SurfaceEExp"){
            TemporalMinorObjectResults.push_back(this -> SurfaceEExp(*it1));
        }
        else{
            cout<<"Can not find Objective function: "<<(*it)<<endl;
            TemporalMinorObjectResults.push_back(0.0);
        }

        it++;
        it1++;

    }

    return TemporalMinorObjectResults;

}

double EvoModel::PointE(list<double> Parameter){
    
   
    VectorXd PointEParameters = VectorXd::Zero((Parameter).size());
    list<double>::iterator it=(Parameter).begin();
    for(int i=0;i<=int((Parameter).size()-1);i++){
        PointEParameters(i) = (*it);
        it++;
    }

    Vector3cd sum=Vector3cd::Zero();
    Vector3cd E_ext=Vector3cd::Zero();
    double x,y,z;
    x=PointEParameters(0);
    y=PointEParameters(1);
    z=PointEParameters(2);
    E_ext(0) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(1) = E0*n_E0(1)*(cos(K*(n_K(0)*y+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
    E_ext(2) = E0*n_E0(2)*(cos(K*(n_K(0)*z+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);   
    cout << E_ext(0) << endl;                                         
    for (int i=0;i<N;i++){
        double rx=x-d*R(3*i);                  //R has no d in it, so needs to time d
        double ry=y-d*R(3*i+1);
        double rz=z-d*R(3*i+2);
        //cout << rx << "," << ry << "," << rz << endl;
        Matrix3cd A=this->A_dic_generator(rx,ry,rz);
        sum(0)+=(A(0,0)*P(3*i)+A(0,1)*P(3*i+1)+A(0,2)*P(3*i+2));
        sum(1)+=(A(1,0)*P(3*i)+A(1,1)*P(3*i+1)+A(1,2)*P(3*i+2));
        sum(2)+=(A(2,0)*P(3*i)+A(2,1)*P(3*i+1)+A(2,2)*P(3*i+2));
    }
    sum(0)=E_ext(0)-sum(0);
    sum(1)=E_ext(1)-sum(1);
    sum(2)=E_ext(2)-sum(2);
    return (sum).norm();
    
}


double EvoModel::SurfaceEExp(list<double> Parameter){
    
    VectorXd SurfaceEExpParameters = VectorXd::Zero((Parameter).size());
    list<double>::iterator it=(Parameter).begin();
    for(int i=0;i<=int((Parameter).size()-1);i++){
        SurfaceEExpParameters(i) = (*it);
        it++;
    }
    double SurfaceEExpRx = SurfaceEExpParameters(0);
    double SurfaceEExpRy = SurfaceEExpParameters(1);
    double SurfaceEExpRz = SurfaceEExpParameters(2);
    double SurfaceEExpL1 = SurfaceEExpParameters(3);
    double SurfaceEExpL2 = SurfaceEExpParameters(4);
    double SurfaceEExpD = SurfaceEExpParameters(5);
    
    int SurfaceEExpNx = int(SurfaceEExpL1/SurfaceEExpD);
    int SurfaceEExpNy = int(SurfaceEExpL2/SurfaceEExpD);
    double x,y,z;
    x=SurfaceEExpRx-SurfaceEExpL1/2;
    y=SurfaceEExpRy-SurfaceEExpL2/2;
    z=SurfaceEExpRz;
    double SurfaceEExpResult = 0.0;

    while(x <= SurfaceEExpRx + SurfaceEExpL1/2){
        while(y <= SurfaceEExpRy + SurfaceEExpL2/2){
            Vector3cd sum=Vector3cd::Zero();
            Vector3cd E_ext=Vector3cd::Zero();
            E_ext(0) = E0*n_E0(0)*(cos(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            E_ext(1) = E0*n_E0(1)*(cos(K*(n_K(0)*y+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);
            E_ext(2) = E0*n_E0(2)*(cos(K*(n_K(0)*z+n_K(1)*y+n_K(2)*z))+sin(K*(n_K(0)*x+n_K(1)*y+n_K(2)*z))*1i);                                            
            for (int i=0;i<N;i++){
                double rx=x-d*R(3*i);                  //R has no d in it, so needs to time d
                double ry=y-d*R(3*i+1);
                double rz=z-d*R(3*i+2);
                Matrix3cd A=this->A_dic_generator(rx,ry,rz);
                sum(0)+=(A(0,0)*P(3*i)+A(0,1)*P(3*i+1)+A(0,2)*P(3*i+2));
                sum(1)+=(A(1,0)*P(3*i)+A(1,1)*P(3*i+1)+A(1,2)*P(3*i+2));
                sum(2)+=(A(2,0)*P(3*i)+A(2,1)*P(3*i+1)+A(2,2)*P(3*i+2));
            }
            //cout<<x<<endl<<y<<endl<<z<<endl;
            sum(0)=E_ext(0)-sum(0);
            sum(1)=E_ext(1)-sum(1);
            sum(2)=E_ext(2)-sum(2);
            
            SurfaceEExpResult = SurfaceEExpResult + exp(((sum).norm())*((sum).norm()));
            x=x + SurfaceEExpD;
            y=y + SurfaceEExpD;
        }
    }
    SurfaceEExpResult = SurfaceEExpResult/(SurfaceEExpNx*SurfaceEExpNy);
    return SurfaceEExpResult;
}

double EvoModel::PointEWithPenalty(list<double> Parameter){
    if(!HavePenalty){
        return this -> PointE(MajorObjectParameters);
    }
    else{
        return this -> PointE(MajorObjectParameters) + L1Norm();
    }
}

double EvoModel::SurfaceEExpWithPenalty(list<double> Parameter){
    if(!HavePenalty){
        return this -> SurfaceEExp(MajorObjectParameters);
    }
    else{
        return this -> SurfaceEExp(MajorObjectParameters) + L1Norm();
    }
}

double EvoModel::L1Norm(){
    double Penalty = 0;
    for (int i=0;i<N;i++){
        Penalty += 0.5-abs(diel_old(3*i)-0.5);
    }
    Penalty = Penalty * PenaltyFactor;
    return Penalty;
}



