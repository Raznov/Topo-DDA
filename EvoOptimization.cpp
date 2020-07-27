#include "definition.h"

void EvoOptimization(Model *model, double epsilon, Vector3d r, int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO){
    
    string save_position="./222210-111105-25-808-01-tight/";
    ofstream convergence;
    convergence.open(save_position+"convergence.txt");
    
    double epsilon_partial=0.001;
    for(int iteration=0;iteration<=MAX_ITERATION_EVO-1;iteration++){
        //solve DDA
        cout<<"###EVO ITERATION "<<iteration<<endl;
        //get object function value
        cout<<"###START ORIGINAL PROBLEM"<<endl;
        
        (*model).bicgstab(MAX_ITERATION, MAX_ERROR);  
        (*model).solve_E(); 
        (*model).output_to_file(save_position+"Model_output/", iteration);
        double obj=(*model).point_E(r);
        convergence << obj << "\n";
        cout<<"objective function at iteration "<<iteration<<" is "<<obj<<endl;
        //get partial derivative of current model
        high_resolution_clock::time_point t1 = high_resolution_clock::now();     
        cout<<"###START PARTIAL DERIVATIVE"<<endl;
        VectorXd devx;
        VectorXcd Adevxp;
        VectorXcd devp;
        tie(devx, Adevxp)=(*model).devx_and_Adevxp(epsilon_partial, r);
        devp=(*model).devp(epsilon_partial, r);
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout<<"------------------------finished in "<<duration/1000<<" s"<<endl;

        //Solving adjoint problem
        cout<<"###START ADJOINT PROBLEM"<<endl;
        (*model).change_E(devp);
        (*model).bicgstab(MAX_ITERATION, MAX_ERROR);
        VectorXcd lambdaT=(*model).get_P();
        (*model).reset_E();                                  //reset E to initial value



        //times lambdaT and Adevxp together
        VectorXcd mult_result;
        list<int> para_nums, para_starts, para_dep_nums, para_dep_starts;
        tie(para_nums, para_starts, para_dep_nums, para_dep_starts)=(*model).get_para_info();
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
        /*
        cout<<"devp1-3"<<devp(1)<<endl<<devp(2)<<endl<<devp(3)<<endl;
        cout<<"lambdaT1-3"<<endl<<lambdaT(1)<<endl<<lambdaT(2)<<endl<<lambdaT(3)<<endl;
        cout<<"Adevxp1-3"<<endl<<Adevxp(1)<<endl<<Adevxp(2)<<endl<<Adevxp(3)<<endl;
        cout<<"devx"<<endl<<devx(1)<<endl<<devx(2)<<endl<<devx(3)<<endl;
        cout<<"mult_result_real"<<endl<<mult_result_real(1)<<endl<<mult_result_real(2)<<endl<<mult_result_real(3)<<endl;
        */
        gradients=devx-mult_result_real;              //What's the legitimacy in here to ignore the imag part?
        //cout<<"gradients1-3"<<endl<<gradients(0)<<endl<<gradients(1)<<endl<<gradients(2)<<endl;
        double step_len = (*model).get_step_length(gradients,epsilon);
        VectorXd step=epsilon*gradients;               //Find the maximum. If -1 find minimum
        cout << "epsilon = " << step_len << endl;
        //cout<<"devp1"<<devp<<endl;
        //cout<<"lambdaT1"<<endl<<lambdaT<<endl;
        //cout<<"Adevxp1"<<endl<<Adevxp<<endl;
        //cout<<"devx"<<endl<<devx<<endl;
        //cout<<"mult_result_real"<<endl<<mult_result_real<<endl;
        //cout<<"gradients"<<endl<<gradients<<endl;
        //cout<<"step"<<endl<<step<<endl;

           
        (*model).change_para_diel(step);
        
    }
    convergence.close();


}