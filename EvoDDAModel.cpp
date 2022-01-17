#include "definition.h"



EvoDDAModel::EvoDDAModel(list<string>* ObjectFunctionNames_, list<list<double>*>* ObjectParameters_, double epsilon_fix_, bool HavePathRecord_, bool HavePenalty_, bool HaveOriginHeritage_, bool HaveAdjointHeritage_, double PenaltyFactor_, string save_position_, CoreStructure* CStr_, list<DDAModel*> ModelList_){
    output_time = 0.0;
    ObjectFunctionNames = ObjectFunctionNames_;
    save_position = save_position_;
    ObjectParameters = ObjectParameters_;
    HavePenalty = HavePenalty_;
    PenaltyFactor = PenaltyFactor_;
    epsilon_fix = epsilon_fix_;
    epsilon_tmp = epsilon_fix;
    HavePathRecord = HavePathRecord_;
    HaveOriginHeritage = HaveOriginHeritage_;
    HaveAdjointHeritage = HaveAdjointHeritage_;
    MaxObj = 0.0;
    Stephold = 0;
    CStr = CStr_;
    ModelList = ModelList_;
    ModelNum = ModelList.size();
    cout << "ModelNum" << ModelNum << endl;
    MaxObjarray = VectorXd::Zero(ModelNum);
    Originarray = VectorXd::Zero(ModelNum);
    PreviousObj = 0.0;
    CutoffHold = 0;

    SpacePara* spacepara = (*CStr).get_spacepara();
    VectorXd* Para = (*spacepara).get_Para();
    int n_para = (*Para).size();                     //Total number of parameters
    gradientsquare = VectorXd::Zero(n_para);


    list<DDAModel*>::iterator it_ModelList = ModelList.begin();
    for (int i = 0; i <= ModelNum - 1; i++) {
        if ((*(*(*(it_ModelList))).get_Core()).get_CStr() != CStr) {
            cout << "The DDAModel number: " << i << " does not share the same CoreStructure" << endl;
        }
        it_ModelList++;
    }

    int N = (*CStr).get_N();
    for (int i = 0; i <= ModelNum - 1; i++) {
        VectorXcd Ptmp = VectorXcd::Zero(N*3);
        PforOrigin.push_back(Ptmp);
        PforAdjoint.push_back(Ptmp);
        PforOriginMax.push_back(Ptmp);
        PforAdjointMax.push_back(Ptmp);

    }

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

    //-----------------generate obj list----------------------------
    it_ModelList = ModelList.begin();
    for (int i = 0; i <= ModelNum - 1; i++) {
        double origin = 0.0;
        ObjectiveDDAModel* objective = ObjectiveFactory(MajorObjectFunctionName, MajorObjectParameters, *it_ModelList);
        ObjList.push_back(objective);
        it_ModelList++;
    }
    

}

tuple<VectorXd, VectorXcd> EvoDDAModel::devx_and_Adevxp_tmp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin){
    int N = (*CurrentModel).get_N();
    SpacePara* spacepara = (*CurrentModel).get_spacepara();
    VectorXi* geometryPara = (*spacepara).get_geometryPara();
    VectorXd* Para = (*spacepara).get_Para();
    VectorXi* Free = (*spacepara).get_Free();

    VectorXd* diel_old = (*CurrentModel).get_diel_old();
    VectorXcd* material = (*CurrentModel).get_material();
    double lam = (*CurrentModel).get_lam();
    double K = (*CurrentModel).get_K();
    double d = (*CurrentModel).get_d();
    Vector3d n_E0 = (*CurrentModel).get_nE0();
    Vector3d n_K = (*CurrentModel).get_nK();
    VectorXcd* al = (*CurrentModel).get_al();
    VectorXcd* P = (*CurrentModel).get_P();

    VectorXcd Adevxp=VectorXcd::Zero(3*N);

    int n_para=(*Free).size();
    int n_para_all = (*Para).size();
    
    VectorXd devx=VectorXd::Zero(n_para);
    vector<list<int>> Paratogeometry(n_para_all);
    //cout << "n_para:  " << n_para << endl;
    //cout << "geometry_Para size: " << (*geometryPara).size() << endl;
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[(*geometryPara)(i)]).push_back(i);
    }



    //Vector3d diel_old_tmp = Vector3d::Zero();
    //Vector3cd diel_tmp = Vector3cd::Zero();

    for(int i = 0; i <= n_para - 1; i++){
        int FreeParaPos = (*Free)(i);

        double diel_old_tmp = (*Para)(FreeParaPos);
        int sign = 0;
        if (diel_old_tmp >= epsilon) {
            sign = -1;
        }
        else {
            sign = 1;
        }
        diel_old_tmp += sign * epsilon;

        

        //because i am changing diel_old_tmp as local variable and this does not influence diel_old, the singleresponse will not respond to this change
        //if (objective->Have_Devx) objective->SingleResponse(position1, true);
        int labelfloor = int(floor((*diel_old)(i)));
        int labelnext = labelfloor + 1;
        if (labelfloor >= 1) {
            labelnext = labelfloor;
        }
        std::complex<double> diel_tmp = (*material)(labelfloor) + ((*diel_old)(i) - double(labelfloor)) * ((*material)(labelnext) - (*material)(labelfloor));

        //if (objective->Have_Devx) objective->SingleResponse(position1, false);
        complex<double> oneoveralpha = (1.0 / Get_Alpha(lam, K, d, diel_tmp, n_E0, n_K));
        
        
        list<int>::iterator it = Paratogeometry[FreeParaPos].begin();
        for (int j = 0; j <= Paratogeometry[FreeParaPos].size() - 1; j++) {
            int position = *it;
            complex<double> change = (oneoveralpha - (*al)(3 * position)) / (sign * epsilon);
            Adevxp(3 * position) = change;
            Adevxp(3 * position+1) = change;
            Adevxp(3 * position+2) = change;
            it++;
        }
                   
        devx(i)=(objective->GroupResponse()-origin)/(sign*epsilon);
        /*
        if(objective->Have_Devx) objective->SingleResponse(position1, true);
                
        if(objective->Have_Devx) objective->SingleResponse(position1, false);

        it2 = (*it1).begin();

        for (int j = 0; j <= (*it1).size()-1; j++) {
            int position2 = (*it2);
                  
            if(objective->Have_Devx) objective->SingleResponse(position2, true);                  
                    
            if(objective->Have_Devx) objective->SingleResponse(position2, false);

            it2++;
                    
        }
        */
            
    }
    for(int i=0;i<=3*N-1;i++){
        Adevxp(i) = Adevxp(i) * ((*P)(i));
    }

    return make_tuple(devx, Adevxp);  
}

tuple<VectorXd, VectorXcd> EvoDDAModel::devx_and_Adevxp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin) {
    int N = (*CurrentModel).get_N();
    SpacePara* spacepara = (*CurrentModel).get_spacepara();
    VectorXi* geometryPara = (*spacepara).get_geometryPara();
    VectorXd* Para = (*spacepara).get_Para();
    VectorXi* Free = (*spacepara).get_Free();
    //vector<list<int>>* Paratogeometry = (*spacepara).get_Paratogeometry();
    
    /*
    list<int> tmplist = (*Paratogeometry)[3000];
    list<int>::iterator tmpit = (*Paratogeometry)[2021].begin();
    for (int i = 0; i <= (*Paratogeometry)[3000].size() - 1; i++) {
        cout << (*tmpit) << endl;
        tmpit++;
    }
    */
    
    VectorXd* diel_old = (*CurrentModel).get_diel_old();
    VectorXcd* material = (*CurrentModel).get_material();
    double lam = (*CurrentModel).get_lam();
    double K = (*CurrentModel).get_K();
    double d = (*CurrentModel).get_d();
    Vector3d n_E0 = (*CurrentModel).get_nE0();
    Vector3d n_K = (*CurrentModel).get_nK();
    VectorXcd* al = (*CurrentModel).get_al();
    VectorXcd* P = (*CurrentModel).get_P();

    VectorXcd Adevxp = VectorXcd::Zero(3 * N);

    int n_para = (*Free).size();
    int n_para_all = (*Para).size();

    VectorXd devx = VectorXd::Zero(n_para);
    
    //cout << "n_para:  " << n_para << endl;
    //cout << "geometry_Para size: " << (*geometryPara).size() << endl;
    



    vector<list<int>> Paratogeometry(n_para_all);
    for (int i = 0; i <= N - 1; i++) {
        (Paratogeometry[(*geometryPara)(i)]).push_back(i);
    }

    if (!(*spacepara).get_Filter()) {
        for (int i = 0; i <= n_para - 1; i++) {
            int FreeParaPos = (*Free)(i);

            double diel_old_origin = (*Para)(FreeParaPos);
            double diel_old_tmp = diel_old_origin;
            int sign = 0;
            if (diel_old_origin >= epsilon) {
                sign = -1;
            }
            else {
                sign = 1;
            }
            diel_old_tmp += sign * epsilon;



            list<int>::iterator it = Paratogeometry[FreeParaPos].begin();
            for (int j = 0; j <= Paratogeometry[FreeParaPos].size() - 1; j++) {
                //cout << (*it) << endl;
                int position = *it;
                complex<double> alphaorigin = (*al)(3 * position);
                if (objective->Have_Devx) objective->SingleResponse(position, true);
                (*CStr).UpdateStrSingle(position, diel_old_tmp);
                (*CurrentModel).UpdateAlphaSingle(position);
                if (objective->Have_Devx) objective->SingleResponse(position, false);
                complex<double> change = ((*al)(3 * position) - alphaorigin) / (sign * epsilon);
                Adevxp(3 * position) = change;
                Adevxp(3 * position + 1) = change;
                Adevxp(3 * position + 2) = change;

                it++;
            }

            devx(i) = (objective->GroupResponse() - origin) / (sign * epsilon);  //If some obj has x dependency but you denote the havepenalty as false, it will still actually be calculated in an efficient way.
            it = Paratogeometry[FreeParaPos].begin();
            for (int j = 0; j <= Paratogeometry[FreeParaPos].size() - 1; j++) {
                int position = *it;
                if (objective->Have_Devx) objective->SingleResponse(position, true);
                (*CStr).UpdateStrSingle(position, diel_old_origin);
                (*CurrentModel).UpdateAlphaSingle(position);
                if (objective->Have_Devx) objective->SingleResponse(position, false);
                it++;
            }

        }
    }
    else {
        for (int i = 0; i <= n_para - 1; i++) {
            int FreeParaPos = (*Free)(i);

            double diel_old_origin = (*Para)(FreeParaPos);
            double diel_old_tmp = diel_old_origin;
            int sign = 0;
            if (diel_old_origin >= epsilon) {
                sign = -1;
            }
            else {
                sign = 1;
            }
            diel_old_tmp += sign * epsilon;



            list<int>::iterator it = Paratogeometry[FreeParaPos].begin();
            for (int j = 0; j <= Paratogeometry[FreeParaPos].size() - 1; j++) {
                //cout << (*it) << endl;
                int position = *it;
                complex<double> alphaorigin = (*al)(3 * position);
                if (objective->Have_Devx) objective->SingleResponse(position, true);
                (*CStr).UpdateStrSingle(position, diel_old_tmp);
                (*CurrentModel).UpdateAlphaSingle(position);
                if (objective->Have_Devx) objective->SingleResponse(position, false);
                complex<double> change = ((*al)(3 * position) - alphaorigin) / (sign * epsilon);
                Adevxp(3 * position) = change;
                Adevxp(3 * position + 1) = change;
                Adevxp(3 * position + 2) = change;

                it++;
            }

            devx(i) = (objective->GroupResponse() - origin) / (sign * epsilon);  //If some obj has x dependency but you denote the havepenalty as false, it will still actually be calculated in an efficient way.
            it = Paratogeometry[FreeParaPos].begin();
            for (int j = 0; j <= Paratogeometry[FreeParaPos].size() - 1; j++) {
                int position = *it;
                if (objective->Have_Devx) objective->SingleResponse(position, true);
                (*CStr).UpdateStrSingle(position, diel_old_origin);
                (*CurrentModel).UpdateAlphaSingle(position);
                if (objective->Have_Devx) objective->SingleResponse(position, false);
                it++;
            }

        }
    }
    
    
    for (int i = 0; i <= 3 * N - 1; i++) {
        Adevxp(i) = Adevxp(i) * ((*P)(i));
    }

    return make_tuple(devx, Adevxp);
}

VectorXcd EvoDDAModel::devp(double epsilon, DDAModel* CurrentModel, ObjectiveDDAModel* objective, double origin){
    //move origin=objective0->GetVal() outside because it is the same for one partial derivative of the entire structure
    VectorXcd* P = (*CurrentModel).get_P();
    VectorXcd result=VectorXcd::Zero((*P).size());
    for(int i=0;i<= (*P).size() -1;i++){
        int position = i/3;
        
        objective->SingleResponse(position, true);
        
        (*P)(i)= (*P)(i)+epsilon;
        
        objective->SingleResponse(position, false);
        
        result(i)+=(objective->GroupResponse()-origin)/epsilon;
        
        objective->SingleResponse(position, true);
        
        (*P)(i)= (*P)(i)-epsilon;
        
        complex<double> tmp=epsilon*1.0i;
        
        (*P)(i)= (*P)(i)+tmp;
        
        objective->SingleResponse(position, false);
        
        complex<double> tmpRes = (objective->GroupResponse()-origin)/tmp;
        result(i)+=tmpRes;
        
        objective->SingleResponse(position, true);
        
        (*P)(i)= (*P)(i)-tmp;
        
        objective->SingleResponse(position, false);
    }
    //cout << "Devp_sum: " << result.sum() << endl;
    return result;
}



void EvoDDAModel::EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, double start_num){
    ofstream convergence;
    ofstream Originiterations;
    ofstream Adjointiterations;
    int TotalOriginIt = 0;
    int TotalAdjointIt = 0;
    
    convergence.open(save_position + "convergence.txt");
    Originiterations.open(save_position + "Originiterations.txt");
    Adjointiterations.open(save_position + "Adjointiterations.txt");
    //Parameters for Adam Optimizer.
    //double beta1 = 0.9;
    //double beta2 = 0.99;
    //---------------new beta for THG test----------------
    double beta1 = 0.9;
    double beta2 = 1 - (1 - beta1) * (1 - beta1);
    VectorXd V;
    VectorXd S;
    

    high_resolution_clock::time_point TotalTime0 = high_resolution_clock::now();

    double epsilon_partial=0.001;
    for(int iteration=0;iteration<=MAX_ITERATION_EVO-1;iteration++){
        //solve DDA
        cout << "######################EVO ITERATION " << iteration << "#######################" << endl;
        //get object function value
        cout<<"-----------------------------START ORIGINAL PROBLEM---------------------------"<<endl;

        double obj;
        VectorXd objarray = VectorXd::Zero(ModelNum);
        list<DDAModel*>::iterator it_ModelList = ModelList.begin();
        list<ObjectiveDDAModel*>::iterator it_ObjList = ObjList.begin();

        high_resolution_clock::time_point out_start = high_resolution_clock::now();
        (*CStr).output_to_file(save_position + "CoreStructure/", iteration + start_num, "simple");
        high_resolution_clock::time_point out_end = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(out_end - out_start).count();
        output_time += duration;

        list<VectorXcd>::iterator it_PforOrigin = PforOrigin.begin();
        for (int i = 0; i <= ModelNum - 1; i++) {
            //cout << (*(it_PforOrigin))(0) << endl;
            (*(*it_ModelList)).InitializeP(*(it_PforOrigin));
            (*(*it_ModelList)).bicgstab(MAX_ITERATION, MAX_ERROR);
            if (HaveOriginHeritage == true) {
                *(it_PforOrigin) = *((*(*it_ModelList)).get_P());
            }
            it_PforOrigin++;
            (*(*it_ModelList)).update_E_in_structure();
            if (iteration == MAX_ITERATION_EVO - 1) {                                    //useless fix, not gonna to use RResultswithc = true feature in the future
                (*(*it_ModelList)).solve_E();
            }

            (*(*it_ModelList)).solve_E();
            out_start = high_resolution_clock::now();
            (*(*it_ModelList)).output_to_file(save_position + "Model_output/", iteration + start_num);
            out_end = high_resolution_clock::now();
            duration = duration_cast<milliseconds>(out_end - out_start).count();
            output_time += duration;

            objarray(i) = (*(*it_ObjList)).GetVal();

            Originiterations << (*(*it_ModelList)).get_ITERATION() << endl;
            TotalOriginIt += (*(*it_ModelList)).get_ITERATION();

            it_ModelList++;
            it_ObjList++;
        }
        obj = objarray.sum()/ModelNum;                              //Take average
        convergence << obj << " ";
        cout << "objective function at iteration " << iteration << " is " << obj << endl;

        high_resolution_clock::time_point t0 = high_resolution_clock::now();

        VectorXd* diel_old = (*CStr).get_diel_old();
        VectorXd* diel_old_max = (*CStr).get_diel_old_max();

        double epsilon = epsilon_fix;
        

        if (HavePathRecord) {

            if ((abs(obj - PreviousObj)) / PreviousObj <= 0.0001 || epsilon <= 0.0001) {
                CutoffHold += 1;
            }
            else {
                if (CutoffHold > 0) {
                    CutoffHold -= 1;
                }
            }
            cout << "CutoffHold" << CutoffHold << endl;
            PreviousObj = obj;
            if (CutoffHold >= 3) {
                cout << "Three times with small change in obj, break the iterations" << endl;
                //break;
            }

            if (obj < MaxObj) {
                epsilon_tmp = epsilon_tmp / 10;
                Stephold = 0;
                (*diel_old) = (*diel_old_max);
                it_ModelList = ModelList.begin();
                it_PforOrigin = PforOrigin.begin();
                list<VectorXcd>::iterator it_PforAdjoint = PforAdjoint.begin();
                list<VectorXcd>::iterator it_PforOriginMax = PforOriginMax.begin();
                list<VectorXcd>::iterator it_PforAdjointMax = PforAdjointMax.begin();
                for (int i = 0; i <= ModelNum - 1; i++) {
                    VectorXcd* P = (*(*it_ModelList)).get_P();
                    VectorXcd* P_max = (*(*it_ModelList)).get_P_max();
                    VectorXcd* al = (*(*it_ModelList)).get_al();
                    VectorXcd* al_max = (*(*it_ModelList)).get_al_max();
                    (*P) = (*P_max);
                    (*al) = (*al_max);
                    objarray(i) = MaxObjarray(i);
                    *(it_PforOrigin) = *(it_PforOriginMax);
                    *(it_PforAdjoint) = *(it_PforAdjointMax);
                    it_ModelList++;
                    it_PforOrigin++;
                    it_PforAdjoint++;
                    it_PforOriginMax++;
                    it_PforAdjointMax++;
                }
  
                obj = MaxObj;
                cout << "New Obj smaller then Old One, back track to previous structure and search with new step size: " << epsilon_tmp << endl;
                /*
                if (obj != objective->GetVal()) {
                    cout << "Reset failed, objective is not equal to MaxObj" << endl;
                }
                */
            }
            else {

                if ((abs(obj - PreviousObj)) / PreviousObj <= 0.0001) {
                    CutoffHold += 1;
                }
                else {
                    if (CutoffHold > 0) {
                        CutoffHold -= 1;
                    }
                }

                (*diel_old_max) = (*diel_old);
                it_ModelList = ModelList.begin();
                it_PforOrigin = PforOrigin.begin();
                list<VectorXcd>::iterator it_PforAdjoint = PforAdjoint.begin();
                list<VectorXcd>::iterator it_PforOriginMax = PforOriginMax.begin();
                list<VectorXcd>::iterator it_PforAdjointMax = PforAdjointMax.begin();
                for (int i = 0; i <= ModelNum - 1; i++) {
                    VectorXcd* P = (*(*it_ModelList)).get_P();
                    VectorXcd* P_max = (*(*it_ModelList)).get_P_max();
                    VectorXcd* al = (*(*it_ModelList)).get_al();
                    VectorXcd* al_max = (*(*it_ModelList)).get_al_max();
                    (*P_max) = (*P);
                    (*al_max) = (*al);
                    MaxObjarray(i) = objarray(i);
                    *(it_PforOriginMax) = *(it_PforOrigin);
                    *(it_PforAdjointMax) = *(it_PforAdjoint);
                    it_ModelList++;
                    it_PforOrigin++;
                    it_PforAdjoint++;
                    it_PforOriginMax++;
                    it_PforAdjointMax++;
                }
                MaxObj = obj;
                Stephold += 1;
                
                if (Stephold >= 2) {
                    epsilon_tmp = epsilon_tmp * 10;
                    cout << "Two times increase with previous step size, try with larger step size: " << epsilon_tmp << endl;
                    Stephold = 0;
                }
                else {
                    cout << "Not smaller obj nor three continuous increase. Current step size is: " << epsilon_tmp << endl;
                }   
            }
            epsilon = epsilon_tmp;
        }

        for (int i = 0; i <= ModelNum - 1; i++) {        //Origins equals to the objective functions of each model of current structure
            Originarray(i) = objarray(i);
        }
           
        
        /*
        if((*ObjectFunctionNames).size()>1){
            list<double> obj_minor =  this->MinorObjective();
            list<double>::iterator it_obj_minor = obj_minor.begin();
            for(int i=0; i<=obj_minor.size()-1; i++){
                convergence << *it_obj_minor << " ";
                it_obj_minor++;
            }
        }
        */
        convergence << "\n";

        
        
        
        SpacePara* spacepara = (*CStr).get_spacepara();
        VectorXi* geometryPara = (*spacepara).get_geometryPara();
        VectorXd* Para = (*spacepara).get_Para();
        VectorXi* Free = (*spacepara).get_Free();
        int n_para_all = (*Para).size();
        int n_para = (*Free).size();
        int N = (*CStr).get_N();


        VectorXd gradients = VectorXd::Zero(n_para);

        it_ModelList = ModelList.begin();
        it_ObjList = ObjList.begin();
        list<VectorXcd>::iterator it_PforAdjoint = PforAdjoint.begin();
        for (int i = 0; i <= ModelNum - 1; i++) {
            
            //----------------------------------------get partial derivative of current model---------------------------
            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            //cout << "---------------------------START PARTIAL DERIVATIVE of Model" << i << " ----------------------" << endl;
            VectorXd devx;
            VectorXcd Adevxp;
            VectorXcd devp;
            tie(devx, Adevxp) = this->devx_and_Adevxp(epsilon_partial, *it_ModelList, *it_ObjList, Originarray(i));
            devp = this->devp(epsilon_partial, *it_ModelList, *it_ObjList, Originarray(i));
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(t2 - t1).count();
            //cout << "------------------------PARTIAL DERIVATIVE finished in " << duration / 1000 << " s-------------------------" << endl;

            //------------------------------------Solving adjoint problem-----------------------------------------
            //cout << "---------------------------START ADJOINT PROBLEM of Model" << i << " ----------------------" << endl;
            (*(*it_ModelList)).change_E(devp);

            (*(*it_ModelList)).InitializeP(*(it_PforAdjoint));

            (*(*it_ModelList)).bicgstab(MAX_ITERATION, MAX_ERROR);

            if (HaveAdjointHeritage == true) {
                *(it_PforAdjoint) = *((*(*it_ModelList)).get_P());
            }
            it_PforAdjoint++;

            VectorXcd* P = (*(*it_ModelList)).get_P();
            VectorXcd lambdaT = (*P);
            (*(*it_ModelList)).reset_E();                                  //reset E to initial value
            Adjointiterations << (*(*it_ModelList)).get_ITERATION() << endl;
            TotalAdjointIt += (*(*it_ModelList)).get_ITERATION();


            //times lambdaT and Adevxp together
            VectorXcd mult_result;
            mult_result = VectorXcd::Zero(n_para);               //multiplication result has the length of parameter
            vector<list<int>> Paratogeometry(n_para_all);
            for (int i = 0; i <= N - 1; i++) {
                (Paratogeometry[(*geometryPara)(i)]).push_back(i);
            }

            for (int i = 0; i <= n_para - 1; i++) {
                int FreeParaPos = (*Free)(i);

                list<int>::iterator it = Paratogeometry[FreeParaPos].begin();
                for (int j = 0; j <= Paratogeometry[FreeParaPos].size() - 1; j++) {
                    int position = *it;
                    mult_result(i) += lambdaT(3 * position) * Adevxp(3 * position);
                    mult_result(i) += lambdaT(3 * position + 1) * Adevxp(3 * position + 1);
                    mult_result(i) += lambdaT(3 * position + 2) * Adevxp(3 * position + 2);
                    it++;
                }
            }


            VectorXd mult_result_real = VectorXd::Zero(n_para);
            for (int i = 0; i <= n_para - 1; i++) {
                complex<double> tmp = mult_result(i);
                mult_result_real(i) = tmp.real();
            }
            gradients += devx - mult_result_real;              //What's the legitimacy in here to ignore the imag part?
            
            it_ModelList++;
            it_ObjList++;
        }

        //The final gradients for this iteration
        gradients = gradients / ModelNum;

        if ((*((*CStr).get_spacepara())).get_Filter() && (iteration >= 1)) {
            //For iteration=0, Para, Para_origin, Para_filtered are all the same. No need for updating gradients.
            //current_it is actually the it in current evo-1 as the str is updated in iteration-1.
            gradients = gradients_filtered(gradients, iteration-1, MAX_ITERATION_EVO-1);  
        }
            
        



        double epsilon_final=epsilon;


        /*
        if ((iteration == 87) || (iteration == 88)) {
            string name1, name2;
            name1 = save_position + to_string(iteration) + "V.txt";
            name2 = save_position + to_string(iteration) + "S.txt";
            
            ofstream fout1(name1), fout2(name2);
            fout1 << V << endl;
            fout2 << S << endl;

            fout1.close();
            fout2.close();

        }
        */
        if(method == "Adam"){
            cout << "Using Adam Optimizer." << endl;
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

            if (iteration <= 3) {
                epsilon_final = 0.1;
            }
            else {
                epsilon_final = epsilon;
            }
        } 

        if (method == "Adamdecay") {
            cout << "Using Adam Optimizer With decay." << endl;
            if (iteration == 0) {
                V = (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
                S = (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));
            }
            else {
                V = beta1 * V + (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
                S = beta2 * S + (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));
            }
            for (int i = 0; i <= n_para - 1; i++) {
                gradients(i) = V(i) / (sqrt(S(i)) + 0.00000001);
            }

            int decaybarrier = 150;

            if (iteration <= 3) {
                epsilon_final = 0.1;
            }
            else if (iteration >= decaybarrier) {
                double expconst = 10;
                epsilon_final = epsilon * exp(-(iteration - decaybarrier) / expconst);
                //epsilon_final = epsilon / exp(-(iteration - decaybarrier));
            }
            else {
                epsilon_final = epsilon;
            }
        }
        
        if (method == "Adagrad") {
            for (int i = 0; i <= n_para - 1; i++) {
                gradientsquare(i) += pow(gradients(i), 2)/100000;
                gradients(i) = gradients(i) / sqrt(gradientsquare(i) + 1);
            }

        }
        cout << "abs(gradients.cwiseAbs().mean() before: "<< abs(gradients.cwiseAbs().mean()) << endl;

        if (abs(gradients.cwiseAbs().mean()) < 0.1) {
            gradients /= abs(gradients.cwiseAbs().mean()) / 0.1;
        }
        
        cout << "abs(gradients.cwiseAbs().mean() after: " << abs(gradients.cwiseAbs().mean()) << endl;
        VectorXd step = epsilon_final * gradients;            //Find the maximum. If -1 find minimum
        //cout << "epsilon = " << epsilon << endl;
        //cout << "step = "<< step.mean() << endl;

           
        (*CStr).UpdateStr(step, iteration, MAX_ITERATION_EVO - 1);
        it_ModelList = ModelList.begin();
        for (int i = 0; i <= ModelNum - 1; i++) {
            (*(*it_ModelList)).UpdateAlpha();                  //Dont forget this, otherwise bicgstab wont change
            it_ModelList++;
        }
        
    }
    
    
    Originiterations << TotalOriginIt << endl;
    Adjointiterations << TotalAdjointIt << endl;

    convergence.close();
    Originiterations.close();
    Adjointiterations.close();

}

void EvoDDAModel::EvoOptimization(int MAX_ITERATION, double MAX_ERROR, int MAX_ITERATION_EVO, string method, VectorXd* V_, VectorXd* S_) {
    ofstream convergence;
    ofstream Originiterations;
    ofstream Adjointiterations;
    int TotalOriginIt = 0;
    int TotalAdjointIt = 0;

    convergence.open(save_position + "convergence.txt");
    Originiterations.open(save_position + "Originiterations.txt");
    Adjointiterations.open(save_position + "Adjointiterations.txt");
    //Parameters for Adam Optimizer.
    //double beta1 = 0.9;
    //double beta2 = 0.99;
    //---------------new beta for THG test----------------
    double beta1 = 0.7;
    double beta2 = 1 - (1 - beta1) * (1 - beta1);
    VectorXd V = *V_;
    VectorXd S = *S_;


    high_resolution_clock::time_point TotalTime0 = high_resolution_clock::now();

    double epsilon_partial = 0.001;
    for (int iteration = 0; iteration <= MAX_ITERATION_EVO - 1; iteration++) {
        //solve DDA
        cout << "######################EVO ITERATION " << iteration << "#######################" << endl;
        //get object function value
        cout << "-----------------------------START ORIGINAL PROBLEM---------------------------" << endl;

        double obj;
        VectorXd objarray = VectorXd::Zero(ModelNum);
        list<DDAModel*>::iterator it_ModelList = ModelList.begin();
        list<ObjectiveDDAModel*>::iterator it_ObjList = ObjList.begin();
        (*CStr).output_to_file(save_position, iteration);
        list<VectorXcd>::iterator it_PforOrigin = PforOrigin.begin();
        for (int i = 0; i <= ModelNum - 1; i++) {
            //cout << (*(it_PforOrigin))(0) << endl;
            (*(*it_ModelList)).InitializeP(*(it_PforOrigin));
            (*(*it_ModelList)).bicgstab(MAX_ITERATION, MAX_ERROR);
            if (HaveOriginHeritage == true) {
                *(it_PforOrigin) = *((*(*it_ModelList)).get_P());
            }
            it_PforOrigin++;
            (*(*it_ModelList)).update_E_in_structure();
            if (iteration == MAX_ITERATION_EVO - 1) {                                    //useless fix, not gonna to use RResultswithc = true feature in the future
                (*(*it_ModelList)).solve_E();
            }
            //(*(*it_ModelList)).output_to_file(save_position + "Model_output\\", iteration, i);
            objarray(i) = (*(*it_ObjList)).GetVal();

            Originiterations << (*(*it_ModelList)).get_ITERATION() << endl;
            TotalOriginIt += (*(*it_ModelList)).get_ITERATION();

            it_ModelList++;
            it_ObjList++;
        }
        obj = objarray.sum() / ModelNum;                              //Take average
        convergence << obj << " ";
        cout << "objective function at iteration " << iteration << " is " << obj << endl;

        high_resolution_clock::time_point t0 = high_resolution_clock::now();

        VectorXd* diel_old = (*CStr).get_diel_old();
        VectorXd* diel_old_max = (*CStr).get_diel_old_max();

        double epsilon = epsilon_fix;


        if (HavePathRecord) {

            if ((abs(obj - PreviousObj)) / PreviousObj <= 0.0001 || epsilon <= 0.0001) {
                CutoffHold += 1;
            }
            else {
                if (CutoffHold > 0) {
                    CutoffHold -= 1;
                }
            }
            cout << "CutoffHold" << CutoffHold << endl;
            PreviousObj = obj;
            if (CutoffHold >= 3) {
                cout << "Three times with small change in obj, break the iterations" << endl;
                //break;
            }

            if (obj < MaxObj) {
                epsilon_tmp = epsilon_tmp / 10;
                Stephold = 0;
                (*diel_old) = (*diel_old_max);
                it_ModelList = ModelList.begin();
                it_PforOrigin = PforOrigin.begin();
                list<VectorXcd>::iterator it_PforAdjoint = PforAdjoint.begin();
                list<VectorXcd>::iterator it_PforOriginMax = PforOriginMax.begin();
                list<VectorXcd>::iterator it_PforAdjointMax = PforAdjointMax.begin();
                for (int i = 0; i <= ModelNum - 1; i++) {
                    VectorXcd* P = (*(*it_ModelList)).get_P();
                    VectorXcd* P_max = (*(*it_ModelList)).get_P_max();
                    VectorXcd* al = (*(*it_ModelList)).get_al();
                    VectorXcd* al_max = (*(*it_ModelList)).get_al_max();
                    (*P) = (*P_max);
                    (*al) = (*al_max);
                    objarray(i) = MaxObjarray(i);
                    *(it_PforOrigin) = *(it_PforOriginMax);
                    *(it_PforAdjoint) = *(it_PforAdjointMax);
                    it_ModelList++;
                    it_PforOrigin++;
                    it_PforAdjoint++;
                    it_PforOriginMax++;
                    it_PforAdjointMax++;
                }

                obj = MaxObj;
                cout << "New Obj smaller then Old One, back track to previous structure and search with new step size: " << epsilon_tmp << endl;
                /*
                if (obj != objective->GetVal()) {
                    cout << "Reset failed, objective is not equal to MaxObj" << endl;
                }
                */
            }
            else {

                if ((abs(obj - PreviousObj)) / PreviousObj <= 0.0001) {
                    CutoffHold += 1;
                }
                else {
                    if (CutoffHold > 0) {
                        CutoffHold -= 1;
                    }
                }

                (*diel_old_max) = (*diel_old);
                it_ModelList = ModelList.begin();
                it_PforOrigin = PforOrigin.begin();
                list<VectorXcd>::iterator it_PforAdjoint = PforAdjoint.begin();
                list<VectorXcd>::iterator it_PforOriginMax = PforOriginMax.begin();
                list<VectorXcd>::iterator it_PforAdjointMax = PforAdjointMax.begin();
                for (int i = 0; i <= ModelNum - 1; i++) {
                    VectorXcd* P = (*(*it_ModelList)).get_P();
                    VectorXcd* P_max = (*(*it_ModelList)).get_P_max();
                    VectorXcd* al = (*(*it_ModelList)).get_al();
                    VectorXcd* al_max = (*(*it_ModelList)).get_al_max();
                    (*P_max) = (*P);
                    (*al_max) = (*al);
                    MaxObjarray(i) = objarray(i);
                    *(it_PforOriginMax) = *(it_PforOrigin);
                    *(it_PforAdjointMax) = *(it_PforAdjoint);
                    it_ModelList++;
                    it_PforOrigin++;
                    it_PforAdjoint++;
                    it_PforOriginMax++;
                    it_PforAdjointMax++;
                }
                MaxObj = obj;
                Stephold += 1;

                if (Stephold >= 2) {
                    epsilon_tmp = epsilon_tmp * 10;
                    cout << "Two times increase with previous step size, try with larger step size: " << epsilon_tmp << endl;
                    Stephold = 0;
                }
                else {
                    cout << "Not smaller obj nor three continuous increase. Current step size is: " << epsilon_tmp << endl;
                }
            }
            epsilon = epsilon_tmp;
        }

        for (int i = 0; i <= ModelNum - 1; i++) {        //Origins equals to the objective functions of each model of current structure
            Originarray(i) = objarray(i);
        }


        /*
        if((*ObjectFunctionNames).size()>1){
            list<double> obj_minor =  this->MinorObjective();
            list<double>::iterator it_obj_minor = obj_minor.begin();
            for(int i=0; i<=obj_minor.size()-1; i++){
                convergence << *it_obj_minor << " ";
                it_obj_minor++;
            }
        }
        */
        convergence << "\n";




        SpacePara* spacepara = (*CStr).get_spacepara();
        VectorXi* geometryPara = (*spacepara).get_geometryPara();
        VectorXd* Para = (*spacepara).get_Para();
        VectorXi* Free = (*spacepara).get_Free();
        int n_para_all = (*Para).size();
        int n_para = (*Free).size();
        int N = (*CStr).get_N();


        VectorXd gradients = VectorXd::Zero(n_para);

        it_ModelList = ModelList.begin();
        it_ObjList = ObjList.begin();
        list<VectorXcd>::iterator it_PforAdjoint = PforAdjoint.begin();
        for (int i = 0; i <= ModelNum - 1; i++) {

            //----------------------------------------get partial derivative of current model---------------------------
            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            //cout << "---------------------------START PARTIAL DERIVATIVE of Model" << i << " ----------------------" << endl;
            VectorXd devx;
            VectorXcd Adevxp;
            VectorXcd devp;
            tie(devx, Adevxp) = this->devx_and_Adevxp(epsilon_partial, *it_ModelList, *it_ObjList, Originarray(i));
            devp = this->devp(epsilon_partial, *it_ModelList, *it_ObjList, Originarray(i));
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(t2 - t1).count();
            //cout << "------------------------PARTIAL DERIVATIVE finished in " << duration / 1000 << " s-------------------------" << endl;

            //------------------------------------Solving adjoint problem-----------------------------------------
            //cout << "---------------------------START ADJOINT PROBLEM of Model" << i << " ----------------------" << endl;
            (*(*it_ModelList)).change_E(devp);

            (*(*it_ModelList)).InitializeP(*(it_PforAdjoint));
            (*(*it_ModelList)).bicgstab(MAX_ITERATION, MAX_ERROR);
            if (HaveAdjointHeritage == true) {
                *(it_PforAdjoint) = *((*(*it_ModelList)).get_P());
            }
            it_PforAdjoint++;

            VectorXcd* P = (*(*it_ModelList)).get_P();
            VectorXcd lambdaT = (*P);
            (*(*it_ModelList)).reset_E();                                  //reset E to initial value
            Adjointiterations << (*(*it_ModelList)).get_ITERATION() << endl;
            TotalAdjointIt += (*(*it_ModelList)).get_ITERATION();


            //times lambdaT and Adevxp together
            VectorXcd mult_result;
            mult_result = VectorXcd::Zero(n_para);               //multiplication result has the length of parameter
            vector<list<int>> Paratogeometry(n_para_all);
            for (int i = 0; i <= N - 1; i++) {
                (Paratogeometry[(*geometryPara)(i)]).push_back(i);
            }

            for (int i = 0; i <= n_para - 1; i++) {
                int FreeParaPos = (*Free)(i);
                list<int>::iterator it = Paratogeometry[FreeParaPos].begin();
                for (int j = 0; j <= Paratogeometry[FreeParaPos].size() - 1; j++) {
                    int position = *it;
                    mult_result(i) += lambdaT(3 * position) * Adevxp(3 * position);
                    mult_result(i) += lambdaT(3 * position + 1) * Adevxp(3 * position + 1);
                    mult_result(i) += lambdaT(3 * position + 2) * Adevxp(3 * position + 2);
                    it++;
                }
            }


            VectorXd mult_result_real = VectorXd::Zero(n_para);
            for (int i = 0; i <= n_para - 1; i++) {
                complex<double> tmp = mult_result(i);
                mult_result_real(i) = tmp.real();
            }
            gradients += devx - mult_result_real;              //What's the legitimacy in here to ignore the imag part?

            it_ModelList++;
            it_ObjList++;
        }

        //The final gradients for this iteration
        gradients = gradients / ModelNum;
        double epsilon_final = epsilon;

        if ((iteration == 87) || (iteration == 88)) {
            string name1, name2;
            name1 = save_position + to_string(iteration) + "V.txt";
            name2 = save_position + to_string(iteration) + "S.txt";

            ofstream fout1(name1), fout2(name2);
            fout1 << V << endl;
            fout2 << S << endl;

            fout1.close();
            fout2.close();

        }

        if (method == "Adam") {
            cout << "Using Adam Optimizer." << endl;

            V = beta1 * V + (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
            S = beta2 * S + (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));

            for (int i = 0; i <= n_para - 1; i++) {
                gradients(i) = V(i) / (sqrt(S(i)) + 0.00000001);
            }

            if (iteration <= 3) {
                epsilon_final = 0.1;
            }
            else {
                epsilon_final = epsilon;
            }
        }

        if (method == "Adamdecay") {
            cout << "Using Adam Optimizer With decay." << endl;

            if (iteration == 0) {
                V = (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
                S = (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));
            }
            else {
                V = beta1 * V + (1 - beta1) * gradients / (1 - pow(beta1, iteration + 1));
                S = beta2 * S + (1 - beta2) * (gradients.array().pow(2).matrix()) / (1 - pow(beta2, iteration + 1));
            }

            for (int i = 0; i <= n_para - 1; i++) {
                gradients(i) = V(i) / (sqrt(S(i)) + 0.00000001);
            }

            int decaybarrier = 100;

            if (iteration <= 3) {
                epsilon_final = 0.1;
            }
            else if (iteration >= decaybarrier) {
                double expconst = 10;
                epsilon_final = epsilon * exp(-(iteration - decaybarrier) / expconst);
                //epsilon_final = epsilon / exp(-(iteration - decaybarrier));
            }
            else {
                epsilon_final = epsilon;
            }
        }

        if (method == "Adagrad") {
            for (int i = 0; i <= n_para - 1; i++) {
                gradientsquare(i) += pow(gradients(i), 2) / 100000;
                gradients(i) = gradients(i) / sqrt(gradientsquare(i) + 1);
            }

        }
        cout << "abs(gradients.cwiseAbs().mean() before: " << abs(gradients.cwiseAbs().mean()) << endl;

        if (abs(gradients.cwiseAbs().mean()) < 0.1) {
            gradients /= abs(gradients.cwiseAbs().mean()) / 0.1;
        }

        cout << "abs(gradients.cwiseAbs().mean() after: " << abs(gradients.cwiseAbs().mean()) << endl;
        VectorXd step = epsilon_final * gradients;            //Find the maximum. If -1 find minimum
        //cout << "epsilon = " << epsilon << endl;
        //cout << "step = "<< step.mean() << endl;


        (*CStr).UpdateStr(step, iteration, MAX_ITERATION_EVO - 1);
        it_ModelList = ModelList.begin();
        for (int i = 0; i <= ModelNum - 1; i++) {
            (*(*it_ModelList)).UpdateAlpha();                  //Dont forget this, otherwise bicgstab wont change
            it_ModelList++;
        }

    }


    Originiterations << TotalOriginIt << endl;
    Adjointiterations << TotalAdjointIt << endl;

    convergence.close();
    Originiterations.close();
    Adjointiterations.close();

}

ObjectiveDDAModel* EvoDDAModel::ObjectiveFactory(string ObjectName, list<double> ObjectParameters, DDAModel* ObjDDAModel){
    if (HavePenalty) {
        cout << "Using L1 Penalty with Penalty Factor " << PenaltyFactor << endl;
    }
    if (MajorObjectFunctionName == "PointE"){
        return new ObjectivePointEDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    }
    if (MajorObjectFunctionName == "PointEList") {
        return new ObjectivePointListEDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    }
    if (MajorObjectFunctionName == "PointI") {
        return new ObjectivePointIDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    }
    if (MajorObjectFunctionName == "IntegratedE") {
        return new ObjectiveIntegratedEDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    }
    if (MajorObjectFunctionName == "MidAvgE") {
        return new ObjectiveMidAvgEDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    }
    if (MajorObjectFunctionName == "scattering0D") {
        return new Objectivescattering0D(ObjectParameters, ObjDDAModel, this, HavePenalty);
    }
    if (MajorObjectFunctionName == "Abs") {
        return new ObjectiveAbs(ObjectParameters, ObjDDAModel, this, HavePenalty);
    }

    /*
    else if (MajorObjectFunctionName == "SurfaceEExp"){
        return new ObjectiveSurfaceEExp(ObjectParameters, this, HavePenalty);
    }
    else if (MajorObjectFunctionName == "ExtSurfaceEExp"){
        return new ObjectiveExtSurfaceEExp(ObjectParameters, this, HavePenalty);
    }
    else if (MajorObjectFunctionName == "ExtSurfaceEExp_CPU") {
        return new ObjectiveExtSurfaceEExp_CPU(ObjectParameters, this, HavePenalty);
    }
    else if (MajorObjectFunctionName == "ExtSurfaceEMax") {
        return new ObjectiveExtSurfaceEMax(ObjectParameters, this, HavePenalty);
    }
    else if (MajorObjectFunctionName == "ExtSurfaceEExp_CPU_Old") {
        return new ObjectiveExtSurfaceEExp_CPU_Old(ObjectParameters, this, HavePenalty);
    }
    else if (MajorObjectFunctionName == "ObjectiveG") {
        return new ObjectiveG(ObjectParameters, this, HavePenalty);
    }*/
    else{
        // NOT FINALIZED. SHOULD RAISE AN EXCEPTION HERE.
        cout << "NOT A LEGIT OBJECTIVE NAME!" << endl;
        return new ObjectivePointEDDAModel(ObjectParameters, ObjDDAModel, this, HavePenalty);
    }

}

double EvoDDAModel::L1Norm(){
    double Penalty = 0;
    int N = (*CStr).get_N();
    VectorXd* diel_old = (*CStr).get_diel_old();
    for (int i=0;i<N;i++){
        Penalty += 0.5-abs((*diel_old)(3*i)-0.5);
    }
    Penalty = Penalty * PenaltyFactor;
    return Penalty;
}

double EvoDDAModel::get_output_time() {
    return output_time / 1000;      //seconds
}

double PtoFderivative(const double input, const double beta, const double ita) {
    double result = 0.0;
    if (input <= ita && input >= 0.0) {
        return beta * exp(-beta * (1 - input / ita)) + exp(-beta);
    }
    else if (input > ita && input <= 1.0) {
        return beta * exp(-beta * (input - ita) / (1 - ita)) + exp(-beta);
    }
    else {
        cout << "ERROR: PtoFderivative(const double input, const double beta, const double ita)--input out of range" << endl;
        throw 1;
    }
}

VectorXd EvoDDAModel::gradients_filtered(VectorXd gradients, int current_it, int Max_it) {
    SpacePara* sp = (*CStr).get_spacepara();
    FilterOption* fo = (*sp).get_Filterstats();
    const VectorXd* Para_filtered = (*sp).get_Para_filtered();
    const VectorXi* Free = (*sp).get_Free();
    (*fo).update_beta(current_it, Max_it);                     //current_it is actually the it in current evo-1 as the str is updated in iteration-1.
    const double gbeta = (*fo).get_beta();
    const double gita = (*fo).get_ita();
    

    int NFpara = gradients.size();
    VectorXd result = VectorXd::Zero(NFpara);
    const vector<vector<WeightPara>>* FreeWeight = (*sp).get_FreeWeight();
    if (NFpara != (*FreeWeight).size()) {
        cout << "ERROR: EvoDDAModel::gradients_filtered--NFpara != (*FreeWeight).size()" << endl;
        throw 1;
    }
    for (int i = 0; i <= NFpara - 1; i++) {
        //int position = (*Free)(i);
        //double pftmp = (*Para_filtered)(position);
        //double ptofd = PtoFderivative(pftmp, gbeta, gita);

        int num_weight = ((*FreeWeight)[i]).size();
        double sum_weight = 0.0;
        for (int j = 0; j <= num_weight - 1; j++) {
            double weight = (*FreeWeight)[i][j].weight;
            int weightpos = (*FreeWeight)[i][j].position;
            sum_weight += weight;
        }
        for (int j = 0; j <= num_weight - 1; j++) {
            double weight = (*FreeWeight)[i][j].weight;
            int weightpos = (*FreeWeight)[i][j].position;
            result(weightpos) += (weight / sum_weight);
        }
    }

    for (int i = 0; i <= NFpara - 1; i++) {
        int position = (*Free)(i);
        double pftmp = (*Para_filtered)(position);
        double ptofd = PtoFderivative(pftmp, gbeta, gita);
        result(position) = result(position) * ptofd * gradients(position);
    }

    return result;

}