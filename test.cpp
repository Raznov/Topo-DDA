#include "definition.h"

#include "definition.h"
#define PI 3.14159265


int main() {

    string save_position;
    Vector3i bind; 
    Vector3d l;
    int MAX_ITERATION_EVO = 200;

    //-------------------------------------------bind4------------------------------------------
    /*
    //thick100
    save_position = "./thick100-diel2d5-phi0theta0-lam500-size1000-focus50-bind4/";
    bind << 4, 4, 4;
    l << 40.0, 40.0, 4.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick100-diel2d5-phi0theta0-lam500-size1500-focus50-bind4/";
    bind << 4, 4, 4;
    l << 60.0, 60.0, 4.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick100-diel2d5-phi0theta0-lam500-size2000-focus50-bind4/";
    bind << 4, 4, 4;
    l << 80.0, 80.0, 4.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick100-diel2d5-phi0theta0-lam500-size2500-focus50-bind4/";
    bind << 4, 4, 4;
    l << 100.0, 100.0, 4.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);
    //thick200
    save_position = "./thick200-diel2d5-phi0theta0-lam500-size1000-focus50-bind4/";
    bind << 4, 4, 4;
    l << 40.0, 40.0, 8.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick200-diel2d5-phi0theta0-lam500-size1500-focus50-bind4/";
    bind << 4, 4, 4;
    l << 60.0, 60.0, 8.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick200-diel2d5-phi0theta0-lam500-size2000-focus50-bind4/";
    bind << 4, 4, 4;
    l << 80.0, 80.0, 8.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick200-diel2d5-phi0theta0-lam500-size2500-focus50-bind4/";
    bind << 4, 4, 4;
    l << 100.0, 100.0, 8.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);
    //thick300
    save_position = "./thick300-diel2d5-phi0theta0-lam500-size1000-focus50-bind4/";
    bind << 4, 4, 4;
    l << 40.0, 40.0, 12.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick300-diel2d5-phi0theta0-lam500-size1500-focus50-bind4/";
    bind << 4, 4, 4;
    l << 60.0, 60.0, 12.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick300-diel2d5-phi0theta0-lam500-size2000-focus50-bind4/";
    bind << 4, 4, 4;
    l << 80.0, 80.0, 12.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);
    */
    save_position = "./thick300-diel2d5-phi0theta0-lam500-size2500-focus50-bind4/";
    bind << 4, 4, 4;
    l << 100.0, 100.0, 12.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);
    //thick400
    save_position = "./thick400-diel2d5-phi0theta0-lam500-size1000-focus50-bind4/";
    bind << 4, 4, 4;
    l << 40.0, 40.0, 16.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick400-diel2d5-phi0theta0-lam500-size1500-focus50-bind4/";
    bind << 4, 4, 4;
    l << 60.0, 60.0, 16.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);

    save_position = "./thick400-diel2d5-phi0theta0-lam500-size2500-focus50-bind4/";
    bind << 4, 4, 4;
    l << 100.0, 100.0, 16.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO);



    return 0;

}



