#include "definition.h"

#include "definition.h"
#define PI 3.14159265


int main() {

    string save_position;
    Vector3i bind;
    Vector3d l;
    int MAX_ITERATION_EVO = 200;

    //-------------------------------------------bind2------------------------------------------
    //thick100

    Vector3d move_focus;
    move_focus << 0.5, 0.5, 0.0;
    save_position = "./thick400-diel2d5-phi0theta0-lam500-size2500-focus50-bind2-trymid/";
    bind << 2, 2, 2;
    l << 100.0, 100.0, 16.0;
    Evo_single(save_position, bind, l, MAX_ITERATION_EVO, move_focus);


    return 0;

}



