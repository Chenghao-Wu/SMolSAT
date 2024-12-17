/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/


#ifndef STATIC_TRAJECTORY_LIST
#define STATIC_TRAJECTORY_LIST

#include <iostream>
#include "trajectory_list.h"
#include "analysis_base.h"
#include "system.h"
#include "trajectory_set.h"

class Static_Trajectory_List: public Trajectory_List, public Analysis_Base
{

  public:
    Static_Trajectory_List();
    Static_Trajectory_List(std::shared_ptr<System> syst, int capacity=0);
    void reset(std::shared_ptr<System> syst, int capacity=0);

    void atomkernel(Trajectory * traj);
    
    void set(std::shared_ptr<System>  syst, Trajectory_Set * trajectory_set);		//initialize trajectory list based on trajectory set
};


#endif
