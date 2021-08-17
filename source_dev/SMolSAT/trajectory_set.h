/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef TRAJECTORY_SET
#define TRAJECTORY_SET

#include "trajectory.h"

//namespace std {

class Trajectory_Set
{
    int n_trajectories;
    Trajectory * trajectories;

  public:
    Trajectory_Set();
    Trajectory_Set(const Trajectory_Set &);
    ~Trajectory_Set();
    Trajectory_Set operator =(const Trajectory_Set &);

    virtual Trajectory * show_trajectory(int index){return &trajectories[index];};
    virtual int show_n_trajectories()const{return n_trajectories;};
};

//}

#endif