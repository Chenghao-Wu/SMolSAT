/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "trajectory_set.h"
#include "trajectory.h"

using namespace std;

Trajectory_Set::Trajectory_Set()
{
    n_trajectories=0;
    trajectories=new Trajectory [n_trajectories];
}

Trajectory_Set::Trajectory_Set(const Trajectory_Set & copy)
{
    int trajii;

    n_trajectories = copy.n_trajectories;
    trajectories=new Trajectory [n_trajectories];

    for(trajii=0;trajii<n_trajectories;trajii++)
    {
        trajectories[trajii]=copy.trajectories[trajii];
    }
}

Trajectory_Set::~Trajectory_Set()
{
    delete [] trajectories;
}


Trajectory_Set Trajectory_Set::operator=(const Trajectory_Set & copy)
{
    int trajii;

    if(this!=&copy)
    {
        delete [] trajectories;

        n_trajectories = copy.n_trajectories;
        trajectories=new Trajectory [n_trajectories];

        for(trajii=0;trajii<n_trajectories;trajii++)
        {
            trajectories[trajii]=copy.trajectories[trajii];
        }
    }
    
    return *this;
}
