/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

//#include "system.h"
#include "trajmath.h"
#include "trajectory.h"
#include <Eigen/Dense>


#ifndef RGTENSOR
#define RGTENSOR

using namespace std;
  
class System;

typedef float threefloat [3];
  
class RgTensor 
{    
  private:
    
    std::shared_ptr<System> system;
    Trajectory * trajectory;		//trajectory for which the rg tensor is calculated
    sixfloat ** rgtensor;		//Rg tensor of trajectory as function of time, determined separately for the trajectory starting from the beginning of each block: rgtensor[blockii][timeii][index]
    int n_blocks; 
    int * blocksize;
    threefloat ** eigenvalues; 		//Eigenvalues of Rg tensor as a function of time, determined separately for the trajectory starting from the beginning of each block: eigenvalues[blockii][index];
    threefloat * mean_eigenvalues;
    int maxtimes;
    int * blocks_per_time;
    
    void calc_tensor();
    void calc_eigenvalues();
    
    void eigkernel(Eigen::Matrix3f rgarray, int blockii,int timeii);
  public:
    RgTensor(std::shared_ptr<System> sys, Trajectory* traj);
    
    void reset(Trajectory*traj);
    
    void show_eigenvalues(int timeii,threefloat eigout){eigout=mean_eigenvalues[timeii];};
    float show_eigenvalues(int timeii,int ii){return mean_eigenvalues[timeii][ii];};
    void write(string);
    void write(ofstream&)const;

    
};


#endif