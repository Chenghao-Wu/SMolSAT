/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/
#ifndef MULTIBODY_SET
#define MULTIBODY_SET

#include <vector>

#include "multibody.h"
#include "trajectory.h"
#include "trajectory_set.h"
//#include "system.h"

class System;

class Multibody_Set: public Trajectory_Set
{
    int n_multibodies;
    Multibody * multibodies;

  public:
    Multibody_Set();					//default constructor
    Multibody_Set(const Multibody_Set &);		//copy constructor
    ~Multibody_Set();					//destructor
    Multibody_Set operator =(const Multibody_Set &);	//equality operator
    Multibody_Set(int multibody_count);			//constructor to set number of multibodies
    


    /*Methods to define one or more multibodies in the array of multibodies*/
    void set_multibody(int multibody_index, int n_bodies, Trajectory** bodies){multibodies[multibody_index].set(n_bodies,bodies);};
    void set_multibody(int multibody_index, const Multibody & multibody){multibodies[multibody_index]=multibody;};
    
    void set(System * system, int multibody_count);		//reset number of multibodies and reinitialize array of multibodies
    void set(vector<Multibody> mbodies);		//define set of multibodies based on vector of multibodies
    Multibody * show_multibody(int index);
    int show_n_multibodies(){return n_multibodies;};
    
    /*Method to compute COM or centroid trajectories for all multibodies in set*/
    void compute_trajectories(bool centertype, int traj_type);

    /*Replace trajectory_set methods to instead use the multibody list*/
    Trajectory * show_trajectory(int index){return &multibodies[index];};
    virtual int show_n_trajectories()const{return n_multibodies;};



};

#endif
