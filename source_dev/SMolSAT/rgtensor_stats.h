/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "analysis_base.h"
#include "rgtensor.h"
#include "trajectory_list.h"
#include "trajectories.h"

#ifndef RGTENSOR_STATS
#define RGTENSOR_STATS

using namespace std;

class RgTensor_Stats:public Analysis_Base
{
  float * mean_gyration_radius;
  threefloat * eigenvalues; // stores mean eigenvalues as a function of time
  threefloat ** single_eigenvalues; // stores single-trajectory eigenvalues as a function of time;
  float ** single_relative_asphericity; // stores single-trajectory relative asphericity as a function of time;
  float * mean_relative_asphericity;
  float ** single_gyration_radius; // stores single-trajectory gyration radius as a function of time;
  float ** rel_asphericity_dist;
  float ** gyration_rad_dist;
  
  int n_times;
  int n_trajectories;
  void atomkernel(int,int,int,int){};
  void displacementkernel(int,int,int,int,int,int,int){};
  
  bool rel_asphericity_dist_calc;
  bool gyration_rad_dist_calc;
  int rel_asphericity_bins;
  int gyration_rad_bins;
  float max_gyr_rad_binned;
  
  
  /*dummy variables*/
  int traj_index;	//variable to keep track of which member of the list we are on
  
  
  public:
    RgTensor_Stats(std::shared_ptr<System>);
    
    Analysis_Type what_are_you(){Analysis_Type type = rgtensor_stats; return type;};		//virtual method to report the type of analysis
    
    void analyze (Trajectory_List*);
    void listkernel(Trajectory*);
    
    void write(string)const;
    void write(ofstream&)const;
    void write_rel_asphericity_dist(string);
    void write_gyration_rad_dist(string);
    
    void calc_rel_asphericity_dist(int nbins);
    void calc_gyration_rad_dist(int,float);  
    
    void run(Trajectories trjs,string listname);

};

void export_RgTensor_Stats(pybind11::module& m);

#endif