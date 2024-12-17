/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef VAN_HOVE_DISTINCT
#define VAN_HOVE_DISTINCT

#include "space-time_correlation_function.h"
#include "trajectory_list_bins.h"
#include "trajectories.h"

using namespace std;

class Van_Hove_Distinct: public Space_Time_Correlation_Function
{ 
    Trajectory_List_Bins * binned_trajectories;
    
    Trajectory_List * currentlist0;
    Trajectory_List * currentlist1;
    
    bool use_binned;
    
    int nx, ny, nz;		//number of bins in the x, y, and z dimensions

  public:
    Van_Hove_Distinct();
    
    Van_Hove_Distinct(std::shared_ptr<System>sys, Trajectory_List_Bins binnedtraj, int bin_count, float value_max=0);
    Van_Hove_Distinct(std::shared_ptr<System>sys, int bin_count, float value_max=0);
    void set(std::shared_ptr<System>sys, int bin_count, float value_max=0);
    
    
    void analyze(Trajectory_List*);
    void analyze(Trajectory_List*, Trajectory_List*);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory*, int, int, int);
    void listkernel2(Trajectory*, Trajectory*, int, int, int);
    
    void write(string filename)const;
    void write(ofstream& output)const;

    void run(Trajectories cl,string listname1,string listname2);

};

void export_Van_Hove_Distinct(pybind11::module& m);

#endif