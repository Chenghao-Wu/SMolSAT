/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef ANALYSIS_ONETIME
#define ANALYSIS_ONETIME

#include "analysis_base.h"

using namespace std;

class Analysis_Onetime: public Analysis_Base
{
  protected:
    int time_scheme;		//determines what times to loop over. If time_scheme = -1, loop over all times. If time_scheme is zero or positive, only use one time per block, with the value setting the time index offset from the beginning of the block
    
    virtual void preprocess(){};
    virtual void preprocess2(){};
    virtual void timekernel (int timeii){cout<<"Error: Trajectory list targets with one list not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed;	//method that is looped over by analyze with single trajectory list
    virtual void timekernel2 (int timeii){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed;	//method that is looped over by analyze with two trajectory lists
    int determine_n_times();
  public:
    
    Analysis_Onetime();
    Analysis_Onetime(const Analysis_Onetime&);
    
    virtual void analyze (Trajectory_List *);
    virtual void analyze (Trajectory_List *, Trajectory_List*);
    virtual void bin_hook(Trajectory_List *,int,int,int){cout<<"Error: Trajectory list bins not implemented for this analysis method.\n";};
    virtual void postprocess_bins(){};
    
    int system_time(int timeindex);	//return index of system time given index of internal time
};

void export_Analysis_Onetime(pybind11::module& m);

#endif