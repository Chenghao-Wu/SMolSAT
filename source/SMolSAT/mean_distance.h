/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef MEANDISTANCE
#define MEANDISTANCE

#include "system.h"
#include <sstream>
#include "analysis_onetime.h"
#include "trajectories.h"

using namespace std;

class Mean_Distance: public Analysis_Base
{
    float max_distance;
    int n_bins;
    float bin_size;
    int n_times;
    float * time_m_sqr_dist;
    int * weighting;
    int * n_atoms_i;
    int * n_atoms_j;
    bool in_mole;
    
  public:
    Mean_Distance();			//default constructor    
    Mean_Distance(const Mean_Distance &);		//copy constructor
    Mean_Distance(std::shared_ptr<System> sys, bool in_mole=0);
    
    Analysis_Type what_are_you(){Analysis_Type type = mean_square_distance; return type;};		//virtual method to report the type of analysis

    void analyze (Trajectory_List*, Trajectory_List*);

    void listkernel(Trajectory* , int timegapii, int thisii, int nextii);
    void listkernel2(Trajectory* , Trajectory* ,int timegapii,int thisii, int nextii);

    void postprocess_list();
    
    void write(string);
    void write(ofstream& output);
    
    void run(Trajectories trjs,string listname1,string listname2);
    
//	bool isThreadSafe(){return true;};
};

void export_Mean_Distance(pybind11::module& m);

#endif
