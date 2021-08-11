/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef RADIAL_DISRIBUTION_FUNCTION
#define RADIAL_DISRIBUTION_FUNCTION

#include "system.h"
#include <sstream>
#include "analysis_onetime.h"
#include "trajectories.h"

using namespace std;

class Radial_Distribution_Function: public Analysis_Onetime
{
    float max_distance;
    int n_bins;
    float bin_size;
    int n_times;
    float ** time_rdf;
    float * mean_rdf;
    int * n_atoms_i;
    int * n_atoms_j;
    
    
  public:
    Radial_Distribution_Function();			//default constructor    
    Radial_Distribution_Function(const Radial_Distribution_Function &);		//copy constructor
    Radial_Distribution_Function(std::shared_ptr<System> sys, int nbins, int timescheme, float maxdistance=0);
    
    Radial_Distribution_Function operator = (const Radial_Distribution_Function &);	//assignment
    
    //Radial_Distribution_Function operator+ (const Radial_Distribution_Function &);
    
    void set(std::shared_ptr<System> sys, int nbins, int timescheme, float maxdistance=0);
    
    Analysis_Type what_are_you(){Analysis_Type type = radial_distribution_function; return type;};		//virtual method to report the type of analysis
    
    void preprocess(){trajectory_list2=trajectory_list;};
    void timekernel(int timeii){timekernel2(timeii);};
    void timekernel2(int timeii);
    void listkernel(Trajectory *, int, int, int);
    void listkernel2(Trajectory *, Trajectory *, int, int, int);
    void postprocess_list();
    void bin(int, float);
    
    void write(string)const;
    void write(ofstream& output)const;
    
    void structure_factor(string,int);
    
    void run(Trajectories trjs,string listname);
    void run(Trajectories trjs,string listname1,string listname2);
    
//	bool isThreadSafe(){return true;};
};

void export_Radial_Distribution_Function(pybind11::module& m);

#endif
