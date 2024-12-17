/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef MEANSQUAREDISTANCE
#define MEANSQUAREDISTANCE

#include "system.h"
#include <sstream>
#include "analysis_onetime.h"
#include "trajectories.h"

using namespace std;

class MeanSquared_Distance: public Analysis_Onetime
{
    float * time_m_sqr_dist;
    int * weighting;

    int n_times;
    float ** time_rdf;
    int * n_atoms_i;
    int * n_atoms_j;
    bool is_inter;
    
  public:
    MeanSquared_Distance();			//default constructor    
    MeanSquared_Distance(const MeanSquared_Distance &);		//copy constructor
    MeanSquared_Distance(std::shared_ptr<System> sys, bool is_inter=0);
    
    MeanSquared_Distance operator = (const MeanSquared_Distance &);	//assignment
    
    //MeanSquared_Distance operator+ (const MeanSquared_Distance &);
    
    void set(std::shared_ptr<System> sys);
    
    Analysis_Type what_are_you(){Analysis_Type type = radial_distribution_function; return type;};		//virtual method to report the type of analysis
    
    void preprocess(){trajectory_list2=trajectory_list;};
    void timekernel(int timeii){timekernel2(timeii);};
    void timekernel2(int timeii);
    void listkernel(Trajectory *, int, int, int);
    void listkernel2(Trajectory *, Trajectory *, int, int, int);
    void postprocess_list();
    void bin(int, float);
    
    void write(string);
    void write(ofstream& output);
    
    void run(Trajectories trjs,string listname);
    void run(Trajectories trjs,string listname1,string listname2);
    
//	bool isThreadSafe(){return true;};
};

void export_MeanSquared_Distance(pybind11::module& m);

#endif
