/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef VAN_HOVE_SELF
#define VAN_HOVE_SELF

#include "space-time_correlation_function.h"
#include "trajectories.h"

using namespace std;
	
class Van_Hove_Self: public Space_Time_Correlation_Function
{
    void initialize(std::shared_ptr<System> sys, int bin_count, float value_max);
    
    //void atom_list(int atomcount, int* species_index, int* molecule_index, int* atom_type, int* atom_index);


  public:
    Van_Hove_Self(); 
    Van_Hove_Self(std::shared_ptr<System> sys, int bin_count, float value_max=0);
    void set(std::shared_ptr<System> sys, int bin_count, float value_max=0);
    
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel (Trajectory* current_trajectory, int timegapii, int thisii, int nextii);

    void write(string filename)const;
    void write(ofstream& output)const;
    
    void run(Trajectories cl,string listname);

};

void export_Van_Hove_Self(pybind11::module& m);

#endif