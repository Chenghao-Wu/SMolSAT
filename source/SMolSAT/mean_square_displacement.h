/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef MEAN_SQUARE_DISPLACEMENT
#define MEAN_SQUARE_DISPLACEMENT

#include "system.h"
#include "trajectories.h"
#include <sstream>

#include "../extern/pybind11/include/pybind11/pybind11.h"

using namespace std;

class PYBIND11_EXPORT Mean_Square_Displacement: public Analysis_Base
{
    int n_times;
    float * msd;
    float * weighting;
    float * timetable;
    void initialize(std::shared_ptr<System> sys);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    float atomcount;
    
    /*internal calculation variables*/
    int currenttime, nexttime, currenttimegap;
    
    
  public:
    Mean_Square_Displacement();			//default constructor
    Mean_Square_Displacement(const Mean_Square_Displacement &);		//copy constructor
    Mean_Square_Displacement(std::shared_ptr<System> sys);
    Mean_Square_Displacement operator = (const Mean_Square_Displacement &);	//assignment
    
    Analysis_Type what_are_you(){Analysis_Type type = mean_square_displacement; return type;};		//virtual method to report the type of analysis
    
    float * normalized();
    void write(string)const;
    void write(ofstream&)const;
    void set(std::shared_ptr<System> sys){initialize(sys);};

    void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *, int, int, int);
    void postprocess_list();
    
    void bin_hook(Trajectory_List*,int,int,int);
    void postprocess_bins();
    
    float show(int t)const{return msd[t];};			//method to return one timestep of msd array
//	bool isThreadSafe(){return true;};

    void run(Trajectories cl,string listname);
};

void export_Mean_Square_Displacement(pybind11::module& m);

#endif
