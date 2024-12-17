/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef MEAN_SQUARE_DISPLACEMENT_2D
#define MEAN_SQUARE_DISPLACEMENT_2D

#include "system.h"
#include "trajectories.h"
#include <sstream>

#include "../extern/pybind11/include/pybind11/pybind11.h"

#include <sstream>
using namespace std;

class PYBIND11_EXPORT Mean_Square_Displacement_2D: public Analysis_Base
{
    int n_times, atomcount;
    float * msd;
    int * weighting;
    float * timetable;
    string plane;
    typedef float (Trajectory::*length)(int,int)const;		
    length distancefun;
       
    /*calculation variables*/
    int currenttime, nexttime, currenttimegap;
    
    void initialize(std::shared_ptr<System> sys,string);
    
  public:
    Mean_Square_Displacement_2D();					// default constructor
    ~Mean_Square_Displacement_2D();					// destructor
    Mean_Square_Displacement_2D(std::shared_ptr<System> sys,string);
    Mean_Square_Displacement_2D(const Mean_Square_Displacement_2D &);	// copy constructor
    Mean_Square_Displacement_2D operator = (const Mean_Square_Displacement_2D &);	// assignment operator
       
  
    Analysis_Type what_are_you(){Analysis_Type type = mean_square_displacement_2d; return type;};		//virtual method to report the type of analysis
    
    void single_atom(int species_index, int molecule_index, int atom_type, int atom_index);
    float * normalized();
    void write(string)const;
    void write(ofstream&)const;
    
    void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
    void postprocess_list();
    
    void bin_hook(Trajectory_List*,int,int,int);
    void postprocess_bins();   
    
    float show(int t){return msd[t];};			//method to return one timestep of msd array

    void run(Trajectories cl,string listname);
};

void export_Mean_Square_Displacement_2D(pybind11::module& m);

#endif