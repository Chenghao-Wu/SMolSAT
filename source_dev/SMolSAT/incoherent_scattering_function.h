/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef INCOHERENT_SCATTERING_FUNCTION
#define INCOHERENT_SCATTERING_FUNCTION

#include "correlation_2d.h"
#include "wave_vectors.h"
#include "trajectories.h"

using namespace std;

class Incoherent_Scattering_Function:public Correlation_2D
{
    /*computational members - just used in calculations; no useful value later on*/
    int currenttime, nexttime, timegap;

    bool fullblock;
    
  public:
    Incoherent_Scattering_Function();		//default constructor
    Incoherent_Scattering_Function(const Incoherent_Scattering_Function &);		//copy constructor
    Incoherent_Scattering_Function operator=(const Incoherent_Scattering_Function &);
    //~Incoherent_Scattering_Function();
    
    Incoherent_Scattering_Function(std::shared_ptr<System> sys, std::shared_ptr<Wave_Vectors> wv, bool fblock=0);
    Incoherent_Scattering_Function(std::shared_ptr<System> sys, std::shared_ptr<Wave_Vectors> wv, int inner, int outer, bool fblock = 0);
    
    Analysis_Type what_are_you(){Analysis_Type type = incoherent_scattering_function; return type;};		//virtual method to report the type of analysis
    
    
     void analyze(Trajectory_List *,Trajectory_List *){cout<<"Error: Trajectory list targets with two lists not implemented for this analysis method.\n";}; //analysis method for when two trajectory lists are needed
    void analyze(Trajectory_List * t_list);
    void list_displacementkernel(int,int,int);
    void listkernel(Trajectory *);
    void listkernel(Trajectory *, int, int, int);
    
    void bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii);
    void postprocess_bins();
    void write(string filename)const;
    void write(ofstream& output)const;

    void run(Trajectories cl,string listname);
    //bool isThreadSafe(){return true;};
};

void export_Incoherent_Scattering_Function(pybind11::module& m);


#endif