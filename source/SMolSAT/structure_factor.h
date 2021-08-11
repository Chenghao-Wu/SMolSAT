/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/
#include <complex>
#include <string>
#include "analysis_base.h"
#include "wave_vectors.h"
#include "trajectory_list.h"
#include "system.h"
#include "analysis_onetime.h"
#include "trajectories.h"

using namespace std;

class Structure_Factor:public Analysis_Onetime
{
    //Trajectory_List * trajlist1;
    //Trajectory_List * trajlist2;
    float * structure_factor;
    int n_wavenumbers;
    
    int n_atoms;
    int currenttime;
    
    complex<float> ** wavedensity1;
    complex<float> ** wavedensity2;
    
    complex<float> ** current_wavedensity;
    
    void atomkernel(int,int,int,int){};
    
    Wave_Vectors wavevectors;
    
    void preprocess();
    void preprocess2();
    void timekernel (int timeii);	//method that is looped over by analyze with single trajectory list
    void timekernel2 (int timeii);	//method that is looped over by analyze with two trajectory lists
    void postprocess_list();
    
    
  public:
    Structure_Factor();
    ~Structure_Factor();
    Structure_Factor(const Structure_Factor &);
    Structure_Factor operator = (const Structure_Factor &);
    Structure_Factor(std::shared_ptr<System> sys, string plane,float max_length_scale, int timescheme = -1);

    void analyze(Trajectory_List * t_list);
    void analyze(Trajectory_List * t_list,Trajectory_List* t_list2);

    void analyze_wave_density(Trajectory_List * t_list);
    void listkernel(Trajectory* current_trajectory);
    void write(string filename)const;
    void write(ofstream& output)const;
    //bool isThreadSafe(){return true;};    

    void run(Trajectories cl,string listname);
    void run(Trajectories cl,string listname1,string listname2);
};
void export_Structure_Factor(pybind11::module& m);