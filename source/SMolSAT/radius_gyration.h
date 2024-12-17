/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef RADIUSGYRATION
#define RADIUSGYRATION

#include <sstream>
#include <string>

#include "multibody_analysis.h"
#include "coordinate.h"
#include "trajectories.h"
#include "trajmath.h"

using namespace std;

class Radius_Gyration: public Multibody_Analysis
{
    int n_times;
    float * rg;
    int * weighting;
    float * timetable;

    vector<Trajectory*> trajectories;
    Coordinate * coordinates;
    
    void initialize(std::shared_ptr<System>);

    sixfloat * gytensor;

    int atomcount;
    
    Coordinate dimensions;
    
    typedef Coordinate(Radius_Gyration::*vector_prep)(Coordinate)const;

    
  public:
    Radius_Gyration();			//default constructor
    Radius_Gyration(const Radius_Gyration &);		//copy constructor
    Radius_Gyration operator = (const Radius_Gyration &);	//assignment
    ~Radius_Gyration();
    
    Radius_Gyration(std::shared_ptr<System>);
    
    
    
    //Analysis_Type what_are_you(){Analysis_Type type = gyration_radius; return type;};		//virtual method to report the type of analysis
    
    void set(std::shared_ptr<System> sys){initialize(sys);};
    
    void analyze(Multibody_List * mblist);
    void listkernel(Multibody *, int, int, int);
    void postprocess_list();
    
    void write(string) const;
    void write(ofstream&) const;

    void write_tensor(string) const;
    void write_tensor(ofstream&) const;
    
    void run(Trajectories trjs,string listname);
    //void bin_hook(Trajectory_List*,int,int,int);
    //void postprocess_bins();

//	bool isThreadSafe(){return true;};
};

void export_Radius_Gyration(pybind11::module& m);

#endif
