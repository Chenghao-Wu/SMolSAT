/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef BOND_AUTOCORRELATION_FUNCTION
#define BOND_AUTOCORRELATION_FUNCTION

#include <sstream>
#include <string>

#include "multibody_analysis.h"
#include "coordinate.h"
#include "trajectories.h"

using namespace std;

class Bond_Autocorrelation_Function: public Multibody_Analysis
{
    int n_times;
    float * baf;
    int * weighting;
    float * timetable;
    void initialize(std::shared_ptr<System>);
    void initialize(std::shared_ptr<System>, string dimselect);  
    void initialize(std::shared_ptr<System>, int l_type, string dimselect); 
    int atomcount;
    
    Coordinate dimensions;
    
    typedef Coordinate(Bond_Autocorrelation_Function::*vector_prep)(Coordinate)const;
    vector_prep vprep;
    
    Coordinate prep_inplane(Coordinate)const;
    Coordinate prep_outofplane(Coordinate)const;
    
    typedef float(Bond_Autocorrelation_Function::*polynomial_choice)(float)const;
    polynomial_choice legendre_p;
    
    float legendre_1(float val)const;
    float legendre_2(float val)const;


    
    
  public:
    Bond_Autocorrelation_Function();			//default constructor
    Bond_Autocorrelation_Function(const Bond_Autocorrelation_Function &);		//copy constructor
    Bond_Autocorrelation_Function operator = (const Bond_Autocorrelation_Function &);	//assignment
    ~Bond_Autocorrelation_Function();
    
    Bond_Autocorrelation_Function(std::shared_ptr<System>);
    Bond_Autocorrelation_Function(std::shared_ptr<System>, string dim);
    
    Bond_Autocorrelation_Function(std::shared_ptr<System>, int l_type, string dim);
    
    
    //Analysis_Type what_are_you(){Analysis_Type type = gyration_radius; return type;};		//virtual method to report the type of analysis
    
    void set(std::shared_ptr<System> sys){initialize(sys);};
    
    void analyze(Multibody_List * mblist);
    void list_displacementkernel(int,int,int);
    void listkernel(Multibody *, int, int, int);
    void postprocess_list();
    
    void write(string) const;
    void write(ofstream&) const;
    
    void run(Trajectories trjs,string listname);
    //void bin_hook(Trajectory_List*,int,int,int);
    //void postprocess_bins();

//	bool isThreadSafe(){return true;};
};

void export_Bond_Autocorrelation_Function(pybind11::module& m);

#endif
