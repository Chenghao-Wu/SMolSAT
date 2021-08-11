/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef MULTIBODY_ANALYSIS
#define MULTIBODY_ANALYSIS

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include "multibody_list.h"
#include "../extern/pybind11/include/pybind11/pybind11.h"

using namespace std;

class Multibody_Analysis
{
  protected:
    std::shared_ptr<System>  system;
    time_t start;			//timer start
    time_t finish;			//timer stop
    Multibody_List * multibody_list;
  
  public:
    Multibody_Analysis();
    Multibody_Analysis(const Multibody_Analysis &);
    Multibody_Analysis operator=(const Multibody_Analysis &);
    ~Multibody_Analysis(){};
    
    
    virtual void analyze(Multibody_List *){cout<<"Error: Multibody list targets not implemented for this analysis method.\n";};
    virtual void list_displacementkernel(int, int, int){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    virtual void listkernel(Multibody*,int,int,int){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    virtual void postprocess(){cout<<"Error: Multibody list targets not fully implemented for this analysis method.\n";};
    
    virtual void write(string)const{cout<<"Error: No standard write method implemented for this multibody_analysis method.\n";};
    virtual void write(ofstream&)const{cout<<"Error: No standard write method implemented for this multibody_analysis method.\n";};
    
    virtual bool isThreadSafe() {return false;};
};

void export_Multibody_Analysis(pybind11::module& m);

#endif