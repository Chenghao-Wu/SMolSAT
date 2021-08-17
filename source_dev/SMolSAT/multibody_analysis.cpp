/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/
#include "multibody_analysis.h"

using namespace std;
namespace py = pybind11;

Multibody_Analysis::Multibody_Analysis()
{
  system = 0;
  multibody_list = 0;
}

Multibody_Analysis::Multibody_Analysis(const Multibody_Analysis & copy)
{
  system = copy.system;
  multibody_list = copy.multibody_list;
}


Multibody_Analysis Multibody_Analysis::operator=(const Multibody_Analysis & copy)
{
  if(this!=&copy)
  {
    system = copy.system;
    multibody_list = copy.multibody_list;
  }
  
  return *this;
}

void export_Multibody_Analysis(py::module& m)
    {
    py::class_<Multibody_Analysis, std::shared_ptr<Multibody_Analysis> >(m,"Multibody_Analysis")
    ;
    }
