/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "value_list.h"

using namespace std;

namespace py = pybind11;


void export_Value_List(py::module& m)
    {
    py::class_<Value_List <float>, std::shared_ptr<Value_List <float>> /* <--- trampoline*/ >(m,"Value_List")
    ;
    }