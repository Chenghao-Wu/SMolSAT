/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/


#include "coordinate.h"
#include "trajectory.h"
#include "atom_trajectory.h"
#include "molecule.h"
#include "multibody.h"
#include "system.h"
#include "analysis_base.h"
#include "analysis_onetime.h"
#include "correlation_2d.h"
#include "multibody_analysis.h"
#include "trajectories.h"
#include "value_list.h"
#include "mean_square_displacement.h"
#include "displacement_list.h"
#include "bond_autocorrelation_function.h"
#include "rgtensor.h"
#include "rgtensor_stats.h"
#include "wave_vectors.h"
#include "structure_factor.h"
#include "radial_distribution_function.h"
#include "non_gaussian_parameter.h"
#include "incoherent_scattering_function.h"
#include "space-time_correlation_function.h"
#include "van_hove_self.h"
#include "van_hove_distinct.h"

#include "../extern/pybind11/include/pybind11/pybind11.h"
#include "../extern/pybind11/include/pybind11/stl_bind.h"

//! Create the python module
/*! each class sets up its own python exports in a function export_ClassName
    create the SMolSAT python module and define the exports here.
*/
using namespace std;

PYBIND11_MODULE(_SMolSAT, m)
{
    export_Coordinate(m);
    export_Trajectory(m);
    export_Atom_Trajectory(m);
    export_Molecule(m);
    export_Multibody(m);
    export_System(m);
    export_Analysis_Base(m);
    export_Analysis_Onetime(m);
    export_Correlation_2D(m);
    export_Multibody_Analysis(m);
    export_Value_List(m);
    export_Trajectories(m);
    export_Displacement_List(m);
    export_Mean_Square_Displacement(m);
    export_Bond_Autocorrelation_Function(m);
    export_RgTensor_Stats(m);
    export_Wave_Vectors(m);
    export_Structure_Factor(m);
    export_Radial_Distribution_Function(m);
    export_Non_Gaussian_Parameter(m);
    export_Incoherent_Scattering_Function(m);
    export_Space_Time_Correlation_Function(m);
    export_Van_Hove_Self(m);
    export_Van_Hove_Distinct(m);
}