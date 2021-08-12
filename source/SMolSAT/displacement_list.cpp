/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "displacement_list.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "static_trajectory_list.h"

#include <omp.h>

using namespace std;
namespace py = pybind11;


Displacement_List::Displacement_List()
{
  timegap=0;
  n_times = 0;
  syst = 0;
  system=0;
  
  included = new Boolean_List[n_times];
  time_conversion = new int [0];
  defined_times = new bool [0];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(syst);
  }


}


Displacement_List::Displacement_List(const Displacement_List & copy)
{
  int timeii;
  
  syst=copy.syst;
  system = copy.system;
  timegap=copy.timegap;
  n_times = copy.n_times;

  included = new Boolean_List [n_times];
  time_conversion = new int [syst->show_n_timesteps()];
  defined_times = new bool [syst->show_n_timesteps()];
  values=copy.values;
  for(timeii=0;timeii<n_times;timeii++)
  {
	  included[timeii] = copy.included[timeii];
  }
  
  for(timeii=0;timeii<syst->show_n_timesteps();timeii++)
  {
    time_conversion[timeii] = copy.time_conversion[timeii];
    defined_times[timeii] = copy.defined_times[timeii];
  }

}



/** **/
Displacement_List::Displacement_List(std::shared_ptr<System> sys)
{
  int timeii;

  system = sys;
  syst=sys;
  n_times = syst->show_n_timesteps();
  timegap=0;
  included = new Boolean_List[n_times];
  values.resize(n_times);
  
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(syst);
    values[timeii].resize(system->show_n_trajectories());
  }
  
  time_conversion = new int [syst->show_n_timesteps()];
  defined_times = new bool [syst->show_n_timesteps()];
  for(int timeii=0;syst->show_n_timesteps();timeii++)
  {
    time_conversion[timeii]=timeii;
    defined_times[timeii]=1;
  }

}


/** **/
Displacement_List::Displacement_List(std::shared_ptr<System> sys,int t_gap)
{
  int timeii;
  timegap = t_gap;
  
  system = sys;
  syst=sys;
  n_times = syst->show_n_timesteps();
  
  included = new Boolean_List[n_times];
  
    values.resize(n_times);
  
  for(int timeii=0;timeii<n_times;timeii++)
  {
    included[timeii].set(syst);
    values[timeii].resize(system->show_n_trajectories());
  }
  
  time_conversion = new int [n_times];
  defined_times = new bool [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    time_conversion[timeii]=int(float(timeii - sys->show_frt())/float(system->show_n_exponential_steps()));
    defined_times[timeii]=1;
  }

}



Displacement_List Displacement_List::operator = (const Displacement_List & copy)
{
  int timeii;

  if(this!=&copy)
  {

    system = copy.system;
    syst=copy.syst;
    n_times = copy.n_times;
    timegap = copy.timegap;

    included = new Boolean_List [n_times];
    time_conversion = new int [syst->show_n_timesteps()];
    defined_times = new bool [syst->show_n_timesteps()];
    values=copy.values;
    for(timeii=0;timeii<n_times;timeii++)
    {
      included[timeii] = copy.included[timeii];
    }
  
    for(timeii=0;timeii<syst->show_n_timesteps();timeii++)
    {
      time_conversion[timeii] = copy.time_conversion[timeii];
      defined_times[timeii] = copy.defined_times[timeii];
    }

  }

  return *this;

}




/*Methods to do analysis using trajectory list*/

void Displacement_List::analyze(Trajectory_List * t_list)
{
  trajectory_list=t_list;
  system->displacement_list(this,timegap,false);
  postprocess_list();
}

void Displacement_List::list_displacementkernel(int timegapii,int thisii, int nextii)
{
  (trajectory_list[0]).listloop(this,timegapii, thisii, nextii);
}



void Displacement_List::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
  included[convert_time(thisii)](current_trajectory->show_trajectory_ID(),1);
  values[convert_time(thisii)][current_trajectory->show_trajectory_ID()]=current_trajectory->distance(thisii,nextii);
}


void Displacement_List::postprocess_list()
{
}



/*Method to write MSD data to file*/

void Displacement_List::write(string filename)const
{
  int timeii;

  cout << "\nWriting msd to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Displacement data created by SMolDAT v." << VERSION << "\n";
  output << "Timegap " << syst->displacement_times(timegap)<< "\n";
  output << "Mean_displacement " << mean() << "\n";
  output << "Mean_square_displacement " << power_mean(2) << "\n";
}


void Displacement_List::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting msd to file.";

  output << "Displacement data created by SMolDAT v." << VERSION << "\n";
  output << "Timegap " << syst->displacement_times(timegap)<< "\n";
  output << "Mean_displacement " << mean() << "\n";
  output << "Mean_square_displacement " << power_mean(2) << "\n";
}


void Displacement_List::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  trajectory_list=t_list;

  list_displacementkernel(timegapii, thisii, nextii);

}



void Displacement_List::postprocess_bins()
{
  postprocess_list();
}

void Displacement_List::run(Trajectories trjs,string displacement_listname,string listname)
{
  Displacement_List*dlist_pointer;
  dlist_pointer = new Displacement_List();

  cout << "Calculating list of displacement scalars.\n";
  start = time(NULL);

  analyze(trjs.trajectories[listname]);
  (*dlist_pointer)= (*this) ;
  trjs.add_value_list(dlist_pointer,displacement_listname);
  finish = time(NULL);
  cout << "\nCalculated list of displacement scalars in " << finish-start<<" seconds.";
}

void export_Displacement_List(py::module& m)
    {
    py::class_<Displacement_List, std::shared_ptr<Displacement_List> >(m,"Displacement_List",py::base<Value_List<float>>(),py::base<Analysis_Base>())
    .def(py::init< std::shared_ptr<System>, int >())
    //.def("analyze", static_cast<void (Displacement_List::*)(Trajectory_List* )> (&Displacement_List::analyze))
    .def("run",&Displacement_List::run)
    .def("write", static_cast<void (Displacement_List::*)(string )const> (&Displacement_List::write))
    .def("write", static_cast<void (Displacement_List::*)(ofstream& )const> (&Displacement_List::write))
    ;
    }