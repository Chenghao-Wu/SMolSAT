/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "analysis_onetime.h"
#include "system.h"

using namespace std;
namespace py=pybind11;

Analysis_Onetime::Analysis_Onetime():Analysis_Base()
{
  time_scheme=-1;
}


Analysis_Onetime::Analysis_Onetime(const Analysis_Onetime& copy):Analysis_Base(copy)
{
  time_scheme=copy.time_scheme;
}





void Analysis_Onetime::analyze(Trajectory_List * t_list)
{
  int timeii;
  
  trajectory_list=t_list;

  preprocess();
  
  if(time_scheme==-1)
  {
    for (timeii=0; timeii<system->show_n_timesteps();timeii++)
    {
      timekernel(timeii);
    }
  }
  else if(time_scheme<-1)
  {
    timekernel(-time_scheme-2);
  }
  else
  {
    for (timeii=time_scheme; timeii<system->show_n_exponentials();timeii+=system->show_n_exponential_steps())
    {
      timekernel2(timeii);
    }
  }
  
  postprocess_list();
}



void Analysis_Onetime::analyze(Trajectory_List * t_list,Trajectory_List* t_list2)
{
  int timeii;
  
  trajectory_list=t_list;
  trajectory_list2=t_list2;
  
  preprocess2();
  
  if(time_scheme==-1)
  {
    for (timeii=0; timeii<system->show_n_timesteps();timeii++)
    {
      timekernel2(timeii);
    }
  }
  else if(time_scheme<-1)
  {
    timekernel2(-time_scheme-2);
  }
  else
  {
    for (timeii=time_scheme; timeii<system->show_n_exponentials();timeii+=system->show_n_exponential_steps())
    {
      timekernel2(timeii);
    }
  }
  
  postprocess_list();
}



int Analysis_Onetime::determine_n_times()
{
  int timeii;
  int timecount = 0;
  
  if(time_scheme==-1)
  {
    for (timeii=0; timeii<system->show_n_timesteps();timeii++)
    {
      timecount++;
    }
  }
  else if(time_scheme<-1)
  {
      timecount=1;
  }
  else
  {
    for (timeii=time_scheme; timeii<system->show_n_exponentials();timeii+=system->show_n_exponential_steps())
    {
      timecount++;
    }
  }
  
  return timecount;
}

float * Analysis_Onetime::determine_times()
{
  float * times;
  int n_times=determine_n_times();
  times=new float [n_times];
  int timeii;
  if(time_scheme==-1)
  {
    times = system->show_timelist();
  }
  else
  {
    //cout<<(1*system->show_n_exponential_steps()+time_scheme)*system->show_time_unit()<<endl;
    for (timeii=time_scheme; timeii<system->show_n_exponentials();timeii+=system->show_n_exponential_steps())
    {
      times[timeii]=system->show_time(timeii*system->show_n_exponential_steps()+time_scheme);
    }
  }
  return times;
}


int Analysis_Onetime::system_time(int timeindex)
{
  if(time_scheme==-1)
  {
    return timeindex;
  }
  else
  {
    return timeindex*system->show_n_exponential_steps()+time_scheme;
  }
}

void export_Analysis_Onetime(py::module& m)
    {
    py::class_<Analysis_Onetime, std::shared_ptr<Analysis_Onetime> >(m,"Analysis_Onetime",py::base<Analysis_Base>())
    ;
    }