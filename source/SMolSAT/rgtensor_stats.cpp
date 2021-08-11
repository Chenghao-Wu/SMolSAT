/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/



/*Note that currently this class can only use the trajectory list run method, and that it will only use the 
list of trajectories at the first time stored in the trajectory list; it is not currently implemented to work
with time variant trajectory lists*/

#include "rgtensor_stats.h"
#include "version.h"
#include <math.h>
#include "system.h"

using namespace std;
namespace py=pybind11;

RgTensor_Stats::RgTensor_Stats(std::shared_ptr<System> sys)
{
  int timeii;
  system=sys;
  
  n_times = system->show_n_timegaps();
  
  eigenvalues = new threefloat [n_times]; 
  mean_gyration_radius = new float [n_times];
  mean_relative_asphericity = new float [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    eigenvalues[timeii][0]=0;
    eigenvalues[timeii][1]=0;
    eigenvalues[timeii][2]=0;
    
    mean_gyration_radius[timeii] = 0;
    mean_relative_asphericity[timeii] = 0;
  }
  
  traj_index = 0;
  
}
 
void RgTensor_Stats::analyze(Trajectory_List* trajlist)
{
  int timeii, trajii;
  trajectory_list=trajlist;
  n_trajectories = trajectory_list->show_n_trajectories(0);
  
  single_eigenvalues = new threefloat * [n_times];
  single_relative_asphericity = new float * [n_times];
  single_gyration_radius = new float * [n_times+1];
  for(timeii=0;timeii<n_times;timeii++)
  {
    single_eigenvalues[timeii] = new threefloat [n_trajectories];
    single_relative_asphericity[timeii] = new float [n_trajectories];
    single_gyration_radius[timeii] = new float [n_trajectories];
    for(trajii=0;trajii<n_trajectories;trajii++)
    {
      single_eigenvalues[timeii][trajii][0]=0;
      single_eigenvalues[timeii][trajii][1]=0;
      single_eigenvalues[timeii][trajii][2]=0;
      
      single_relative_asphericity[timeii][trajii] = 0;
      single_gyration_radius[timeii][trajii] = 0;
    }
  }
  
  trajectory_list->listloop(this,0);
  
}

void RgTensor_Stats::listkernel(Trajectory* trajectory)
{
  int timeii;
  float temp1, temp2, temp3;		//temporary storage of eigenvalues
  float gyration_radius_sq;
  
  RgTensor rgtensor(system,trajectory);
  
  for(timeii=1;timeii<n_times;timeii++)
  {   
    temp1=rgtensor.show_eigenvalues(timeii,0);
    temp2=rgtensor.show_eigenvalues(timeii,1);
    temp3=rgtensor.show_eigenvalues(timeii,2);
    
    eigenvalues[timeii][0]+=temp1/trajectory_list->show_n_trajectories(0);
    eigenvalues[timeii][1]+=temp2/trajectory_list->show_n_trajectories(0);
    eigenvalues[timeii][2]+=temp3/trajectory_list->show_n_trajectories(0);
    
    gyration_radius_sq=temp1+temp2+temp3;
    
    single_eigenvalues[timeii][traj_index][0]=temp1;
    single_eigenvalues[timeii][traj_index][1]=temp2;
    single_eigenvalues[timeii][traj_index][2]=temp3;
    
    single_gyration_radius[timeii][traj_index]=pow(double(gyration_radius_sq),double(0.5));
    single_relative_asphericity[timeii][traj_index]=(temp3-.5*temp2-.5*temp1)/gyration_radius_sq;
    mean_gyration_radius[timeii]+=single_gyration_radius[timeii][traj_index]/n_trajectories;
    
  }
  traj_index++;
}


void RgTensor_Stats::calc_rel_asphericity_dist(int nbins)
{
  
  int timeii,trajii,binii;
  
  rel_asphericity_bins = nbins;
  
  rel_asphericity_dist = new float * [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    rel_asphericity_dist[timeii] = new float [rel_asphericity_bins];
    
    for(binii=0;binii<rel_asphericity_bins;binii++)
    {
      rel_asphericity_dist[timeii][binii]=0;
    }
    #ifdef never
    for(trajii=0;trajii<n_times;trajii++)
    {
      rel_asphericity_dist[timeii][int(single_relative_asphericity[timeii][trajii]*float(rel_asphericity_bins))]+=1/float(n_trajectories);	//increment appropriate bin count
    }
    
    #endif
  }
  
}


void RgTensor_Stats::calc_gyration_rad_dist(int nbins, float max_rad)
{
  
  int timeii,trajii,binii;
  
  gyration_rad_bins = nbins;
  max_gyr_rad_binned=max_rad;
  
  gyration_rad_dist = new float * [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    gyration_rad_dist[timeii] = new float [gyration_rad_bins];
    
    for(binii=0;binii<gyration_rad_bins;binii++)
    {
      gyration_rad_dist[timeii][binii]=0;
    }
    #ifdef never
    for(trajii=0;trajii<n_times;trajii++)
    {
      binii=int(single_gyration_radius[timeii][trajii]*float(gyration_rad_bins)/max_gyr_rad_binned);
      if(binii>gyration_rad_bins){binii=gyration_rad_bins;};
      gyration_rad_dist[timeii][binii]+=1/float(n_trajectories);	//increment appropriate bin count
    }
    #endif
  }
}


void RgTensor_Stats::write_rel_asphericity_dist(string filename)
{
  int timeii, binii;
  float * times;
  
  ofstream output(filename.c_str());
  cout << "Writing gyration tensor relative asphericity distribution to file " << filename << ".\n";                                                                                                      
  output << "Gyration tensor data created by AMDAT v." << VERSION << "\n";  
  
  times = system->displacement_times();
  
  for(binii=0;binii<rel_asphericity_bins;binii++)
  {
    output << "\t" << (float(binii))/float(rel_asphericity_bins);
  }
  
  output << "\t" << 1 << "\n";
  
  for(timeii=1;timeii<n_times;timeii++)
  {
    output << times[timeii];
    for(binii=0;binii<rel_asphericity_bins;binii++)
    {
      output << "\t" << rel_asphericity_dist[timeii][binii];
    }
    output << "\n";
  }
  
}

void RgTensor_Stats::write_gyration_rad_dist(string filename)
{
  int timeii, binii;
  float * times;
  
  ofstream output(filename.c_str());
  cout << "Writing trajectory gyration radius distribution to file " << filename << ".\n";                                                                                                      
  output << "Gyration tensor data created by AMDAT v." << VERSION << "\n";  
  
  times = system->displacement_times();
  
  for(binii=0;binii<=gyration_rad_bins;binii++)
  {
    output << "\t" << (float(binii))*max_gyr_rad_binned/float(gyration_rad_bins);
  }
  
  output << "\toverflow\n";
  
  for(timeii=1;timeii<n_times;timeii++)
  {
    output << times[timeii];
    for(binii=0;binii<=gyration_rad_bins;binii++)
    {
      output << "\t" << gyration_rad_dist[timeii][binii];
    }
    output << "\n";
  }
  
}



void RgTensor_Stats::write(string filename)const
{
  int timeii;
  float * times;
  
  ofstream output(filename.c_str());
  cout << "Writing gyration tensor eigenvalues to file " << filename << ".\n";
  output << "Gyration tensor data created by AMDAT v." << VERSION << "\n";  
  
  times = system->displacement_times();
  
  for(timeii=1;timeii<n_times;timeii++)
  {
    output << times[timeii] << "\t" << eigenvalues[timeii][0] << "\t" << eigenvalues[timeii][1] << "\t" << eigenvalues[timeii][2] << "\n";
  }
}


void RgTensor_Stats::write(ofstream& output)const
{
  int timeii;
  float * times;
  
  cout << "Writing gyration tensor eigenvalues to file.\n";
  output << "Gyration tensor data created by AMDAT v." << VERSION << "\n";  
  
  times = system->displacement_times();
  
  for(timeii=1;timeii<n_times;timeii++)
  {
    output << times[timeii] << "\t" << eigenvalues[timeii][0] << "\t" << eigenvalues[timeii][1] << "\t" << eigenvalues[timeii][2] << "\n";
  }
}

void RgTensor_Stats::run(Trajectories trjs,string listname)
{
  cout << "\nCalculating Rg statistics.\n";cout.flush();
  start = time(NULL);
  analyze(trjs.trajectories[listname]);
  finish = time(NULL);
  cout << "\nCalculated Rg statistics in " << finish-start<<" seconds."<<endl;
}



void export_RgTensor_Stats(py::module& m)
    {
    py::class_<RgTensor_Stats, std::shared_ptr<RgTensor_Stats> >(m,"RgTensor_Stats",py::base<Analysis_Base>())
    .def(py::init< std::shared_ptr<System> >())
    //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
    .def("run",&RgTensor_Stats::run)
    .def("write", static_cast<void (RgTensor_Stats::*)(string )const> (&RgTensor_Stats::write))
    .def("write", static_cast<void (RgTensor_Stats::*)(ofstream& )const> (&RgTensor_Stats::write))
    .def("calc_rel_asphericity_dist",&RgTensor_Stats::calc_rel_asphericity_dist)
    .def("write_rel_asphericity_dist",&RgTensor_Stats::write_rel_asphericity_dist)
    .def("calc_gyration_rad_dist",&RgTensor_Stats::calc_gyration_rad_dist)
    .def("write_gyration_rad_dist",&RgTensor_Stats::write_gyration_rad_dist)
    ;
    }
