/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include <math.h>
#include "mean_distance.h"
#include "version.h"

#define PI 3.1459265

using namespace std;
namespace py=pybind11;

Mean_Distance::Mean_Distance():Analysis_Base()
{
  max_distance=0;
  n_bins=0;
  bin_size=0;
  n_times=0;
  
  time_m_sqr_dist= new float [n_times];
  weighting = new int [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
}



Mean_Distance::Mean_Distance(std::shared_ptr<System> sys, bool inmole)
{
  int timeii, binii;
  float minboxsize;

  system=sys;
  in_mole=inmole;
  
  
  n_times = system->show_n_timegaps();
  
  time_m_sqr_dist= new float [n_times];
  weighting = new int [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_m_sqr_dist[timeii]=0;
    weighting[timeii]=0;
    n_atoms_i[timeii]=0;
    n_atoms_j[timeii]=0;
    
  }
}


Mean_Distance::Mean_Distance(const Mean_Distance & copy):Analysis_Base(copy)
{
  int timeii, binii;
  
  system=copy.system;
  
  n_times=copy.n_times;
  
  time_m_sqr_dist= new float [n_times];
  weighting = new int [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_m_sqr_dist[timeii]=copy.time_m_sqr_dist[timeii];
    weighting[timeii]=copy.weighting[timeii];
    n_atoms_i[timeii]=copy.n_atoms_i[timeii];
    n_atoms_j[timeii]=copy.n_atoms_j[timeii];
    
  }
}

void Mean_Distance::analyze(Trajectory_List* trajlist1, Trajectory_List* trajlist2)
{
  //n_atoms_i[timeii]=trajectory_list->show_n_trajectories(0);
  //n_atoms_j[timeii]=trajectory_list2->show_n_trajectories(0);
  int timeii;
  trajectory_list=trajlist1;
  trajectory_list2=trajlist2;
  for(timeii=0;timeii<n_times;timeii++)
  {
    trajectory_list->listloop(this,0, timeii, 0);
  }
  postprocess_list();
}

void Mean_Distance::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{ 
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}


void Mean_Distance::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii)
{
  float distance;
  //if(traj1!=traj2)
  if(in_mole)
  {
    if(traj1->show_moleculeID()==traj2->show_moleculeID())
    {
    if(traj1!=traj2)
    {
    distance=(traj2->show_unwrapped(thisii)-(traj1->show_unwrapped(thisii))).length();	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    time_m_sqr_dist[thisii]+=distance;
    weighting[thisii]+=1;
    }
    }
  }else{
    if(traj1!=traj2)
    {
    distance=(traj2->show_unwrapped(thisii)-(traj1->show_unwrapped(thisii))).length();	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    time_m_sqr_dist[thisii]+=distance;
    weighting[thisii]+=1;
    }
  }
}
 
 
 void Mean_Distance::postprocess_list()
 {
   int timeii;
   for(timeii=0;timeii<n_times;timeii++)
    {
      time_m_sqr_dist[timeii]/=weighting[timeii];
    }
 }
 
 void Mean_Distance::write(string filename)
 {
  int timeii;
  float * times;
  ofstream output (filename.c_str());

  cout << "\nWriting to file " <<filename<<".";cout.flush();

  /*Write first row - list of bin numbers*/
  output << "Mean squared distance data created by SMolDAT v." << VERSION << "\n";
  
  times = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    output << times[timeii]<<"\t"<<time_m_sqr_dist[timeii]<<"\n";
  }

  output.close();
 }
 
 
  void Mean_Distance::write(ofstream& output)
 {
  int timeii;
  float * times;
  cout << "\nWriting to file.";

  /*Write first row - list of bin numbers*/
  output << "Mean squared distance data created by SMolDAT v." << VERSION << "\n";

  times = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    output << times[timeii]<<"\t"<<time_m_sqr_dist[timeii]<<"\n";
  }

 }

void Mean_Distance::run(Trajectories cl,string listname, string listname2)
{
  

  cout << "\nCalculating mean squared distance.\n";cout.flush();
  start = time(NULL);
  analyze(cl.trajectories[listname],cl.trajectories[listname2]);
  finish = time(NULL);
  cout << "\nCalculated mean squared distance in " << finish-start<<" seconds.\n";

}

void export_Mean_Distance(py::module& m)
  {
  py::class_<Mean_Distance, std::shared_ptr<Mean_Distance> >(m,"Mean_Distance",py::base<Analysis_Base>())
  .def(py::init< std::shared_ptr<System>, bool>())
  //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
  .def("run",static_cast<void (Mean_Distance::*)(Trajectories, string, string )> (&Mean_Distance::run))
  .def("write", static_cast<void (Mean_Distance::*)(string )> (&Mean_Distance::write))
  .def("write", static_cast<void (Mean_Distance::*)(ofstream& )> (&Mean_Distance::write))
  ;
  }