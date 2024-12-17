/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include <math.h>
#include "mean_square_distance.h"
#include "version.h"
#include <pybind11/stl.h>

#define PI 3.1459265

using namespace std;
namespace py=pybind11;

MeanSquared_Distance::MeanSquared_Distance():Analysis_Onetime()
{
  n_times=0;
  
  time_m_sqr_dist= new float [n_times];
  weighting = new int [n_times];

  time_rdf= new float * [n_times];
  n_atoms_i=new int [n_times];
  n_atoms_j = new int [n_times];
  
}



MeanSquared_Distance::MeanSquared_Distance(std::shared_ptr<System> sys,  bool isinter)
{
  int timeii, binii;
  
  is_inter=isinter;

  system=sys;
  
  time_scheme = -1;
  n_times=determine_n_times();

  time_m_sqr_dist= new float [n_times];
  weighting = new int [n_times];

  time_rdf= new float * [n_times];
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


MeanSquared_Distance::MeanSquared_Distance(const MeanSquared_Distance & copy):Analysis_Onetime(copy)
{
  int timeii, binii;
  
  system=copy.system;

  time_scheme = copy.time_scheme;
  n_times=copy.n_times;
  
  time_m_sqr_dist= new float [n_times];
  weighting = new int [n_times];

  time_rdf= new float * [n_times];
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

MeanSquared_Distance MeanSquared_Distance::operator=(const MeanSquared_Distance & copy)
{
  int timeii, binii;
  
  if(this!=&copy)
  {
  
  delete [] n_atoms_i;
  delete [] n_atoms_j;
  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] time_rdf[timeii];
  }
  
  delete [] time_rdf;
  delete [] time_m_sqr_dist;
  delete [] weighting;

  system=copy.system;
  
  time_scheme = copy.time_scheme;
  n_times=copy.n_times;
  
  time_m_sqr_dist= new float [n_times];
  weighting = new int [n_times];

  time_rdf= new float * [n_times];
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
  return *this;
}

void MeanSquared_Distance::set(std::shared_ptr<System> sys)
{
  int timeii, binii;
  
  delete [] n_atoms_i;
  delete [] n_atoms_j;
  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] time_rdf[timeii];
  }
  delete [] time_rdf;
  delete [] time_m_sqr_dist;
  delete [] weighting;

  system=sys;

  

  time_scheme = -1;
  n_times=determine_n_times();
  
  time_m_sqr_dist= new float [n_times];
  weighting = new int [n_times];

  time_rdf= new float * [n_times];
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


void MeanSquared_Distance::timekernel2(int timeii)
{
   n_atoms_i[timeii]=trajectory_list->show_n_trajectories(system_time(timeii));
   n_atoms_j[timeii]=trajectory_list2->show_n_trajectories(system_time(timeii));
   trajectory_list->listloop(this,0, timeii, 0);
}


void MeanSquared_Distance::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
  trajectory_list2->listloop2(this, current_trajectory, 0, thisii, 0);
}


void MeanSquared_Distance::listkernel2(Trajectory* traj1, Trajectory* traj2,int timegapii,int thisii, int nextii)
{
  float distance;
  //cout<<traj1->show_moleculeID()<<" "<<traj2->show_moleculeID()<<endl;
  if(is_inter)
  {
  if(traj1->show_moleculeID()==traj2->show_moleculeID())
  {
    if(traj1!=traj2)
    {
    distance=pow((traj2->show_unwrapped(thisii)-(traj1->show_unwrapped(thisii))).length(),2);	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    //cout<<distance<<endl;
    time_m_sqr_dist[thisii]+=distance;
    weighting[thisii]+=1;
    }
  }
  }
  else{
  if(traj1!=traj2)
  {
    distance=pow((traj2->show_unwrapped(thisii)-(traj1->show_unwrapped(thisii))).length(),2);	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    time_m_sqr_dist[thisii]+=distance;
    weighting[thisii]+=1;
  }
  }
  
}


 
void MeanSquared_Distance::postprocess_list()
{
  int timeii;
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_m_sqr_dist[timeii]/=weighting[timeii];
  }
}
 
 
 void MeanSquared_Distance::write(string filename)
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
 
 
  void MeanSquared_Distance::write(ofstream& output)
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

std::vector<double> MeanSquared_Distance::get_time_m_sqr_dist_vD() 
{
  int timeii;
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    time_m_sqr_dist_vD.push_back(time_m_sqr_dist[timeii]);
  }
  return time_m_sqr_dist_vD;
}


void MeanSquared_Distance::run(Trajectories cl,string listname1,string listname2)
{
  cout << "\nCalculating mean squared distance.\n";cout.flush();
  start = time(NULL);
  analyze(cl.find_trajectorylist(listname1),cl.find_trajectorylist(listname2));
  finish = time(NULL);
  cout << "\nCalculated mean squared distance in " << finish-start<<" seconds.\n";

}

 void export_MeanSquared_Distance(py::module& m)
    {
    py::class_<MeanSquared_Distance, std::shared_ptr<MeanSquared_Distance> >(m,"MeanSquared_Distance",py::base<Analysis_Onetime>())
    .def(py::init< std::shared_ptr<System>, bool >())
    //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
    .def("run",static_cast<void (MeanSquared_Distance::*)(Trajectories, string, string )> (&MeanSquared_Distance::run))
    .def("write", static_cast<void (MeanSquared_Distance::*)(string )> (&MeanSquared_Distance::write))
    .def("write", static_cast<void (MeanSquared_Distance::*)(ofstream& )> (&MeanSquared_Distance::write))
    .def("get", &MeanSquared_Distance::get_time_m_sqr_dist_vD)  // Remove incorrect static_cast
    ;
    }