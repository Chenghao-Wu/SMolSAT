/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "van_hove_distinct.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;
namespace py=pybind11;

/*----------------------------------------------------------------------------------------------*/

/*Constructor*/

Van_Hove_Distinct::Van_Hove_Distinct()
{
  system = 0;
  n_bins = 0;
  max_value = 0;
  n_times = 0;
  bin_size = 0;
  timetable = 0;
  correlation = new float * [1];
  weighting = new int [0];
  correlation[0]=new float [0];
}



/*----------------------------------------------------------------------------------------------*/

Van_Hove_Distinct::Van_Hove_Distinct(std::shared_ptr<System>sys, Trajectory_List_Bins binnedtraj, int bin_count, float value_max)
{
  int timeii, binii;
  
  system=sys;
  n_bins=bin_count;
  
  if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
  else max_value=value_max;
  
  bin_size = (max_value)/float(n_bins);
  
  n_times = system->show_n_timegaps();;

  timetable = system->displacement_times();
  
  use_binned = true;
  
  trajectory_list = new Trajectory_List [2];
//  currentlists = new Trajectory_List [2];
  
  correlation = new float * [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
    correlation[timeii]=new float [n_bins];
    for(binii=0;binii<n_bins;binii++)
    {
      correlation[timeii][binii]=0;
    }
  }
}


Van_Hove_Distinct::Van_Hove_Distinct(std::shared_ptr<System>sys, int bin_count, float value_max)
{
  int timeii, binii;
  
  system=sys;
  n_bins=bin_count;
  
  if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
  else max_value=value_max;
  
  bin_size = (max_value)/float(n_bins);
  
  n_times = system->show_n_timegaps();;

  timetable = system->displacement_times();
  
  use_binned = false;
  
  trajectory_list = new Trajectory_List [2];
//  currentlists = new Trajectory_List [2];
 
  correlation = new float * [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
    correlation[timeii]=new float [n_bins];
    for(binii=0;binii<n_bins;binii++)
    {
      correlation[timeii][binii]=0;
    }
  }
  
}



void Van_Hove_Distinct::set(std::shared_ptr<System>sys, int bin_count, float value_max)
{
  int timeii, binii;
  
  delete [] trajectory_list;
//  delete [] currentlists;
  
  for(timeii=0;timeii<n_times;timeii++)
  {
    delete [] correlation[timeii];
  }
  
  delete [] correlation;
  delete [] weighting;
  
  system=sys;
  n_bins=bin_count;
  
  if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
  else max_value=value_max;
  
  bin_size = (max_value)/float(n_bins);
  
  n_times = system->show_n_timegaps();;

  timetable = system->displacement_times();
  
  use_binned = false;
  
  
  
  trajectory_list = new Trajectory_List [2];
//  currentlists = new Trajectory_List [2];
  
  correlation = new float * [n_times];
  weighting = new int [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]=0;
    correlation[timeii]=new float [n_bins];
    for(binii=0;binii<n_bins;binii++)
    {
      correlation[timeii][binii]=0;
    }
  }
 
  
}


void Van_Hove_Distinct::analyze(Trajectory_List * t_list1)
{

  trajectory_list=t_list1;
  trajectory_list2=t_list1;
 
  system->displacement_list(this);
  postprocess_list();
}


void Van_Hove_Distinct::analyze(Trajectory_List * t_list1, Trajectory_List * t_list2)
{
  trajectory_list=t_list1;
  trajectory_list2=t_list2;
  
  system->displacement_list(this);
  postprocess_list();
  
}




void Van_Hove_Distinct::list_displacementkernel(int timegapii, int thisii, int nextii)
{
 
    weighting[timegapii]+=trajectory_list->show_n_trajectories(thisii);
//    current_list1->listloop(this,timegapii, thisii, nextii);
    trajectory_list->listloop(this,timegapii, thisii, nextii);
}


void Van_Hove_Distinct::listkernel(Trajectory* traj1, int timegapii, int thisii, int nextii)
{
  trajectory_list2->listloop2(this, traj1, timegapii, thisii, nextii);
}




void Van_Hove_Distinct::listkernel2(Trajectory* traj1, Trajectory* traj2, int timegapii, int thisii, int nextii)
{
  float distance;
  if(traj1!=traj2)
  {
    distance=(traj2->show_coordinate(nextii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size());	//calculate shortest distance between two coordinates, taking into account periodic boundaries
    bin(timegapii,distance);
    //cout<<(traj2->show_coordinate(nextii)-(traj1->show_coordinate(thisii))).length_unwrapped(system->size())<<endl;
  }
  else
  {
    weighting[timegapii]--;
  }

}
void Van_Hove_Distinct::write(string filename)const
{
  int timeii;
  int binii;
  
  ofstream output (filename.c_str());		//open correlation file
  
  output << "Van Hove Self Function created by SMolDAT v." << VERSION << "\n"; 
  output << n_bins << " bins\n";
  output << n_times << " times\n\n";
  
  cout << "\nWriting Van Hove Self Function to file " <<filename<<"." ;
  
  output << "\t";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << bin_size/2+float(binii)*bin_size << "\t";		//write bins at this time to file
  }
  output << "\n";
  
  for(timeii=0;timeii<n_times;timeii++)
  {
   output << timetable[timeii] << "\t";
    for(binii=0;binii<n_bins-1;binii++)
    {
      output << correlation[timeii][binii] << "\t";		//write bins at this time to file
    }
    output << correlation[timeii][n_bins-1] << "\n";		//write last bin at this time to file
  }
  output.close();					//close file  exit(0);
}


void Van_Hove_Distinct::write(ofstream& output)const
{
  int timeii;
  int binii;
  
  output << "Van Hove Self Function created by SMolDAT v." << VERSION << "\n"; 
  output << n_bins << " bins\n";
  output << n_times << " times\n\n";
  
  cout << "\nWriting Van Hove Self Function to file." ;
  
  output << "\t";
  
  for(binii=0;binii<n_bins;binii++)
  {
    output << bin_size/2+float(binii)*bin_size << "\t";		//write bins at this time to file
  }
  output << "\n";
  
  for(timeii=0;timeii<n_times;timeii++)
  {
   output << timetable[timeii] << "\t";
    for(binii=0;binii<n_bins-1;binii++)
    {
      output << correlation[timeii][binii] << "\t";		//write bins at this time to file
    }
    output << correlation[timeii][n_bins-1] << "\n";		//write last bin at this time to file
  }
}

void Van_Hove_Distinct::run(Trajectories cl,string listname1,string listname2)
{
  system->boxify();

  cout << "\nCalculating distinct part of Van Hove correlation function.\n";cout.flush();
  start = time(NULL);
  analyze(cl.trajectories[listname1],cl.trajectories[listname2]);
  finish = time(NULL);
  cout << "\nCalculated distinct Van Hove in in " << finish-start<<" seconds.\n";
}


void export_Van_Hove_Distinct(py::module& m)
    {
    py::class_<Van_Hove_Distinct, std::shared_ptr<Van_Hove_Distinct> >(m,"Van_Hove_Distinct",py::base<Space_Time_Correlation_Function>())
    .def(py::init< std::shared_ptr<System>,int , float >())
    //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
    .def("run",&Van_Hove_Distinct::run)
    .def("write", static_cast<void (Van_Hove_Distinct::*)(string )const> (&Van_Hove_Distinct::write))
    .def("write", static_cast<void (Van_Hove_Distinct::*)(ofstream& )const> (&Van_Hove_Distinct::write))
    ;
    }