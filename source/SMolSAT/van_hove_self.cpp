/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "van_hove_self.h"
#include <stdlib.h>
#include <iostream>
#include "version.h"

using namespace std;
namespace py=pybind11;

Van_Hove_Self::Van_Hove_Self(std::shared_ptr<System> sys, int bin_count, float value_max)
{
	initialize(sys, bin_count, value_max);
}



Van_Hove_Self::Van_Hove_Self()
{
	system = 0;
	n_bins = 0;
	max_value=0;
	bin_size=0;
	n_times=0;
	timetable=0;

	correlation = new float * [1];
	weighting = new int [0];
	correlation[0]=new float [0];
}



void Van_Hove_Self::set(std::shared_ptr<System> sys, int bin_count, float value_max)
{
	clear_memory();
	initialize(sys, bin_count, value_max);
}


void Van_Hove_Self::initialize(std::shared_ptr<System> sys, int bin_count, float value_max)
{
	int timeii;
	int binii;

	n_bins = bin_count;
	system = sys;

	if(value_max==0) max_value = (system->size().min())/2;	//if no max range given, set it to be half the minimum dimension of the box at the initial time.
	else max_value=value_max;

	bin_size = (max_value)/float(n_bins);

	n_times = system->show_n_timegaps();


	timetable = system->displacement_times();

	 //allocate memory for van hove self-correlation function and weighting and initialize to zero
	correlation = new float * [n_times];
	weighting = new int [n_times];
	for(timeii=0;timeii<n_times;timeii++)
	{
		correlation[timeii]=new float [n_bins];
		weighting[timeii] = 0;
		for(binii=0;binii<n_bins;binii++)
		{
			correlation[timeii][binii]=0;
		}
	}
}






void Van_Hove_Self::analyze(Trajectory_List * t_list)
{

	trajectory_list=t_list;

	system->displacement_list(this);
	postprocess_list();
}

void Van_Hove_Self::list_displacementkernel(int timegapii, int thisii, int nextii)
{

	//weighting[timegapii]+=trajectory_list[0].show_n_trajectories(thisii);
	(trajectory_list[0]).listloop(this, timegapii, thisii, nextii);
}

void Van_Hove_Self::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
	float distance = (current_trajectory->show_unwrapped(nextii)-current_trajectory->show_unwrapped(thisii)).length();
	bin(timegapii,distance);
	weighting[timegapii]++;


}

void Van_Hove_Self::write(string filename)const
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


void Van_Hove_Self::write(ofstream& output)const
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

void Van_Hove_Self::run(Trajectories cl,string listname)
{
  cout << "\nCalculating self part of Van Hove correlation function.\n";cout.flush();
  start = time(NULL);
  analyze(cl.trajectories[listname]);
  finish = time(NULL);
  cout << "\nCalculated self Van Hove in in " << finish-start<<" seconds.\n";
}

void export_Van_Hove_Self(py::module& m)
    {
    py::class_<Van_Hove_Self, std::shared_ptr<Van_Hove_Self> >(m,"Van_Hove_Self",py::base<Space_Time_Correlation_Function>())
    .def(py::init< std::shared_ptr<System>,int , float >())
    //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
    .def("run",&Van_Hove_Self::run)
    .def("write", static_cast<void (Van_Hove_Self::*)(string )const> (&Van_Hove_Self::write))
    .def("write", static_cast<void (Van_Hove_Self::*)(ofstream& )const> (&Van_Hove_Self::write))
    ;
    }