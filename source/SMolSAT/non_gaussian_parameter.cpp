/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "non_gaussian_parameter.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"

using namespace std;
namespace py=pybind11;

Non_Gaussian_Parameter::Non_Gaussian_Parameter()
{
  system = 0;
  n_times = 0;
  
  ngp = new float [n_times];
  weighting = new long int [n_times];
  timetable = 0;
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
}



Non_Gaussian_Parameter::Non_Gaussian_Parameter(const Non_Gaussian_Parameter & copy)
{
  system = copy.system;
  n_times = copy.n_times;
  msd = copy.msd;
  
  ngp = new float [n_times];
  weighting = new long int [n_times];
  timetable = system -> displacement_times();
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]=copy.ngp[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
  atomcount = copy.atomcount;
}


Non_Gaussian_Parameter Non_Gaussian_Parameter::operator =(const Non_Gaussian_Parameter & copy)
{
 if(this!=&copy)
 {
  system = copy.system;
  n_times = copy.n_times;
  msd = copy.msd;
  
  ngp = new float [n_times];
  weighting = new long int [n_times];
  timetable = system -> displacement_times();
  for(int timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]=copy.ngp[timeii];
    weighting[timeii]=0;;
  }
  atomcount = copy.atomcount;
 }
 return *this;
}


Non_Gaussian_Parameter::Non_Gaussian_Parameter(std::shared_ptr<System> sys)
{
  int timeii;

  system = sys;

  Mean_Square_Displacement m(sys);

  msd=m;

  n_times = system->show_n_timegaps();

  //allocate memory for ngp data and msd data
  ngp = new float [n_times];
  weighting = new long int [n_times];
  timetable = system -> displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    ngp[timeii]=0;
    weighting[timeii]=0;;
  }
  atomcount = 0;
}


void Non_Gaussian_Parameter::analyze(Trajectory_List * t_list)
{

	trajectory_list=t_list;
	system->displacement_list(this,false);
	postprocess_list();
}



void Non_Gaussian_Parameter::list_displacementkernel(int timegapii,int thisii, int nextii)
{
	currenttime=thisii;
	nexttime=nextii;
	currenttimegap=timegapii;
	#pragma omp atomic
	weighting[timegapii]+=trajectory_list->show_n_trajectories(thisii);
	(trajectory_list[0]).listloop(this,timegapii,thisii,nextii);
}



void Non_Gaussian_Parameter::listkernel(Trajectory* current_trajectory)
{
//	#pragma omp atomic
//	ngp[currenttimegap]+=pow(current_trajectory->distance(currenttime,nexttime),4);
}


void Non_Gaussian_Parameter::listkernel(Trajectory* current_trajectory, int timegapii,int thisii, int nextii)
{
	#pragma omp atomic
	ngp[timegapii]+=pow(current_trajectory->distance(thisii,nextii),4);
}


void Non_Gaussian_Parameter::postprocess_list()
{
	int timeii;

	for(timeii=0;timeii<n_times;timeii++)
	{
		ngp[timeii] *= (3.0/(float(weighting[timeii])))/(5*pow((msd.show(timeii)),2.0));
		ngp[timeii] -= 1.0;
	}
}



void Non_Gaussian_Parameter::write(string filename)const
{
  int timeii;

  cout << "\nWriting non-Gaussian parameter to file " << filename <<".";

  ofstream output(filename.c_str());

  output << "Non-Gaussian parameter data created by AMDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<ngp[timeii]<<"\n";
  }
}

void Non_Gaussian_Parameter::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting non-Gaussian parameter to file.";

  output << "Non-Gaussian parameter data created by AMDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<ngp[timeii]<<"\n";
  }
}


int Non_Gaussian_Parameter::max()const
{
  int maxtime = -1; 
  int timeii;
  float maxvalue=0;

  for(timeii=1;timeii<n_times;timeii++)
  {
    if(ngp[timeii]>maxvalue)
    {
      maxvalue = ngp[timeii];
      maxtime = timeii;
    }
  }

  return maxtime;
}


void Non_Gaussian_Parameter::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  trajectory_list=t_list;

  list_displacementkernel(timegapii, thisii, nextii);
  
}


void Non_Gaussian_Parameter::postprocess_bins()
{
  postprocess_list();
}

void Non_Gaussian_Parameter::run(Trajectories trjs,string listname)
{

  msd.analyze(trjs.trajectories[listname]);

  cout << "\nCalculating mean square displacement.\n";cout.flush();
  start = time(NULL);
  analyze(trjs.trajectories[listname]);
  finish = time(NULL);
  cout << "\nCalculated mean square displacement in " << finish-start<<" seconds."<<endl;
}


void export_Non_Gaussian_Parameter(py::module& m)
    {
    py::class_<Non_Gaussian_Parameter, std::shared_ptr<Non_Gaussian_Parameter> >(m,"Non_Gaussian_Parameter",py::base<Analysis_Base>())
    .def(py::init< std::shared_ptr<System> >())
    //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
    .def("run",&Non_Gaussian_Parameter::run)
    .def("write", static_cast<void (Non_Gaussian_Parameter::*)(string )const> (&Non_Gaussian_Parameter::write))
    .def("write", static_cast<void (Non_Gaussian_Parameter::*)(ofstream& )const> (&Non_Gaussian_Parameter::write))
    ;
    }
