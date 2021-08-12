/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include <math.h>
#include <stdlib.h>
#include <vector>
#include <omp.h>

#include "incoherent_scattering_function.h"
#include "system.h"
#include "version.h"

using namespace std;
namespace py=pybind11;

Incoherent_Scattering_Function::Incoherent_Scattering_Function()
{

  system = 0;
  wavevectors = 0;
  fullblock = 0;
  
  bin_size=0;
  
  first_bin_index = 0;
  last_bin_index = 0;
  n_spacebins = 0;
  
  n_times = 0;
  firsttime = 0;
  lasttime = -1;
  
  timetable = 0;
  timegap_weighting=0;
  
  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = 0;
    correlation[timeii] = new float[n_spacebins];
    for(int wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }
}



Incoherent_Scattering_Function::Incoherent_Scattering_Function(const Incoherent_Scattering_Function & copy)
{
  system=copy.system;
  wavevectors = copy.wavevectors;
  bin_size = copy.bin_size;
  first_bin_index = copy.first_bin_index;
  last_bin_index = copy.last_bin_index;
  n_spacebins = copy.n_spacebins;
  fullblock = copy.fullblock;
  
  n_times=copy.n_times;
  firsttime=copy.firsttime;
  lasttime=copy.lasttime;
  timetable = copy.timetable;
  timegap_weighting = system->timegap_weighting(fullblock);
  
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = copy.n_atoms[timeii];
    correlation[timeii] = new float[n_spacebins];
    for(int wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=copy.correlation[timeii][wavenumberii];
    }
  }
}

#ifdef NEVER

Incoherent_Scattering_Function::~Incoherent_Scattering_Function()
{
  for(int timeii=0;timeii<n_times;timeii++)
  {
    delete [] correlation[timeii];
  }
  delete [] correlation;
  delete [] n_atoms;
}


#endif


Incoherent_Scattering_Function Incoherent_Scattering_Function::operator =(const Incoherent_Scattering_Function& copy)
{
  if(this!=&copy)
  {
  for(int timeii=0;timeii<n_times;timeii++)
  {
    delete [] correlation[timeii];
  }
  delete [] correlation;
  delete [] n_atoms;
  system=copy.system;
  wavevectors = copy.wavevectors;
  bin_size = copy.bin_size;
  first_bin_index = copy.first_bin_index;
  last_bin_index = copy.last_bin_index;
  n_spacebins = copy.n_spacebins;
  fullblock = copy.fullblock;
  
  n_times=copy.n_times;
  firsttime=copy.firsttime;
  lasttime=copy.lasttime;
  timetable = copy.timetable;
  timegap_weighting = system->timegap_weighting(fullblock);
  
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
  for(int timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = copy.n_atoms[timeii];
    correlation[timeii] = new float[n_spacebins];
    for(int wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=copy.correlation[timeii][wavenumberii];
    }
  }
  }
  return *this;
}


Incoherent_Scattering_Function::Incoherent_Scattering_Function(std::shared_ptr<System> sys, std::shared_ptr<Wave_Vectors> wv, bool fblock)
{
  int timeii, wavenumberii;

  system = sys;

  wavevectors = wv;
  fullblock = fblock;
  bin_size = wavevectors->show_delta_wavenumber();
  timegap_weighting = system->timegap_weighting(fullblock);
  n_atoms_represented = 0;
  
  
  first_bin_index = 0;
  last_bin_index = wavevectors->show_n_wavenumbers()-1;
  n_spacebins = last_bin_index-first_bin_index+1;

  n_times = system->show_n_timegaps();
  firsttime = 0;
  lasttime = n_times-1;
  timetable = system->displacement_times();

  /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = 0;
    correlation[timeii] = new float[n_spacebins];
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }
}



Incoherent_Scattering_Function::Incoherent_Scattering_Function(std::shared_ptr<System> sys, std::shared_ptr<Wave_Vectors> wv, int inner, int outer, bool fblock)
{

  int timeii, wavenumberii;

  system = sys;

  wavevectors = wv;

  bin_size = wavevectors->show_delta_wavenumber();
  fullblock = fblock;
  
  timegap_weighting = system->timegap_weighting(fullblock);
  n_atoms_represented = 0;

  first_bin_index = inner;
  last_bin_index = outer;
  n_spacebins = last_bin_index-first_bin_index+1;

  n_times = system->show_n_timegaps();
  firsttime = 0;
  lasttime = n_times-1;
  timetable = system->displacement_times();

   /*allocate memory for intermediate scattering function and zero out values*/
  correlation = new float*[n_times];
  n_atoms = new float [n_times];
  for(timeii=0;timeii<n_times;timeii++)
  {
    n_atoms[timeii] = 0;
    correlation[timeii] = new float[n_spacebins];
    for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
    {
      correlation[timeii][wavenumberii]=0;
    }
  }
}



/*-----------------------------------------------------*/
/*---Methods to use trajectory list loops over atoms---*/
/*-----------------------------------------------------*/

void Incoherent_Scattering_Function::analyze(Trajectory_List * t_list)
{
	trajectory_list=t_list;

	system->displacement_list(this,fullblock);
	postprocess_list();
}

void Incoherent_Scattering_Function::list_displacementkernel(int timegapii,int thisii, int nextii)
{
	currenttime=thisii;
	nexttime=nextii;
	timegap=timegapii;
//	//cout<<trajectory_list->show_n_trajectories(currenttime)<<"\t"<<timegap_weighting[timegap]<<"\t";
//	//n_atoms[timegap]+=float(trajectory_list->show_n_trajectories(currenttime))/float(timegap_weighting[timegap]);
//	//n_atoms[timegap]+=float(trajectory_list->show_n_trajectories(currenttime));
//	//cout << n_atoms[timegap]<<"\n";
//	trajectory_list->listloop(this,currenttime);
	trajectory_list->listloop(this,timegapii,thisii,nextii);

}

/* This is deprecated in favor of the multithreaded version below */
void Incoherent_Scattering_Function::listkernel(Trajectory* current_trajectory)
{
	int wavenumberii;
	int wavevectorii;
	vector<Coordinate>vectorlist;
	int vectorcount;
	Coordinate coordinate1;
	Coordinate coordinate2;
	
	double tempcorrelation=0;

	n_atoms[timegap]++;
	/*increment fourier bins*/
	for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
	{
		tempcorrelation=0;
		vectorlist = wavevectors->vectorlist(wavenumberii+first_bin_index);
		vectorcount = wavevectors->vectorcount(wavenumberii+first_bin_index);
		for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
		{
			coordinate1 = current_trajectory->show_unwrapped(currenttime);
			coordinate2 = current_trajectory->show_unwrapped(nexttime);
			tempcorrelation += double(cos(vectorlist[wavevectorii]&(coordinate2-coordinate1))) / double(vectorcount);
		}
		correlation[timegap][wavenumberii] += float(tempcorrelation);
	}	
}

/* This version is for multihreading */
void Incoherent_Scattering_Function::listkernel(Trajectory* current_trajectory, int timegapii, int thisii, int nextii)
{
	int wavenumberii;
	int wavevectorii;
	vector<Coordinate>vectorlist;
	int vectorcount;
	Coordinate coordinate1;
	Coordinate coordinate2;
	
	double tempcorrelation=0;

	n_atoms[timegapii]++;
	/*increment fourier bins*/
	for(wavenumberii=0;wavenumberii<n_spacebins;wavenumberii++)
	{
		tempcorrelation=0;
		vectorlist = wavevectors->vectorlist(wavenumberii+first_bin_index); //first_bin_index is a global variable! Watch out when threading!
		vectorcount = wavevectors->vectorcount(wavenumberii+first_bin_index); //first_bin_index is a global variable! Watch out when threading!
		for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
		{
			coordinate1 = current_trajectory->show_unwrapped(thisii);
			coordinate2 = current_trajectory->show_unwrapped(nextii);
			tempcorrelation += double(cos(vectorlist[wavevectorii]&(coordinate2-coordinate1))) / double(vectorcount);
		}
		#pragma omp atomic
		correlation[timegapii][wavenumberii] += float(tempcorrelation);
	}	
}

/*-----------------------------------------------------*/
/*---------Methods to use trajectory list bins---------*/
/*-----------------------------------------------------*/


void Incoherent_Scattering_Function::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  
  trajectory_list=t_list;
  list_displacementkernel(timegapii, thisii, nextii);

}

void Incoherent_Scattering_Function::postprocess_bins()
{
  postprocess_list();
}

/*Write correlation object to file*/
void Incoherent_Scattering_Function::write(string filename)const
{
  int timeii, binii;

  ofstream output (filename.c_str());

  cout << "\nWriting to file " <<filename<<"."<<endl;cout.flush();

  /*Write first row - list of bin numbers*/
  output << "Incoherent Scattering Function created by SMolDAT v." << VERSION << "\n";
  output << "\t";
  for(binii=first_bin_index;binii<=last_bin_index;binii++)
  {
    output << wavevectors->show_approx_wavenumber(binii) << "\t";  //error here
  }
  output << "\n";

  /*Write correlation function to file*/
  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    output << timetable[timeii] << "\t"; 	//first column is time data
    for(binii=0;binii<n_spacebins;binii++)
    {
      output << correlation[timeii][binii] << "\t";
//      output << correlation[timeii][binii] << "\t"; //print unnormalized data to file
    }
    output << "\n";
  }
  output.close();
}

void Incoherent_Scattering_Function::write(ofstream& output)const
{
  int timeii, binii;

  cout << "\nWriting to file."<<endl;

  /*Write first row - list of bin numbers*/
  output << "Incoherent Scattering Function created by SMolDAT v." << VERSION << "\n";
  output << "\t";
  for(binii=first_bin_index;binii<=last_bin_index;binii++)
  {
    output << wavevectors->show_approx_wavenumber(binii) << "\t";  //error here
  }
  output << "\n";

  /*Write correlation function to file*/
  for(timeii=firsttime;timeii<=lasttime;timeii++)
  {
    output << timetable[timeii] << "\t"; 	//first column is time data
    for(binii=0;binii<n_spacebins;binii++)
    {
      output << correlation[timeii][binii] << "\t";
//      output << correlation[timeii][binii] << "\t"; //print unnormalized data to file
    }
    output << "\n";
  }
}

void Incoherent_Scattering_Function::run(Trajectories cl,string listname)
{
  cout << "\nCalculating Self Intermediate Scattering Function calculated.\n";cout.flush();
  start = time(NULL);
  analyze(cl.trajectories[listname]);
  finish = time(NULL);
  cout << "\nSelf Intermediate Scattering Function calculated in " << finish-start<<" seconds.\n";
}

void export_Incoherent_Scattering_Function(py::module& m)
    {
    py::class_<Incoherent_Scattering_Function, std::shared_ptr<Incoherent_Scattering_Function> >(m,"Incoherent_Scattering_Function",py::base<Correlation_2D>())
    .def(py::init< std::shared_ptr<System>, std::shared_ptr<Wave_Vectors>, int,int , bool >())
    //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
    .def("run",&Incoherent_Scattering_Function::run)
    .def("write", static_cast<void (Incoherent_Scattering_Function::*)(string )const> (&Incoherent_Scattering_Function::write))
    .def("write", static_cast<void (Incoherent_Scattering_Function::*)(ofstream& )const> (&Incoherent_Scattering_Function::write))
    ;
    }