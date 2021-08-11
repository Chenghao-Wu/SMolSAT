/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include "version.h"
#include "structure_factor.h"
#include <omp.h>


using namespace std;
namespace py=pybind11;

Structure_Factor::Structure_Factor()
{
  int vectorcount;
  trajectory_list=0;
  trajectory_list2=0;
  n_wavenumbers=0;
  current_wavedensity=0;
  
  n_wavenumbers = 0;
  n_atoms=0;
  currenttime=0;
  time_scheme=-1;

  structure_factor = new float [n_wavenumbers];
  wavedensity1 = new complex<float> * [n_wavenumbers];
  wavedensity2 = new complex<float> * [n_wavenumbers];
  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    vectorcount = wavevectors.vectorcount(wavenumberii);
    wavedensity1[wavenumberii]= new complex<float> [vectorcount];
    wavedensity2[wavenumberii]= new complex<float> [vectorcount];
    for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)		//loop over wavevectors for this wavenumber
    {
      wavedensity1[wavenumberii][wavevectorii].real(0);
      wavedensity1[wavenumberii][wavevectorii].imag(0);
      wavedensity2[wavenumberii][wavevectorii].real(0);
      wavedensity2[wavenumberii][wavevectorii].imag(0);
    }
  }
}

// Standard Constructor
Structure_Factor::Structure_Factor(std::shared_ptr<System> sys, string plane ,float max_length_scale, int timescheme)
{
  system = sys;

  if(max_length_scale<0&&timescheme<-1)
  {
    max_length_scale=sys->size(-timescheme-2).min();
  }

  Wave_Vectors wv(sys,plane,max_length_scale);

  int vectorcount;
  trajectory_list=0;
  trajectory_list2=0;
  current_wavedensity=0;

  wavevectors = wv;

  time_scheme=timescheme;
  n_wavenumbers = wavevectors.show_n_wavenumbers();

  structure_factor = new float [n_wavenumbers];
  wavedensity1 = new complex<float> * [n_wavenumbers];
  wavedensity2 = new complex<float> * [n_wavenumbers];

  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    vectorcount = wavevectors.vectorcount(wavenumberii);
    wavedensity1[wavenumberii]= new complex<float> [vectorcount];
    wavedensity2[wavenumberii]= new complex<float> [vectorcount];
    for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)		//loop over wavevectors for this wavenumber
    {
      wavedensity1[wavenumberii][wavevectorii].real(0);
      wavedensity1[wavenumberii][wavevectorii].imag(0);
      wavedensity2[wavenumberii][wavevectorii].real(0);
      wavedensity2[wavenumberii][wavevectorii].imag(0);
    }
  }
}

Structure_Factor::~Structure_Factor()
{
  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
   {
     delete [] wavedensity1[wavenumberii];
     delete [] wavedensity2[wavenumberii];
  }
  delete [] structure_factor;
  delete [] wavedensity1;
  delete [] wavedensity2;
}
/*---------------------------------------------------------------------------------------------*/



Structure_Factor::Structure_Factor(const Structure_Factor & copy)
:Analysis_Onetime(copy)
{
  int vectorcount;
    time_scheme = copy.time_scheme;
    system = copy.system;
    trajectory_list=copy.trajectory_list;
    trajectory_list2=copy.trajectory_list2;
    n_wavenumbers=copy.n_wavenumbers;
    wavevectors=copy.wavevectors;
    n_atoms=copy.n_atoms;

    structure_factor = new float [n_wavenumbers];
    wavedensity1 = new complex<float> * [n_wavenumbers];
    wavedensity2 = new complex<float> * [n_wavenumbers];

    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors.vectorcount(wavenumberii);
      wavedensity1[wavenumberii]=new complex<float> [vectorcount];
      wavedensity2[wavenumberii]=new complex<float> [vectorcount];
      for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
      {
	wavedensity1[wavenumberii][wavevectorii] = copy.wavedensity1[wavenumberii][wavevectorii];
	wavedensity2[wavenumberii][wavevectorii] = copy.wavedensity2[wavenumberii][wavevectorii];
      }
    }
}


Structure_Factor Structure_Factor::operator=(const Structure_Factor & copy)
{
  if(this!=&copy)
  {
    int vectorcount;
    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      delete [] wavedensity1[wavenumberii];
      delete [] wavedensity2[wavenumberii];
    }
    delete [] structure_factor;
    delete [] wavedensity1;
    delete [] wavedensity2;

    time_scheme = copy.time_scheme;
    system = copy.system;
    trajectory_list=copy.trajectory_list;
    trajectory_list2=copy.trajectory_list2;
    n_wavenumbers=copy.n_wavenumbers;
    wavevectors=copy.wavevectors;
    n_atoms=copy.n_atoms;

    structure_factor = new float [n_wavenumbers];
    wavedensity1 = new complex<float> * [n_wavenumbers];
    wavedensity2 = new complex<float> * [n_wavenumbers];

    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors.vectorcount(wavenumberii);
      wavedensity1[wavenumberii]=new complex<float> [vectorcount];
      wavedensity2[wavenumberii]=new complex<float> [vectorcount];
      for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
      {
	wavedensity1[wavenumberii][wavevectorii] = copy.wavedensity1[wavenumberii][wavevectorii];
	wavedensity2[wavenumberii][wavevectorii] = copy.wavedensity2[wavenumberii][wavevectorii];
      }
    }
  }
  return *this;
}


void Structure_Factor::preprocess()
{
  int wavenumberii;
  n_atoms=0;
  
  /*allocate memory for wavedensity and structure factor and zero out structure factor*/
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    
    structure_factor[wavenumberii] = 0;
    wavedensity1[wavenumberii]=new complex<float> [wavevectors.vectorcount(wavenumberii)];
    

  }
}

void Structure_Factor::preprocess2()
{
  int wavenumberii;
  
  n_atoms=0;
  
  /*allocate memory for wavedensity and structure factor and zero out structure factor*/
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
    structure_factor[wavenumberii] = 0;
    wavedensity1[wavenumberii]=new complex<float> [wavevectors.vectorcount(wavenumberii)];
    wavedensity2[wavenumberii]=new complex<float> [wavevectors.vectorcount(wavenumberii)];
  }
}

void Structure_Factor::analyze(Trajectory_List * t_list)
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
      timekernel(timeii);
    }
  }
  
  postprocess_list();
}

void Structure_Factor::analyze(Trajectory_List * t_list,Trajectory_List* t_list2)
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

/*Calculates symmetric structure factor*/
void Structure_Factor::timekernel(int timeii)
{
  int vectorcount;
  
   //Zero out wave densities
    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors.vectorcount(wavenumberii);
      for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)		//loop over wavevectors for this wavenumber
      {
        wavedensity1[wavenumberii][wavevectorii].real(0);
        wavedensity1[wavenumberii][wavevectorii].imag(0);
      }
    }

    
    /*calculate wave densities*/
    currenttime=timeii;
    current_wavedensity = wavedensity1;
    analyze_wave_density(trajectory_list); 
    //Calculate structure factor
    for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors.vectorcount(wavenumberii);
      for(int wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
      {
        structure_factor[wavenumberii] += (wavedensity1[wavenumberii][wavevectorii].real()*wavedensity1[wavenumberii][wavevectorii].real() + wavedensity1[wavenumberii][wavevectorii].imag()*wavedensity1[wavenumberii][wavevectorii].imag())/(float(vectorcount));
      }
    }
    n_atoms+=trajectory_list->show_n_trajectories(timeii);
}


//*Calculated asymmetric structure factor*//
void Structure_Factor::timekernel2(int timeii)
{
  int wavenumberii, wavevectorii;
  int vectorcount;
  
  //Zero out wave densities
    for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors.vectorcount(wavenumberii);
      for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)		//loop over wavevectors for this wavenumber
      {
        wavedensity1[wavenumberii][wavevectorii].real(0);
        wavedensity1[wavenumberii][wavevectorii].imag(0);
        wavedensity2[wavenumberii][wavevectorii].real(0);
        wavedensity2[wavenumberii][wavevectorii].imag(0);
      }
    }

    /*calculate wave densities*/
    currenttime=timeii;
    current_wavedensity = wavedensity1;
    analyze_wave_density(trajectory_list);
    current_wavedensity = wavedensity2;
    analyze_wave_density(trajectory_list2);

    //Calculate structure factor
    for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
    {
      vectorcount = wavevectors.vectorcount(wavenumberii);
      for(wavevectorii=0;wavevectorii<vectorcount;wavevectorii++)
      {
        structure_factor[wavenumberii] += (wavedensity1[wavenumberii][wavevectorii].real()*wavedensity2[wavenumberii][wavevectorii].real() + wavedensity1[wavenumberii][wavevectorii].imag()*wavedensity2[wavenumberii][wavevectorii].imag())/(float(vectorcount));
      }
    }

    n_atoms+=trajectory_list->show_n_trajectories(timeii);
}

void Structure_Factor::postprocess_list()
{
  for(int wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)
  {
      if (float(n_atoms)!=0)
      {
         structure_factor[wavenumberii]/=float(n_atoms);
      }
    else
    {
        structure_factor[wavenumberii]=0;
    }
  }
}


void Structure_Factor::analyze_wave_density(Trajectory_List * t_list)
{
  t_list->listloop(this,currenttime);
}

void Structure_Factor::listkernel(Trajectory* current_trajectory)
{
  int wavenumberii, vectorii;
  vector<Coordinate>vectorlist;
  Coordinate coordinate;
  int vectorcount;
  float k_dot_r;
  complex<double> tempcomplex;

  coordinate = current_trajectory->show_coordinate(currenttime);	//get atom coordinate at give time // the currenttime reference is... actually probably ok.
  for(wavenumberii=0;wavenumberii<n_wavenumbers;wavenumberii++)		//loop over wavenumbers for which wavedensity will be calculated
  {
    vectorlist = wavevectors.vectorlist(wavenumberii);	//call up wave vector list for this wavenumber //The wavevectors reference may be what's breaking parallelization
    vectorcount = wavevectors.vectorcount(wavenumberii);	//get count of wave vectors for this wavenumber
    for(vectorii=0;vectorii<vectorcount;vectorii++)		//loop over wavevectors for this wavenumber
    {
      k_dot_r = vectorlist[vectorii]&coordinate;		//calculate dot product of wave vector and present atomecoordinate
      tempcomplex.real(cos(k_dot_r));
      tempcomplex.imag(sin(k_dot_r));
//       #pragma omp atomic
       
       current_wavedensity[wavenumberii][vectorii] += tempcomplex;
      //current_wavedensity[wavenumberii][vectorii].real() += cos(k_dot_r);	//add contribution to real part of wave density // This COULD break paralleliation. Might need to omp flush the variable?
//      #pragma omp atomic
      //current_wavedensity[wavenumberii][vectorii].imag() += sin(k_dot_r);	//add contribution to imaginary part of wave density
    }
  }
}

/*Write correlation object to file*/
void Structure_Factor::write(string filename)const
{
  int binii;
  Coordinate mean_wavevector;

  ofstream output (filename.c_str());

  cout << "\nWriting to file " <<filename<<".";cout.flush();

  /*Write first row - list of bin numbers*/
  output << "Structure factor data created by AMDAT v." << VERSION << "\n";
  output << "Nominal_q\tMean_q\tStddev_q\tMean_qx\tMean_qy\tMean_qz\tS(q)\n";
  for(binii=0;binii<n_wavenumbers;binii++)
  {
    mean_wavevector = wavevectors.show_mean_wavevector(binii);
    output << wavevectors.show_approx_wavenumber(binii) << "\t" << wavevectors.show_mean_wavenumber(binii) << "\t" << wavevectors.show_stdev_wavenumber(binii) << "\t"<< mean_wavevector.show_x() << "\t" << mean_wavevector.show_y() << "\t" << mean_wavevector.show_z() << "\t" << structure_factor[binii]  <<  endl;
  }
  output << "\n";

  output.close();
}


void Structure_Factor::write(ofstream& output)const
{
  int binii;
  Coordinate mean_wavevector;


  cout << "\nWriting to file.";

  /*Write first row - list of bin numbers*/
  output << "Structure factor data created by AMDAT v." << VERSION << "\n";
  output << "Nominal_q\tMean_q\tStddev_q\tMean_qx\tMean_qy\tMean_qz\tS(q)\n";
  for(binii=0;binii<n_wavenumbers;binii++)
  {
    mean_wavevector = wavevectors.show_mean_wavevector(binii);
    output << wavevectors.show_approx_wavenumber(binii) << "\t" << wavevectors.show_mean_wavenumber(binii) << "\t" << wavevectors.show_stdev_wavenumber(binii) << "\t"<< mean_wavevector.show_x() << "\t" << mean_wavevector.show_y() << "\t" << mean_wavevector.show_z() << "\t" << structure_factor[binii]  <<  endl;
  }
  output << "\n";

}

void Structure_Factor::run(Trajectories cl,string listname)
{
  /*Reading wave vectors from file*/
  cout << "\nReading wave vectors from file.\n";cout.flush();

  cout <<"\n"<< system->show_n_trajectories()<<"\t"<<wavevectors.show_n_wavenumbers()<<"\t"<<cl.find_trajectorylist(listname)->show_n_trajectories(0)<<"\t"<<time_scheme<<"\t";cout.flush();


  cout << "\nCalculating structure factor.\n";cout.flush();
  start = time(NULL);
  analyze(cl.trajectories[listname]);
  finish = time(NULL);
  cout << "\nCalculated structure factor in " << finish-start<<" seconds.\n";

}

void Structure_Factor::run(Trajectories cl,string listname1,string listname2)
{
  /*Reading wave vectors from file*/
  cout << "\nReading wave vectors from file.\n";cout.flush();

  cout <<"\n"<< system->show_n_trajectories()<<"\t"<<wavevectors.show_n_wavenumbers()<<"\t"<<cl.find_trajectorylist(listname1)->show_n_trajectories(0)<<"\t"<<time_scheme<<"\t";cout.flush();

  cout << "\nCalculating structure factor.\n";cout.flush();
  start = time(NULL);
  analyze(cl.find_trajectorylist(listname1),cl.find_trajectorylist(listname1));
  finish = time(NULL);
  cout << "\nCalculated structure factor in " << finish-start<<" seconds.\n";

}

void export_Structure_Factor(py::module& m)
    {
    py::class_<Structure_Factor, std::shared_ptr<Structure_Factor> >(m,"Structure_Factor",py::base<Analysis_Onetime>())
    .def(py::init< std::shared_ptr<System>, string, float, int >())
    //.def("analyze", static_cast<void (Mean_Square_Displacement::*)(Trajectory_List* )> (&Mean_Square_Displacement::analyze))
    .def("run",static_cast<void (Structure_Factor::*)(Trajectories, string )> (&Structure_Factor::run))
    .def("run",static_cast<void (Structure_Factor::*)(Trajectories, string, string )> (&Structure_Factor::run))
    .def("write", static_cast<void (Structure_Factor::*)(string )const> (&Structure_Factor::write))
    .def("write", static_cast<void (Structure_Factor::*)(ofstream& )const> (&Structure_Factor::write))
    ;
    }