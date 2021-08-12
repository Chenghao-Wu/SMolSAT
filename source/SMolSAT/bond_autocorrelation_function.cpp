/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "bond_autocorrelation_function.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "multibody_list.h"
#include <omp.h>
#include "system.h"

using namespace std;
namespace py = pybind11;

Bond_Autocorrelation_Function::Bond_Autocorrelation_Function()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  baf = new float [n_times];
  weighting = new int [n_times];

  atomcount = 0;
  dimensions.set(1,1,1);
  
  vprep=&Bond_Autocorrelation_Function::prep_inplane;
  legendre_p=&Bond_Autocorrelation_Function::legendre_2;
}


Bond_Autocorrelation_Function::Bond_Autocorrelation_Function(const Bond_Autocorrelation_Function & copy)
{
  int timeii;

  system = copy.system;
  multibody_list = copy.multibody_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=copy.baf[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
  dimensions=copy.dimensions;
  
  vprep=copy.vprep;
  legendre_p=copy.legendre_p;
}

Bond_Autocorrelation_Function::~Bond_Autocorrelation_Function()
{
  delete [] baf;
  delete [] weighting;
  delete [] timetable;
}


/** **/
Bond_Autocorrelation_Function::Bond_Autocorrelation_Function(std::shared_ptr<System> sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
  dimensions.set(1,1,1);
  vprep=&Bond_Autocorrelation_Function::prep_inplane;
  legendre_p=&Bond_Autocorrelation_Function::legendre_2;

}


Bond_Autocorrelation_Function::Bond_Autocorrelation_Function(std::shared_ptr<System> sys, string dimselect)
{
  int timeii;

  system = sys;

  Coordinate dim;
  if(dimselect=="xyz")
    {
      dim.set(1,1,1);
    }
    else if(dimselect=="xy")
    {
      dim.set(1,1,0);
    }
    else if(dimselect=="xz")
    {
      dim.set(1,0,1);
    }
    else if(dimselect=="yz")
    {
      dim.set(0,1,1);
    }
    else if(dimselect=="x")
    {
      dim.set(1,0,0);
    }
    else if(dimselect=="y")
    {
      dim.set(0,1,0);
    }
    else if(dimselect=="z")
    {
      dim.set(0,0,1);
    }

  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
  dimensions=dim;
  
  if(dimensions.sum()==2||dimensions.sum()==3)
  {
    vprep=&Bond_Autocorrelation_Function::prep_inplane;
  }
  else if(dimensions.sum()==1)
  {
    vprep=&Bond_Autocorrelation_Function::prep_outofplane;
  }
  
  legendre_p=&Bond_Autocorrelation_Function::legendre_2;

}


Bond_Autocorrelation_Function::Bond_Autocorrelation_Function(std::shared_ptr<System> sys,int l_type, string dimselect)
{
  int timeii;

  system = sys;

  Coordinate dim;
  if(dimselect=="xyz")
    {
      dim.set(1,1,1);
    }
    else if(dimselect=="xy")
    {
      dim.set(1,1,0);
    }
    else if(dimselect=="xz")
    {
      dim.set(1,0,1);
    }
    else if(dimselect=="yz")
    {
      dim.set(0,1,1);
    }
    else if(dimselect=="x")
    {
      dim.set(1,0,0);
    }
    else if(dimselect=="y")
    {
      dim.set(0,1,0);
    }
    else if(dimselect=="z")
    {
      dim.set(0,0,1);
    }

  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
  dimensions=dim;
  
  if(dimensions.sum()==2||dimensions.sum()==3)
  {
    vprep=&Bond_Autocorrelation_Function::prep_inplane;
  }
  else if(dimensions.sum()==1)
  {
    vprep=&Bond_Autocorrelation_Function::prep_outofplane;
  }
  
  if(l_type==2)
  {
    legendre_p=&Bond_Autocorrelation_Function::legendre_2;
  }
  else if(l_type==1)
  {
    legendre_p=&Bond_Autocorrelation_Function::legendre_1;
  }

}



Bond_Autocorrelation_Function Bond_Autocorrelation_Function::operator = (const Bond_Autocorrelation_Function & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  multibody_list = copy.multibody_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  delete [] baf;
  delete [] weighting;

  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=copy.baf[timeii];
    weighting[timeii]=copy.weighting[timeii];
  }
  dimensions=copy.dimensions;
  legendre_p=copy.legendre_p;

  
  }

  return *this;

}


void Bond_Autocorrelation_Function::initialize(std::shared_ptr<System>  sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  delete [] baf;
  delete [] weighting;

  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
  dimensions.set(1,1,1);
  legendre_p=&Bond_Autocorrelation_Function::legendre_2;
}


void Bond_Autocorrelation_Function::initialize(std::shared_ptr<System>  sys, string dimselect)
{
  int timeii;

  system = sys;
  
  Coordinate dim;
  if(dimselect=="xyz")
    {
      dim.set(1,1,1);
    }
    else if(dimselect=="xy")
    {
      dim.set(1,1,0);
    }
    else if(dimselect=="xz")
    {
      dim.set(1,0,1);
    }
    else if(dimselect=="yz")
    {
      dim.set(0,1,1);
    }
    else if(dimselect=="x")
    {
      dim.set(1,0,0);
    }
    else if(dimselect=="y")
    {
      dim.set(0,1,0);
    }
    else if(dimselect=="z")
    {
      dim.set(0,0,1);
    }

  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  delete [] baf;
  delete [] weighting;

  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
  dimensions=dim;
  
  
  if(dimensions.sum()==2||dimensions.sum()==3)
  {
    vprep=&Bond_Autocorrelation_Function::prep_inplane;
  }
  else if(dimensions.sum()==1)
  {
    vprep=&Bond_Autocorrelation_Function::prep_outofplane;
  }
  
  legendre_p=&Bond_Autocorrelation_Function::legendre_2;
}



void Bond_Autocorrelation_Function::initialize(std::shared_ptr<System>  sys, int l_type, string dimselect)
{
  int timeii;

  system = sys;

  Coordinate dim;
  if(dimselect=="xyz")
    {
      dim.set(1,1,1);
    }
    else if(dimselect=="xy")
    {
      dim.set(1,1,0);
    }
    else if(dimselect=="xz")
    {
      dim.set(1,0,1);
    }
    else if(dimselect=="yz")
    {
      dim.set(0,1,1);
    }
    else if(dimselect=="x")
    {
      dim.set(1,0,0);
    }
    else if(dimselect=="y")
    {
      dim.set(0,1,0);
    }
    else if(dimselect=="z")
    {
      dim.set(0,0,1);
    }

  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  delete [] baf;
  delete [] weighting;

  baf = new float [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    baf[timeii]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
  dimensions=dim;
  
  
  if(dimensions.sum()==2||dimensions.sum()==3)
  {
    vprep=&Bond_Autocorrelation_Function::prep_inplane;
  }
  else if(dimensions.sum()==1)
  {
    vprep=&Bond_Autocorrelation_Function::prep_outofplane;
  }
  
  if(l_type==2)
  {
    legendre_p=&Bond_Autocorrelation_Function::legendre_2;
  }
  else if(l_type==1)
  {
    legendre_p=&Bond_Autocorrelation_Function::legendre_1;
  }
}


/*Methods to do analysis using trajectory list*/

void Bond_Autocorrelation_Function::analyze(Multibody_List * t_list)
{
  multibody_list=t_list;
  system->displacement_list(this);
  postprocess_list();
}

void Bond_Autocorrelation_Function::list_displacementkernel(int timegapii,int thisii, int nextii)
{
  weighting[timegapii]+=multibody_list->show_n_multibodies(thisii);
  multibody_list->listloop(this,timegapii, thisii, nextii);
}



void Bond_Autocorrelation_Function::listkernel(Multibody* current_multibody, int timegapii,int thisii, int nextii)
{
  float dotproduct = (this->*vprep)((*current_multibody)(1)->show_unwrapped(thisii)-(*current_multibody)(0)->show_unwrapped(thisii))&(this->*vprep)((*current_multibody)(1)->show_unwrapped(nextii)-(*current_multibody)(0)->show_unwrapped(nextii));	//compute dot product between unit vectors at initial and later times

  
  baf[timegapii]+=(this->*legendre_p)(dotproduct);//increment baf by chosen legendre polynomial of dot product above
}


void Bond_Autocorrelation_Function::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        baf[timeii] /= float(weighting[timeii]);

  }
}



/*Method to write MSD data to file*/

void Bond_Autocorrelation_Function::write(string filename)const
{
  int timeii;

  cout << "\nWriting baf to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Bond autocorrelation function data created bys SMolDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<baf[timeii]<<"\n";
  }
}

void Bond_Autocorrelation_Function::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting baf to file.";

  output << "Bond autocorrelation function data created bys SMolDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<baf[timeii]<<"\n";
  }
}



Coordinate Bond_Autocorrelation_Function::prep_inplane(Coordinate coord)const
{
  return (coord*dimensions).unit_vector();
}


Coordinate Bond_Autocorrelation_Function::prep_outofplane(Coordinate coord)const
{
  Coordinate addend(1,1,1);
  Coordinate newcoord;
  newcoord.set((coord*dimensions).length(),((dimensions*(-1)+addend)*coord).length(),0);
  return newcoord.unit_vector();
}


#ifdef NEVER
void Bond_Autocorrelation_Function::bin_hook(Trajectory_List * t_list, int timegapii, int thisii, int nextii)
{
  multibody_list=t_list;

  list_displacementkernel(timegapii, thisii, nextii);

}



void Bond_Autocorrelation_Function::postprocess_bins()
{
  postprocess_list();
}

#endif

float Bond_Autocorrelation_Function::legendre_1(float val)const
{
  return val;
}

float Bond_Autocorrelation_Function::legendre_2(float val)const
{
  return 0.5*(3.0*val*val - 1.0);
}

void Bond_Autocorrelation_Function::run(Trajectories trjs,string multibody_list_name)
{
  
  Multibody_List * multibodylist;
  
  multibodylist = trjs.find_multibody_list(multibody_list_name);

  cout << "\nCalculating bond autocorrelation function.\n";cout.flush();
  start = time(NULL);
  analyze(multibodylist);
  finish = time(NULL);
  cout << "\nCalculated bond autocorrelation function in " << finish-start<<" seconds."<<endl;
}


void export_Bond_Autocorrelation_Function(py::module& m)
    {
    py::class_<Bond_Autocorrelation_Function, std::shared_ptr<Bond_Autocorrelation_Function> >(m,"Bond_Autocorrelation_Function",py::base<Multibody_Analysis>())
    .def(py::init< std::shared_ptr<System>, string >())
    //.def("analyze", static_cast<void (Bond_Autocorrelation_Function::*)(Trajectory_List* )> (&Bond_Autocorrelation_Function::analyze))
    .def("run",&Bond_Autocorrelation_Function::run)
    .def("write", static_cast<void (Bond_Autocorrelation_Function::*)(string )const> (&Bond_Autocorrelation_Function::write))
    .def("write", static_cast<void (Bond_Autocorrelation_Function::*)(ofstream& )const> (&Bond_Autocorrelation_Function::write))
    ;
    }