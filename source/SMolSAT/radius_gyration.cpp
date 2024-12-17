/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "radius_gyration.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "version.h"
#include "multibody_list.h"
#include "system.h"

using namespace std;
namespace py = pybind11;

Radius_Gyration::Radius_Gyration()
{
  n_times = 0;

   //allocate memory for mean square displacement data
  rg = new float [n_times];
  gytensor = new sixfloat [n_times];

  weighting = new int [n_times];

  atomcount = 0;
}


Radius_Gyration::Radius_Gyration(const Radius_Gyration & copy)
{
  int timeii;

  system = copy.system;
  multibody_list = copy.multibody_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  rg = new float [n_times];
  gytensor = new sixfloat [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    rg[timeii]=copy.rg[timeii];
    gytensor[timeii][0] = copy.gytensor[timeii][0];
    gytensor[timeii][1] = copy.gytensor[timeii][1];
    gytensor[timeii][2] = copy.gytensor[timeii][2];
    gytensor[timeii][3] = copy.gytensor[timeii][3];
    gytensor[timeii][4] = copy.gytensor[timeii][4];
    gytensor[timeii][5] = copy.gytensor[timeii][5];
    weighting[timeii]=copy.weighting[timeii];
  }
}

Radius_Gyration::~Radius_Gyration()
{
  delete [] rg;
  delete [] gytensor;
  delete [] weighting;
  delete [] timetable;
}


/** **/
Radius_Gyration::Radius_Gyration(std::shared_ptr<System> sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data
  rg = new float [n_times];
  gytensor = new sixfloat [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    rg[timeii]=0;
    gytensor[timeii][0]=0;
    gytensor[timeii][1]=0;
    gytensor[timeii][2]=0;
    gytensor[timeii][3]=0;
    gytensor[timeii][4]=0;
    gytensor[timeii][5]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;

}


Radius_Gyration Radius_Gyration::operator = (const Radius_Gyration & copy)
{
  int timeii;

  if(this!=&copy)
  {

  system = copy.system;
  multibody_list = copy.multibody_list;

  n_times = copy.n_times;
  atomcount = copy.atomcount;

  delete [] rg;
  delete [] gytensor;
  delete [] weighting;

  rg = new float [n_times];
  gytensor = new sixfloat [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();

  for(timeii=0;timeii<n_times;timeii++)
  {
    rg[timeii]=copy.rg[timeii];
    gytensor[timeii][0] = copy.gytensor[timeii][0];
    gytensor[timeii][1] = copy.gytensor[timeii][1];
    gytensor[timeii][2] = copy.gytensor[timeii][2];
    gytensor[timeii][3] = copy.gytensor[timeii][3];
    gytensor[timeii][4] = copy.gytensor[timeii][4];
    gytensor[timeii][5] = copy.gytensor[timeii][5];
    weighting[timeii]=copy.weighting[timeii];
  }

  }

  return *this;

}


void Radius_Gyration::initialize(std::shared_ptr<System>  sys)
{
  int timeii;

  system = sys;
  n_times = system->show_n_timegaps();

   //allocate memory for mean square displacement data

  delete [] rg;
  delete [] gytensor;
  delete [] weighting;

  rg = new float [n_times];
  gytensor = new sixfloat [n_times];
  weighting = new int [n_times];

  timetable = system->displacement_times();
  for(timeii=0;timeii<n_times;timeii++)
  {
    rg[timeii]=0;
    gytensor[timeii][0]=0;
    gytensor[timeii][1]=0;
    gytensor[timeii][2]=0;
    gytensor[timeii][3]=0;
    gytensor[timeii][4]=0;
    gytensor[timeii][5]=0;
    weighting[timeii]=0;
  }
  atomcount = 0;
}



/*Methods to do analysis using trajectory list*/

void Radius_Gyration::analyze(Multibody_List * t_list)
{
  int timeii;
  multibody_list=t_list;
  for(timeii=0;timeii<n_times;timeii++)
  {
    weighting[timeii]+=multibody_list->show_n_multibodies(timeii);
    multibody_list->listloop(this,0, timeii, 0);
  }
  postprocess_list();
}

void Radius_Gyration::listkernel(Multibody* current_multibody, int timegapii,int thisii, int nextii)
{

  sixfloat gytensor_this;

  rg[thisii]+=current_multibody->square_gyration_radius(thisii);

  trajectories=current_multibody->show_bodies();
  int n_trajectories=trajectories.size();
  coordinates = new Coordinate [n_trajectories];

    for(int trajii=0;trajii<n_trajectories;trajii++)
    {
      coordinates[trajii]=trajectories[trajii]->show_unwrapped(thisii);
	  //coordinates[trajii]=consistent_position(trajii,timeii);
    }
    gyration_tensor(coordinates, n_trajectories, gytensor_this);

  gytensor[thisii][0]+=gytensor_this[0];
  gytensor[thisii][1]+=gytensor_this[1];
  gytensor[thisii][2]+=gytensor_this[2];
  gytensor[thisii][3]+=gytensor_this[3];
  gytensor[thisii][4]+=gytensor_this[4];
  gytensor[thisii][5]+=gytensor_this[5];
  // this is to calculate the centroid position of the multibody!!! This is not that correct for the atomic system
}

void Radius_Gyration::postprocess_list()
{

   for(int timeii=0;timeii<n_times;timeii++)
  {

        rg[timeii] /= float(weighting[timeii]);
        gytensor[timeii][0]/=float(weighting[timeii]);
        gytensor[timeii][1]/=float(weighting[timeii]);
        gytensor[timeii][2]/=float(weighting[timeii]);
        gytensor[timeii][3]/=float(weighting[timeii]);
        gytensor[timeii][4]/=float(weighting[timeii]);
        gytensor[timeii][5]/=float(weighting[timeii]);

  }
}

void Radius_Gyration::write(string filename)const
{
  int timeii;

  cout << "\nWriting radius of gyration to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Radius of gyration data created bys SMolDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<rg[timeii]<<"\n";
  }
}

void Radius_Gyration::write_tensor(string filename)const
{
  int timeii;

  cout << "\nWriting gyration tensor to file "<<filename<<".";

  ofstream output(filename.c_str());

  output << "Gyration tensor data (xx xy xz yy yz zz) created bys SMolDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<gytensor[timeii][0]<<"\t"<<gytensor[timeii][1]<<"\t"<<gytensor[timeii][2]<<"\t"<<gytensor[timeii][3]<<"\t"<<gytensor[timeii][4]<<"\t"<<gytensor[timeii][5]<<"\n";
  }
}

void Radius_Gyration::write(ofstream& output)const
{
  int timeii;

  cout << "\nWriting radius of gyration to file.";

  output << "Radius of gyration data created bys SMolDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<rg[timeii]<<"\n";
  }
}

void Radius_Gyration::write_tensor(ofstream& output)const
{
  int timeii;

  cout << "\nWriting gyration tensor  to file "<<".";

  output << "Gyration tensor data (xx xy xz yy yz zz) created bys SMolDAT v." << VERSION << "\n";
  for(timeii=0;timeii<n_times;timeii++)
  {
    output << timetable[timeii]<<"\t"<<gytensor[timeii][0]<<"\t"<<gytensor[timeii][1]<<"\t"<<gytensor[timeii][2]<<"\t"<<gytensor[timeii][3]<<"\t"<<gytensor[timeii][4]<<"\t"<<gytensor[timeii][5]<<"\n";
  }
}

void Radius_Gyration::run(Trajectories trjs,string multibody_list_name)
{
  
  Multibody_List * multibodylist;
  
  multibodylist = trjs.find_multibody_list(multibody_list_name);

  cout << "\nCalculating radius of gyration.\n";cout.flush();
  start = time(NULL);
  analyze(multibodylist);
  finish = time(NULL);
  cout << "\nCalculated radius of gyration in " << finish-start<<" seconds."<<endl;
}


void export_Radius_Gyration(py::module& m)
    {
    py::class_<Radius_Gyration, std::shared_ptr<Radius_Gyration> >(m,"Radius_Gyration",py::base<Multibody_Analysis>())
    .def(py::init< std::shared_ptr<System> >())
    //.def("analyze", static_cast<void (Radius_Gyration::*)(Trajectory_List* )> (&Radius_Gyration::analyze))
    .def("run",&Radius_Gyration::run)
    .def("write", static_cast<void (Radius_Gyration::*)(string )const> (&Radius_Gyration::write))
    .def("write", static_cast<void (Radius_Gyration::*)(ofstream& )const> (&Radius_Gyration::write))
    .def("write_tensor", static_cast<void (Radius_Gyration::*)(string )const> (&Radius_Gyration::write_tensor))
    .def("write_tensor", static_cast<void (Radius_Gyration::*)(ofstream& )const> (&Radius_Gyration::write_tensor))
    ;
    }