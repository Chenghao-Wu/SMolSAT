/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/
#include "trajectories.h"
#include "mean_square_displacement.h"


#include "../extern/pybind11/include/pybind11/stl.h"

using namespace std;

namespace py = pybind11;


Trajectories::Trajectories(std::shared_ptr<System> sys)
{
  analyte = sys;
}

void Trajectories::initialize_lists()
 {
   //list of particle lists
   n_gaussian_comparisons = 0;
   vhs_defined = 0;
   n_trajectory_list_bins = 0;

 }

  /*finds trajectorylist object by custom name*/
Value_List<float>* Trajectories::find_value_list(string listname, bool allow_nofind)const
{
  
  Value_List<float>* v_list;

  try
  {
    v_list = value_lists.at(listname);
  }
  catch(out_of_range & sa)
  {
    if(allow_nofind)
    {
      v_list=0;
    }
    else
    {
      cout << "\nError: value_list " << listname << " does not exist.\n";
      exit(0);
    }
  }
  
  return v_list;

}



 /*--------------------------------------------------------------------------------*/



void Trajectories::add_value_list(Value_List<float>* av_list, string listname)
{
  
  bool result;

  result=(value_lists.insert(listname,av_list));

  if(!result)
  {
    cout << "\nWarning: neighbor_list "<< listname<<" not created because a neighbor_list with this name already exists. Replacement of a neighbor_list requires that you first delete the existing list with the same name.\n";
  }


}


 /*--------------------------------------------------------------------------------*/



void Trajectories::delete_value_list(string listname)
{
  
  bool result;

  Value_List<float> * vlist;
  
  //check if the specified neighbor_list exists
  vlist = find_value_list(listname, 1);
  if(vlist==0)
  {
    cout << "\nWarning: neighbor_list "<< listname<<" not deleted because it does not exist.\n";
  }
  else
  {
    value_lists.erase(listname);	//remove this neighbor_list from list of neighbor_lists
    delete vlist;			//deallocate memory for this neighbor_list
  }


}


  /*finds trajectorylist object by custom name*/
Trajectory_List* Trajectories::find_trajectorylist(string listname, bool allow_nofind)const
{
    Trajectory_List * trajectory_list;

  try
  {
    trajectory_list = trajectories.at(listname);
  }
  catch(out_of_range & sa)
  {
    if(allow_nofind)
    {
      trajectory_list=0;
    }
    else
    {
      cout << "\nError: trajectory_list " << listname << " does not exist.\n";
      exit(0);
    }
  }

  return trajectory_list;

}

void Trajectories::add_trajectorylist(Trajectory_List * t_list, string listname)
{
 bool result;

  result=(trajectories.insert(listname,t_list));

  if(!result)
  {
    cout << "\nWarning: trajectory_list "<< listname<<" not created because a trajectory_list with this name already exists. Replacement of a trajectory_list requires that you first delete the existing list with the same name.\n";
  }


}

void Trajectories::combine_trajectories(string newlistname, vector<string> listnames)
{
    int argii;

    int n_args=listnames.size();

    Trajectory_List * new_trajectory_list;
    new_trajectory_list = new Trajectory_List;

    //newlistname=args[1];

    (*new_trajectory_list)=(*find_trajectorylist(listnames[0]));

    for(argii=1;argii<n_args;argii++)
    {
        (*new_trajectory_list)=(*new_trajectory_list)||(*find_trajectorylist(listnames[argii]));
    }

    add_trajectorylist(new_trajectory_list,newlistname);

}

void Trajectories::delete_trajectory_list(string listname)
{
  
  Trajectory_List * trajectory_list;

  trajectory_list = find_trajectorylist(listname,1);
  if(trajectory_list==0)
  {
    cout << "\nWarning: trajectory_list " << listname << " does not exist and therefore cannot be deleted.";
  }
  else
  {
    trajectories.erase(listname);
    delete trajectory_list;
  }
}


Multibody_List* Trajectories::find_multibody_list(string listname,bool allow_nofind)const
{
  Multibody_List * multibody_list;

  try
  {
    multibody_list = multibody_lists.at(listname);
  }
  catch(out_of_range & sa)
  {
    if(allow_nofind)
    {
      multibody_list=0;
    }
    else
    {
      cout << "\nError: multibody_list " << listname << " does not exist.\n";
      exit(0);
    }
  }

  return multibody_list;
}


void Trajectories::add_multibody_list(Multibody_List* multibody_list,string multibody_list_name)
{
  bool result;

  result=(multibody_lists.insert(multibody_list_name,multibody_list));

  if(!result)
  {
    cout << "\nWarning:multibody_list "<< multibody_list_name<<" not created because a multibody_list with this name already exists. Replacement of a multibody_list requires that you first delete the existing list with the same name.\n";
  }
}


void Trajectories::delete_multibody_list(string listname)
{
  
  Multibody_List * multibody_list;

  multibody_list = find_multibody_list(listname,1);
  if(multibody_list==0)
  {
    cout << "\nWarning: multibody_list " << listname << " does not exist and therefore cannot be deleted.";
  }
  else
  {
    multibody_lists.erase(listname);
    delete [] multibody_list;
  }
}

void Trajectories::combine_multibody_lists(string newlistname, vector<string> listnames)
{
    int argii;

    int n_args=listnames.size();

    Multibody_List * new_multibody_list;
    new_multibody_list = new Multibody_List;

    (*new_multibody_list)=(*find_multibody_list(listnames[0]));

    for(argii=1;argii<n_args;argii++)
    {
        (*new_multibody_list)=(*new_multibody_list)+(*find_multibody_list(listnames[argii]));
    }

    cout<<"\nA combined motibody list "<<newlistname<<" created.";
    add_multibody_list(new_multibody_list,newlistname);

}

void Trajectories::delete_multibodies(string multibody_set_name)
{  
   analyte->delete_multibody_set(multibody_set_name);
}

void Trajectories::add_trajectorylist_bins(Trajectory_List_Bins * traj_list_bins, string listname)
{
  /** Adds new Binned Trajectory List object to stored values
  * @author Mark Mackura
  * @date 4/11/2012
  **/
  int trajnum;

  trajnum = find_trajectorylist_bins(listname);
  if(trajnum==-1)
  {
    trajnum = n_trajectory_list_bins;
    binned_trajectories.push_back(traj_list_bins);
    n_trajectory_list_bins++;
  }
  else
  {
    binned_trajectories[trajnum] = traj_list_bins;
  }
  trajectory_list_bin_names[trajnum] = listname;
}


int Trajectories::find_trajectorylist_bins(string listname)
{
  /** Finds trajectory_list_bins object by custom name
  * @author Mark Mackura
  * @date 4/11/2012
  **/
    int listii;
    for(listii=0;listii<n_trajectory_list_bins;listii++)
    {
      if(listname==trajectory_list_bin_names[listii])
      {
	return listii;
      }
    }
    return -1;
}


void Trajectories::remove_bin_list(string listname)
{
  cout << "\nRemoving binned trajectory list " << listname << endl;

  int listii = find_trajectorylist_bins(listname);
  if (listii<0)
  {
      cout << "\nBinned list " << listname << " not found! Cannot remove it!" << endl;
  }
  else
  {
      delete binned_trajectories.at(listii);
      binned_trajectories.erase(binned_trajectories.begin() + listii);
      for (int i=listii; i<n_trajectory_list_bins; i++)
      {
          trajectory_list_bin_names[i]=trajectory_list_bin_names[i+1];
      }
      n_trajectory_list_bins--;
  }
}

void Trajectories::create_list(string listname, string args)
{
    Static_Trajectory_List* trajectory;
    Trajectory_List * trajpointer;

    trajectory = new Static_Trajectory_List;
    trajectory->reset(analyte);
    trajpointer=(Trajectory_List*)trajectory;
    trajectory->analyze(args);
    add_trajectorylist(trajpointer, listname);	//add trajectory list to array
    cout<<"\nTrajectory list "<<listname<<" created with "<<trajpointer->show_n_trajectories(0)<< " trajectories."<<endl;
}

void Trajectories::create_multibodies(string multibody_list_name, string trajectory_type_name, string centertypename, string args)
{
    //int n_args = tokenize(args).count()+4;
    string trajectory_list_name;
    int centertype;

    Multibody_Set* multibody_set_pointer;
    Multibody_List* new_multibody_list;
    Trajectory_Set * trajectory_set_pointer;
    Static_Trajectory_List * new_trajectory_list;

    new_trajectory_list = new Static_Trajectory_List;
    new_multibody_list=new Multibody_List;

    trajectory_list_name = multibody_list_name;

    multibody_set_pointer = analyte->create_multibody_set(multibody_list_name, args);    //create multibody set with name that is the same as the multibody list. This is where the multibodies are created.

    new_multibody_list->set(analyte,multibody_set_pointer);
    add_multibody_list(new_multibody_list,multibody_list_name);

    if(centertypename == "centroid")
    {
      centertype = 0;
    }
    else if(centertypename == "com")
    {
      centertype = 1;
    }
    else
    {
      cout << "\n Type of multibody center '" << centertypename << "' not recognized. Allowable options are 'centroid' and 'com'";
      exit (0);
    }

    trajectory_set_pointer = analyte->create_trajectory_set(trajectory_list_name,multibody_list_name,trajectory_type_name, centertype);

    new_trajectory_list->set(analyte,trajectory_set_pointer);
    add_trajectorylist(new_trajectory_list, trajectory_list_name);
}


void export_Trajectories(pybind11::module& m)
{
    py::class_<Trajectories, std::shared_ptr<Trajectories> >(m,"Trajectories")
    .def(py::init< std::shared_ptr<System> >(),py::arg("system"))
    .def("create_list", &Trajectories::create_list,py::arg("name"),py::arg("args"))
    .def("create_multibodies", &Trajectories::create_multibodies,py::arg("name"),py::arg("trj_list_name"),py::arg("center_type"),py::arg("args"))
    .def("combine_multibody_lists", &Trajectories::combine_multibody_lists,py::arg("name"),py::arg("multibodies"))
    .def("combine_trajectories", &Trajectories::combine_trajectories,py::arg("name"),py::arg("trjs"))
    .def("delete_trajectory_list", &Trajectories::delete_trajectory_list)
    .def("delete_multibodies", &Trajectories::delete_multibodies)
    .def("delete_multibody_list", &Trajectories::delete_multibody_list)
    ;
}
