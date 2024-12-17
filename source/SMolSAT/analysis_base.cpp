/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include <iostream>
#include <stdlib.h>
#include "system.h"
#include "analysis_base.h"
#define ARGMAX 10000

using namespace std;
namespace py = pybind11;

/*----------------Constructors and assignment----------------------*/

/*Default constructor*/
Analysis_Base::Analysis_Base()
{
  system = 0;
  trajectory_list=0;
}


/*Copy constructor*/
Analysis_Base::Analysis_Base(const Analysis_Base & copy)
{
  system = copy.system;
  trajectory_list = copy.trajectory_list;
}
  

/*Assignment operator*/
Analysis_Base Analysis_Base::operator =(const Analysis_Base & copy)
{
  if(this!=&copy)
  {
    system = copy.system;
    trajectory_list = copy.trajectory_list;
  }
  
  return *this;
}


/*Methods for employing the loops implemented in class System over various sets of atoms*/

void Analysis_Base::atom_species(int species_index, int atom_type, int atom_index)
{preprocess();system->loop_atom_species(this, species_index, atom_type, atom_index);postprocess();}

void Analysis_Base::atom_species(string species_name, string atomtype_name, int atom_index)
{
  int species_index, atom_type;
  species_index = system->show_species_index(species_name);
  atom_type = system->show_atomtype_index(atomtype_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  if(atom_type==-1){cout<<"Error: atom type "<<atomtype_name<<" not found.\n";exit(1);}
  system->loop_atom_species(this, species_index, atom_type, atom_index);
}

void Analysis_Base::type_molecule(int species_index, int molecule_index, int atom_type)
{preprocess();system->loop_type_molecule(this, species_index, molecule_index, atom_type);postprocess();}

void Analysis_Base::type_molecule(string species_name, int molecule_index, string atomtype_name)
{
  int species_index, atom_type;
  species_index = system->show_species_index(species_name);
  atom_type = system->show_atomtype_index(atomtype_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  if(atom_type==-1){cout<<"Error: atom type "<<atomtype_name<<" not found.\n";exit(1);}
  system->loop_type_molecule(this, species_index, molecule_index, atom_type);
}
		
void Analysis_Base::molecule(int species_index, int molecule_index)
{preprocess();system->loop_molecule(this, species_index, molecule_index);postprocess();}

void Analysis_Base::molecule(string species_name, int molecule_index)
{
  int species_index;
  species_index = system->show_species_index(species_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  system->loop_molecule(this, species_index, molecule_index);
}

void Analysis_Base::species(int species_index)
{preprocess();system->loop_species(this,species_index);postprocess();}


void Analysis_Base::species(string species_name)
{
  int species_index;
  species_index = system->show_species_index(species_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  system->loop_species(this,species_index);
}


void Analysis_Base::type_species(int species_index, int atom_type)
{preprocess();system->loop_type_species(this,species_index,atom_type);postprocess();}

void Analysis_Base::type_species(string species_name, string atomtype_name)
{
  int species_index, atom_type;
  species_index = system->show_species_index(species_name);
  atom_type = system->show_atomtype_index(atomtype_name);
  if(species_index==-1){cout<<"Error: species "<<species_name<<" not found.\n";exit(1);}
  if(atom_type==-1){cout<<"Error: atom type "<<atomtype_name<<" not found.\n";exit(1);}
  system->loop_type_species(this,species_index,atom_type);
}


void Analysis_Base::type_system(int atom_type)
{preprocess();system->loop_type_system(this,atom_type);postprocess();}


void Analysis_Base::type_system(string atomtype_name)
{
  int atom_type;
  atom_type = system->show_atomtype_index(atomtype_name);
  if(atom_type==-1){cout<<"Error: atom type "<<atomtype_name<<" not found.\n";exit(1);}
  system->loop_type_system(this,atom_type);
}


void Analysis_Base::all()
{system->loop_system(this);}


void Analysis_Base::analyze(string runline)
{
  Tokenize tokenize;
  
  int n_args;
  int args_used=0;
  int args_needed;
  string args [ARGMAX];
  string command;
  n_args = tokenize(runline,args);
  
  if(n_args==0){cout << "Error: No atom set command found.";exit(1);}
  
  preprocess();
  
  while(args_used<n_args)
  {
    command = args[args_used];
    args_used++;
    if(command == "atom_species")
    {
      args_needed = atom_species();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
      atom_species(args[args_used],args[args_used+1],atoi(args[args_used+2].c_str()));
      args_used+=args_needed;
    }
    else if(command == "type_molecule")
    {
      args_needed = type_molecule();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
      
      type_molecule(args[args_used],atoi(args[args_used+1].c_str()),args[args_used+2]);
      args_used+=args_needed;
    }
    else if(command == "molecule")
    {
      args_needed = molecule();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
      molecule(args[args_used],atoi(args[args_used+1].c_str()));
      args_used+=args_needed;
    }
    else if(command == "species")
    {
      args_needed = species();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
    
      species(args[args_used]);
      args_used+=args_needed;
    }
    else if(command == "type_species")
    {
      args_needed = type_species();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
    
      type_species(args[args_used],args[args_used+1]);
      args_used+=args_needed;
    }
    else if(command == "type_system")
    {
      args_needed = type_system();
      if(n_args-args_used<args_needed){cout<<"Error: Insufficient arguements for command " << command << ".";exit(1);}
    
      type_system(args[args_used]);
      args_used+=args_needed;
    }
    else if(command == "all")
    {
      args_needed = 0;
    
      all();
    }
    else
    {
      cout << "\nAtom set command "<<command<< " for analysis not understood.\n";
      exit(1);
    }
  }
  
  postprocess();
}

void export_Analysis_Base(py::module& m)
    {
    py::class_<Analysis_Base, std::shared_ptr<Analysis_Base> >(m,"Analysis_Base")
    ;
    }
