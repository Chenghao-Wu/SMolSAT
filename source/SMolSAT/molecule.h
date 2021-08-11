/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#ifndef MOLECULE
#define MOLECULE


#include "atom_trajectory.h"
#include "coordinate.h"
#include "trajectory.h"
#include "multibody.h"

#include "../extern/pybind11/include/pybind11/pybind11.h"

//namespace std {

class Molecule
{
    int total_atoms;
    int n_atomtypes;
    int * n_atoms;
    
    int moleculeID;
    int species;		//stores species of molecule
    int n_timesteps;
    
    bool atoms_unwrapped;
    
    Atom_Trajectory** atoms;
    void clear_memory();			//method to clear memory allocated to arrays so that they may be recreated
    
  public:
    Molecule(int types=0,int ts=0);
    Molecule(int,int*,int);
    ~Molecule();
    void set(int, int *,int);
    int typecount(int);				//method to return number of atoms of a given type
    int atomcount();				//method to return number of atoms of a given type
    float * msd(int, int, int);		//method to pass on msd of a particular atom
    Atom_Trajectory show_atom(int,int);	//method to return atom_trajectory object in molecule
    int unwrap_atoms(const Coordinate &);			//unwrap all atom trajectories (only works w/ constant volume)
    int wrap_atoms(const Coordinate *, Coordinate **);				//wrap all atom trajectories				
    
    Coordinate show_unwrapped(int atom_type, int atom_index, int timestep)const;	//show unwrapped coordinate of a given atom at a given time
    
    Atom_Trajectory* show_atom_trajectory(int, int);
    float atom_distance(int,int,int,int);
    
    void set_coordinate(int, int, const Coordinate &, int);
    void set_unwrapped(int, int, const Coordinate &, int);
    void set_velocity(int, int, const Coordinate &, int);
    void set_moleculeID(int ID){moleculeID=ID;};	//set unique molecule ID
    void set_species(int spec){species=spec;};		//set species
    int show_species()const{return species;};		//return species index
    
    void ID_to_atoms();		//Pass unique ID of molecule down to constituent atoms
    
    Multibody create_multibody(System * sys, Coordinate boxsize) const;	//return multibody consisting of all atoms in molecule, requires box size at t=0
    Multibody create_multibody(System * sys, int typeii, Coordinate boxsize) const;	//return multibody consisting of all atoms of type typeii in molecule, requires box size at t=0
    Multibody create_multibody(System * sys, int n_bodies, int * typeii, int * index, Coordinate boxsize) const; //return multibody consisting of atoms specified by a list of types and indices, requires box size at t=0
};
void export_Molecule(pybind11::module& m);
//}

#endif