/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include <vector>
#include <iostream>
#include "system.h"
#ifndef BOOLEAN_LIST
#define BOOLEAN_LIST

class Boolean_List
{
  std::shared_ptr<System> system;
  //bool * included;		//array holding a boolean for each trajectory in the system
  std::vector <bool> included;	//array holding a boolean for each trajectory in the system
  void update_size();
  
  public:
  Boolean_List();
  Boolean_List(std::shared_ptr<System> sys);
  Boolean_List(std::shared_ptr<System> sys, int * inc, int n_included);
  Boolean_List(const Boolean_List &); //MEM - copy constructor
  ~Boolean_List();
  Boolean_List operator= (const Boolean_List &);

  void set(std::shared_ptr<System> sys);
  void set(std::shared_ptr<System> sys, int * inc, int n_included);

  bool operator() (int);
  void operator() (int,bool);
  Boolean_List operator&& (Boolean_List &);
  Boolean_List operator|| (Boolean_List &);
  bool operator== (Boolean_List &);
  bool operator!= (Boolean_List &);

  int show_n_included()const;
  int show_trajectory_ids(int, int *);
  int first_included()const;
  int show_size()const;

};
//}
#endif
