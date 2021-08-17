/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include "coordinate.h"
#include "system.h"
#include <string>
#include <iostream>

#ifndef WAVE_VECTORS
#define WAVE_VECTORS

using namespace std;

class Wave_Vectors
{
  protected:
    std::shared_ptr<System> system;

    float wavegrid_spacing;		//spacing of grid of wave-vectors; by default equal to 2pi/L, where L is the minimum system dimension
    float maxrange;			//distance wavegrid extends (symmetrically) from origin along each axis
    float delta_wavenumber;		//thickness of spherical wavenumber shells
    string plane;
    int n_atoms_looped;
    int n_wavenumbers;
   
    
    //int * n_wavevectors;	//number of wavevectors for each wavenumber
    //float * approx_wavenumber;
    //int * binsize;		//stores allocation size of each 'bin' of wavevectors
    //Coordinate ** wavevector;	//list of wavevectors by wavenumber bin
    
    vector<float> approx_wavenumber;
    vector<vector<Coordinate>> wavevector; 


    void calculate (std::shared_ptr<System> sys, int shellcount);
    //void calculate (Coordinate boxsize, float deltak, float kmax, int maxvectors);

    virtual void bin (int, int, int){};
    void read_vectors();
    void read_vectors_2d(string);
    void read_vectors_1d(string);


  public:
    Wave_Vectors();
    ~Wave_Vectors(); //destructor
    Wave_Vectors(const Wave_Vectors &); //copy constructor
    Wave_Vectors(std::shared_ptr<System> sys);
    Wave_Vectors(std::shared_ptr<System> sys, int shellcount=300);
    Wave_Vectors(std::shared_ptr<System> sys, string plane, float max_length_scale=0);	//used for structure factor calculation
  
    Wave_Vectors operator = (const Wave_Vectors &); //assignment operator
    
    vector<Coordinate> vectorlist(int index)const{return wavevector[index];};
    int vectorcount(int index)const{return wavevector[index].size();};
    float show_delta_wavenumber()const{return delta_wavenumber;};
    int show_n_wavenumbers()const{return wavevector.size();};
    float show_approx_wavenumber(int index)const{return approx_wavenumber[index];}
    float show_mean_wavenumber(int index)const;		//return average wavenumber of vectors corresponding to wavenumber
    float show_stdev_wavenumber(int index)const;
    Coordinate show_mean_wavevector(int index)const;
    
};

void export_Wave_Vectors(pybind11::module& m);

#endif

/*Specify system length and the corresponding wavenumber gives both the minimum wavenumber and the wavenumber increment for the grid.  The wavenumber increment for spherical shells is then half of this.  We'll use a symetrical wavenumber grid even if the system is asymmetric, using the minimum dimension as the relevant system dimension.*/
