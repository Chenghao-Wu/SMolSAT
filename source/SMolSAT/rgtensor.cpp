/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

/*need to incorporate center of mass part*/

#include "rgtensor.h"
#include "math.h"
#include <stdlib.h>
#include "version.h"
#include "system.h"

using namespace std;
using namespace Eigen;

RgTensor::RgTensor(std::shared_ptr<System> sys,Trajectory*traj)
{
  int timeii;
  
  system=sys;
  trajectory=traj;
  
  n_blocks = system->show_n_exponentials();
  
  rgtensor = new sixfloat * [n_blocks];
  maxtimes = system->show_n_timegaps();
  
  blocks_per_time = new int [maxtimes];
  mean_eigenvalues = new threefloat [maxtimes];
  
  blocksize = new int [n_blocks];
  for(timeii=0;timeii<maxtimes;timeii++)
  {
    blocks_per_time[timeii]=0;
    mean_eigenvalues[timeii][0]=0;
    mean_eigenvalues[timeii][1]=0;
    mean_eigenvalues[timeii][2]=0;
  }
  
  
  calc_tensor();
  calc_eigenvalues();
}


void RgTensor::calc_tensor()
{
  int blockii;
  int * timelist;
  float * times;
  
  sixfloat * temptensor;
  
  Coordinate * coordinatelist;
  
  for(blockii=0;blockii<n_blocks;blockii++)
  {
    
    blocksize[blockii] = system->big_blocksize(blockii);
    timelist = new int [blocksize[blockii]];
    times = new float [blocksize[blockii]];
    system->big_block(blockii,timelist);	//list of times for block (counting allowed inter-block timegaps)
    system->show_times(blocksize[blockii],timelist,times);
    
    temptensor = new sixfloat [blocksize[blockii]];
    coordinatelist=trajectory->show_coordinates(timelist,blocksize[blockii]);	//get list of coordinates	
    gyration_tensor(coordinatelist,times,blocksize[blockii],temptensor);		//calculate gyration tensor as function of time
    rgtensor[blockii]=temptensor;			//copy gyration tensor into long-term memory
  }
  
  
}

void RgTensor::calc_eigenvalues()
{
  int blockii,timeii;
  Matrix3f rgarray;
  eigenvalues = new threefloat * [n_blocks];
  for(blockii=0;blockii<n_blocks;blockii++)
  {
    eigenvalues[blockii] = new threefloat [blocksize[blockii]];
    for(timeii=0;timeii<blocksize[blockii];timeii++)
    {
      /*create Eigen Matrix corresponding to current gyration tensor*/
      rgarray(0,0)=rgtensor[blockii][timeii][0];
      rgarray(1,0)=rgtensor[blockii][timeii][1];
      rgarray(2,0)=rgtensor[blockii][timeii][2];
      rgarray(0,1)=rgtensor[blockii][timeii][1];
      rgarray(1,1)=rgtensor[blockii][timeii][3];
      rgarray(2,1)=rgtensor[blockii][timeii][4];
      rgarray(0,2)=rgtensor[blockii][timeii][2];
      rgarray(1,2)=rgtensor[blockii][timeii][4];
      rgarray(2,2)=rgtensor[blockii][timeii][5];
      
      eigkernel(rgarray,blockii,timeii);
      
      mean_eigenvalues[timeii][0]+=eigenvalues[blockii][timeii][0];
      mean_eigenvalues[timeii][1]+=eigenvalues[blockii][timeii][1];
      mean_eigenvalues[timeii][2]+=eigenvalues[blockii][timeii][2];
      
      blocks_per_time[timeii]++;
    }
  }
  
  for(timeii=0;timeii<maxtimes;timeii++)
  {
    mean_eigenvalues[timeii][0]/=blocks_per_time[timeii];
    mean_eigenvalues[timeii][1]/=blocks_per_time[timeii];
    mean_eigenvalues[timeii][2]/=blocks_per_time[timeii];
  }
}

void RgTensor::eigkernel(Matrix3f rgarray, int blockii, int timeii)
{
  //Vector1f eigvals;
  
  SelfAdjointEigenSolver<Matrix3f> eigensolver(rgarray);
  //JAMA::Eigenvalue <float> eigmat(rgarray);
  //eigmat.getRealEigenvalues(eigvals);
      
  eigenvalues[blockii][timeii][0]=eigensolver.eigenvalues()(0);
  eigenvalues[blockii][timeii][1]=eigensolver.eigenvalues()(1);
  eigenvalues[blockii][timeii][2]=eigensolver.eigenvalues()(2);
  
}


void RgTensor::write(string filename)
{
  int timeii;
  float * times;
  
  ofstream output(filename.c_str());
  output << "Single particle gyration tensor data created by SMolDAT v." << VERSION << "\n"; 
  
  times = system->displacement_times();
  
  for(timeii=0;timeii<maxtimes;timeii++)
  {
    //output << times[timeii] << "\t" << mean_eigenvalues[timeii][0] << "\t" << mean_eigenvalues[timeii][1] << "\t" << mean_eigenvalues[timeii][2] << "\n";
    output << times[timeii] << "\t" << eigenvalues[0][timeii][0] << "\t" << eigenvalues[0][timeii][1] << "\t" << eigenvalues[0][timeii][2] << "\n";
  }
}


void RgTensor::write(ofstream& output)const
{
  int timeii;
  float * times;
  
  output << "Single particle gyration tensor data created by SMolDAT v." << VERSION << "\n"; 
  
  times = system->displacement_times();
  
  for(timeii=0;timeii<maxtimes;timeii++)
  {
    //output << times[timeii] << "\t" << mean_eigenvalues[timeii][0] << "\t" << mean_eigenvalues[timeii][1] << "\t" << mean_eigenvalues[timeii][2] << "\n";
    output << times[timeii] << "\t" << eigenvalues[0][timeii][0] << "\t" << eigenvalues[0][timeii][1] << "\t" << eigenvalues[0][timeii][2] << "\n";
  }
}