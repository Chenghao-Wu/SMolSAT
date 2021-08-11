/*Soft-Matter Molecular Simulation Analysis Toolkit (SMolSAT)*/
/*Written by Zhenghao Wu, modified from AMDAT written by David S. Simmons*/

#include <iostream>
#include <stdio.h>
#include "version.h"

using namespace std;

#include "progress.h"

void print_progress(int done,int total)
{
	float percent;
	//if(!fileoutput)
	//{
	  percent = int(100*float(done)/float(total));
	  cout.width(3);
	  cout << "\b\b\b\b" << percent << "%";
	  cout.flush();
	//}
}