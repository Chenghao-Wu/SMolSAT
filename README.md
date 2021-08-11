# SMolSAT.py: Soft-Matter Molecular Simulation Analysis Toolkit

### Install Requirement
* cmake 3.12
* g++ 7.0 

### How to install?
```bash
cd source
mkdir build
cd build
cmake ../.
make
make install
```

### How to use SMolSAT.py?
First, we need to include the package in the python PATH
```python
import sys
sys.path.append('install_path/source')
import SMolSAT
```

Then, you can try SMolSAT.py
```python
# define your system
sys=SMolSAT.System(ensemble='nv') 
# define timetype of your trajectory: linear or exponential
sys.set_exponential_timetype(n_blocks=1 ,block_size=62 ,exp_base=1.2,time_unit=0.0002)
# define all atom types appear in the system
sys.atomtype_list=["1","2","3"] 
# define the species of your interest
sys.add_species(name="name",number=1,atoms=[1, 0, 0])
# read the trajectory: only support lammps custom type now
sys.read_trajectory(type='custom',file='')

# define particle trajectories for analysis
particle_trjs=SMolSAT.Trajectories(system=sys)
# support choose all, type, species_type, atom_type, atom_index 
particle_trjs.create_list(name="name", args="type_system 1") 

# start analysis
msd=SMolSAT.msd(system=sys,trajs=particle_trjs,listname="name",out="msd.dat")
#msd.plot()
```
### Analysis Methods
1. gyration tensor
2. end-to-end distance
3. mean square displacement
4. non-gaussian parameter
5. radial distribution function
6. structure factor
7. incoherent scattering function
8. bond autocorrelation function
9. van hove correlation function
   * self part
   * distinct part