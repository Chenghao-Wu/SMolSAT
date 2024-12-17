# SMolSAT.py: Soft-Matter Molecular Simulation Analysis Toolkit

## [Documentation](https://zhenghaowu.notion.site/SMolSAT-py-ba41b0ef45404db6aa1ef6bb518af701)

### Install Requirement
* cmake 3.12
* g++ 7.0 
* python 3.x+
   * numpy
  
  
### How to install?
```bash
cd source
mkdir build
cd build
cmake -DPYTHON_EXECUTABLE=$(python3 -c "import sys; print(sys.executable)") ../.
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

Then, you can try SMolSAT.py for the **example** bead-spring polymer system
```python
import sys
sys.path.append('/home/zwu/pyAMDAT/source/')
import SMolSAT

# System definition
ss=SMolSAT.System(ensemble='nv')
# set exponential time scheme
ss.set_exponential_timetype(n_blocks=1 ,block_size=43 ,exp_base=1.2,time_unit=0.01)

# There are 500 Kremer-Grest bead-spring polymer chains
# Each chain possesses 20 beads
# The trajectory type is LAMMPS CUSTOM

num_chains=500
num_bead=20

# set atom type list
ss.atomtype_list=["1"]
# add species
ss.add_species(name="polymer",number=500,atoms=[20])
    
# read trajectory file
ss.read_trajectory(type='custom',file='trajectory_KG_bulk_T1500.prd.custom')

list_=SMolSAT.Trajectories(system=ss)
# create particle trajectory for all atoms
list_.create_list(name="all",    args="all")

# Analysis methods

# inter-molecular radial distribution function
rdf=SMolSAT.rdf(system=ss,nbins=500,max_length_scale=10,timescheme=0,trajs=list_,listname="all",out="rdf_interMole.dat",is_inter=True)

# inter-segmental radial distribution function
rdf=SMolSAT.rdf(system=ss,nbins=500,max_length_scale=10,timescheme=0,trajs=list_,listname="all",out="rdf.dat",is_inter=False)

# radius of gyration
list_.create_multibodies(name="chains",trj_list_name="chain",center_type="centroid", args="all_molecule")
rg=SMolSAT.rg2(system=ss,trajs=list_,listname="chains",out="rg.dat",tensor=True) # argmument tensor enables calculations of gyration tensor out file is out+"tensor"

# structure factor
strfac=SMolSAT.structure_factor(system=ss,plane="xyz",max_length_scale=0,timescheme=0,trajs=list_,listname="all",out="strfac.dat")

# chain end-to-end distance
list_.create_list(name="head_atom", args="atom_species polymer 1 0")
list_.create_list(name="end_atom", args="atom_species polymer 1 19")
end2end_distance2=SMolSAT.mean_squared_distance(system=ss,trajs=list_,listname1="head_atom",listname2="end_atom",out="end2end2.dat",in_mole=True)
end2end_distance=SMolSAT.mean_distance(system=ss,trajs=list_,listname1="head_atom",listname2="end_atom",out="end2end.dat",in_mole=True)

# mean square displacement of segments
msd=SMolSAT.msd(system=ss,trajs=list_,listname="all",out="msd_all.dat")
msd_xy=SMolSAT.msd2d(system=ss,trajs=list_,listname="all",out="msd_xy.dat",plane="xy")
msd_yz=SMolSAT.msd2d(system=ss,trajs=list_,listname="all",out="msd_yz.dat",plane="yz")
msd_xz=SMolSAT.msd2d(system=ss,trajs=list_,listname="all",out="msd_xz.dat",plane="xz")

# mean square displacement of chain center-of-mass
# create particle trajectory for all monomers by combining atoms using geometric centroid method
msd=SMolSAT.msd(system=ss,trajs=list_,listname="chains",out="msd_chains.dat")

# self-part intermediate scattering function
# q value equals the first peak position of the structure factor
isfs=SMolSAT.isfs(system=ss,plane="xyz",max_length_scale=12.4385,index1=26,index2=26,timescheme=0,trajs=list_,listname="all",out="isfs.dat")

# non-gassian parameter
ngp=SMolSAT.ngp(system=ss,trajs=list_,listname="all",out="ngp.dat")

# bond vector autocorrelation function
# create bond vectors list
for bondii in range(19):
    list_.create_multibodies(name="bond_"+str(bondii),trj_list_name="bond",center_type="centroid", args="species_atomlist polymer 1 "+str(bondii)+" 1 "+str(bondii+1))
multibodies=[]
for bondii in range(19):
    multibodies.append("bond_"+str(bondii))
list_.combine_multibody_lists(name="bonds",multibodies=multibodies)
list_.combine_trajectories(name="bonds",trjs=multibodies)

baf=SMolSAT.baf(system=ss,plane="xyz",trajs=list_,listname="bonds",out="baf.dat")

# end-to-end vector autocorrelation function
# the end-to-end vector is defined as the vector between the first and end beads of a single polymer chain
list_.create_multibodies(name="e2e",trj_list_name="e2e",center_type="centroid", args="species_atomlist polymer 1 0 1 19")
e2e_acf=SMolSAT.baf(system=ss,plane="xyz",trajs=list_,listname="e2e",out="e2e_acf.dat")
```
### Analysis Methods
1. gyration tensor
2. end-to-end distance
3. mean square internal distance
4. mean square displacement
5. non-gaussian parameter
6. radial distribution function (rdf)
7. structure factor (strfac)
8. incoherent scattering function (isfs)
9. bond autocorrelation function (baf)
10. van hove correlation function (vhc)
   * self part
   * distinct part
