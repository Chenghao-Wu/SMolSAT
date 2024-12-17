import sys
sys.path.append('/home/ubuntu/packages/SMolSAT/source')
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
ss.read_trajectory(type='custom',file='../example/bead_spring_polymer/trajectory_KG_bulk_T1500.prd.custom')

list_=SMolSAT.Trajectories(system=ss)
# create particle trajectory for all atoms
list_.create_list(name="all",    args="all")

msid=SMolSAT.mean_square_internal_distance(system=ss,
                        trajs=list_,
                        species='polymer',
                        beads=[1,0,1,1,1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,1,10,1,11,1,12,1,13,1,14,1,15,1,16,1,17,1,18,1,19],
                        out='test_msid.dat',
                        in_mole=True)

#msid.write("test_msid.dat")