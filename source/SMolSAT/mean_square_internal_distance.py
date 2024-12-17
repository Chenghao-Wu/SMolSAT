#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np

class mean_square_internal_distance:
    def __init__(self,  system=None,
                        trajs=None,
                        species=None,
                        beads=None,
                        out=None,
                        in_mole=False,
                        timescheme=0):
        self.system=system
        self.trajs=trajs
        self.out=out
        self.in_mole=in_mole
        self.timescheme=timescheme
        for i in range(len(beads)//2):
            trajs.create_list(name=f"bead_{i}", args=f"atom_species {species} {beads[2*i]} {beads[2*i+1]}")
        data_all=[]
        for bead_spacing_i in range(1,len(beads)//2):
            data=[]
            for i in range(0,len(beads)//2-bead_spacing_i):
                self.analysis=MeanSquared_Distance(system,in_mole)
                print(f"bead_{i}",f"bead_{i+bead_spacing_i}")
                self.analysis.run(trajs,f"bead_{i}",f"bead_{i+bead_spacing_i}")
                _data=np.array(self.analysis.get())
                data.append(np.mean(_data))
            data_all.append(np.mean(data))
        self.data=np.array(data_all)
        self.write(out)
        
    def get(self):
        return self.data[['dist']]
    
    def get_t(self):
        return self.data["t"]

    def write(self, filename):
        with open(filename, 'w') as f:
            f.write("bead_spacing\tmean_square_internal_distance\n")  # Write header
            for i in range(len(self.data)):
                f.write(f"{i+1}\t{self.data[i]}\n")
                
