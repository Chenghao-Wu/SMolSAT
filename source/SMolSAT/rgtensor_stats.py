#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np


class rg_stats:
    def __init__(self,  system=None,
                        trajs=None,
                        listname=None,
                        out=None):
        self.system=system
        self.trajs=trajs
        self.listname=listname
        self.out=out

        self.analysis=RgTensor_Stats(system)
        self.analysis.run(trajs,listname)
        self.analysis.write(out)
        self.read()

    def get(self):
        return self.data[['xx','yy','zz']]
    
    def get_t(self):
        return self.data["t"]

    def read(self):
        with open(self.out) as f:
            lines = (line.strip() for line in f if len(line.strip().split())==4)
            self.data = np.loadtxt(lines, delimiter='\t',dtype={'names': ('t', 'xx','yy','zz'),'formats': ('f4','f4','f4','f4')})

    def calc_gyration_rad_dist(self,out=None,n_bin=None,r_cut=None):
        self.analysis.calc_gyration_rad_dist(n_bin,r_cut)
        self.write_gyration_rad_dist(out)
    
    def calc_rel_asphericity_dist(self,out=None,n_bin=None):
        self.analysis.calc_rel_asphericity_dist(n_bin)
        self.write_rel_asphericity_dist(out)