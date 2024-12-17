#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np

class mean_squared_distance:
    def __init__(self,  system=None,
                        trajs=None,
                        listname1=None,
                        listname2=None,
                        out=None,
                        in_mole=False,
                        timescheme=0):
        self.system=system
        self.trajs=trajs
        self.listname1=listname1
        self.listname2=listname2
        self.out=out
        self.in_mole=in_mole
        self.timescheme=timescheme
        self.analysis=MeanSquared_Distance(system,in_mole)
        self.analysis.run(trajs,listname1,listname2)
        self.analysis.write(out)
        self.read()

    def get(self):
        return self.data[['dist']]
    
    def get_t(self):
        return self.data["t"]

    def read(self):
        with open(self.out) as f:
            lines = (line.strip() for line in f if len(line.strip().split())==2)
            self.data = np.loadtxt(lines, delimiter='\t',dtype={'names': ('t', 'dist'),'formats': ('f4','f4')})
