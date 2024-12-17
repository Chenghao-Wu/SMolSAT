#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np

class rg2:
    def __init__(self,  system=None,
                        trajs=None,
                        listname=None,
                        out=None,tensor=False):
        self.system=system
        self.trajs=trajs
        self.listname=listname
        self.out=out

        self.analysis=Radius_Gyration(system)
        self.analysis.run(trajs,listname)
        self.analysis.write(out)
        
        self.read()
        if tensor:
            self.analysis.write_tensor(out+"_tensor")
            self.read_tensor()

    def get(self):
        return self.data[['rg']]
    
    def get_t(self):
        return self.data["t"]

    def read(self):
        with open(self.out) as f:
            lines = (line.strip() for line in f if len(line.strip().split())==2)
            self.data = np.loadtxt(lines, delimiter='\t',dtype={'names': ('t', 'rg'),'formats': ('f4','f4')})

    def read_tensor(self):
        with open(self.out+"_tensor") as f:
            lines = (line.strip() for line in f if len(line.strip().split())==7)
            self.tensor = np.loadtxt(lines, delimiter='\t',dtype={'names': ('t', 'xx', 'xy', 'xz', 'yy', 'yz', 'zz'),'formats': ('f4','f4','f4','f4','f4','f4','f4')})

