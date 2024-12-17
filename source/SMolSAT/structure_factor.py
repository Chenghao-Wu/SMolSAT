#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np


class structure_factor:
    def __init__(self,system=None,plane=None,max_length_scale=None,timescheme=-1,trajs=None,listname=None,listname1=None,listname2=None,out=None):
        self.system=system
        self.trajs=trajs
        self.listname=listname
        self.out=out

        self.analysis=Structure_Factor(system,plane,max_length_scale,timescheme)
        if listname!=None:
            self.analysis.run(trajs,listname)
        else:
            self.analysis.run(trajs,listname1,listname2)
        self.analysis.write(out)
        self.read()

    def get(self):
        return self.data["S(q)"]
    
    def get_q(self):
        return self.data["Mean_q"]

    def read(self):
        with open(self.out) as f:
            lines = (line.strip() for line in f if len(line.strip().split())==7)
            self.data = np.loadtxt(lines, delimiter='\t',skiprows=2,dtype={'names': ('Nominal_q','Mean_q','Stddev_q','Mean_qx','Mean_qy','Mean_qz','S(q)'),'formats': ('f4','f4','f4','f4','f4','f4','f4')})
