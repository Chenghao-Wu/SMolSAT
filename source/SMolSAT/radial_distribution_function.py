#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np


class rdf:
    def __init__(self,system=None,nbins=None,max_length_scale=None,timescheme=-1,trajs=None,listname=None,listname1=None,listname2=None,out=None,is_inter=False):
        self.system=system
        self.trajs=trajs
        self.listname=listname
        self.out=out
        
        self.analysis=Radial_Distribution_Function(system,nbins,timescheme,max_length_scale,is_inter)
        if listname!=None:
            self.analysis.run(trajs,listname)
        else:
            self.analysis.run(trajs,listname1,listname2)
        self.analysis.write(out)
        self.read()

    def get(self):
        return self.data["g(r)"]
    
    def get_q(self):
        return self.data["bin"]

    def read(self):
        with open(self.out) as f:
            lines = (line.strip() for line in f if len(line.strip().split())==2)
            self.data = np.loadtxt(lines, delimiter='\t',skiprows=2,dtype={'names': ('bin','g(r)'),'formats': ('f4','f4')})