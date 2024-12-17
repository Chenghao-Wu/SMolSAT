#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np

class vhs:
    def __init__(self,system=None,trajs=None,listname=None,out=None,nbins=None,r_cut=None):
        self.system=system
        self.trajs=trajs
        self.listname=listname
        self.out=out
        self.nbins=nbins
        self.r_cut=r_cut

        self.analysis=Van_Hove_Self(system,nbins,r_cut)
        self.analysis.run(trajs,listname)
        self.analysis.write(out)
        self.read()

    def get(self):
        return self.data["ngp"]
    
    def get_t(self):
        return self.data[:,0]

    def read(self):
        with open(self.out) as f:
            line_header = (line.strip() for line in f if len(line.strip().split())==self.nbins)
            self.r_bins=np.loadtxt(line_header, delimiter='\t')
            
        with open(self.out) as f:
            lines = (line.strip() for line in f if len(line.strip().split())==self.nbins+1)
            self.data = np.loadtxt(lines, delimiter='\t')