#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

from ._SMolSAT import *

import numpy as np
import matplotlib.pyplot as plt
import mpltex


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

    def plot(self,file=None):
        
        linestyle=mpltex.linestyles(hollow_styles=[True],lines=["-"],markers=[])
        #create the plot object
        fig, ax = plt.subplots(nrows=1)

        ax.plot(self.data["t"],self.data["dist"],**next(linestyle),label=None,zorder=3)
        
        #set the plot parameters
        ax.legend(loc="best",ncol=1)
        #ax.set_xlim(90,1100)
        #ax.set_ylim(1e-7,1e-3)
        ax.set_xlabel(r"time")
        ax.set_ylabel(r"mean squared distance")
        fig.tight_layout(pad=0.1) 
        if file!=None:
            fig.savefig(file)
