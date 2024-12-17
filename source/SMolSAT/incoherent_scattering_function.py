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


class isfs:
    def __init__(self,system=None,trajs=None,listname=None,out=None,plane=None,index1=None,index2=None,max_length_scale=None,timescheme=-1):
        self.system=system
        self.trajs=trajs
        self.listname=listname
        self.out=out
        self.wv=Wave_Vectors(system,plane,max_length_scale)
        self.analysis=Incoherent_Scattering_Function(system,self.wv,index1,index2,timescheme)
        self.analysis.run(trajs,listname)
        self.analysis.write(out)
        self.read()

    def get(self):
        return self.data["isfs"]
    
    def get_q(self):
        return self.data["t"]

    def read(self):
        with open(self.out) as f:
            lines = (line.strip() for line in f if len(line.strip().split())==2)
            self.data = np.loadtxt(lines, delimiter='\t',skiprows=2,dtype={'names': ('t','isfs'),'formats': ('f4','f4')})

    def plot(self,file=None):
        
        linestyle=mpltex.linestyles(hollow_styles=[True],lines=["-"],markers=[])
        #create the plot object
        fig, ax = plt.subplots(nrows=1)

        ax.plot(self.data["t"],self.data["isfs"],**next(linestyle),label=None,zorder=3)
        
        #set the plot parameters
        ax.legend(loc="best",ncol=1)
        #ax.set_xlim(90,1100)
        #ax.set_ylim(1e-7,1e-3)
        ax.set_xlabel(r"t")
        ax.set_ylabel(r"isfs")
        fig.tight_layout(pad=0.1) 
        if file!=None:
            fig.savefig(file)
