#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import numpy as np
import mpltex
import sys


#@mpltex.aps_decorator
def plot_msd(t=None,msd=None,save=False):
    linestyle=mpltex.linestyles(hollow_styles=[True],lines=["-"],markers=[])
    #create the plot object
    fig, ax = plt.subplots(nrows=1)

    ax.plot(t,msd,**next(linestyle),label=None,zorder=3)
    
    #set the plot parameters
    ax.legend(loc="best",ncol=1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    #ax.set_xlim(90,1100)
    #ax.set_ylim(1e-7,1e-3)
    ax.set_xlabel(r"time")
    ax.set_ylabel(r"msd")
    fig.tight_layout(pad=0.1) 
    if save:
        fig.savefig("msd.png")
