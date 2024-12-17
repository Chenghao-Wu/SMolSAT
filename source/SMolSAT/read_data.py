#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:58:21 2018

@author: bruce
"""

import numpy as np

def read_msd(msd_file):
    with open(msd_file) as f:
        lines = (line.strip() for line in f if len(line.strip().split())==2)
        data = np.loadtxt(lines, delimiter='\t',dtype={'names': ('t', 'msd'),'formats': ('f4','f4')})
    return data