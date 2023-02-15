#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 16:18:47 2021

@author: linnea
"""

 ## code to run for tglv
import pytest
import numpy as np
import matplotlib.pyplot as plt
from scipy import io
import pyssage
from pyssage import quadvar


#runfile('/Users/linnea/Google Drive/Doktorand!/CML/pyssage-master/pyssage/quadvar.py', wdir='/Users/linnea/Google Drive/Doktorand!/CML/pyssage-master/pyssage')

list_s = [1, 2, 3, 5, 10, 50]

for svals in list_s:

    path = "/Users/linnea/Documents/MATLAB/D101_r20_grids_s{}.mat".format(svals)
    
    print(path)
    
    gridmat = io.loadmat(path)
    
    grid = gridmat['grids']
    
    k=100
    
    tlqv9_D101_r20_s1 = np.zeros((k, 50, 2))
    
    j=0
    for i in range(400, 400+k):
        tlqv9_D101_r20_s1[j] = quadvar.four_tlqv(grid[i], 1,  0, 1, 1)
        j=j+1
        print(j)
    
    savefile = "tlqv4_D101_s{}_r20.mat".format(svals)
    matlabel = "tlqv4_D101_s{}_r20".format(svals)
    print(matlabel)
    
    io.savemat(savefile, mdict={matlabel: tlqv9_D101_r20_s1})

#plt.plot(stat_4gridarray[0,:,0],stat_4gridarray[0,:,1], '*')

