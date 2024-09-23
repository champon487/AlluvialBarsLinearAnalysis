# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 09:17:43 2019

@author: river603
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import csv
import os
import glob
import re
import copy

def zerowaveremove(b,lr,li):
        
    element_to_remove = 0.0 
    
    indexes_to_remove = [i for i, x in enumerate(li) if x == element_to_remove]

    # 要素を削除する
    for index in sorted(indexes_to_remove, reverse=True):
        del li[index]
        del lr[index]
        del b[index]
        
    return b,li,lr
    

path_io = 'C:/Users/Toshiki Iwasaki/北海道大学大学院工学研究院 Dropbox/iwasaki toshiki/temp/linearstability/adami/'
case = "Test"
  

with open(path_io+'ForcedLam.dat','r') as fpin2:
    lines = fpin2.readlines()
        
    b2   = []
    lr1 = []
    li1 = []
    lr2 = []
    li2 = []
    lr3 = []
    li3 = []
    lr4 = []
    li4 = []
        
    for line in lines:
        value = line.replace('\n', '').split()
        b2.append(float(value[0]))
        lr1.append(float(value[1]))
        li1.append(float(value[2]))
        lr2.append(float(value[3]))
        li2.append(float(value[4]))
        lr3.append(float(value[5]))
        li3.append(float(value[6]))
        lr4.append(float(value[7]))
        li4.append(float(value[8]))
         
fpin2.close()
       

plt.figure(figsize=(6,4),tight_layout=True)

ax1 = plt.subplot2grid((1,2),(0,0),colspan=1)
ax2 = plt.subplot2grid((1,2),(0,1),colspan=1)

ax1.set_xlabel('Spatial growth rate, $\lambda_r$')
ax1.set_ylabel('Aspect ratio, $\\beta$')
ax1.set_xlim([-1,1])
ax1.set_ylim([0,30])

ax2.set_xlabel('Nondimensional wavenumber, $\lambda_s$')
#ax2.set_ylabel('Aspect ratio, $\\beta$')
ax2.set_xlim([-1,1])
ax2.set_ylim([0,30])

ax1.plot([0,0],[0,30],'--',color='k',linewidth=0.5)
ax1.plot(lr2,b2,'-',color='g')
ax1.plot(lr4,b2,'-',color='r')
ax1.plot(lr1,b2,'--',color='b')
ax1.plot(lr3,b2,'--',color='k')
ax2.plot(li2,b2,'-',color='g')
ax2.plot(li4,b2,'-',color='r')
ax2.plot(li1,b2,'--',color='b')
ax2.plot(li3,b2,'--',color='k')

'''
plt.plot(lr1,bb1,'.',color='b',markersize=2)
plt.plot(lr2,bb2,'.',color='b',markersize=2)
plt.plot(lr3,bb3,'.',color='b',markersize=2)
plt.plot(lr4,bb4,'.',color='b',markersize=2)
plt.plot(li1,bb1,'.',color='r',markersize=2)
plt.plot(li2,bb2,'.',color='r',markersize=2)
plt.plot(li3,bb3,'.',color='r',markersize=2)
plt.plot(li4,bb4,'.',color='r',markersize=2)
'''
#plt.plot(lr1,bb1,color='b')
#plt.plot(lr2,bb2,color='b')
#plt.plot(lr3,bb3,color='b')
#plt.plot(lr4,bb4,color='b')
#plt.plot(li1,bb1,color='r')
#plt.plot(li2,bb2,color='r')
#plt.plot(li3,bb3,color='r')
#plt.plot(li4,bb4,color='r')

plt.savefig(path_io+case+ "_ForcedBar_Diagram.jpg", dpi=300)
        
plt.close()
    