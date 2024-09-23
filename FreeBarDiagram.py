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

path_io = 'C:/Users/Toshiki Iwasaki/北海道大学大学院工学研究院 Dropbox/iwasaki toshiki/temp/linearstability/adami/'
case = "Test"

with open(path_io+'contor.dat','r') as fpin:
    
    line = fpin.readline()
    head = line.replace('\n', '').split()
    
    beta = float(head[0])
    
    line = fpin.readline()
    head = line.replace('\n', '').split()
    
    nx = int(head[0])-1
    ny = int(head[1])-1
    
    print(nx,ny)
    
    lines = fpin.readlines()
        
    l   = []
    b   = []
    ome   = []
    some  = []
        
    for line in lines:
        value = line.replace('\n', '').split()
        l.append(float(value[0]))
        b.append(float(value[1]))
        ome.append(float(value[2]))
        some.append(float(value[3]))
        
fpin.close()
        
x   = [[0 for i in range(nx+1)] for j in range(ny+1)]
y   = [[0 for i in range(nx+1)] for j in range(ny+1)]
Omega = [[0 for i in range(nx+1)] for j in range(ny+1)]
omega = [[0 for i in range(nx+1)] for j in range(ny+1)]
    
ij = 0

for j in range(ny+1):
     for i in range(nx+1):
         x[j][i]   = l[ij]
         y[j][i]   = b[ij]
         Omega[j][i] = ome[ij]
         omega[j][i] = some[ij]
                
         ij = ij+1         

dOme = [0 for i in range(nx)]

b_omax = []
l_omax = []

for j in range(0,ny+1,4):
    for i in range(nx):
        dOme[i] = abs(float(Omega[j][i+1])-float(Omega[j][i]))
        
    iome = dOme.index(min(dOme))
    
    if float(Omega[j][iome])>0:
        l_omax.append(x[j][iome])
        b_omax.append(y[j][iome])
        

plt.figure(figsize=(4.5,4),tight_layout=True)
plt.xlabel('Nondimensional wavenumber, $\lambda$')
plt.ylabel('Aspect ratio, $\\beta$')
#plt.subplots_adjust(bottom=0.15)
plt.axis([0,1.75,0.,30])

levels = np.linspace(-2,2,31)
labels = np.linspace(-2,2,5)
         
plt.contour(x, y, Omega, colors=['black'], levels=[0],lw=5)
plt.contour(x, y, omega, colors=['grey'], levels=[0],lw=5)
plt.contourf(x, y, Omega, levels, cmap=cm.bwr, extend="both")
        
plt.colorbar(ticks=labels, label='Growth rate, $\Omega$', fraction=0.03)

plt.plot(l_omax,b_omax,color='k',linestyle="dashed")

plt.plot([0,1.75],[beta,beta],'-',color='yellow')


plt.savefig(path_io+case+ "_FreeBar_Diagram.jpg", dpi=300)
        
plt.close()


with open(path_io+'Omega-lam.dat','r') as fpin:
    
    line = fpin.readline()
    head = line.replace('\n', '').split()
    
    beta = float(head[0])
    
    line = fpin.readline()
    head = line.replace('\n', '').split()
    
    nx = int(head[0])-1
        
    lines = fpin.readlines()
        
    l   = []
    ome   = []
    some  = []
        
    for line in lines:
        value = line.replace('\n', '').split()
        l.append(float(value[0]))
        ome.append(float(value[1]))
        some.append(float(value[2]))
        
fpin.close()

lam = np.array(l)
growth = np.array(ome)
celerity = np.array(some)

plt.figure(figsize=(7.,2.5),tight_layout=True)
plt.xlabel('Nondimensional wavenumber, $\lambda$')
plt.ylabel('Growth rate, $\Omega$')
#plt.subplots_adjust(bottom=0.15)
plt.xlim(0,1.75)

plt.plot(lam,growth,'-',color='k')

plt.plot([0,1.75],[0,0],'--',color='k',lw=0.5)

plt.savefig(path_io+case+"GrowthRate.jpg", dpi=300)
        
plt.close()
    