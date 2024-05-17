# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
S0=75 #mm
k=0.001 # /d
t0=0
tmax=200
steps=200
dt=1
a=1
time_series= np.arange(t0,tmax+1,dt)

def deri(S):
    dSdt=(-k)*S**a
    return dSdt

S=np.ones(steps+1)
S[0]=S0
for t in range(steps):
    S[t+1]=S[t]+deri(S[t])
print(S)
ax1=plt.subplot(211)
ax1.plot(time_series,S)
plt.title('linear reservor without external forces')
plt.xlabel('Time Series  /day')
plt.ylabel('Storage')

#------------------------------------------#
pre=np.loadtxt('C:/Users/lzhi/data.txt')
pre=pd.DataFrame(pre)
preci=pre.iloc[:,3]
eva=pre.iloc[:,4]
steps=668
t0=0
t=np.arange(t0,669,1)
S2=np.zeros(steps+2)
dSdt=np.zeros(steps+1)
S2[0]=S0
for i in t:
    dSdt[i]=(-k)*S2[i]**a+preci[i]-eva[i]
    S2[i+1]=S2[i]+dSdt[i]
ax2=plt.subplot(212)
ax2.plot(t,S2[:-1])
plt.title('linear reservor with external forces')
plt.xlabel('Time Series /day')
plt.ylabel('Storage')
plt.show()

real_S=np.ones(steps+1)
for j in range(len(time_series)):
    real_S[j]=np.exp(-k*t[j]+np.log(S0))
ax1.plot(time_series,real_S[:201])
plt.figure
plt.plot(time_series,S,time_series,real_S[:201])
plt.title('Comparison of analytical solution and numerical solution without forces')
plt.show()

plt.figure
plt.plot(t,S2[:-1],t,real_S)
plt.title('Comparison of analytical solution and numerical solution without forces')
plt.show()