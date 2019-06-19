# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 16:40:38 2018

@author: allen
"""
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from pylab import * 
# Parameters.
x = np.linspace(4,100,2000)/10**(6)
# question a
ep = 0.95
def L(lab):
    h = 6.63 * (10**(-34))
    c = 3*(10**8)
    k = 1.38 *(10 **(-23))
    T = 291
    sigma = 5.67 * (10**(-8))
    x = h*c/k/lab/T
    L = 2 * h *(c**2)/(lab**5)/(exp(x)-1)
    return L

P = L(x)
plt.figure
plt.plot(x*10**6,P)
plt.xlabel('Wavelength  , micrometer')
plt.ylabel('Spectral radiation W/(m^2/μm)')
plt.title('Spectral Radiation in TUDelft')
plt.show()

# question b
x = np.linspace(11,15,1000)/10**6
R = np.zeros(len(x))
import scipy.integrate as integrate
for i in range(len(x)):
    R[i] = integrate.quad(L , 11*10**(-6) , x[i])[0]*ep*pi

plt.figure
plt.plot(x*10**6,R)
plt.xlabel('Wavelength  , μm')
plt.ylabel('Total Radiation W/m^2')
plt.title('Total Radiation in Delft with Wavelength variance')
plt.show()
