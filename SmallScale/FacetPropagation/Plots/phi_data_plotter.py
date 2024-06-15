import numpy as np
from matplotlib import pyplot as plt
import math
import pandas as pd
from scipy.interpolate import make_interp_spline
import csv

#Object that holds the phi data for a given simulation. Saves a HELL of a lot
#of time on importing multiple sets of data.
class Radiation:
    def __init__ (self, data):
        self.intensity = np.empty([2, data.shape[0]])
        self.real = np.empty(len(data), dtype=float)
        self.im = np.empty(len(data), dtype=float)

        for i in range(data.shape[0]):
	        self.real[i] = data[i,5] + data[i,7]
	        self.im[i] = data[i,6] + data[i,8]
	        self.intensity[0,i] = math.sqrt(self.real[i]**2 + self.im[i]**2)
        
#Load in your phi_XYZ.out file and declare it as a radiation object.
#You might need to declare a new phi directory.
data = np.loadtxt('phis/phi_XYZ.out')
ask_radiation = Radiation(data)

#Did your plots?
plt.plot(data[:,0], ask_radiation.intensity[0,:], label='Gremlins')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=3)
plt.xlabel("Viewing Angle " r"$\theta$ ($\degree$)")
plt.ylabel("Electric Field Magnitude (V/m)")	
plt.yscale('log')
plt.ylim(1e-20, 1e-15)
plt.title("I Really Need to Read Through Code Presets")
plt.savefig("Askaryan_Radiation_Test.png")
plt.show()
