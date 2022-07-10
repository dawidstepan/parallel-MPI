#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Dawid Stepanovic
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_1():
	
	data = np.genfromtxt('~/MPI/log.txt', delimiter=',', dtype=None, encoding='UTF-8') # change directory
	data = pd.DataFrame(data)
		
	plt.plot(data.iloc[:,1],data.iloc[0,2]/data.iloc[:,2],label='speedup', color='r', marker='.',lw=0.5)	    
	plt.xlabel('p - number of processes')
	plt.ylabel(r'speedup, $S_{p} = \frac{T_{1}}{T_{p}}$')
	plt.minorticks_on()
	plt.grid()
	plt.legend()
	plt.savefig('speedup.png')
	plt.show() 

def plot_2():
	
	data = np.genfromtxt('~/MPI/log.txt', delimiter=',', dtype=None, encoding='UTF-8') # change directory
	data = pd.DataFrame(data)
		
	plt.plot(data.iloc[:,1],data.iloc[:,2],label='runtime', color='b', marker='.',lw=0.5)	    
	plt.xlabel('p - number of processes')
	plt.ylabel('runtime in seconds [sec]')
	plt.minorticks_on()
	plt.grid()
	plt.legend()
	plt.savefig('runtime.png')
	plt.show() 
	
plot_1()
plot_2()