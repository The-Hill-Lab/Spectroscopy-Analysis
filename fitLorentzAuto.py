#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version 1.1 8/1/2019

Caleb M. Hill
Assistant Professor
Department of Chemistry
University of Wyoming
caleb.hill@uwyo.edu

This program fits a series of spectra to Lorentzian functions. The input data is expected to be a two-dimensional array with Ns rows and Nf columns. The output is a two-dimensional array with Ns rows and 6 columns.
-Ns: Number of spectra
-Nf: Number of points in each spectra
-Output Column 1: x0 / px (center frequency/wavelength)
-Output Column 2: w / px (FWHM)
-Output Column 3: A / px*a.u. (Area under curve)
-Output Column 4: Error in x0 / px
-Output Column 5: Error in w / px
-Output Column 6: Error in A / px*a.u.
In these expressions, "px" represents the wavelength/frequency units in pixels, and a.u. represents the arbitrary units of intensity.
"""

import tkinter as tk
from tkinter import filedialog
import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as opt

#Get file path via GUI
root = tk.Tk()
root.withdraw()
filePath = filedialog.askopenfilename()

#Import data
rawData = sp.genfromtxt(filePath,delimiter="\t") #Why?
totSpec = sp.sum(rawData,axis=1)/rawData.shape[1]
plt.plot(totSpec)
plt.ylabel('Intensity')
plt.xlabel('pixel')
plt.show()

#Trim fitting region
xi = int(input("Where should fitting start?"))
xf = int(input("Where should fitting end?"))
fitData = rawData[xi:xf,:]

#Create initial guesses
x0g = 0.5*(xf-xi)
wg = x0g
Ag = totSpec.max()*sp.pi*wg/2
xfit = list(range(xf-xi))

#Define fitting functions
def lorentzian(x,x0,w,A):
	return A/sp.pi*(0.5*w)/((x-x0)**2+(0.5*w)**2)

def lorFit(data):
	pg = sp.array([x0g,wg,Ag])
	popt, pcov = opt.curve_fit(lorentzian,xfit,data,p0=pg)
	x0 = popt[0] + xi + 1
	w = popt[1]
	A = popt[2]
	perr = sp.sqrt(sp.diag(pcov))
	x0err = perr[0]
	werr = perr[1]
	Aerr = perr[2]
	return [x0,w,A,x0err,werr,Aerr]

#Main loop
N = rawData.shape[1]
fitOut = sp.zeros((N,6))
for n in range(N):
	fitOut[n,:] = lorFit(fitData[:,n])

#Save file	
path = filePath[0:-4] + '_lorFits.txt'
sp.savetxt(path,fitOut,delimiter='\t')