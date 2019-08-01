#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version 1.1 8/1/2019

Caleb M. Hill
Assistant Professor
Department of Chemistry
University of Wyoming
caleb.hill@uwyo.edu

This program fits a collection of spectra to Lorentzian functions. The input data is expected to exist as individual files within a directory. Each spectra is sequentially loaded and fit through a GUI. The resulting fits are output in a single Ns x 6 text file after all spectra have been fit.
-Ns: Number of spectra (files in directory)
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
from os import listdir
from os.path import isfile, join
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import scipy.optimize as opt

#Get all filepaths in a directory
filePath = tk.filedialog.askdirectory()
files = [f for f in listdir(filePath) if isfile(join(filePath,f))]
Nf = len(files)

#Create the output array
fitParamArray = sp.zeros((Nf,2))
n = 0

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

#Function for fitting
def fitData(data,x1,x2):
    x0g = 0.5*(x2-x1)
    wg = x0g
    Ag = sp.amax(data)*sp.pi*wg/2
    xfit = list(range(x1,x2,1))
    pg = sp.array([x0g,wg,Ag])
    popt, pcov = opt.curve_fit(lorentzian,xfit,data[x1:x2,0],p0=pg)
    perr = sp.sqrt(sp.diag(pcov))
    ax.plot(xfit,lorentzian(xfit,*popt))
    plt.pause(0.5)
    return popt

#Function for plotting
def plotData(n):
    data = pd.read_csv(join(filePath,files[n]),dtype=sp.float64,sep='\t',header=None).as_matrix()
    ax.cla()
    ax.plot(data)
    plt.pause(0.01)
    return data

#Stuff for the region selection
def line_select_callback(eclick, erelease):
    global data, n, fitParamArray
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    print(n,files[n])
    p = fitData(data,int(x1),int(x2))
    fitParamArray[n,:] = [n, p[0]]
    if n+1==Nf:
        plt.close()
        sp.savetxt(join(filePath,'fitParameters.txt'),fitParamArray,delimiter='\t')
    else:
        n += 1
        data = plotData(n)

def toggle_selector(event):
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)

#Create and initialize figure
fig, ax = plt.subplots(1,1)
toggle_selector.RS = RectangleSelector(ax, line_select_callback, drawtype='box', useblit=True, button=[1, 3], minspanx=5, minspany=5, spancoords='pixels', interactive=True)
plt.connect('key_press_event', toggle_selector)

data = plotData(0)
plt.show()