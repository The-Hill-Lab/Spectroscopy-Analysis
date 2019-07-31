#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version 1.3 4/12/2019

Caleb M. Hill
Assistant Professor
Department of Chemistry
University of Wyoming
caleb.hill@uwyo.edu

This program is designed to extract a series of spectra from a data file consisting of a series of CCD images. The expected format for the input data is a two-dimensional array, with Nt*Ny points in the first dimension and Nx + 1 points in the second. 
Nt: # of successive CCD images acquired (# of time points)
Ny: # of spatial pixels in the y-direction
Nx: # of spatial pixels in the x-direction

Spectra are extracted from desired regions of the CCD through a GUI employing a rectangular selection (note that only the vertical boundaries of the selection affect the output). A simple background subtraction is carried out based on the regions immediately above and below the region of interest. The size of this region is controlled using the "bgWindow" variable below.
"""

import tkinter as tk
from tkinter import filedialog
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import pandas as pd
import scipy.optimize as opt
		
#Get file path via GUI
root = tk.Tk()
root.withdraw()
filePath = filedialog.askopenfilename()

#Import data
Ny = int(input("What is the image height (in px)?"))
rawData = pd.read_csv(filePath,dtype=sp.int32).as_matrix()
bgWindow = 5
Nt = int(rawData.shape[0] / Ny)
Nx = rawData.shape[1] - 1

dataCube = sp.zeros((Ny,Nx,Nt),dtype=sp.int32)
for nx in range(Nx):
	for ny in range(Ny):
		for nz in range(Nt):
			dataCube[ny,nx,nz] = rawData[nz*Ny+ny,nx+1]
	
totInt = sp.sum(dataCube,axis=2)

#File to output spectra of a single region
ns = 1
def exportSpectra(x1,x2,y1,y2):
	global ns
	Pspec = y2-y1+1
	Pbg = 2*bgWindow
	subCube = dataCube[y1:(y2+1),:,:]
	bgCube = dataCube[(y1-bgWindow):(y2+bgWindow+1),:,:]
	spectra = (1+Pspec/Pbg)*sp.sum(subCube,axis=(0)) - Pspec/Pbg*sp.sum(bgCube,axis=(0))
	path = filePath + '_' + str(ns) + '_x' + str(int((x1+x2)/2)) + '_y' + str(int((y1+y2)/2)) + '.txt'
	sp.savetxt(path,spectra,delimiter='\t')
	ns += 1

#Stuff for the region selection - carried out by fitting the y-profile to a gaussian
def gaussian(y,y0,w,A,bg):
	return bg + A*sp.exp(-(y-y0)**2/(2*w**2))

def gaussFit(data,y1,y2):
	y0g = 0.5*(y2+y1)
	wg = 0.25*(y2-y1)
	Ag = max(data) - min(data)
	bgg = min(data)
	pg = sp.array([y0g,wg,Ag,bgg])
	yfit = list(range(y1,y2+1))
	popt, pcov = opt.curve_fit(gaussian,yfit,data,p0=pg)
	y1d = int(round(popt[0] - popt[1]))
	y2d = int(round(popt[0] + popt[1]))
	return y1d, y2d

#This is the GUI function for region selection
def line_select_callback(eclick, erelease):
	'eclick and erelease are the press and release events'
	x1, y1 = eclick.xdata, eclick.ydata
	x2, y2 = erelease.xdata, erelease.ydata
	x1 = int(round(x1))
	x2 = int(round(x2))
	y1 = int(round(y1))
	y2 = int(round(y2))
	yProf = sp.sum(totInt[y1:(y2+1),:],axis=(1))
	y1d, y2d = gaussFit(yProf,y1,y2)
	exportSpectra(x1,x2,y1d,y2d)
	print("Selected Region: (%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
	print("Fit Region: (%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1d, x2, y2d))

#This can be used to deactivate the region selection, if desired. Can be helpful in resizing the window and such.
def toggle_selector(event):
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)

#Creates the figure and binds the GUI functionality
fig, current_ax = plt.subplots()  
totIntPlot = plt.imshow(totInt)
toggle_selector.RS = RectangleSelector(current_ax, line_select_callback, drawtype='box', useblit=True, button=[1, 3], minspanx=5, minspany=5, spancoords='pixels', interactive=True)
plt.connect('key_press_event', toggle_selector)
plt.show()