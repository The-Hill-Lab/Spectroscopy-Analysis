#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version 1.2 7/31/2019

Caleb M. Hill
Assistant Professor
Department of Chemistry
University of Wyoming

This program is designed to extract spectra from "hyperspectral imaging" (HSI) data files. The expected input files are 2-dimensional arrays, with Nx*Ny points in the first dimension and Nz points in the second. Nx and Ny correspond to the number of spatial pixels in the x and y directions, with Nz corresponds to the number of wavelength values (i.e., the number of pixels along the long dimension of the employed CCD).

Data is extracted in a "point-and-click" fashion through a GUI. Two files are generated for each point: (1) the raw data containing the integrated spectra and (2) a background corrected spectrum. The size of the region of interest is set through an input prompt. The background correction is a simple subtraction based on the average spectrum of the pixels just outside the region of interest.
"""

import tkinter as tk
from tkinter import filedialog
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
		
#Get file path via GUI
root = tk.Tk()
root.withdraw()
filePath = filedialog.askopenfilename()
	
#Import data
rawData = pd.read_csv(filePath,dtype=sp.float64,sep=',',header=None).as_matrix()
x_index = rawData[:,0]
Nx = int(max(x_index))
Ny = int(rawData.shape[0]/Nx)
Nz = rawData.shape[1] - 1

#Establish analysis window
w = int(input("What's the analysis window size? (e.g., x by x)"))
d = int((w-1)/2)

#Change image scaling (useful if there is a bright defect preventing visualization of the features of interest)
imgScale = float(input("How should the image max be scaled? (use 1 if not sure)"))

dataCube = sp.zeros((Ny,Nx,Nz))
for nx in range(Nx):
	for ny in range(Ny):
		for nz in range(Nz):
			dataCube[ny,nx,nz] = rawData[nx*Ny+ny,nz+1]
			
#Save integrated intensity image
path = filePath + '_img.txt'
image = sp.sum(dataCube, axis=(2))
sp.savetxt(path, image, delimiter='\t')			
			
#Function to output spectra of a single region
ns = 1
def exportSpectra(x,y):
	global ns
	x1 = max(x-d,0)
	x2 = min(x+d,Nx)
	y1 = max(y-d,0)
	y2 = min(y+d,Ny)
	Np = (x2-x1+1)*(y2-y1+1)
	xb1 = max(x1-1,0)
	xb2 = min(x2+1,Nx)
	yb1 = max(y1-1,0)
	yb2 = min(y2+1,Ny)
	Nb = (xb2-xb1+1)*(yb2-yb1+1)-Np
	subCube = dataCube[y1:(y2+1),x1:(x2+1),:]
	bgCube = dataCube[yb1:(yb2+1),xb1:(xb2+1),:]
	rawspectra = sp.sum(subCube,axis=(0,1))
	bgcorrspectra = rawspectra - (sp.sum(bgCube,axis=(0,1))-rawspectra)*(Np/Nb)
	path = filePath + '_' + str(ns) + '_x' + str(int((x1+x2)/2)) + '_y' + str(int((y1+y2)/2)) + '.txt'
	sp.savetxt(path,rawspectra,delimiter='\t')
	bgpath = filePath + '_' + str(ns) + '_x' + str(int((x1+x2)/2)) + '_y' + str(int((y1+y2)/2)) + '_bgCorr.txt'
	sp.savetxt(bgpath,bgcorrspectra,delimiter='\t')
	ns += 1

#GUI function for region selection
def onclick(event):
	x = int(round(event.xdata))
	y = int(round(event.ydata))
	exportSpectra(x,y)
	print(x,y)

#Create the figure and bind the GUI function
totInt = sp.sum(dataCube,axis=2)
totIntMax = sp.amax(totInt)
fig, current_ax = plt.subplots()  
totIntPlot = plt.imshow(totInt,interpolation='none',vmax=imgScale*totIntMax)
plt.connect('button_press_event', onclick)
plt.show()
