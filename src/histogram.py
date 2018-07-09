#
# histogram: Builds a histogram from a float xdata vector
#
# Usage: 
#       import histogram
#       x, y = histogram.histogram(xdata)
#
# Options: min [minimum value to be considered] 
#          max [maximum value to be considered] 
#          step [step of summation]
#          nsteps [number of steps]
# 
#          Normalization (default: int=0):
#          int=0 : normalize to unity the integral of the histogram. 
#          int=1 : normalize the maximum probability to unity
#
# All parameters except the input data (xdata) are optional. Default values are:
#
# minimum value: The mininum value from data.
# maximum value: The maximum value from data.
# step: (max - min) / nsteps
# number of steps: 100
#
# L. Martinez, IFSC-USP, Sep 6, 2011.
# Translated to python: IQ/UNICAMP on July 8, 2018
#

import numpy as np

def histogram(x,\
              min=None,max=None,nsteps=None,step=None,int=None) :

  if min == None : min = np.min(x)
  if max == None : max = np.max(x)
  if nsteps == None : nsteps = 100
  if step == None : step = ( max - min ) / nsteps
  if int == None : int = 0
  dstep = ( max - min ) / nsteps

  xhist = np.zeros(nsteps+1,dtype=float)
  yhist = np.zeros(nsteps+1,dtype=float)

  # Assigning the value of the x coordinate of each histogram step 

  for i in range(0,nsteps+1) :
    xhist[i] = min + (i-1)*dstep

  # Computing the histogram

  for val in x : 
    for i in range(0,nsteps+1) :
      if ( val >= ( xhist[i] - step ) and \
           val <= ( xhist[i] + step ) ) :
        yhist[i] = yhist[i] + 1

  # Normalizing the histogram 
  
  ndata = len(x)
  integral = 0.
  pmax = 0.
  for i in range(0,nsteps+1) :
    yhist[i] = yhist[i] / ndata
    if int == 0 : integral = integral + yhist[i]*dstep
    if int == 1 : pmax = np.max([pmax,yhist[i]])

  if int == 0 :
    for i in range(0,nsteps+1) :
      yhist[i] = yhist[i] / integral

  if int == 1 :
    for i in range(0,nsteps+1) :
      yhist[i] = yhist[i] / pmax

  return xhist, yhist














