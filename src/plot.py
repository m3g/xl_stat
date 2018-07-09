#!/usr/bin/python

import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

#matplotlib.rc('text', usetex=True)     
#protein_type = sys.argv[1]

data_file = open( 'xl_descriptor.tab' )

score_simxl, spec_count, species_count, consistency = \
   np.loadtxt(data_file,usecols=(0,1,2,3),unpack=True,comments="#",dtype=float)

# filter by specific score_simxl

score_choice = 4.0

print list(sorted(set(score_simxl)))

ndata = np.size(score_simxl) 
filter = np.zeros(ndata,dtype=bool)
i = 0
for value in score_simxl : 
  if score_simxl[i] == score_choice :
    filter[i] = True
  i = i + 1

#x = consistency[filter]
x = species_count[filter]
y = spec_count[filter]

plt.xlabel('Species count')
plt.ylabel('Spectrum count')

plt.xlim(-1,6)
plt.ylim(-1,6)

plt.plot(x,y,'o')
plt.show()


#plt.subplots_adjust(left=0.14, 
#                    bottom=0.10, 
#                    right=0.95, 
#                    top=0.90, 
#                    wspace=0.3, 
#                    hspace=0.3)
#plt.gcf().set_size_inches(6,6)
#
#plt.savefig('~/Dropbox/temp/plot.pdf')













