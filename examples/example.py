#!/usr/bin/python

import sys
# include path for xl_stat.py file
sys.path.append('../src/')

import xl_stat
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import histogram

#
#  Set system to be analyzed
#

plot_dir = "./"

#protein_name='SalbIII'
#xml_file = 'salbiii_hitsDetail.dat'
#topolink_log = 'salbiii_topolink.log'
#domain = [2,134]

protein_name='ALB'
xml_file = 'alb_hitsDetail.dat'
topolink_log = 'alb_topolink.log'
domain = [1,201] ; protein_name=protein_name+'-D1'
#domain = [202,390] ; protein_name=protein_name+'-D2' 
#domain = [391,584] ; protein_name=protein_name+'-D3' 

# topolink input file containing the linktype definitions:

topolink_input = "topolink.inp"

#
# Read all files to get link data
#

nlinks, links = xl_stat.read_all(xml_file=xml_file,\
                                 topolink_log=topolink_log,\
                                 topolink_input=topolink_input,\
                                 domain=domain)

# indicators available: 

scores = [ 'Average Score1', 'Average Score2', 
           'Number of Scans', 'Number of Species', \
           'Maximum Score1', 'Maximum Score2' ]

#
# Remove one of the links from the list, if you want
#
#links = xl_stat.remove(links,'S133-K99')

#
# Plot one score as a function of the other, for a given tolerance relative to dmax
#
#x, y = xl_stat.setplot(links,x='Consistency',y='Number of Scans',tol=5.)
#plt.plot(x,y,'o')
#plt.xlim(-0.5,1.5) # Uncomment if x is consistency for a nice plot
#plt.show()
#sys.exit()

#
# Plot the point-biserial correlation as function of the tolerance, for one score
# tol=[-3.,15.,0.5] is the minimum, maximum and step of the tolerance relative do dmax.
#
#x, pbs = xl_stat.pbs_vs_tol(links,score='Consistency',tol=[-3., 20., 0.5])
#plt.plot(x,pbs)
#plt.show()
#sys.exit()

#
# For all scores, plot the point-biserial correlation as a function
# of the deviation fro dmax
#

iplot = 0
for score in scores :
  iplot=iplot+1
  x, pbs = xl_stat.pbs_vs_tol(links,score=score,tol=[-3., 20., 0.5]) 
  plt.subplot(3,2,iplot)
  plt.plot(x,pbs,color='black')
  plt.title(score,size=12)
  plt.xlabel('d$_{max}$ deviation',size=12)
  plt.ylabel('correlation',size=12)


plt.subplots_adjust(left=0.14, 
                    bottom=0.10, 
                    right=0.95, 
                    top=0.90, 
                    wspace=0.4, 
                    hspace=0.6)
plt.gcf().set_size_inches(6,8)
plt.savefig(plot_dir+'/'+protein_name+'_pbs_correlations.pdf')
plt.close()

#
# Plot the actual data from which the PBS correlations are computed
# (the scores as a function of the consistency with the model, for
#  a given tolerance)
#

tol = 5.
iplot = 0
for score in scores : 
  x, y = xl_stat.setplot(links,x='Consistency',y=score,tol=tol)
  pbs = xl_stat.point_biserial(x,y)
  iplot=iplot+1
  plt.subplot(3,2,iplot)
  plt.plot(x,y,'o',alpha=0.3,color='black')
  plt.xlabel('Consistency',size=12)
  plt.ylabel(score,size=12)
  plt.xlim(-0.5,1.5)
  plt.annotate('R='+'{:3.2f}'.format(pbs),xy=(0.02,0.9),size=12,xycoords='axes fraction')
 
plt.subplots_adjust(left=0.14, 
                    bottom=0.10, 
                    right=0.95, 
                    top=0.90, 
                    wspace=0.4, 
                    hspace=0.6)
plt.gcf().set_size_inches(6,8)
plt.savefig(plot_dir+'/'+protein_name+'_scores_consistency.pdf')
plt.close()


sys.exit()


#
# Testing 22   17 2.75 1.47 1.90 0.77
#2.25 0.09 1.80
#
#2.02 0.14 1.00 13.2

#voltar
test = True
if test :
  score = [ 4.40, 2.53, 9.10, 13.2 ]
  score = [ 3.01, 2.53, 2, 24 ]
  score = [ 3.7, 2.53, 4, 33 ]
              
  nget = 0
  nc = 0
  for link in links :
    #if link.avgscore1 >= score[0] or \
    #   link.avgscore2 >= score[1] or \
    #   link.nspecies >= score[2] or \
    #   link.nscans >= score[3] :
    if link.maxscore1 >= score[0] or \
       link.nspecies >= score[2] or \
       link.nscans >= score[3] :
      nget = nget + 1
      xl_stat.write(link)
      if link.consistency : nc = nc + 1
  print xml_file, domain
  print '{:4} {:4} {:3.2f} {:3.2f} {:3.2f} {:3} {:3}'.format(nget, nc, score[0], score[1], score[2], score[3], float(nc)/nget)

  sys.exit()

#
# Search best set of scores
#


else :
   
  # Number of residues in domain

  n = domain[1]-domain[0]

  # We want a set of about 2N/10 constraints, where N is the number of
  # residues of the protein
  
  n_min = int((n-0.05*n)*(2./10.))
  n_max = int((n+0.05*n)*(2./10.))
  print n_min, n_max
  
  # We will search for all possible combinations of three of the scores to
  # get the best set (the one with the greatest number of true positives)
  
  nsteps = 10
  minscores = np.array([ min(y[0]), min(y[1]), min(y[2]), min(y[5]) ])
  maxscores = np.array([ max(y[0]), max(y[1]), max(y[2]), max(y[5]) ])
  step = (maxscores - minscores)/nsteps
  
  print 'Nget  Nc  avsc1   avsc2  nspec  nscans nc/nget  totnc'
  ncmax = 0
  for i in range(0,nsteps) :
    #for j in range(0,nsteps) :
       for k in range(0,nsteps) :
         for l in range(0,nsteps) : 
  
           score0 = minscores[0] + i*step[0]
           #score1 = minscores[1] + j*step[1]
           score2 = minscores[2] + k*step[2]
           score3 = minscores[3] + l*step[3]
  
           nget = 0
           nc = 0
           for link in links : 
              #if link.avgscore1 >= score0 or \
              #   link.avgscore2 >= score1 or \
              #   link.nspecies >= score2 or \
              #   link.nscans >= score3 :
              if link.maxscore1 >= score0 or \
                 link.nspecies >= score2 or \
                 link.nscans >= score3 :
                nget = nget + 1
                if link.consistency : nc = nc + 1
           if nget >= n_min and nget <= n_max : 
             #print '{:4} {:3} {:4.2f} {:3.2f} {:3.2f} {:3} {:3.2f} {:3}'.format(nget, nc, score0, score1, score2, score3, float(nc)/nget, ncons)
             print '{:4} {:3} {:4.2f} {:3.2f} {:3} {:3.2f} {:3}'.format(nget, nc, score0, score2, score3, float(nc)/nget, ncons)


sys.exit()

