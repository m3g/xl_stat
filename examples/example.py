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

protein_name='SalbIII'
xml_file = 'salbiii_hitsDetail.dat'
topolink_log = 'salbiii_topolink.log'
domain = [2,134]
xic_file_name='salbiii_xic.dat'

#protein_name='ALB'
#xml_file = 'alb_hitsDetail.dat'
#topolink_log = 'alb_topolink.log'
#domain = [1,201] ; protein_name=protein_name+'-D1'
#domain = [202,390] ; protein_name=protein_name+'-D2' 
#domain = [391,584] ; protein_name=protein_name+'-D3' 
#xic_file_name=None

# topolink input file containing the linktype definitions:

topolink_input = "topolink.inp"

#
# Read all files to get link data
#

nlinks, links = xl_stat.read_all(xml_file=xml_file,\
                                 topolink_log=topolink_log,\
                                 topolink_input=topolink_input,\
                                 xic_file_name=xic_file_name,\
                                 domain=domain)

# indicators available: 

scores = [ 'Average Score1', 'Average Score2', 
           'Number of Scans', 'Number of Species', \
           'Sum of log(XIC)', 'Maximum log(XIC)' ]

#
# Remove one of the links from the list, if you want
#
#links = xl_stat.remove(links,'S133-K99')

#
# Filter the links to get only those with xic data
#

tol=5.
nc = 0
for link in links :
  if xl_stat.setconsistency(link,tol=tol) : nc = nc + 1
print ' Total number of links: ', len(links), ' nc = ', nc

nlinks = len(links)-1
for ilink in range(nlinks,0,-1) :
  if not links[ilink].hasxic :
    links = xl_stat.remove(links,links[ilink].name)

nc = 0
for link in links :
  if xl_stat.setconsistency(link,tol=tol) : nc = nc + 1
print ' Number of links with XIC data: ', len(links), ' nc = ', nc

#
# Plot one score as a function of the other, for a given tolerance relative to dmax
#
#x, y = xl_stat.setplot(links,x='Average Score1',y='Sum of log(XIC)',tol=5.)
#plt.plot(x,y,'o')
##plt.xlim(-0.5,1.5) # Uncomment if x is consistency for a nice plot
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
plt.savefig(plot_dir+'/'+protein_name+'_pbs_correlations.png')
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
plt.savefig(plot_dir+'/'+protein_name+'_scores_consistency.png')
plt.close()

Write = False
if Write : 
  #
  # Output the result of using a set of parameters to filter links
  #
  
  print '# protein: ',protein_name, domain
  
             #  Name                 value
  scores = [ [ 'Average Score1'     , 2.00 ] ,\
             [ 'Average Score2'     , 2.00 ] ,\
             [ 'Number of Species'  ,   4  ] ,\
             [ 'Number of Scans'    ,  13  ] ]
  
  filter_type='and'
  #filter_type='or'
  
  filtered_links = xl_stat.filter(links,scores,filter_type)

  # Print links selected
  tol=5.0
  nc=0
  for link in filtered_links : 
    xl_stat.write(link,tol=5.0)
    if xl_stat.setconsistency(link,tol=5.0) :
      nc=nc+1
  print '# Total: ',len(filtered_links), ' NC = ',nc

Search = False
if Search :

  #
  # Search best set of scores
  #
  
  # the scores list contains the range of scores to be search, with steps
               # Name                 min    max  step
  scores = [ [ 'Average Score1'   ,   0.,  4.00,  0.5 ] , \
             [ 'Average Score2'   ,   0.,    1.,  0.5 ] , \
             [ 'Number of Species',    0,     1,    1 ] , \
             [ 'Number of Scans'  ,    0,     1,    1 ] , \
             [ 'Sum of log(XIC)'  ,   0.,    60,   5. ] ]
  
  # Range of filtered selection sizes to be considered (here about 20% of domains size)
  
  n = domain[1] - domain[0]         # domain size
  n_min = int((n-0.05*n)*(2./10.))  # 5% less than 20% of domain size
  n_max = int((n+0.05*n)*(2./10.))  # 5% more than 20% of domain size
  nfilter = [ n_min, n_max ]
  tol=5.
  filter_type='and'
   
  xl_stat.search_filters(links,scores=scores,filter_type=filter_type,nfilter=nfilter,tol=tol) 





