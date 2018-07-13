#
# Compute scores from data read
# Functions to parse sim-xl xml outputs and topolink log files
#
# L. Martinez, Institute of Chemistry - University of Campinas
# Jul 10, 2018 
# http://m3g.iqm.unicamp.br/topolink
#

#
# Function that starts everything
#

from sys import exit

def read_all(xml_file=None,\
             topolink_log=None,\
             topolink_input=None,\
             xic_file_name=None,\
             domain=None) :

  from sys import exit

  error = False
  if xml_file == None : 
    print "ERROR: You need to provide the xml_file."
    error = True
  if topolink_log == None : 
    print "ERROR: You need to provide the topolink log file."
    error = True
  if topolink_input == None : 
    print "ERROR: You need to provide the topolink input file containing linktype information."
    error = True

  if error : exit()

  if domain == None : domain = [0,100000]

  # Read xml file from SIM-XL

  nlinks, links = readxml(xml_file,xic_file_name,domain)

  # Read topolink input to get the length of the linkers

  for link in links :
    link.dmax = getdmax(topolink_input,link)

  # Read topolink log file to get the euclidean and topological distances

  for link in links :
    link.deuc, link.dtop = readlog(topolink_log,link)

  # Set consistency with 0. tolerance:

  for link in links :
    link.consistency = setconsistency(link,tol=0.)

  return nlinks, links


#
# Function and classes for reading links
#

class Link :
  def __init__(self,name,nscans) :
    self.name = name
    self.nscans = int(nscans)
    self.consistency = False
    self.deuc = 0.
    self.dtop = 0.

  def init_scans(self) :
    if self.nscans == 0 : 
      print ' ERROR: nscans = 0 ', self.name
    self.score1 = [ 0. for i in range(0,self.nscans) ]
    self.score2 = [ 0. for i in range(0,self.nscans) ]
    self.mplush = [ 0. for i in range(0,self.nscans) ]
    self.iscan = [ 0 for i in range(0,self.nscans) ]
    self.xic = [ 0. for i in range(0,self.nscans) ]
    self.logxic = [ -100. for i in range(0,self.nscans) ]

  def init_index(self) :
    self.str = self.name.split('-')
    self.residue1_index = int(self.str[0][1:])
    self.residue2_index = int(self.str[1][1:])

  def mplush_repeat(self,mplush,unique) :
    for data in unique :
      if abs( mplush - data ) < 1.0 : 
        return True
    return False

  def set_scores(self) :

    from numpy import log10

    unique = [self.mplush[0]]
    for i in range(1,self.nscans) :
      if not self.mplush_repeat(self.mplush[i],unique) :
        unique.append(self.mplush[i])
    self.nspecies = len(unique)

    self.avgscore1 = sum(self.score1)/len(self.score1)
    self.avgscore2 = sum(self.score2)/len(self.score2)
    self.maxscore1 = max(self.score1)
    self.maxscore2 = max(self.score2)
    self.sumscore1 = sum(self.score1)
    self.sumscore2 = sum(self.score2)

    # The xic score considers only scans for which proper xic data 
    # was reported (that is, that was not set to -1.00)

    self.hasxic = False

    self.sumxic = 0.
    self.sumlogxic = 0.
    nxic = 0
    for xic in self.xic :
      if xic >= 0. :  
        self.hasxic = True
        nxic = nxic + 1
        self.sumxic = self.sumxic + xic
        self.sumlogxic = self.sumlogxic + log10(xic)

    if self.hasxic : 
      self.avgxic = self.sumxic / float(nxic)
      self.avglogxic = self.sumlogxic / float(nxic)
      self.maxxic = max(self.xic)
      self.maxlogxic = log10(self.maxxic)

  def getscore(self,score) :
    if score == 'Consistency' : return self.consistency
    if score == 'Average Score1' : return self.avgscore1
    if score == 'Average Score2' : return self.avgscore2 
    if score == 'Maximum Score1' : return self.maxscore1 
    if score == 'Maximum Score2' : return self.maxscore2 
    if score == 'Sum of Score1' : return self.sumscore1 
    if score == 'Sum of Score2' : return self.sumscore2 
    if score == 'Number of Scans' : return self.nscans 
    if score == 'Number of Species' : return self.nspecies 
    if score == 'Average XIC' : return self.avgxic
    if score == 'Sum of XIC' : return self.sumxic
    if score == 'Maximum XIC' : return self.sumxic
    if score == 'Average log(XIC)' : return self.avglogxic
    if score == 'Sum of log(XIC)' : return self.sumlogxic
    if score == 'Maximum log(XIC)' : return self.maxlogxic
    return None


#
# Function that compares two links
#

def is_link(name1,name2) :
  r1 = name1.split('-')
  r2 = name2.split('-')
  if ( r1[0] == r2[0] and r1[1] == r2[1] ) or \
     ( r1[1] == r2[0] and r1[0] == r2[1] ) :
    return True
  return False

# Function that checks if the link belongs to the chosen domain

def in_domain(name,domain) :
  r = name.split('-')
  # Remove wrong assignments for which the two residues are the same
  if r[0] == r[1] : return False
  # Check if pair is in domain
  if ( int(r[0][1:]) >= int(domain[0]) and \
       int(r[0][1:]) <= int(domain[1]) ) and \
     ( int(r[1][1:]) >= int(domain[0]) and \
       int(r[1][1:]) <= int(domain[1]) ) :
    return True
  return False

# Function that reads the XIC file

def readxic(xic_file_name,link) :

  from numpy import zeros

  xic_file = open(xic_file_name)
  xic = zeros(link.nscans,dtype=float)
  xic = xic - 1. # xic=-1. means that no xic was reported

  iread = -1
  for line in xic_file :
    data = line.split()
    if iread != -1 :
      xic_read = float(data[len(data)-1])
      if xic_read > 0. : xic[iread] = xic_read
      iread = -1
      continue
    for i in range(0,link.nscans) :
      try : 
        test = int(data[0]) 
        if data[0] == link.iscan[i] :
          iread = i
          break
      except :
        break

  return xic

# Function that reads link log and sets link data

def readlog(logfile_name,link) :

  logfile = open(logfile_name)
  for line in logfile :
    if line[2:7] == "LINK:" : 
      residue1_name = oneletter(line[8:12].strip())
      residue1_index = line[15:19].strip()
      residue2_name = oneletter(line[25:29].strip())
      residue2_index = line[32:36].strip()
      name = residue1_name+residue1_index+'-'+residue2_name+residue2_index
      deuc = float(line[43:51])
      dtop = float(line[52:61].replace(">",""))
      if is_link(name,link.name) :
        logfile.close()
        return deuc, dtop
  logfile.close()
  print " Warning: Could not find LINK data in TopoLink log for link: "+link.name
  return -1., -1.

# Function that reads the linker length from a topolink input file with linktype
# definitions

def getdmax(linktype_file,link) :

  linktype_file = open(linktype_file)
  for line in linktype_file :
    if not comment(line) :
      if "linktype" in line : 
        data = line.split()
        residue1_name = oneletter(data[1])
        residue2_name = oneletter(data[5])
        dmax = float(data[9])
        name = link.name.split('-')
        if ( name[0][0] == residue1_name and name[1][0] == residue2_name ) or \
           ( name[0][0] == residue2_name and name[1][0] == residue1_name ) : 
          linktype_file.close()
          return dmax
  linktype_file.close()
  return -1.

#
# Print link in the topolink format
#

def write(link,tol=None) :

  if tol == None : tol = 0.

  data = link.name.split("-")
  res1_name = threeletter(data[0][0])
  res2_name = threeletter(data[1][0])
  res1_index = int(data[0][1:])
  res2_index = int(data[1][1:])
  consistency = setconsistency(link,tol=tol)
  print '{} {:3} {} {:3} {} {:3.2f} {:3.2f}'.format(res1_name,res1_index,res2_name,res2_index,consistency,link.dtop,link.dmax)

#
# Set consistency, given a tolerance
#

def setconsistency(link,tol=None) :

  if tol == None : tol = 0.
  if link.dtop >= 0. and link.dtop <= link.dmax + tol :
    consistency = True
  else :
    consistency = False

  return consistency

# 
# Compute point-biserial correlation (x assumes 0 or 1 values, y is continuous)
# x is boolean (True or False for each group)
#

def point_biserial(x,y) :

  import numpy as np

  ndata = np.size(x)
  g1 = np.zeros(ndata,dtype=bool)
  g2 = np.zeros(ndata,dtype=bool)
  i = 0
  for value in x :
    if x[i] :
      g1[i] = True
    else :
      g2[i] = True
    i = i + 1

  data1 = y[g1]
  data2 = y[g2]

  n1 = len(data1)
  n2 = len(data2)
  if n1 == 0 or n2 == 0 : return 0.

  avg1 = float(sum(data1))/n1
  avg2 = float(sum(data2))/n2
  sd = np.std(y)

  pbs = ( ( avg1 - avg2 ) / sd ) * np.sqrt( float(n1*n2) / (ndata**2) )

  return pbs

#
# Function that reads the xml file generated by sim-xl
#

def readxml(data_file_name,xic_file_name,domain) :

  # Count the number of links listed

  data_file = open(data_file_name)
  nlinks = 0
  for line in data_file :
    if not comment(line) : 
      if newlink(line) != None :
        name = newlink(line) 
        if in_domain(name,domain) : nlinks = nlinks + 1
  data_file.seek(0)

  # Create list of links

  links = [ Link("None",0) for i in range(nlinks) ]  

  # Count the number of scans per link

  ilink = -1
  for line in data_file :
    if not comment(line) : 
      if newlink(line) != None :
        name = newlink(line) 
        if in_domain(name,domain) :
          ilink = ilink + 1
          links[ilink].name = name
      if "Scan" in line : 
        if in_domain(name,domain) : 
          links[ilink].nscans = links[ilink].nscans + 1
  data_file.seek(0)

  # Initialize the scans per link

  for link in links :
    link.init_scans()
    link.init_index()

  # Read the score parameters for each scan

  ilink = -1
  for line in data_file :
    
    if not comment(line) : 
      if newlink(line) != None :
        name = newlink(line) 
        if in_domain(name,domain) :
          ilink = ilink + 1
          iscan = -1
  
      if "Scan" in line : 
        if in_domain(name,domain) : 
          iscan = iscan + 1
    
          line = line.replace("Scan:"," ")
          line = line.replace("Secondary Score:"," ")
          line = line.replace("Score:"," ")
          line = line.replace("Experimental M+H:"," ")
          data = line.split()
    
          links[ilink].score1[iscan] =  float(data[1])
          links[ilink].score2[iscan] =  float(data[2])
          links[ilink].mplush[iscan] =  float(data[3])
          links[ilink].iscan[iscan] = data[0]
    
  # Read XIC data

  if xic_file_name != None :
    for link in links :
      link.xic = readxic(xic_file_name,link)

  # Compute scores from data read

  for link in links :
    link.set_scores()

  return nlinks, links

#
# Remove a specific link from the list
#

def remove(links,name) :
  i=-1
  for link in links :
    i=i+1
    if link.name == name : 
      del links[i]
      return links
  print ' Error in remove function : link ', link.name, ' not found in link list. ' 
  exit()

#
# Check if a line is commented
#

def comment(line) :
  if len(line.strip()) == 0 : return True
  if line.strip()[0] == "#" :
    return True
  else :
    return False

#
# Check if a line in a xml file is a new link line, and return
# its name if so
#
    
def newlink(line) : 

  residues = { 'C', 'D', 'S', 'Q', 'K',
               'I', 'P', 'T', 'F', 'N', 
               'G', 'H', 'L', 'R', 'W', 
               'A', 'V', 'E', 'Y', 'M' }

  if "(" in line : 
    line = line.replace(" (","",2)
    line = line.replace(") - ","-")
    line = line.replace(")","")
    name = line.strip()
    data = name.split('-')
    if len(data) == 2 : 
      if ( data[0][0] in residues ) and ( data[1][0] in residues ) :
        return name

  return None

#
# Function that returns data to be ploted
#

def setplot(links,x,y,tol=None) :

  from sys import exit
  import numpy as np

  if tol == None : tol = 0.

  if links[0].getscore(x) == None : 
    print ' ERROR: x must be one of the available indicators. '
    exit()

  if links[0].getscore(y) == None : 
    print ' ERROR: y must be one of the available indicators. '
    exit()

  # Set xplot and yplot vectors that will be returned

  nlinks = len(links) 
  xplot = np.zeros(nlinks)
  yplot = np.zeros(nlinks)

  i=-1
  for link in links : 
    i=i+1

    if x == 'Consistency' : 
      xplot[i] = setconsistency(link,tol) # Check consistency according to given tol
    else : 
      xplot[i] = link.getscore(x)

    yplot[i] = link.getscore(y)

  return xplot, yplot
  
#
# Function that sets the data for plotting the correlation of a score
# as a function of the tolerance
#
# tol has three elements: minimum, maximum and step: ex: [-3., 15., 0.5]
# and corresponds to the deviations from the dmax set for each link
#

def pbs_vs_tol(links,score,tol=None) :

  from sys import exit
  import numpy as np

  if tol == None :
    tol = [ -3., 15., 0.5 ]

  if links[0].getscore(score) == None :
    print ' ERROR: score must be one of the available indicators. '
    exit()

  ntol = int((tol[1] - tol[0])/tol[2])
  x = np.zeros(ntol)
  i=-1
  for xtol in x :
    i=i+1
    x[i] = tol[0] + tol[2]*i
  pbs = np.zeros(ntol)
  
  nlinks = len(links)
  y = np.zeros(nlinks,dtype=float)
  consistency = np.zeros(nlinks,dtype=bool)

  i=-1
  for xtol in x :
    i=i+1

    j=-1
    for link in links :
      j=j+1
      consistency[j] = setconsistency(link,xtol)
      y[j] = link.getscore(score)

    # Compute pbs
    pbs[i] = point_biserial(consistency,y)

  return x, pbs

#
# Filter links according to spectral parameters
#

def filter(links,scores,filter_type) :

  filtered = []

  if filter_type == 'or' :
    for link in links :
      select = False
      iscore=-1
      for score in scores :
        iscore=iscore+1
        if link.getscore(scores[iscore][0]) >= scores[iscore][1] :
          select = True
      if select : filtered.append(link)

  if filter_type == 'and' :
    for link in links :
      select = True
      iscore=-1
      for score in scores :
        iscore=iscore+1
        if link.getscore(scores[iscore][0]) < scores[iscore][1] :
          select = False
      if select : filtered.append(link)

  return filtered

#
# Function to map all possible filter options
#

def search_filters(links,scores=None,filter_type=None,nfilter=None,tol=None) :
  
  import numpy as np
  from sys import exit

  if scores == None : 
    print 'ERROR: The scores search list must be set. '
    exit()

  for score in scores :
    if links[0].getscore(score[0]) == None :
      print 'ERROR: Some score incorrectly set in score list. '
      exit()

  if filter_type == None :  
    print "ERROR: Please set filter_type to 'or' or 'and'"
    exit()
  
  if nfilter == None :  
    nfilter = [0,-1] 

  if tol == None :  
    tol = 0.

  step = [ np.arange(scores[i][1],scores[i][2]+scores[i][3],scores[i][3]) \
          for i in range(0,len(scores)) ]

  print '# N(selected)  N(consistent) f(correct) '+\
        '{} {} {} {}'.format(scores[0][0],scores[1][0],scores[2][0],scores[3][0])

  for s0 in step[0] : 
    for s1 in step[1] :
      for s2 in step[2] :
        for s3 in step[3] :

          nselected = 0
          nc = 0
          for link in links :

            if filter_type == 'or' :
              if link.getscore(scores[0][0]) >= s0 or \
                 link.getscore(scores[1][0]) >= s1 or \
                 link.getscore(scores[2][0]) >= s2 or \
                 link.getscore(scores[3][0]) >= s3 :
                nselected = nselected + 1
                if setconsistency(link,tol=tol) :
                  nc = nc + 1
            if filter_type == 'and' :
              if link.getscore(scores[0][0]) >= s0 and \
                 link.getscore(scores[1][0]) >= s1 and \
                 link.getscore(scores[2][0]) >= s2 and \
                 link.getscore(scores[3][0]) >= s3 :
                nselected = nselected + 1
                if setconsistency(link,tol=tol) :
                  nc = nc + 1
  
          if nselected >= nfilter[0] and \
             ( nselected <= nfilter[1] or nfilter[1] < 0 ) :
            print '{:^13}{:^16}{:^10.2f}{:^16.2f}{:^16.2f}{:^15}{:^20}'\
                  .format(nselected,nc,float(nc)/nselected, s0, s1, s2, s3)

#
# Convert one letter and three-letter amino acid codes
#

def oneletter(x):
  d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
       'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
       'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
       'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
  return d[x]

def threeletter(x):
  d = { 'C':'CYS', 'D':'ASP', 'S':'SER', 'Q':'GLN', 'K':'LYS',
        'I':'ILE', 'P':'PRO', 'T':'THR', 'F':'PHE', 'N':'ASN', 
        'G':'GLY', 'H':'HIS', 'L':'LEU', 'R':'ARG', 'W':'TRP', 
        'A':'ALA', 'V':'VAL', 'E':'GLU', 'Y':'TYR', 'M':'MET' }
  return d[x]



  
