#
# Functions to parse sim-xl xml outputs and topolink log files
#
# L. Martinez, Institute of Chemistry - University of Campinas
# Jul 10, 2018 
# http://m3g.iqm.unicamp.br/topolink
#

#
# Function that starts everything
#

def read_all(xml_file=None,\
             topolink_log=None,\
             topolink_input=None,\
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

  nlinks, links = readxml(xml_file,domain)

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

  def init_index(self) :
    self.str = self.name.split('-')
    self.residue1_index = int(self.str[0][1:])
    self.residue2_index = int(self.str[1][1:])

  def set_scores(self) :
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

  def mplush_repeat(self,mplush,unique) :
    for data in unique :
      if abs( mplush - data ) < 1.0 : 
        return True
    return False

# Function that compares two links

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

def write(link) :
  data = link.name.split("-")
  res1_name = threeletter(data[0][0])
  res2_name = threeletter(data[1][0])
  res1_index = int(data[0][1:])
  res2_index = int(data[1][1:])
  print '{} {:3} {} {:3} {} {:3.2f} {:3.2f}'.format(res1_name,res1_index,res2_name,res2_index,link.consistency,link.dtop,link.dmax)

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

def readxml(data_file_name,domain) :

  # Count the number of links listed

  data_file = open(data_file_name)
  nlinks = 0
  for line in data_file :
    if not comment(line) : 
      if "(" in line : 
        line = line.replace(" (","",2)
        line = line.replace(") - ","-")
        line = line.replace(")","")
        name = line.strip()
        if in_domain(name,domain) : nlinks = nlinks + 1
  data_file.seek(0)

  # Create list of links

  links = [ Link("None",0) for i in range(nlinks) ]  

  # Count the number of scans per link

  ilink = -1
  for line in data_file :
    if not comment(line) : 
      if "(" in line : 
        line = line.replace(" (","",2)
        line = line.replace(") - ","-")
        line = line.replace(")","")
        name = line.strip()
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
      if "(" in line : 
        line = line.replace(" (","",2)
        line = line.replace(") - ","-")
        line = line.replace(")","")
        name = line.strip()
        if in_domain(name,domain) :
          ilink = ilink + 1
          iscan = 0
  
      if "Scan" in line : 
        if in_domain(name,domain) : 
    
          line = line.replace("Scan:"," ")
          line = line.replace("Secondary Score:"," ")
          line = line.replace("Score:"," ")
          line = line.replace("Experimental M+H:"," ")
          data = line.split()
    
          links[ilink].score1[iscan] =  float(data[1])
          links[ilink].score2[iscan] =  float(data[2])
          links[ilink].mplush[iscan] =  float(data[3])
    
          iscan = iscan + 1
  
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
# Function that returns data to be ploted
#

def setplot(links,x,y,tol=None) :

  import numpy as np

  if tol == None : tol = 0.

  data_types = ['Consistency',\
                'Average Score1',\
                'Average Score2',\
                'Maximum Score1',\
                'Maximum Score2',\
                'Sum of Score1',\
                'Sum of Score2',\
                'Number of Scans',\
                'Number of Species' ]

  if not x in data_types : 
    print ' ERROR: x must be one of the following data types: '
    for type in data_types :
      print type

  if not y in data_types : 
    print ' ERROR: x must be one of the following data types: '
    for type in data_types :
      print type

  # Set xplot and yplot vectors that will be returned

  nlinks = len(links) 
  xplot = np.zeros(nlinks)
  yplot = np.zeros(nlinks)

  i=-1
  for link in links : 
    i=i+1

    # Check consistency according to given tol
    consistency = setconsistency(link,tol)

    if x == 'Consistency' : xplot[i] = consistency
    if x == 'Average Score1' : xplot[i] = link.avgscore1
    if x == 'Average Score2' : xplot[i] = link.avgscore2 
    if x == 'Maximum Score1' : xplot[i] = link.maxscore1 
    if x == 'Maximum Score2' : xplot[i] = link.maxscore2 
    if x == 'Sum of Score1' : xplot[i] = link.sumscore1 
    if x == 'Sum of Score2' : xplot[i] = link.sumscore2 
    if x == 'Number of Scans' : xplot[i] = link.nscans 
    if x == 'Number of Species' : xplot[i] = link.nspecies 
    
    if y == 'Consistency' : yplot[i] = link.consistency
    if y == 'Average Score1' : yplot[i] = link.avgscore1
    if y == 'Average Score2' : yplot[i] = link.avgscore2 
    if y == 'Maximum Score1' : yplot[i] = link.maxscore1 
    if y == 'Maximum Score2' : yplot[i] = link.maxscore2 
    if y == 'Sum of Score1' : yplot[i] = link.sumscore1 
    if y == 'Sum of Score2' : yplot[i] = link.sumscore2 
    if y == 'Number of Scans' : yplot[i] = link.nscans 
    if y == 'Number of Species' : yplot[i] = link.nspecies 

  return xplot, yplot
  
#
# Function that sets the data for plotting the correlation of a score
# as a function of the tolerance
#
# tol has three elements: minimum, maximum and step: ex: [-3., 15., 0.5]
# and corresponds to the deviations from the dmax set for each link
#

def pbs_vs_tol(links,score,tol=None) :

  import numpy as np

  if tol == None :
    tol = [ -3., 15., 0.5 ]

  data_types = ['Consistency',\
                'Average Score1',\
                'Average Score2',\
                'Maximum Score1',\
                'Maximum Score2',\
                'Sum of Score1',\
                'Sum of Score2',\
                'Number of Scans',\
                'Number of Species' ]

  if not score in data_types : 
    print ' ERROR: score must be one of the following data types: '
    for type in data_types :
      print type

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

      if score == 'Consistency' : y[j] = link.consistency
      if score == 'Average Score1' : y[j] = link.avgscore1
      if score == 'Average Score2' : y[j] = link.avgscore2 
      if score == 'Maximum Score1' : y[j] = link.maxscore1 
      if score == 'Maximum Score2' : y[j] = link.maxscore2 
      if score == 'Sum of Score1' : y[j] = link.sumscore1 
      if score == 'Sum of Score2' : y[j] = link.sumscore2 
      if score == 'Number of Scans' : y[j] = link.nscans 
      if score == 'Number of Species' : y[j] = link.nspecies 

    # Compute pbs
    pbs[i] = point_biserial(consistency,y)

  return x, pbs
    
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

  













