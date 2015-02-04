#!/usr/bin/env python
# ***********************************************************************************
# *    Copyright (C) by Ran Novitsky Nof                                            *
# *                                                                                 *
# *                                                                                 *
# *    This code is free software: you can redistribute it and/or modify                *
# *    it under the terms of the GNU Lesser General Public License as published by  *
# *    the Free Software Foundation, either version 3 of the License, or            *
# *    (at your option) any later version.                                          *
# *                                                                                 *
# *    This program is distributed in the hope that it will be useful,              *
# *    but WITHOUT ANY WARRANTY; without even the implied warranty of               *
# *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
# *    GNU Lesser General Public License for more details.                          *
# *                                                                                 *
# *    You should have received a copy of the GNU Lesser General Public License     *
# *    along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
# ***********************************************************************************
# by Ran Novitsky Nof (ran.nof@gmail.com) @ BSL, 2015

# reads a spyring output and prints a list of stations and latency mean and std
# spyring should be run as:
# spyring -f -d 1 [RING_NAME] > filename

import sys,os
from numpy import array,append

def get_latencies(filename):
  'open a spyring output file and gather latencies for each station'
  D = {}
  f = open(filename,'r')
  lines = f.readlines()
  f.close()
  for line in lines:
    if not 'TRACEBUF2' in line: continue #ignore non-data packets
    d = dict([l.split("=") for l in line.split() if '=' in l])
    if not d['msg'] in D: # d['msg'] is station id
      D[d['msg']]=array([]) # init a list for the station if it appears for the first time
    D[d['msg']]=append(D[d['msg']],float(d['latency'])) # add latency to the list.
  return D

if __name__=='__main__':
  if len(sys.argv)<2:
    sys.exit('Usage: EWlatency.py FILENAME.')
  filename=sys.argv[1]
  if not os.path.isfile(filename):
    sys.exit('%s is not a file.',filename)
  D = get_latencies(filename)
  print('%-15s %-8s %-8s'%('Station','Mean','STD'))
  for k,v in D.items():
    print ('%-15s %-8.3f %-8.3f'%(k,v.mean(),v.std()))
