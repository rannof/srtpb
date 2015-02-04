#!/usr/bin/env python

# ***********************************************************************************
# *    Copyright (C) by Ran Novitsky Nof                                            *
# *                                                                                 *
# *    This file is part of srtpb                                                   *
# *                                                                                 *
# *    srtpb is free software: you can redistribute it and/or modify                *
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


# By Ran Novitsky Nof @ BSL, 2015
# ran.nof@gmail.com

# srtpblib is a C library wrapped by swig for python use. must be compiled before use.
from srtpblib import srtpb_record_handler_slink as daliwrite
from srtpblib import *
from obspy.mseed.core import *
from obspy import read
from StringIO import StringIO
import numpy as np
import threading
from threading import Timer
import argparse,os,sys,textwrap

_verbose=False

def getLatency(s):
  '''a argparse function for latency commandline variable
  checking for a number of a file name, and return a number
  or a dictionary based on the file data.
  file should hold a list of NET.STA=latency for stations.
  stations on on the list will have 0 latency.
  '''
  try: # check if its a number like parameter
    l = float(s)
    return l
  except:
    pass
  if not os.path.exists(s): # check if its a file
    msg = "Can't find latency file %r" % s
    raise argparse.ArgumentTypeError(msg)
  try:
    f = open(s,'r')
    lines = [l.strip().split('=') for l in f.readlines()] #reads the lines,split on '=' and create a dictionary
    f.close()
    l = dict([(id,float(latency)) for id,latency in lines])
    return l
  except:
    msg = "File Format Error in %r." % s
    raise argparse.ArgumentTypeError(msg)

class getStartEndTimes(argparse.Action):
  '''A argparse action subclass.
     Called if a start-end time parameters are entered at command-line
  '''
  def __call__(self, parser, namespace, values, option_string=None):
    if len(values)>2: #make sure only two parameters are given
      raise argparse.ArgumentError(self,'Need one or two parameters.')
    start = values[0] # get start time
    if len(values)==2:
      end = values[1] # get end time
    else:
      end=None
    setattr(namespace, 't', [start,end]) # assign them to parameter t

class getReadOptions(argparse.Action):
  '''A argparse action subclass.
     Called if options parameters are entered at command-line
  '''
  def __call__(self, parser, namespace, values, option_string=None):
    try:
      d = dict([v.split('=') for v in values]) # create a dictionary from options
    except ValueError:
      raise argparse.ArgumentError(self,'Each option should be in the form "key=value" (e.g. format=MSEED')
    setattr(namespace, 'o', d) # assign dictionary to parameter o


# set the argument parser
parser = argparse.ArgumentParser(
         formatter_class=argparse.RawDescriptionHelpFormatter,
         description='''Seedlink Reat-Time Playback''',
         epilog=textwrap.dedent('''\
             The code will read seismic waveforms and send it to a seedlink server
             over datalink protocol, simulating a real time behavior, by adjusting
             waveform time to current time.

             Created by Ran Novitsky Nof (ran.nof@gmail.com) @ BSL, 2015'''))

parser.add_argument('p',metavar='PATH',nargs='+',help='Path of data (use quotes for "PATH" if using regex)')
parser.add_argument('-t',metavar='Start',nargs='+',type=UTCDateTime,help='Start and end times. See obspy.UTCDateTime for acceptable string formats. Usually YYYYMMDDhhmmss.f.\nDefault - None',default=[None,None],action=getStartEndTimes)
parser.add_argument('-o',metavar='options',nargs='+',help='additional optional parameters for read function. See obspy.read for help',default={},action=getReadOptions)
parser.add_argument('-s',metavar='seconds',type=float,help='Packet time interval in seconds. Default - 1 sec.',default=1.0)
parser.add_argument('-l',metavar='Latency',type=getLatency,help='Latency of data. Can be one argument for all stations in seconds or a file name \
containing a list of net.sta=X format. Default - 0.',default=0)
parser.add_argument('-H',metavar='Host',help='Host Address (127.0.0.1)',default="127.0.0.1")
parser.add_argument('-P',metavar='Port',help='Port Number (16000)',default="16000")
parser.add_argument('-v',action='store_true',help='Verbose (False).',default=False)
parser.add_argument('--test',action='store_true',help='Testing mode - no data will be sent. (False)',default=False)

def getData(filelist,starttime,endtime,**kwargs):
  '''
  reads the data from the filelist and extract data between start and end.
  kwargs can be any parameter of obspy read function.
  '''
  S = Stream()
  for f in filelist:
    S += read(f,starttime=starttime,endtime=endtime,**kwargs)
  stations = list(set([t.stats.network+'.'+t.stats.station for t in S])) # get a list of stations in NET.STA format
  if len(S)<1:
    raise ValueError,"No data could be found."
  return stations,S

def setLatency(S,latency):
  '''Sets the latency for each trace in stream S
     according to the latency.
     latency can be similar to all stations if latency s float
     or separately assigned according to a dictionary of the form {NET.STA:latency}.
  '''
  if type(latency)==dict:
    [t.stats.__setattr__('latency',latency[t.stats.network+'.'+t.stats.station]) for t in S if latency.has_key(t.stats.network+'.'+t.stats.station)]
  elif isinstance(latency,(int,float)):
    [t.stats.__setattr__('latency',latency) for t in S]
  else:
    return

def sendStream(S,PacketTimeInterval=1.0,latency=0,T0=None,Tpb=None,test=False):
  '''Send a stream to the server.
  S - Stream
  PacketTimeInterval - time span of each packet to be sent
                       each packet is later compressed and
                       sent as multiple 512byte mseed packets
  latency - latency to assign for each trace with no latency parameter
  T0 - the desired start time of the data.
       will be calculated from data itself if None
  Tpb - Time of playback. will be assigned as now if None
  test - setting test to true will not send the data.

  return:
  timers - timer threads that send the data
  Tpb - the Time of playback start
  starttime - the time of first sample
  endtime -  the time of last sample
  '''
  [t.stats.__setattr__('latency',latency) for t in S if not 'latency' in t.stats] # assign latencies for missing stations
  endtime = S.sort(['endtime'])[-1].stats.endtime # get start time of data
  starttime = S.sort(['starttime'])[0].stats.starttime # get end time of data
  NTI = int(np.ceil((endtime-starttime)/float(PacketTimeInterval))) # number of time intervals
  s = Stream() # start a new stream
  if not T0: T0 = starttime # get T0
  if not Tpb: Tpb = UTCDateTime.utcnow() # get Tpb
  if _verbose: print >> sys.stderr,'Sending (%s):\n%-27s %-15s %-27s (%s)'%(Tpb,'System time','ID','StartTime','latency')
  timers=[]
  for i in xrange(NTI): # for each time interval
    s=S.slice(starttime+i*PacketTimeInterval,starttime+(i+1)*PacketTimeInterval) # slice the data at time interval
    timers += [sendTrace(t,T0,Tpb,t.stats.latency+PacketTimeInterval,test=test) for t in s] # send each trace through timer thread
  return timers,Tpb,starttime,endtime

def slinkOpen(IP="127.0.0.1:16000",v=4):
  'open seedlink data connection. Uses srtpblib c-python library'
  if not srtpbdaliinit(IP,"ERTPB",v): sys.exit('Cannot connect to server [%s].\nExit.'%IP)

def slinkClose():
  'close seedlink data connection. Uses srtpblib c-python library'
  srtpbdaliquit()

def sendTrace(t,T0,Tpb,latency=0,test=False):
  '''
  Sends a trace to seedlink server.
  t- Trace
  T0 - time of first sample of the whole data
  Tpb - time of playback start
  latency - the latency of trace
            (should include the packet time interval
             since in real life data is not sent
             until all data is available)
  test - setting to True will not send the data.
  '''
  dt = t.stats.starttime-T0 # time since start of playback
  t.stats.starttime = Tpb+dt # assign starttime relative to time of playback
  if _verbose: print >> sys.stderr,'%-27s %-15s %-27s (%f)'%(UTCDateTime.utcnow(),t.id,t.stats.starttime,t.stats.latency)
  t.data=t.data.astype(np.int32) # make sure data is in integers
  starttime=util._convertDatetimeToMSTime(t.stats.starttime) # get time stamp as seedlink likes it
  endtime=util._convertDatetimeToMSTime(t.stats.endtime)# get time stamp as seedlink likes it
  buff = StringIO()
  writeMSEED(Stream(t), buff,11,512) # write 512byte packets of mseed at STEIM2 compression (from obspy.mseed.core) to buffer
  buff.seek(0) # rewinde buffer
  if not test:
    # create a threaded timer calling delayedwrite function, with buffer,starttime and endtime at starttime+letancy relative to now
    T = Timer((Tpb+dt+latency)-UTCDateTime.utcnow(),delayedwrite,[buff,starttime,endtime])
    # start the thread (will start is time but might have a few hundred miliseconds delay, depending on the system)
    T.start()
  else:
    T=None
  return T

def delayedwrite(buff,starttime,endtime):
  '''
  Sends 512kb segments of the buffer to seedlink server.
  uses daliwrite (srtpb_record_handler_slink in srtpblib)
  '''
  for i in xrange(buff.len/512):
    daliwrite(buff.read(512),512,int(starttime),int(endtime))

def main(args):
  slinkOpen(':'.join([args.H,args.P]),args.v) # open seedlink server for data connection
  filelist = args.p # get the file list
  starttime,endtime = args.t # get time interval
  stations,S = getData(filelist,starttime,endtime,**args.o) # read the data
  #S = S.slice(starttime=S.sort(['starttime'])[0].stats.starttime+60*8)
  setLatency(S,args.l) # set latencies
  timers,Tpb,starttime0,endtime0 = sendStream(S,args.s,test=args.test) # send the stream
  if _verbose: print >> sys.stderr,'Waiting for data to be sent...'
  if not args.test: [t.join() for t in timers] # wait for threads to finish
  slinkClose() # close connection
  print "Playback started at: %s\nData span: %s - %s\nDelta (sec): %s"%(Tpb,starttime0,endtime0,Tpb-starttime0)

if __name__=='__main__':
    # parse the arguments
  args = parser.parse_args(sys.argv[1:])
  if args.v:
    _verbose=True
    args.v=4
  else:
    args.v=0
  main(args)
