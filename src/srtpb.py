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
import multiprocessing as mp
from multiprocessing import Pool,Manager,Queue,Lock,Process
import argparse,os,sys,textwrap,time

_verbose=False
_test=False
MAXTHREADS=300 # should be higher than the number of station.channels so at each point of time all data from all stations will be sent.
MAXProcs = mp.cpu_count()-1 # for multiprocessing when slicing a Stream
if not MAXProcs: MAXProcs = 1 # make sure at least one processor is available

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

class procNum(argparse.Action):
  '''A argparse action subclass.
     Called if options parameters are entered at command-line
  '''
  def __call__(self, parser, namespace, value, option_string=None):
    if value>MAXProcs:
      raise argparse.ArgumentError(self,"Can't use more than the system's processors number (max - %d)"% MAXProcs)
    setattr(namespace, 'n', value) # assign dictionary to parameter o

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
parser.add_argument('-m',action='store_true',help='Modify time to current time (relative to playback start time)',default=False)
parser.add_argument('-f',action='store_true',help='Fast replay - No delay (Not a "real time" playback)',default=False)
parser.add_argument('-H',metavar='Host',help='Host Address (127.0.0.1)',default="127.0.0.1")
parser.add_argument('-P',metavar='Port',help='Port Number (16000)',default="16000")
parser.add_argument('-v',action='store_true',help='Verbose (False).',default=False)
parser.add_argument('-d',metavar='delay',type=float,help='Directly stream ordered file to server with delay ( in seconds) every 250 packets.',default=-1)
parser.add_argument('-n',metavar='processors',type=int,help='Number of processors for fast replay',default=1,action=procNum)
parser.add_argument('--test',action='store_true',help='Testing mode - no data will be sent. (False)',default=False)


def slinkOpen(IP="127.0.0.1:16000",v=4):
  'open seedlink data connection. Uses srtpblib c-python library'
  if not srtpbdaliinit(IP,"SRTPB",v): sys.exit('Cannot connect to server [%s].\nExit.'%IP)

def slinkClose():
  'close seedlink data connection. Uses srtpblib c-python library'
  srtpbdaliquit()


def getData(filelist,starttime,endtime,**kwargs):
  '''
  reads the data from the filelist and extract data between start and end.
  kwargs can be any parameter of obspy read function.
  '''
  S = Stream()
  for f in filelist:
    try:
      S += read(f,starttime=starttime,endtime=endtime,**kwargs)
    except Exception,msg:
      print >> sys.stderr,msg.args[0]
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

def getSlice(lst,Qin,Qout,L):
  run = True
  while run:
    try:
      starttime,PacketTimeInterval,i,m,T0,Tpb = Qin.get()
    except:
      run=False
      continue
    s=lst[0].slice(starttime+i*PacketTimeInterval,starttime+(i+1)*PacketTimeInterval)
    for t in s:
      t.stats.starttime = t.stats.starttime-t.stats.latency # shift back traces times according to latencies
      if len(t.data)>1: t.data = t.data[:-1] # avoids overlaps of 1 sample. However, the last sample of a trace in the dataset will be lost.
      if m:
        dt = t.stats.starttime-T0 # time since start of data
        t.stats.starttime = Tpb+dt # assign starttime relative to time of playback
    Qout.put(s)
  L.acquire()
  if _verbose: print >> sys.stderr,'Terminating process',os.getpid()
  L.release()

def putSlice(Qin,L,):
  slinkOpen(':'.join([args.H,args.P]),args.v) # open seedlink server for data connection
  run = True
  while run:
    s = Qin.get()
    if s=='done.':
      run=False
      continue
    buff = StringIO()
    writeMSEED(s, buff,11,512) # write 512byte packets of mseed at STEIM2 compression (from obspy.mseed.core) to buffer
    buff.seek(0) # rewinde buffer
    if not _test:
      for i in xrange(buff.len/512):
        if i%250==0: time.sleep(0.1) # must slow down for ElarmS fast playback mode
        daliwrite(buff.buf[i*512:(i+1)*512],512,0,1)
  slinkClose()


def sendStreamFast(S,PacketTimeInterval=1.0,latency=0,T0=None,Tpb=None,m=False):
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
  m - alter time to current time relative to playback start time.


  return:
  timers - timer threads that send the data
  Tpb - the Time of playback start
  starttime - the time of first sample
  endtime -  the time of last sample
  '''
  tic = UTCDateTime.utcnow()
  if _verbose: print >> sys.stderr,'Sorting started at: %s'%tic
  if _verbose: print >> sys.stderr,'Sorting data in %d original packets...'%len(S)
  S.merge()
  S = S.split()
  [t.stats.__setattr__('latency',latency) for t in S if not 'latency' in t.stats] # assign latencies for missing stations
  [t.stats.__setattr__('starttime',t.stats.starttime+t.stats.latency) for t in S] # shift traces times according to latencies
  for t in S:
    t.data=t.data.astype(np.int32) # make sure data is in integers
  endtime = S.sort(['endtime'])[-1].stats.endtime # get end time of data
  elatency = S[-1].stats.latency
  starttime = S.sort(['starttime'])[0].stats.starttime # get start time of data
  slatency = S[0].stats.latency
  NTI = int(np.ceil((endtime-starttime)/float(PacketTimeInterval))) # number of time intervals
  if _verbose: print >> sys.stderr,'Slicing %d Traces into ~%d packets...'%(len(S),NTI*len(S))
  if not T0: T0 = starttime-slatency # get T0
  if not Tpb: Tpb = UTCDateTime.utcnow() # get Tpb
  # slice the stream to time intervals - use multithreading for speed
  lst.append(S)
  tic = UTCDateTime.utcnow()
  [Qin.put([starttime,PacketTimeInterval,i,m,T0,Tpb]) for i in xrange(NTI)]
  [Qin.put(['Done.']) for i in xrange(args.n)]
  [p.join() for p in SlicerPool]
  toc = UTCDateTime.utcnow()
  if _verbose: print >> sys.stderr,'Data sent to slicing at: %s (%lf)'%(toc,toc-tic)
  [Qout.put('done.') for i in xrange(3)]
  [writer.join() for writer in writers]
  toc = UTCDateTime.utcnow()
  if _verbose: print >> sys.stderr,'Data sent to slink server at: %s (%lf)'%(toc,toc-tic)
  return [],Tpb,starttime,endtime

def sendStream(S,PacketTimeInterval=1.0,latency=0,T0=None,Tpb=None,m=False):
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
  m - alter time to current time relative to playback start time.

  return:
  timers - timer threads that send the data
  Tpb - the Time of playback start
  starttime - the time of first sample
  endtime -  the time of last sample
  '''
  [t.stats.__setattr__('latency',latency) for t in S if not 'latency' in t.stats] # assign latencies for missing stations
  endtime = S.sort(['endtime'])[-1].stats.endtime # get end time of data
  starttime = S.sort(['starttime'])[0].stats.starttime # get start time of data
  NTI = int(np.ceil((endtime-starttime)/float(PacketTimeInterval))) # number of time intervals
  s = Stream() # start a new stream
  if not T0: T0 = starttime # get T0
  if not Tpb: Tpb = UTCDateTime.utcnow() # get Tpb
  if _verbose: print >> sys.stderr,'Sending (%s):\n%-27s %-15s %-27s (%s)'%(Tpb,'System time','ID','StartTime','latency')
  timers=[]
  for i in xrange(NTI): # for each time interval
    s=S.slice(starttime+i*PacketTimeInterval,starttime+(i+1)*PacketTimeInterval) # slice the data at time interval
    timers += [sendTrace(t,T0,Tpb,t.stats.latency+PacketTimeInterval,m=m) for t in s] # send each trace through timer thread
    while threading.activeCount()>MAXTHREADS:
      time.sleep(0.01)
  return timers,Tpb,starttime,endtime

def sendTrace(t,T0,Tpb,latency=0,m=False):
  '''
  Sends a trace to seedlink server.
  t- Trace
  T0 - time of first sample of the whole data
  Tpb - time of playback start
  latency - the latency of trace
            (should include the packet time interval
             since in real life data is not sent
             until all data is available)
  m - alter time to current time relative to playback start time.
  '''
  dt = t.stats.starttime-T0 # time since start of data
  if m: t.stats.starttime = Tpb+dt # assign starttime relative to time of playback
  if len(t.data)>1: t.data = t.data[:-1] # avoids overlaps of 1 sample. However, the last sample of a trace in the dataset will be lost.
  if _verbose: print >> sys.stderr,'%-27s %-15s %-27s (%f)'%(UTCDateTime.utcnow(),t.id,t.stats.starttime,t.stats.latency)
  t.data=t.data.astype(np.int32) # make sure data is in integers
  starttime=util._convertDatetimeToMSTime(t.stats.starttime) # get time stamp as seedlink likes it
  endtime=util._convertDatetimeToMSTime(t.stats.endtime)# get time stamp as seedlink likes it
  buff = StringIO()
  writeMSEED(Stream(t), buff,11,512) # write 512byte packets of mseed at STEIM2 compression (from obspy.mseed.core) to buffer
  buff.seek(0) # rewinde buffer
  if not _test:
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
  filelist = args.p # get the file list
  starttime,endtime = args.t # get time interval
  stations,S = getData(filelist,starttime,endtime,**args.o) # read the data
  setLatency(S,args.l) # set latencies
  if args.f:
    if _verbose: print >> sys.stderr,'Running in Fast mode...'
    timers,Tpb,starttime0,endtime0 = sendStreamFast(S,args.s,m=args.m) # send the stream as fast as possible
  else:
    slinkOpen(':'.join([args.H,args.P]),args.v) # open seedlink server for data connection
    timers,Tpb,starttime0,endtime0 = sendStream(S,args.s,m=args.m) # send the stream
  if _verbose: print >> sys.stderr,'Waiting for data to be sent...'
  if not args.test: [t.join() for t in timers] # wait for threads to finish
  if not args.f: slinkClose() # close connection
  print "Playback started at: %s\nData span: %s - %s\nDelta (sec): %s"%(Tpb,starttime0,endtime0,Tpb-starttime0)


def direct(args):
  assert args.t==[None,None],'Error: Direct Mode cannot use time value.'
  slinkOpen(':'.join([args.H,args.P]),args.v) # open seedlink server for data connection
  filelist = args.p # get the file list
  for file in filelist:
    try:
      buffer = open(file,'r').read()
    except Exception as E:
      print >> sys.stderr,'Error with file %s (%s)'%(file,str(E))
      continue
    for i in xrange(len(buffer)/512):
      if i%250==0: time.sleep(args.d)
      daliwrite(buffer[i*512:(i+1)*512],512,0,1)

if __name__=='__main__':
  startRun = UTCDateTime.now()
  print "Start: %s"%startRun
    # parse the arguments
  args = parser.parse_args(sys.argv[1:])
  if args.v:
    _verbose=True
    args.v=4
  else:
    args.v=0
  if args.test: _test=True
  if args.f:
    # multiprocessing stuff
    if _verbose: print >> sys.stderr,'Staring %d subprocesses for fast slicing.'%args.n
    Qin = Queue(args.n+1)
    Qout = Queue(args.n+1)
    L = Lock()
    lst = Manager().list()
    SlicerPool = [Process(target=getSlice,args=(lst,Qin,Qout,L)) for i in xrange(args.n)]
    [p.start() for p in SlicerPool]
    writers = [Process(target=putSlice,args=(Qout,L)) for i in xrange(3)]
    [writer.start() for writer in writers]
  if not args.d>=0:
    main(args)
  else:
    direct(args)
  endRun = UTCDateTime.now()
  print "End: %s (%lf)"%(endRun,endRun-startRun)

