# SRTPB: Seedlink Real-Time Playback
## A tool for simulating a real time blayback of seismic data.
  
  This code aim is to mimic a real life scenario of seismic data stream  
  by resending archived seismic data to a seedlink server as in real time.  
  The code reads the seismic data, slice it to small time intervals, and sends  
  each slice to a seedlink server in real time, optionally altering the start time of the  
  slice with respect to the time of the playback such that packets arrived to  
  the server will show current real time.  
  The code allows to set a parameter of latency for each station, for realistic  
  situations of telemertry issues.  
  Other options include sending the data in accelerated time mode or streaming
  to a file.  
  
  
 DEPENDENCIES:  
  swig - http://www.swig.org  
  python (tested on 2.7) modules:  
   numpy - http://www.numpy.org  
   obspy - http://www.obspy.org  
  C external software (included as tar files):  
   code sould inclide software and libraries from IRIS (www.iris.edu):  
    ringserve - https://seiscode.iris.washington.edu/projects/ringserver   
    libdali - http://ds.iris.edu/pub/programs/ringserver/dalitool-2013.280.tar.gz  
    slinktool - http://ds.iris.edu/ds/nodes/dmc/software/downloads/slinktool    
   
 INSTALL:  
  on terminal at main directory, run:  
   INSTALL  
  final codes need will be located in bin directory wich can than be moved.  
   
 RUNNING:  
   once you have you waveform files ready, run:  
     ringserver ring.conf  
   this will start a seedlink server to connect to for getting the packets.  
   Start the playback:  
     srtpb.py [commandline options]  
   on a separate terminal you can use slinktool to see the packets arrive:  
     slinktool -p localhost  
   and program that needs to get the waveforms can connect to the ringserver.  
   see ringserver documentations for additional help.  
     
 USAGE:  
   srtpb.py [-h] [-t Start [END]] [-o options [options ...]]  
                [-s seconds] [-l Latency] [-H Host] [-P Port] [-v] [--test]  
                PATH [PATH ...]  
  
   positional arguments:  
    PATH                  Path of data (use quotes for "PATH" if using regex)  
  
   optional arguments:  
    -h, --help            show this help message and exit  
    -t Start [END]        Start and end times. See obspy.UTCDateTime for  
                          acceptable string formats. Usually YYYYMMDDhhmmss.f.  
                          Default - None  
    -o options [options ...]  
                          additional optional parameters for read function. See  
                          obspy.read for help  
    -s seconds            Packet time interval in seconds. Default - 1 sec.  
    -l Latency            Latency of data. Can be one argument for all stations  
                          in seconds or a file name containing a list of  
                          net.sta=X format. Default - 0.  
    -H Host               Host Address (127.0.0.1)  
    -P Port               Port Number (16000)  
    -v                    Verbose (False).  
    --test                Testing mode - no data will be sent. (False)  
  
  
 LICENSE:  
  srtpb is free software: you can redistribute it and/or modify                  
  it under the terms of the GNU Lesser General Public License as published by    
  the Free Software Foundation, either version 3 of the License, or              
  (at your option) any later version.                                            
                                                                                  
  This program is distributed in the hope that it will be useful,               
  but WITHOUT ANY WARRANTY; without even the implied warranty of                 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
  GNU Lesser General Public License for more details.                          
                                                                                 
  You should have received a copy of the GNU Lesser General Public License     
  along with this program.  If not, see <http://www.gnu.org/licenses/>.        

# MSRTSLICE
## A tool for preparing archived data for replays.

  This C code will read mseed files and send them to a mseed output file or a datalink server.  
  The code will slice the data to a 1 sec slices and reorgenize the data according to start time of each slice.
  The Code can handle large amounts of data and is a pure C code. Use -h for help. It is part of the srtpb package. 
