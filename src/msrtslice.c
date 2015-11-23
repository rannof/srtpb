
/***********************************************************************************/
/*    Copyright (C) by Ran Novitsky Nof                                            */
/*                                                                                 */
/*    This Code is free software: you can redistribute it and/or modify            */
/*    it under the terms of the GNU Lesser General Public License as published by  */
/*    the Free Software Foundation, either version 3 of the License, or            */
/*    (at your option) any later version.                                          */
/*                                                                                 */
/*    This program is distributed in the hope that it will be useful,              */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of               */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                */
/*    GNU Lesser General Public License for more details.                          */
/*                                                                                 */
/*    You should have received a copy of the GNU Lesser General Public License     */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>.        */
/***********************************************************************************/

// msrtslice reads mseed files and slice them to a 1 sec slices (or less) records, ordered by start time.
// The result is one miniseed file containing all data in 1 sec 512 records.
// If 1 sec of data is less than 512 size record, the record will be partly filled.
// Gaps and overlaps should be same as the original data.

// INPUT: a list of files to read

// OUTPUT: a mseedfile with 1 sec slices of the data

// COMPILE:
//gcc -O2 -Wall -L../ringserver/libmseed -I../ringserver/libmseed -lmseed -lpthread msrtslice.c ../ringserver/libmseed/libmseed.a -o msrtslice -lm

//rm mseedtest; gcc -O2 -Wall -L../ringserver/libmseed -I../ringserver/libmseed -lmseed -lpthread mseedtest.c ../ringserver/libmseed/libmseed.a -o mseedtest -lm ; rm output.mseed ; time mseedtest /work/work/GIIRTDATA/201505241233.mseed

#include "libmseed.h"
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h> // for usleep function
#include <stdio.h>
#ifndef WIN32
  #include <signal.h>
  static void term_handler (int sig); // handler function for exit signals
#endif

#define LEN(x) (sizeof(x)/sizeof(x[0]))
#define MAX( a, b ) ( ( a > b) ? a : b )
#define MIN( a, b ) ( ( a < b) ? a : b )
#define esc1 "\x1b[1;31m\x1b[10000;1f\x1b[1A\x1b[K" // write to the last line of the terminal (unless more than 10000 lines)
                                    // also makes text bold
#define esc2 "\x1b[0m" // get things to normal
#define VERSION "[libmseed " LIBMSEED_VERSION "]"
#define PACKAGE "msrtslice"


typedef struct MSDATAINDEX{
  hptime_t starttimehp;
  hptime_t endtimehp;
  char *filename;
  off_t pos;
  char ID[20];
  double sps;
  int nsmpls;
  struct MSDATAINDEX *prev;
  struct MSDATAINDEX *next;
} MSDataIndex;

// FUNCTIONS
int ms_file2index(); // read file into data index
MSDataIndex *sortindex(MSDataIndex *lut); // sort data index
int printindexnode(MSDataIndex *I); // print index node
int freeindexnode(MSDataIndex *I); // Free index node
int addtrace(char *trace);//add a trace ID to the traces list
MSDataIndex *merge_sort_list(MSDataIndex *list,int (*compare)(MSDataIndex *one,MSDataIndex *two)); // merge sort no recursion
int compare(MSDataIndex *one,MSDataIndex *two);//compare two nodes merge sort util func
int splitindex(); // read data according to index and split it to 1sec slices
MSRecord *slice_msr(MSRecord *msr,hptime_t slicestart,hptime_t sliceend); // slice mseed record according to start and end times
static int parameter_proc (int argcount, char **argvec); // process user params

// VARIABLES
int verbose = 0; // libmseed verbose level
void log_print (char *message); // print lOG messages
void diag_print (char *message); // print ERR messages
char **filelist; // A list of input files
char *outputfile=NULL; // output file name
int numoffiles=0; // counter for input files
hptime_t SEC=HPTMODULUS; // Time slices HPTMODULE should be 1s in us units
hptime_t mintimehp=LLONG_MAX,maxtimehp=LLONG_MIN; // minimum and maximum times of input
hptime_t userstart=HPTERROR,userend=HPTERROR; // user defined start and stop times
MSDataIndex *lut=NULL; // Data index look up table
MSDataIndex *lastrec=NULL; // last record of lut
char **traces=NULL; // list of traces IDs
int lentraces=0; // number of traces IDs
flag overwrite=0; // 1 overwrite output file
MSRecord *msr=NULL; // mseed record
MSFileParam *msfp=NULL; // mseed file parameter object
char starttimeiso[40],endtimeiso[40];

/***************************************************************************
 * usage:
 * Print the usage message and exit.
 ***************************************************************************/
static void
usage (void)
{
  fprintf (stderr, "%s version: %s\n\n", PACKAGE, VERSION);
  fprintf (stderr, "Usage: %s [options] -o outfile infile(s)\n\n", PACKAGE);
  fprintf (stderr,
     " ## Options ##\n"
     " -V             Report program version\n"
     " -h             Show this usage message\n"
     " -v             Be more verbose, multiple flags can be used\n"
     " -tw [begin]~[end]\n"
     "        specify a time window in year,month,day,hour,min,sec format\n"
     "        example: -tw 2002,08,05,14,00,00~2002,08,05,14,15,00\n"
     "                 -tw 2002/08/05 14:00:00~\n"
     "                 -tw ~2002-08-05T14:00:00\n"
     "        the start or end time are optional, but the tilde must be present\n"
     " -a Append to output file, default is to overwrite\n"
     " -o outfile     Specify the output file, required\n"
     "\n"
     " infile         Input Mini-SEED file(s)\n"
     "\n");
}  /* End of usage() */

int main(int argc,char **argv){
  int nrecs=0,rrecs=0;
  MSDataIndex *I=NULL;

#ifndef WIN32
  /* Signal handling, use POSIX calls with standardized semantics */
  struct sigaction sa;

  sa.sa_flags   = SA_RESTART;
  sigemptyset (&sa.sa_mask);

  sa.sa_handler = term_handler;
  sigaction (SIGINT, &sa, NULL);
  sigaction (SIGQUIT, &sa, NULL);
  sigaction (SIGTERM, &sa, NULL);

  sa.sa_handler = SIG_IGN;
  sigaction (SIGHUP, &sa, NULL);
  sigaction (SIGPIPE, &sa, NULL);
  /*****************************************************************/
#endif

  // Process given parameters (command line and parameter file)
  if (parameter_proc (argc, argv) < 0) return -1;
//  ms_loginit (&log_print, "LOG: ", &diag_print, "ERR: ");// init logging
  if ((nrecs = ms_file2index())<=0){ // scan files and create index of records
    exit(1);
  }
  if (userstart!=HPTERROR && mintimehp<userstart) mintimehp=userstart;
  if (userend!=HPTERROR && maxtimehp>userend) maxtimehp=userend;
  mintimehp = ((hptime_t)(mintimehp/HPTMODULUS))*HPTMODULUS; // round to seconds
  maxtimehp = ((hptime_t)(maxtimehp/HPTMODULUS+0.5))*HPTMODULUS; // round up to seconds
  if (maxtimehp<=mintimehp){
    ms_log(2,"Data range not in user requested time range\n");
    exit(1);
  }
  if (verbose) ms_log(1,"Sorting %d records...\n",nrecs);
  lut = merge_sort_list(lut,compare); // sort index by start time
  if (verbose) ms_log(1,"Total number of traces: %d\n",lentraces);
  if (verbose) ms_log(1,"Start: %s\nEnd: %s\n",ms_hptime2isotimestr(mintimehp,starttimeiso,0),ms_hptime2isotimestr(maxtimehp,endtimeiso,0));
  rrecs = splitindex(); // split the records
  if (verbose) ms_log(1,"Efficiency: %0.2lf%%\n",nrecs/(double)rrecs*100); // ratio between input records and records read times. 100% means each record was read once.
  // free index list nodes
  if (lut){
    while (lut){
      I=lut;
      lut=lut->next;
      //printindexnode(I);
      freeindexnode(I);
    }
  }
  return 0;
}

int ms_file2index(){
  // Read file into data index
  MSDataIndex *currec=NULL;
  int nrecs=0; // total number of records
  char ID[20];
  off_t fpos=0;
  int last,currentfile=0;
  int retcode=0;
  char *msfile=NULL;

  for(currentfile=0;currentfile<numoffiles;currentfile++){
    msfile = filelist[currentfile];
    if (verbose) ms_log(1,"Reading %d) %s...\n",currentfile+1,msfile);
    fpos=0;last=0;retcode=0;
    while ( (retcode = ms_readmsr_r (&msfp, &msr, msfile, 0, &fpos, &last,1, 0, 0)) == MS_NOERROR )//read record from file
    {
      // check if record end time in user time frame
      if ( ! (userstart==HPTERROR) && msr_endtime(msr)+HPTMODULUS/msr->samprate <= userstart ) continue;
      if ( ! (userend==HPTERROR) && msr->starttime >= userend ) continue;
      msr_srcname(msr,ID,0); // get ID
      addtrace(ID); // add trace to traces list
      if ( ! (currec = (MSDataIndex *) malloc (sizeof(MSDataIndex) )))//allocate index node
      {
        ms_log (2, "Error, cannot allocate buffer for Index table\n");
        return -1;
      }
      // populate index node
      if ( ! (currec->filename = (char*)malloc(sizeof(char)*strlen(msfile))))// allocate for file name
      {
        ms_log (2, "Error, cannot allocate buffer for file name\n");
        return -1;
      }
        strcpy(currec->filename,msfile);// get file name
        currec->pos=fpos; // get position on file
        strcpy(currec->ID,ID); // get trace ID
        currec->sps=msr->samprate;
        currec->nsmpls=msr->samplecnt;
        currec->starttimehp=msr->starttime; // start time as hp
        currec->endtimehp=msr_endtime(msr)+HPTMODULUS/msr->samprate; // end time - this is the actual end of record not the last sample time
        currec->next=NULL;
        nrecs++; // increase number of records
        mintimehp = MIN(mintimehp,currec->starttimehp);// get min time
        maxtimehp = MAX(maxtimehp,currec->endtimehp);// get max time
        if (  lut == NULL ){ // in case we start a new list
          lut=currec;} // list starts at first record
        else
          {lastrec->next=currec;} // or we add the record at the end of the last one
	      lastrec = currec; // now last record is current one
	      //printindexnode(currec);
    };
    if ( retcode != MS_ENDOFFILE ) // if we didn't read until end of file for some reason
      ms_log (2, "Error reading input file %s: %s\n", msfile, ms_errorstr(retcode));
    /* Cleanup memory and close file */
    ms_readmsr_r (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, verbose);
  }
  return nrecs;
};

// Add a trace id to the traces list
int addtrace(char *trace){
  int i;
  // make sure trace id is not in the list
  for(i=0;i<lentraces;i++) if (!strcmp(trace,traces[i])) return 0;
  // if trace id not in traces
  traces = (char **)realloc(traces,(lentraces+1)*sizeof(char*)); // resize the list
  traces[lentraces] = (char*)malloc((strlen(trace)+1)*sizeof(char)); // allocate memory for new trace ID
  strcpy(traces[lentraces],trace); // add ID to list
  lentraces=1+lentraces; // increase number of traces
  return 1;
}

// splitindex reads the mseed records,
// slice them  according to slices of 1 sec using slice_msr
// and group them to traces.
// then, each 1 sec trace group slice is written to the output.
int splitindex(){
  int retcode=0,nrecs=0,orecs=0,i=0,totalsecDt=0;
  hptime_t currslicestarttime=0,currsliceendtime=0; // current slice end time
  char starttimeiso[40]="";
  MSDataIndex *currnode=NULL,*laststartnode=NULL; // index aux pointers
  MSTraceGroup *mstg=NULL; // mseed trace group for each slice
  off_t fpos=0; // position in file
  int last=0; // is this the last record in file?
  time_t tick=NULL,tock=NULL; // processing time measurement variables
  double DT=0,tot=0,Hz=0; // processing time measurement variables
  totalsecDt=round((maxtimehp-mintimehp)/SEC+0.5);//total time span in seconds i.e number of slices
  currslicestarttime = mintimehp-SEC; // start from minimal time available
  currsliceendtime = currslicestarttime+SEC; // start from minimal time available
  laststartnode=lut; // init to start of index list
  tock = time(NULL); // init timestamp for slice processing start time
  tick=tock; // init timestamp for slice processing end time
  fprintf(stderr,"\n");
  for(i=0;i<=totalsecDt;i++){ // for each slice
    tock=time(NULL); // timestamp end of last process step
    DT = difftime(tock,tick); // processing time per slice (sec)
    tot += DT; // total time processing
    tick = tock; // end time is the new start time
    currslicestarttime = currsliceendtime; // get the starttime of new slice
	  currsliceendtime = currslicestarttime+SEC; // This is the end time of the segment needed (not the time of last sample since different traces has different sampling rate
	  currnode=laststartnode; // start at last starting point to save some time to find the correct record needed for slice.
    Hz = tot>1e-6 ? (double)i/tot : 0; // rate of processing time per slice
	  fprintf(stderr,"%sProcessing: (%s) %0.1lf%% %0.1lf s (@%0.2lfHz) Remain: %0.1lf s%s\n",esc1,ms_hptime2isotimestr(currslicestarttime,starttimeiso,0),100.0*i/totalsecDt,tot,Hz,(totalsecDt-i)/Hz,esc2); // update process status
	  while (currnode && currnode->endtimehp<=currslicestarttime){
	    currnode=currnode->next;// find the correct starting point
	  }
	  if (!currnode){
	    ms_log(1,"Got to the end of index list.\n");// todo: check if currnode is null and exit
	    break;
	  }
	  laststartnode=currnode; // update the new starting point
	  mstg = mst_initgroup(mstg); // start a new group for slice;
	  while (currnode && currnode->starttimehp<currsliceendtime){ // make sure we are not after slice end time frame
	    if(currnode->endtimehp-currslicestarttime < 0.5*HPTMODULUS/currnode->sps || // make sure we are not before slice start time
	       currsliceendtime-currnode->starttimehp < 0.5*HPTMODULUS/currnode->sps ){ // or after slice end time
        	      currnode=currnode->next;
	      continue;
	    }
	    fpos = -currnode->pos; // get record position in file. negative number will read from this position
      // read msr inc. data from file
	    if (!((retcode = ms_readmsr_r(&msfp, &msr, currnode->filename, 0, &fpos, &last,0, 1, 0)) == MS_NOERROR)){
	      ms_log(2,"Can't read file %s @ %d : %s\n",currnode->filename,fpos,ms_errorstr(retcode));
	      // free msr
	      ms_readmsr_r (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, verbose);
	      continue;
	    };
	    nrecs++; // increase number of input records
	    // extract only the needed parts of the record data to put to slice
	    msr = slice_msr(msr,currslicestarttime,currsliceendtime);
      // add msr to tracelist
	    mst_addmsrtogroup(mstg,msr,0,0,0);
	    // free msr
	    ms_readmsr_r (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, verbose);
      // get next record
      currnode = currnode->next;
	  }
	  mst_groupheal(mstg,-1,-1); // merge traces to a continues stream
	  if (verbose>1) mst_printtracelist(mstg,1,1,1); // some info about the trace
    // Send data to output
    orecs+=mst_writemseedgroup(mstg,outputfile,0,512,DE_STEIM2,1,0); // big endiens!
	  // free trace list
	  mst_freegroup(&mstg);
  };
  if (verbose) ms_log(1,"Number of records read: %d write: %d\n",nrecs,orecs);
  return nrecs;
};

// slice mseed record to time limits. sliceend is the true end of the slice.
// if record starts later than the slice start time, no slicing from start will occure.
// if record ends before slice end time, no slicing from end will occure.
// return value is a new record with adjusted time span and data
MSRecord *slice_msr(MSRecord *msr,hptime_t slicestart,hptime_t sliceend){
  hptime_t starttime,endtime;
  int slicefromstart=0,slicefromend=0; // number of data points to slice from start and end
  int samplesize = 0;
  int64_t numsamples=0;
  void *datasamples=NULL;
  char srcname[50]="";
  if ( msr->samplecnt <= 0 || msr->samprate <= 0.0 )
    return msr; // no data no work!
  numsamples = msr->numsamples;
  if ( (samplesize = ms_samplesize(msr->sampletype)) == 0 )
    {
      ms_log (2, "slice_msr(): Unrecognized sample type: '%c'\n\n",msr->sampletype);
      // free msr
      ms_readmsr_r (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, verbose);
      return NULL;
    }
  starttime = msr->starttime;
  endtime = msr_endtime(msr)+HPTMODULUS/msr->samprate;
  if ( endtime == HPTERROR )
    {
      ms_log (2, "slice_msr(): Error calculating record end time\n");
      // free msr
      ms_readmsr_r (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, verbose);
      return NULL;
    }
  //msr_print(msr,1);
  msr_srcname(msr,srcname,0);
  slicefromstart = (int)(((double)(slicestart-starttime)/HPTMODULUS) * msr->samprate+0.5);
  slicefromend = (int)(((double)(endtime-sliceend)/HPTMODULUS) * msr->samprate+0.5);
  if ((int)round((((double)(endtime-sliceend)/HPTMODULUS) * msr->samprate-slicefromend)*10)==-5){
    slicefromend--;
  }
  if ((int)round((((double)(slicestart-starttime)/HPTMODULUS) * msr->samprate-slicefromstart)*10)==5){
    slicefromstart++;
  }
  slicefromstart = (slicefromstart>0)? slicefromstart : 0 ;
  slicefromend = (slicefromend>0)? slicefromend : 0 ;
  numsamples = (slicefromstart)? numsamples-slicefromstart : numsamples ;
  numsamples = (slicefromend)? numsamples-slicefromend : numsamples ;
  //ms_log (1, "*** %s (%d,%d)| from s: %d (%ld-%ld) from e: %d (%ld-%ld)\n\n",msr_srcname(msr,srcname,0),msr->numsamples,msr->samplecnt,slicefromstart,slicestart,starttime,slicefromend,endtime,sliceend);
  if(numsamples<=0){
    ms_log (1, "slice_msr(): No data left in timespan\n\t %s (%d,%d)| from s: %d (%ld-%ld) from e: %d (%ld-%ld)\n\n",msr_srcname(msr,srcname,0),msr->numsamples,msr->samplecnt,slicefromstart,slicestart,starttime,slicefromend,endtime,sliceend);
    // free msr
    ms_readmsr_r (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, verbose);
    return NULL;
  }
  if(numsamples>msr->numsamples) ms_log (2, "too much, dude\n\n");
  if ((datasamples = (char*)malloc(numsamples*samplesize)) == NULL){
    ms_log (2, "slice_msr(): Can't allocate memory\n\n");
    // free msr
    ms_readmsr_r (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, verbose);
    return NULL;
  }
  memcpy(datasamples,(char*)msr->datasamples+slicefromstart*samplesize,numsamples*samplesize);
  free(msr->datasamples);
  msr->datasamples=datasamples; // update data
  msr->numsamples=numsamples;
  msr->samplecnt=numsamples;
  msr->starttime=starttime+(hptime_t)(slicefromstart/msr->samprate*HPTMODULUS+0.5);
  return msr;
}

void log_print (char *message) {
  /* Send message to external log message facility */
  //send_log(message);
  printf("%s",message);
}

void diag_print (char *message) {
  /* Send message to external error message facility */
  //send_error(message);
  printf("%s",message);
}

// print a node in the index list
// node holds info of the mseed record, file, position in file, time span etc.
int printindexnode(MSDataIndex *I){
 char starttimeiso[40],endtimeiso[40];
 ms_log(1,"%s %d %s %s %s %.1lf %d\n",I->filename,I->pos,I->ID,ms_hptime2isotimestr(I->starttimehp,starttimeiso,1),ms_hptime2isotimestr(I->endtimehp,endtimeiso,1),I->sps,I->nsmpls);
 return 1;
}

// free an index node
int freeindexnode(MSDataIndex *I){
  free(I->filename);
  free(I);
  return 1;
}
// compare two list members order by start time
int compare(MSDataIndex *one,MSDataIndex *two){
  if(one->starttimehp < two->starttimehp) return -1;
  else return 0;
}

// sort the index list by start time
MSDataIndex *merge_sort_list(MSDataIndex *list,int (*compare)(MSDataIndex *one,MSDataIndex *two))
{
    int listSize=1,numMerges,leftSize,rightSize;
    MSDataIndex *tail,*left,*right,*next;
    if (!list || !list->next) return list;  // Trivial case

    do { // For each power of two<=list length
        numMerges=0,left=list;tail=list=0; // Start at the start

        while (left) { // Do this list_len/listSize times:
            numMerges++,right=left,leftSize=0,rightSize=listSize;
            // Cut list into two halves (but don't overrun)
            while (right && leftSize<listSize) leftSize++,right=right->next;
            // Run through the lists appending onto what we have so far.
            while (leftSize>0 || (rightSize>0 && right)) {
                // Left empty, take right OR Right empty, take left, OR compare.
                if (!leftSize)                  {next=right;right=right->next;rightSize--;}
                else if (!rightSize || !right)  {next=left;left=left->next;leftSize--;}
                else if (compare(left,right)<0) {next=left;left=left->next;leftSize--;}
                else                            {next=right;right=right->next;rightSize--;}
                // Update pointers to keep track of where we are:
                if (tail) tail->next=next;  else list=next;
                // Sort prev pointer
                //next->prev=tail; // Optional.
                tail=next;
            }
            // Right is now AFTER the list we just sorted, so start the next sort there.
            left=right;
        }
        // Terminate the list, double the list-sort size.
        tail->next=0,listSize<<=1;
    } while (numMerges>1); // If we only did one merge, then we just sorted the whole list.
    return list;
}

/***************************************************************************
 * parameter_proc:
 *
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
parameter_proc (int argcount, char **argvec)
{
  int optind;
  char openflag[3] = "wb";
  FILE *outfile=NULL;
  char *timewin     = 0;
  filelist = malloc((size_t)(argcount-1)*sizeof(char*));//save room for file names assuming all parameters are filenames

  /* Process all command line arguments */
  for (optind = 1; optind < argcount; optind++){
      if (strcmp (argvec[optind], "-V") == 0){
        ms_log (1, "%s version: %s\n", PACKAGE, VERSION);
        exit (0);
      }
      else if (strcmp (argvec[optind], "-h") == 0){
        usage();
        exit (0);
      }
      else if (strncmp (argvec[optind], "-v", 2) == 0){
        verbose += strspn (&argvec[optind][1], "v");
      }
      else if (strcmp (argvec[optind], "-a") == 0){
        strcpy(openflag,"ab");
      }
      else if (strcmp (argvec[optind], "-tw") == 0){
        timewin = argvec[++optind];
        char *timeptr;
        if (strchr (timewin, '~') == NULL){
          ms_log (2,"time window not in begin~[end] format\n");
          exit (1);
        }
        if ((timeptr = strtok(timewin,"~"))==NULL){
          ms_log (2,"time window must specify a begin time\n");
          exit (1);
        }
        if (timewin[0]!='~'){
          if ((userstart = ms_timestr2hptime(timeptr))==HPTERROR){
            ms_log (2,"malformed start time window specification\n");
            exit (1);
          };
        }
        else{
          if ((userend = ms_timestr2hptime(timeptr))==HPTERROR){
            ms_log (2,"malformed end time window specification\n");
            exit (1);
          };
        };
        if ((timeptr = strtok(NULL,"~"))==NULL) continue;
        if ((userend = ms_timestr2hptime(timeptr))==HPTERROR){
          ms_log (2,"malformed end time window specification\n");
          exit (1);
        }
        if (userstart>=userend){
          ms_log (2,"start time window should be less than end time window\n");
          exit (1);
        }
      }
      else if (strcmp (argvec[optind], "-o") == 0){
        outputfile = argvec[++optind];
      }
      else if (strncmp (argvec[optind], "-", 1) == 0 && strlen (argvec[optind]) > 1 ){
        ms_log (2, "Unknown option: %s\n", argvec[optind]);
        exit (1);
      }
      else{
        filelist[numoffiles] = argvec[optind];
        numoffiles++;
      }
  }
  /* Make sure an inputfile was specified */
  if ( ! filelist ){
      ms_log (2, "No input file(s) specified\n\n");
      ms_log (1, "%s version %s\n\n", PACKAGE, VERSION);
      ms_log (1, "Try %s -h for usage\n", PACKAGE);
      exit (1);
  }
  /* Make sure an outputfile was specified */
  if ( ! outputfile ){
    ms_log (2, "No output file was specified\n\n");
    ms_log (1, "Try %s -h for usage\n", PACKAGE);
    exit (1);
  }
  else if ( (outfile = fopen(outputfile, openflag)) == NULL ){
    ms_log (2, "Error opening output file: %s\n", outputfile);
    exit (1);
  }
  else{
    fclose(outfile);
  }
  /* Report the program version and other parameters*/
  if ( verbose ){
    ms_log (1, "%s version: %s\n", PACKAGE, VERSION);
    ms_log (1, "Output file: %s\n", outputfile);
    ms_log (1, "Input files: %d\n", numoffiles);
    if (timewin && userstart!=HPTERROR) ms_log (1, "User Start: %s\n",ms_hptime2isotimestr(userstart,starttimeiso,0));
    if (timewin && userend!=HPTERROR) ms_log (1, "User End: %s\n",ms_hptime2isotimestr(userend,endtimeiso,0));
  }
  return 0;
}  /* End of parameter_proc() */


#ifndef WIN32
// handle tremination signal
static void term_handler (int sig)
{
 fprintf(stderr,"\n\tTerminated by user.\n\tExit.\n");
 exit(0);
}
#endif
