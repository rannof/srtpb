/**********************************************************************************
*    Copyright (C) by Ran Novitsky Nof                                            *
*                                                                                 *
*    This file is part of srtpb                                                   *
*                                                                                 *
*    srtpb is free software: you can redistribute it and/or modify                *
*    it under the terms of the GNU Lesser General Public License as published by  *
*    the Free Software Foundation, either version 3 of the License, or            *
*    (at your option) any later version.                                          *
*                                                                                 *
*    This program is distributed in the hope that it will be useful,              *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of               *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
*    GNU Lesser General Public License for more details.                          *
*                                                                                 *
*    You should have received a copy of the GNU Lesser General Public License     *
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
***********************************************************************************/

#include "libmseed.h"
#include "libdali.h"
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h> // for usleep function
#include <stdio.h>

DLCP *dlconn; // dali connection - data protocol for seedlink server
int keep_going=1;

int srtpbdaliinit(char *ip,char *progid,int v); // connect to seedlink server on data protocol
// ip - the seedlink server ip:port as string
// progid - name of program connecting
// v - verbosity (0-4)
// return 1 on connection 0 on failed connection

int srtpbdaliquit(); // disconnect from seedlink server
// return 0 on disconnecting

void srtpb_record_handler_slink (char *record, int reclen,long starttime,long endtime); // send packet to server
// return 0 on disconnecting


/*****************************************************************/

int srtpbdaliinit(char *ip,char *progid,int v)
{// connect to seedlink server
  dlconn = dl_newdlcp (ip, progid);
  if (v)
    dl_loginit (v, NULL, NULL, NULL, NULL); // initialize verbousity
  /* Connect to server */
  if ( dl_connect (dlconn) < 0 )
    {
      //fprintf (stderr, "Error connecting to server %s\n",ip);
      return 0;
    }
  else
    return 1;
}

int srtpbdaliquit()
{// Make sure seedlink connection is shut down
  if ( dlconn->link != -1 )
    dl_disconnect (dlconn);
  if ( dlconn )
    dl_freedlcp (dlconn);
  return 0;
}

void srtpb_record_handler_slink (char *record, int reclen,long starttime,long endtime)
{// send a packet to server
  char streamID[20]; // hold stream ID

  if (keep_going && dl_write(dlconn,record,reclen,ms_recsrcname(record,streamID,0),starttime,endtime,0)<0)
  {// in case of a write error - This part should be changed I think. not sure it works at all.
    dl_disconnect (dlconn); // disconnect from server
    dl_log_r (dlconn, 1, 0, "[%s] Disconnecting and trying to reconnect\n",dlconn->addr);
    while (keep_going && dl_connect (dlconn)<0) // keep trying to reconnect
    {
      dl_log_r (dlconn, 1, 0, "[%s] Reconnecting in 1 min...\n",dlconn->addr);
      sleep(60); // every minute
    }
    if (keep_going) dl_log_r (dlconn, 1, 0, "[%s] Server reconnected!\n",dlconn->addr);
    else dl_log_r (dlconn, 0, 0, "[%s] USER TERM SIGNAL!\n",dlconn->addr);
  }
}
