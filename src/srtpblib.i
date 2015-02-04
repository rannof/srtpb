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

%module srtpblib


%{
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "libmseed.h"
#include "libdali.h"
%}

%{
extern int srtpbdaliinit(char *ip,char *progid,int v);
extern int srtpbdaliquit();
extern void srtpb_record_handler_slink(char *record, int reclen,long starttime,long endtime);
%}
int srtpbdaliinit(char *ip,char *progid,int v);
int srtpbdaliquit();
void srtpb_record_handler_slink(char *record, int reclen,long starttime,long endtime);


