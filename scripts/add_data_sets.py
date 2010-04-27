#Copyright 2010 Thomas A Caswell
#tcaswell@uchicago.edu
#http://jfi.uchicago.edu/~tcaswell
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or (at
#your option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses>.

import sqlite3 as sq
import lib.util  as util
import lib.general as gen
import trackpy.cpp_wrapper as cw
from datetime import datetime
from datetime import date
import re

_base_path = '/home/tcaswell/colloids/data/polyNIPAM_batch_12/'

def check_existance(conn,fname):
    '''Checks the connection to see if there is already an entry for
    this filename.  Returns bool'''
    res = sq.execute("select key from dsets where fname")

    return True

def guess_temp(fname):
    """Parses the file name and attempts to guess the temperature"""
    pass

def ask_temp(fname):
    """asks the user for the value of the temperature"""
    pass
def guess_dtype(fname):
    """Parse the file name and attempt to guess the type"""
    pass

def ask_dtype(fname):
    """Asks the user for the dtype"""
    pass

def guess_sname(fname):
    """Parse the file name and attempt to guess the sample name"""
    fsplit = fname.split('/')
    if len(fsplit) > 2:
        return fsplit[1]
    else:
        ask_sname(fname)

def ask_sname(fname):
    """Asks the user for the dtype"""
    print "fname: " + fname
    return raw_input("enter sample name: ")

def guess_date(fname):
    """Tries to guess the date from the name"""
    fsplit = fname.split('/')
    if re.match('\d{8}',fsplit[0]):
        return datetime.strptime(fsplit[0],'%Y%m%d').date().isoformat()
    if re.match('\d{4}-\d{2}-\d{2}',fsplit[0]):
        return datetime.strptime(fsplit[0],'%Y-%m-%d').date().isoformat()
        
    return ask_date(fname)

def ask_date(fname):
    """Asks the user what date to use"""
    print "fname: " + fname
    
    while True:
        idate = raw_input("enter the date in ISO format: ")
        if re.match('\d{4}-\d{1,2}-\d{1,2}',idate):
            try:
                return date(*map(int,idate.split('-'))).isoformat()
            except:
                pass
        
    

def display_meta_data(fname,dtype,sname,temp,date):
    """Displays a nicely formatted output of meta data for confirmation"""
    pass

def strip_base(fname):
    """Takes off a base path the path names"""
    return fname.replace(_base_path,'')
    
def guess_meta_data(fname):
    """Calls the specific functions above"""
    return guess_dtype(fname),guess_sname(fname),guess_temp(fname),guess_date

