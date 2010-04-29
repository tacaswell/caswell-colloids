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
import os.path
import re

_base_path = '/home/tcaswell/colloids/data/polyNIPAM_batch_12/'

def check_existance(conn,fname):
    '''Checks the connection to see if there is already an entry for
    this filename.  Returns bool'''
    res = conn.execute("select key from dsets where fname")
    if res.rowcount >0:
        return True
    
    return False

def guess_temp(fname):
    """Parses the file name and attempts to guess the temperature"""
    fname = fname.split('/')[-1]
    match = re.findall('\d\d-\d',fname)
    if len(match) >0:
        return min(map(lambda x: float(x.replace('-','.')),match))

    
    match = re.findall('\d\d_\d',fname)
    if len(match) >0:    
        return min(map(lambda x: float(x.replace('_','.')),match))


    match = re.findall('\d\d-',fname)
    if len(match) >0:
        return min(map(lambda x: float(x.replace('-','')),match))

    match = re.findall('\d\d_',fname)
    if len(match) >0:
        return min(map(lambda x: float(x.replace('_','')),match))

    return None


def ask_temp(fname):
    """asks the user for the value of the temperature"""
    print "fname: " + fname
    temp = raw_input("enter temperature or None: ")
    try:
        return float(temp)
    except:
        return None

def guess_dtype(fname):
    """Parse the file name and attempt to guess the type"""
    if fname.find('ztl') != -1:
        return 'ztl'
    elif fname.find('z')!=-1:
        return 'z'
    elif fname.find('warming')!=-1:
        return 'ramp'
    else:
        return 't'
    
    pass

def ask_dtype(fname):
    """Asks the user for the dtype"""
    print "fname: " +  fname
    menu = { 1: 't',
             2: 'z',
             3: 'ramp',
             4: 'ztl',
             5: 'None'}
    

    for m in menu:
        print m,menu[m]
    while True:
        dt = raw_input("enter type: ")
        try:
            dt = int(dt)
            if menu.has_key(dt):
                return menu[dt]
        except:
            pass
    pass

def guess_sname(fname):
    """Parse the file name and attempt to guess the sample name"""
    fsplit = fname.split('/')
    if len(fsplit) > 2:
        sname = fsplit[1]
        if fsplit[-1].find('thin') != -1:
            sname = sname + '-thin'

    else:
        sname = ask_sname(fname)

    return sname

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
    print "fname: " + fname
    print 'dtype: '+str(dtype)
    print 'sname: '+str(sname)
    print 'temp: ' +str( temp)
    print 'date: ' +str( date)

    

def strip_base(fname):
    """Takes off a base path the path names"""
    return os.path.realpath(fname).replace(_base_path,'')
    
def guess_meta_data(fname):
    """Calls the specific functions above"""
    return guess_dtype(fname),guess_sname(fname),guess_temp(fname),guess_date(fname)

def ask_meta_data(fname,dtype,sname,temp,date):
    """asks if the values are correct and queries for replacements"""
    print "fname: " + fname
    dtype = query_fun(fname,'dtype',dtype,ask_dtype)
    sname = query_fun(fname,'sname',sname,ask_sname)
    temp = query_fun(fname,'temp',temp,ask_temp)
    date = query_fun(fname,'date',date,ask_date)

    return dtype,sname,temp,date
    
def query_fun(fname,key,value,func):
    """asks if the value is correct for the key, if not, calls func to fix it"""
    print key +": "+ str(value)
    resp = util.query_yes_no('is this correct')
    if resp == 'no':
        value = func(fname)

    return value

def extract_md(fname):
    """takes a guess at parsing, asks for validation, returns tuple"""
    fname = strip_base(fname)
    md = guess_meta_data(fname)
    display_meta_data(fname, *md)
    resp = util.query_yes_no('are these values correct')
    if resp == 'yes':
        return md
    elif resp == 'no':
        return ask_meta_data(fname,*md)

def process_fname(conn,fname):
    """Take in a connection and a database, parse the fname and enter into database"""
    fname = os.path.realpath(fname)
    if check_existance(conn,fname):
        raise "already exists in database"
    md = extract_md(fname)
    print (fname,) + md
    conn.execute('insert into dsets (fname,dtype,sname,temp,date) values (?,?,?,?,?)',(fname,) + md)
    conn.commit()

def visit(conn,dirname,names):
    """Function for walk """

    for f in names:
        fname = dirname + '/'+f
        if os.path.isfile(fname):
            (base,ext) = os.path.splitext(fname)
            if ext == '.tif':
                multi_f = re.findall('file\d{3}',f)
                if len(multi_f) ==0:
                    process_fname(conn,fname)
            
def dir_loop(path,conn):
    """Path is the top level directory to look in, uses walk """
    os.path.walk(path,visit,conn)
