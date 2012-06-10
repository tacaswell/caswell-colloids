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
import os.path
import re
from datetime import datetime, date

import lib.util  as util




def check_existance(conn,fname):
    '''Checks the connection to see if there is already an entry for
    this filename.  Returns bool'''
    res = conn.execute("select dset_key from dsets where fname = ?",(fname,)).fetchall()
    if len(res) >0:
        return True
    
    return False

def guess_temp(fname):
    """Parses the file name and attempts to guess the temperature"""
    fname = fname.split('/')[-1]
    match = re.findall('\d\d-\d',fname)
    if len(match) >0:
        if len(match) == 3:
            # this picks out the reference temp as reported by the 
            return float(match[1].replace('-','.'))
        else:
            return min(map(lambda x: float(x.replace('-','.')),match))

    # commented out code to deal with badly formed file names
    
    ## match = re.findall('\d\d_\d',fname)
    ## if len(match) >0:    
    ##     return min(map(lambda x: float(x.replace('_','.')),match))


    ## match = re.findall('\d\d-',fname)
    ## if len(match) >0:
    ##     return min(map(lambda x: float(x.replace('-','')),match))

    ## match = re.findall('\d\d_',fname)
    ## if len(match) >0:
    ##     return min(map(lambda x: float(x.replace('_','')),match))

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
    elif fname.find('_z')!=-1:
        return 'z'
    elif fname.find('warming')!=-1 or fname.find('cooling')!=-1 or fname.find('melt')!=-1 or fname.find('heating')!=-1 or fname.find('ramp')!=-1:
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
    # see if there is a file naming the sample
    
    path_split = os.path.dirname(fname).split('/')
    for j in range(len(path_split),6,-1):
        f = '/'.join(path_split[0:(j+1)]) + '/sname.txt'
        if os.path.isfile(f):
            sf = open(f)
            sname = sf.readline().strip()
            sf.close()
            return sname

    sname = ask_sname(fname)
    return sname

def ask_sname(fname):
    """Asks the user for the dtype"""
    print "fname: " + fname
    return raw_input("enter sample name: ")

def guess_date(fname):
    """Tries to guess the date from the name"""
    fsplit = fname.split('/')
    print fname
    for f in fsplit:
        if len(f) == 10 and re.match('\d{4}-\d{2}-\d{2}',f):
            print f
            return datetime.strptime(f,'%Y-%m-%d').date().isoformat()
        if re.match('\d{8}',f):
            return datetime.strptime(f,'%Y%m%d').date().isoformat()
        
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
        
    

def display_meta_data(fname,dtype,sname,date):
    """Displays a nicely formatted output of meta data for confirmation"""
    print "fname: " + fname
    print 'dtype: '+str(dtype)
    print 'sname: '+str(sname)
    #print 'temp: ' +str( temp)
    print 'date: ' +str( date)

    
    
def guess_meta_data(fname):
    """Calls the specific functions above"""
    return guess_dtype(fname),guess_sname(fname),guess_date(fname)

def ask_meta_data(fname,dtype,sname,date):
    """asks if the values are correct and queries for replacements"""
    print "fname: " + fname
    dtype = query_fun(fname,'dtype',dtype,ask_dtype)
    sname = query_fun(fname,'sname',sname,ask_sname)
    #temp = query_fun(fname,'temp',temp,ask_temp)
    date = query_fun(fname,'date',date,ask_date)

    return dtype,sname,date
    
def query_fun(fname,key,value,func):
    """asks if the value is correct for the key, if not, calls func to fix it"""
    print key +": "+ str(value)
    resp = util.query_yes_no('is this correct')
    if resp == 'no':
        value = func(fname)

    return value

def extract_md(fname):
    """takes a guess at parsing, asks for validation, returns tuple"""

    md = guess_meta_data(fname)
    display_meta_data(fname, *md)
    resp = util.query_yes_no('are these values correct')
    if resp == 'yes':
        return md
    elif resp == 'no':
        return ask_meta_data(fname,*md)

def process_fname(conn,fname,ftype):
    """Take in a connection and a database, parse the fname and enter into database"""
    fname = os.path.realpath(fname)
    if check_existance(conn,fname):
        print 'already in data set'
        return
    md = extract_md(fname)
    md = (fname,) + md + (ftype,)
    conn.execute('insert into dsets (fname,dtype,sname,ddate,ftype) values (?,?,?,?,?)',md)
    conn.commit()


def process_fname_dry(conn,fname,ftype):
    """Take in a connection and a database, parse the fname and enter into database"""
    fname = os.path.realpath(fname)
    if check_existance(conn,fname):
        return
    md = extract_md(fname)
    print (fname,) + md
    

def _visit(args,dirname,names):
    """Function for walk """

    (fun,conn) = args
    
    # way to skip directories I don't like
    if dirname.find('stupid') != -1:
        return
    # if there are less than 100 names, assume we are in a directory
    # of folders and stacks rather than a folder with a series in it
    if len(names)<50:
        # loop over names
        for f in names:
            # assemble fqn
            fname = dirname + '/'+f
            # if it is a file, not a directory try to process it
            if os.path.isfile(fname):
                (base,ext) = os.path.splitext(fname)
                # if the file is a tiff
                if ext == '.tif':
                    # and not part a mulit part file
                    multi_f = re.findall('file\d{3}',f)
                    if len(multi_f) ==0:
                        fun(conn,fname,'1')
    else:
        
        # get file names that are not directories
        files_lst = [f for f in names if os.path.isfile(dirname + '/' + f)]
        # split file name
        splt_names = [re.search('(.*?)([0-9]*).tif',f).groups() for f in files_lst
                      if f[-4:] == '.tif']
                
        # get unique base names
        seen = []
        unique_base = [c[0] for c in splt_names if not (c[0] in seen or seen.append(c[0]))]
        del seen
        
        for f in unique_base:
            fun(conn,dirname + '/' + f,'2')
        
def f_print(conn,fname,ftype):
    print fname

def print_loop(path,conn):
    """Path is the top level directory to look in, uses walk """
    os.path.walk(path,_visit,(f_print,conn))

    
def check_loop(path,conn):
    """Path is the top level directory to look in, uses walk """
    os.path.walk(path,visit2,conn)

            
def add_loop(path,conn):
    """Path is the top level directory to look in, uses walk """
    os.path.walk(path,_visit,(process_fname,conn))

            
def dry_loop(path,conn):
    """Path is the top level directory to look in, uses walk """
    os.path.walk(path,_visit,(process_fname_dry,conn))
    
    
