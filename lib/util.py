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

import sys

class Cords3D:
    """A class for carrying around sets of 3D data"""
    def __init__(self):
        self.x = None
        self.y = None
        self.z = None
        self.I = None
        
    def __iter__(self):
        return _quad_objects(self.x,self.y,self.z,self.I)

class Cords2D:
    """A class for carrying around sets of 2D data"""
    def __init__(self):
        self.x = None
        self.y = None
        self.I = None
        
    def __iter__(self):
        return _triple_objects(self.x,self.y,self.I)

    
def _triple_objects(a,b,c):
    """ generator for triples"""
    for i in range(len(a)):
        yield (a[i],b[i],c[i])
def _quad_objects(a,b,c,d):
    """ generator for quads"""
    for i in range(len(a)):
        yield (a[i],b[i],c[i],d[i])

class dbase_error(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class cord_pairs:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def __iter__(self):
        return itertool.izip(x,y)

## {{{ http://code.activestate.com/recipes/577058/ (r2)
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.
    
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":"yes",   "y":"yes",  "ye":"yes",
             "no":"no",     "n":"no"}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return default
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")
## end of http://code.activestate.com/recipes/577058/ }}}
