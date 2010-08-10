#Copyright 2009 Thomas A Caswell
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


import xml.dom.minidom
import subprocess
import datetime
import os.path
import sys
writeSet = ['z-position','acquisition-time-local']
"""
list of meta data, put which ever values you want into writeSet

'Exposure'
'Binning'
'Region'
'Subtract'
'Shading'
'Sensor Mode'
'Digitizer'
'Gain'
'Camera Shutter'
'Clear Count'
'Clear Mode'
'Frames to Average'
'Trigger Mode'
'Temperature'
'MetaDataVersion'
'ApplicationName'
'ApplicationVersion'
'plane-type'
'pixel-size-x'
'pixel-size-y'
'bits-per-pixel'
'autoscale-state'
'autoscale-min-percent'
'autoscale-max-percent'
'scale-min'
'scale-max'
'spatial-calibration-state'
'spatial-calibration-x'
'spatial-calibration-y'
'spatial-calibration-units'
'image-name'
'threshold-state'
'threshold-low'
'threshold-high'
'threshold-color'
'zoom-percent'
'gamma'
'look-up-table-type'
'look-up-table-name'
'photonegative-mode'
'gray-calibration-curve-fit-algorithm'
'gray-calibration-values'
'gray-calibration-min'
'gray-calibration-max'
'gray-calibration-units'
'plane-guid'
'acquisition-time-local'
'modification-time-local'
'stage-position-x'
'stage-position-y'
'stage-label'
'z-position'
'wavelength'
'camera-binning-x'
'camera-binning-y'
'camera-chip-offset-x'
'camera-chip-offset-y'
'_IllumSetting_'
'_MagNA_'
'_MagRI_'
'_MagSetting_'
'number-of-planes'
'dtime'
"""

class txtFile:
    def __init__(self,fname):
        self.f = open(fout,'w')
    def add_entry(self,key,val):
        if key in writeSet:
            self.f.write(key + ": " + str(val) + "\n")
    def start_frame(self,n):
        self.f.write("Frame " + str(n) + "\n")
    def end_frame(self):
        self.f.write("\n\n")

    def close(self):
        self.f.close()
    def __del__(self):
        self.close()

class dict_vecs:
    def __init__(self):
        self.d = {}
        for k in writeSet:
            self.d[k] = []
    def add_entry(self,key,val):
        if key in writeSet:
            self.d[key].append(val)
    def start_frame(self,n):
        pass
    def end_frame(self):
        pass
    def close(self):
        pass
    

        
def _write(file,key,val):
    if key == "acquisition-time-local" or key == "modification-time-local":
        tmp = int(val[18:])
        val = val[:18] + "%(#)03d"%{"#":tmp}
    file.add_entry(key,val)
        


def _start_group(f,n):
    f.start_frame(n)


def _end_group(f):
    f.end_frame()
    
    
def _parse_attr(file_obj,dom_obj):
    if dom_obj.getAttribute("id") =="Description":
        _parse_des(file_obj,dom_obj)
    elif dom_obj.getAttribute("type") =="int":
        _write(file_obj,dom_obj.getAttribute("id"),int(dom_obj.getAttribute("value")))
    elif  dom_obj.getAttribute("type") =="float":
        _write(file_obj,dom_obj.getAttribute("id"),float(dom_obj.getAttribute("value")))
    else: 
        _write(file_obj,dom_obj.getAttribute("id"), dom_obj.getAttribute("value").encode('ascii'))

def _parse_des(file_obj,des_obj):
    des_string = des_obj.getAttribute("value")
    des_split = des_string.split("&#13;&#10;")

    for x in des_split:
        tmp_split = x.split(":")
        if len(tmp_split) ==2:
            _write(file_obj,tmp_split[0],tmp_split[1].encode('ascii'))

def _parse_params(file_obj,fname):
    f = open(fname)
    for line in f:
        _parse_params_attr(file_obj,line)
    f.close()
    n_fname  = fname[:(len(fname)-5)] + ".done"
    os.rename(fname,n_fname)

def parse_file(fin,f):
    

    
    print fin
    # make sure the files exist
    if not (os.path.exists(fin) ):
        print "file does not exist"
        exit()


        
    




    #changed to deal with 2.5v2.6
    #g = f.create_group("frame{0:06d}".format(0))

    _start_group(f,0)

    a = subprocess.Popen(["tiffinfo","-0",fin] ,stdout=subprocess.PIPE)
    tiff_string = (a.stdout).readlines();
    tiff_string = "".join(tiff_string)

    xml_str = tiff_string[(tiff_string.find("<MetaData>")):(10 + tiff_string.rfind("MetaData>"))]
    dom = xml.dom.minidom.parseString(xml_str)

    props = dom.getElementsByTagName("prop")
    
    for p in props:
        if p.parentNode.nodeName == "PlaneInfo":
            _parse_attr(f,p)
            if p.getAttribute("id") == "acquisition-time-local":
                tmp = p.getAttribute("value")
                initial_t  = datetime.datetime.strptime(tmp[:17],"%Y%m%d %H:%M:%S") 
                initial_t.replace(microsecond=int(tmp[18:])*1000)


        elif p.parentNode.nodeName == "MetaData":
            _parse_attr(f,p)
        elif p.parentNode.nodeName == "SetInfo":
            _parse_attr(f,p)
            if p.getAttribute("id") == "number-of-planes":
                frame_count = int(p.getAttribute("value"))

    _write(f,"dtime",0.0)

    _end_group(f)
    
    

    

    for frame in range(1,frame_count):
        _start_group(f,frame)
        a = subprocess.Popen(["tiffinfo","-"+str(frame),fin] ,stdout=subprocess.PIPE)
        tiff_string = (a.stdout).readlines();
        tiff_string = "".join(tiff_string)
        xml_str = tiff_string[(tiff_string.find("<MetaData>")):(10 + tiff_string.rfind("MetaData>"))]
        dom = xml.dom.minidom.parseString(xml_str)
        props = dom.getElementsByTagName("prop")



        for p in props:
            if p.parentNode.nodeName == "PlaneInfo":
                _parse_attr(f,p)
                if p.getAttribute("id") == "acquisition-time-local":
                    tmp = p.getAttribute("value")
                    current_t = datetime.datetime.strptime(tmp[:17],"%Y%m%d %H:%M:%S") 
                    current_t = current_t.replace(microsecond=int(tmp[18:])*1000)
        dt = current_t - initial_t
        _write(f,"dtime",dt.seconds + dt.microseconds/(pow(10.,6)))
        initial_t = current_t
        _end_group(f)
    f.close()

    return f

if __name__ == "__main__":
    if len(sys.argv) !=3:
        print "provide a input file and an output file"
        exit()
    else:
        fin = sys.argv[1]
        fout =  sys.argv[2]
    parse_file(fin,txtFile(fout))
