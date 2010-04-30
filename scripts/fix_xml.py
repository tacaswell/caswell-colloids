
 
import os.path
import xml.dom.minidom

conv_mapping =  { 'top_cut' :'p_rad',
 'p_rad' : 'd_rad',
 'd_rad': 'mask_rad'}

def convert_bad_xml(fname):
    doc = xml.dom.minidom.parse(fname)

    stanza = doc.getElementsByTagName('iden').item(0)

    prams = stanza.getElementsByTagName("param")



    for p in prams:
        if p.getAttribute('type') == 'int' and p.getAttribute('key') == 'top_cut': 
            p.setAttribute('key','mask_rad')
            p.setAttribute('value','4')

    f = open(fname,'w')
    doc.writexml(f,addindent='   ',newl='\n')
    f.close()

    

def visit(jnk,dirname,names):
    """Function for walk """
    
    for f in names:
        fname = dirname + '/'+f
        if os.path.isfile(fname):
            (base,ext) = os.path.splitext(fname)
            if ext == '.xml':
                convert_bad_xml(fname)


if __name__ == '__main__':
    path = '/home/tcaswell/colloids/data/polyNIPAM_batch_12/20100426'
    os.path.walk(path,visit,None)

                                               
