#!/apps/python/bin/python

#    Automated Auger skimming submition script
#    Takes from two different skims
#        
#    Author: Josh Pond
#    Email: jp4cj@virginia.edu, jpond@jlab.org


import os, time
from subprocess import Popen

# Run jsub command on jsub_file
def submit(jsub_file):
    cmd = 'jsub '+jsub_file
    proc = Popen(cmd, shell = True, executable = os.environ.get('SHELL', '/bin/tcsh'), env = os.environ)    
    time.sleep(.5) # wait for .5 sec to avoid hangs in the jsub program

# Generate jsub file 
# jsub_filename: file name for jsub file
# inFile: File for skimming
# i: skim number
# z: output index
def gen(jsub_filename,inFile,i,z):    
    with open(jsub_filename,'w+') as jsub_file:
        jsub_file.write('''\
PROJECT: g12
TRACK: analysis
JOBNAME: nStar_skim
OS: centos7
MEMORY: 1 GB
COMMAND: source /group/clas/builds/centos7/environment.csh; /home/jpond/build/bin/n_starAnalyzer -rfile inFile
INPUT_FILES: '''+inFile+'''
INPUT_DATA: inFile
OUTPUT_DATA: file
OUTPUT_TEMPLATE: /volatile/clas/clasg12/jpond/nStar4Pi_skim_out/skim'''+str(i)+'''/nStar_skimmed_'''+str(z)+'''.root

''')
    return jsub_filename

if __name__ == '__main__': 
    count = 0 # used as output index
    for root, dirs, files in os.walk("/mss/clas/g12/production/pass1/bos/1-1ckaon1ctrk/"): #Walk over skim one (recommended to check module doc)
        for name in files:
            if not os.path.isfile(os.path.join("/volatile/clas/clasg12/jpond/","nStar4Pi_skim_out/","skim1","nStar_skimmed_"+str(count)+".root")): # check if already skimmed
                submit(gen(os.path.join("/volatile/clas/clasg12/jpond/","nStar4Pi_skim_jsub/","skim1","nStar_skim_"+str(count)+".jsub"),os.path.join(root,name),1,count)) #submit
            count+=1

    count = 0
    for root, dirs, files in os.walk("/mss/clas/g12/production/pass1/bos/2-2pos1neg_not_1ckaon1ctrk/"): #Walk over skim two
        for name in files:
            if not os.path.isfile(os.path.join("/volatile/clas/clasg12/jpond/","nStar4Pi_skim_out/","skim2","nStar_skimmed_"+str(count)+".root")):
                submit(gen(os.path.join("/volatile/clas/clasg12/jpond/","nStar4Pi_skim_jsub/","skim2","nStar_skim_"+str(count)+".jsub"),os.path.join(root,name),2,count))
            count+=1
