#!/apps/python/bin/python
import os, time
from subprocess import Popen

#    Automated Auger merging submition script
#    Takes from two different skims
#        
#    Author: Josh Pond
#    Email: jp4cj@virginia.edu, jpond@jlab.org

# Run jsub command on jsub_file
def submit(jsub_file):
    cmd = 'jsub '+jsub_file
    proc = Popen(cmd, shell = True, executable = os.environ.get('SHELL', '/bin/tcsh'), env = os.environ)    
    time.sleep(.5) # wait for .5 sec to avoid hangs in the jsub program

# Generate jsub file 
# jsub_filename: file name for jsub file
# i: skim number
# z: merge index first digit
# v: merge index second digit
def gen(jsub_filename,i,z,v):    
    with open(jsub_filename,'w+') as jsub_file:
        jsub_file.write('''\
PROJECT: g12
TRACK: analysis
JOBNAME: nStar_merge
OS: centos7
MEMORY: 1 GB
COMMAND: source /group/clas/builds/centos7/environment.csh; /home/jpond/build/bin/neutron_skim_merge /volatile/clas/clasg12/jpond/nStar4Pi_skim_out/skim'''+str(i)+'''/nStar_skimmed_'''+str(z)+str(v)+'''*.root
OUTPUT_DATA: g12_N_pmN_merged.root
OUTPUT_TEMPLATE: /work/clas/clasg12/jpond/nStar4Pi_skim_merged/nStar_merged_'''+str(i)+'_'+str(z)+str(v)+'''.root

''')
    return jsub_filename

if __name__ == '__main__':
    count = 0 #output index
    for z in range(1,10,1):
        for v in range(10):# merge all skimmed files with index matching 10*, 11*, 12*, ... 56*, 57* ... etc
            submit(gen(os.path.join("/volatile/clas/clasg12/jpond/","nStar4Pi_merge_jsub/","skim1","nStar_merge_"+str(count)+".jsub"),1,z,v))
            count += 1
    count = 0
    for z in range(1,10,1):
        for v in range(10):
            submit(gen(os.path.join("/volatile/clas/clasg12/jpond/","nStar4Pi_merge_jsub/","skim2","nStar_merge_"+str(count)+".jsub"),2,z,v))
            count += 1
