#! /usr/bin/env python

import subprocess
import numpy
import os

partition_info=['128s',16] # = [partition,ncores]
# partition_info=['debug',16] # = [partition,ncores]
time_str='4-00:00:00'
project_name=os.getcwd().split('/')[-3]
myemail=os.environ["MYEMAIL"]

#Log the current submission
logstr='''3_force:

***write me***

Codebase: 
Cluster: ganymede
gittag: '''

#Do the git-ing
cmd='git commit -a -m "Commit before run 3_force" > /dev/null'
print(cmd)
subprocess.call(cmd,shell=True)

cmd='gittag.py > temp.temp'
subprocess.call(cmd,shell=True)
logstr+=open('temp.temp','r').readline()
subprocess.call('rm temp.temp',shell=True)

open('3_force.log','w').write(logstr)

#Setup the versionmap and qsub files
vmap_file=open('versionmap.dat','w')
vmap_file.write('vnum\tL\n')

task_file=open('3_force.task','w')
template_file='3_force.template'
template_contents=open(template_file,'r').read()

vnum=0

for L in xrange(32):
	qsub_file=template_file.replace('.template','_'+str(vnum)+'.qsub')
	fout=open(qsub_file,'w')

	contents=template_contents.replace('###',str(vnum))
        contents=contents.replace('*project*',project_name)
	contents=contents.replace('*111*',str(L))
	vmap_file.write(str(vnum)+'\t'+str(L)+'\n')
	task_file.write('bash 3_force_'+str(vnum)+'.qsub\n')
	fout.write(contents)
	fout.close()
	
	vnum+=1
	

n_nodes=int(numpy.ceil(float(vnum)/partition_info[1]))

# Pad to an even multiple of cores per node
n_cores=n_nodes*partition_info[1]
for j in range(vnum,n_cores):
        task_file.write('echo "Fake run"\n')

# Finally output sbatch file
contents=open('3_force.sbatch.template','r').read()
contents=contents.replace('*nnn*',str(n_cores)) # The total number of processors
contents=contents.replace('*NNN*',str(n_nodes)) # The total number of nodes
contents=contents.replace('*ttt*',time_str) # The wall clock time per processor
contents=contents.replace('*partition*',partition_info[0]) # Partition
contents=contents.replace('*myemail*',myemail) # Email
contents=contents.replace('*project*',project_name) # Project name
open('3_force.sbatch','w').write(contents)
