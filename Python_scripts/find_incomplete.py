# this file runs a single simulation, and takes care of all the parameters except "noise"
# which is given.
import fileinput
import os
from os.path import isdir, isfile, join
import time
import sys
sys.path.insert(0, r'F:\Alon\Code\Scripts\Utility')
from utils import *
import math
import datetime
import shutil

# main
if len(sys.argv) < 2:
    simdir = '.'
else:
    simdir = sys.argv[1]
sims = 0
complete = 0
incomplete = 0
percent_done = []
last_times = []
print 'COMPLETE:'
for output_dir in os.listdir(simdir):
    full_dir = join(simdir,output_dir)
    if isdir(full_dir) and output_dir.startswith('output'):
        sims += 1
        suffix = output_dir[len('output'):]
#       print join(simdir, 'input' + suffix)
        try:
            if isfile(join(simdir, 'input' + suffix)):
                input = read_input(join(simdir, 'input' + suffix))
                shutil.copyfile(join(simdir, 'input' + suffix), join(full_dir, 'input' + suffix))
            else:
                input = read_input_from_dir(full_dir)
        except:
            continue
        last_time = last_time_in_direc(full_dir)#, input=input)
        last_times.append(last_time)
        if last_time < input.simtime:
        #if last_time < 5000:
            incomplete += 1
            percent_done.append(float(last_time / input.simtime))
            newname = join(r'G:\Alon\Important\1000particles_vicsek\incomplete', output_dir)
            if last_time < 500.0:
                print full_dir, 'is REEEALLLY slow'
            print last_time,input.simtime, suffix
            #shutil.move(full_dir, newname)
        else:
            complete += 1
            #print last_time,input.simtime, ' -    ',suffix
            
print 'Total simulations: {0}'.format(sims)
print '                          Complete: {0} ({1:.2f}%)'.format(complete, 100.0*float(complete)/sims)
print '                        Incomplete: {0} ({1:.2f}%)'.format(incomplete, 100.0*float(incomplete)/sims)
print 'Average % done for incomplete sims: {0:.2f}%'.format(100.0*np.mean(percent_done))
print '       Average % done for all sims: {0:.2f}%'.format(100.0*np.mean(percent_done + [1.0]*complete))
print '      Average last time (all sims): {0:.2f}'.format(np.mean(last_times))
print '                  Slowest sim time: {0:.2f}'.format(np.min(last_times))