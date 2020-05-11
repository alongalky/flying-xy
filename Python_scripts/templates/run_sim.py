import concurrent.futures
import os
import glob
import time
import sys
from utils import *
import subprocess
from multiprocessing import Pool
import multiprocessing
from run_single import run_single
import itertools
import argparse
import psutil, os
# parameters
#densityVals = [0.01, 0.05, 0.1, 1.0, 2.0]
#densityVals = [0.1, 0.2, 0.35, 0.5,  0.7, 1.5, 1.0, 2.0]
#densityVals = [0.5, 0.5, 0.5, 0.5]
densityVals = [0.5]

#noiseVals = [[0.005, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]]
"""
noiseVals = [[0.03, 0.05, 0.07, 0.1, 0.11, 0.13, 0.14, 0.16],
             [0.03, 0.05, 0.07, 0.1, 0.11, 0.13, 0.14, 0.16],
             [0.03, 0.05, 0.07, 0.1, 0.11, 0.13, 0.14, 0.16],
             [0.03, 0.05, 0.07, 0.1, 0.11, 0.13, 0.14, 0.16]]
"""
#noiseVals = [[0.01, 0.1, 0.5, 1.0, 1.4, 1.7, 1.8, 1.9,
#              2.0, 2.1, 2.2, 2.5, 2.7, 2.9, 3.0, 3.5]]
noiseVals = [[0.01, 0.1, 0.5, 1.0, 1.4, 1.7, 1.8, 1.9]]
#npartVals = [2048, 4096, 8192, 16384]
npartVals = [8192]
dt = 0.05
norm = 'false'
duration = 1000000.0
minutes_per_sim = 500000

#noiseValsFlyingXY = [3.8, 4.0, 4.2, 4.4, 4.6, 5.0, 5.3, 6.0]
#noiseValsFlyingXY = [3.0]

class runParams():
  def __init__(self,noise,norm, dt, density, npart, duration, maxruntime=2000000):
     self.noise = noise
     self.norm = norm
     self.dt = dt
     self.density = density
     self.npart = npart
     self.duration = duration
     self.maxruntime = maxruntime

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make an .mp4 movie from a simulation output directory')
    parser.add_argument('directory', metavar='dir', type=str, nargs='?',default='.',
                       help='directory with simulation data')                   
    parser.add_argument('-p','--proc', metavar='processors', type=int,# action='store_const', nargs=1,
                   default=multiprocessing.cpu_count(), help='number of CPUs to use')
    parser.add_argument('-o','--override', action='store_true', 
                   default=False, help='Override existing simulation data')

    args = parser.parse_args()
    
    start_time = datetime.datetime.now()
    # main program
    if args.override:
        print printtime() + "Deleting old output files"
        # remove old data
        remove_output_dirs('output')
        remove_inputs_and_errors()
        print printtime() + "Finished removing old output files"
    
    # prepare inputs
    inputs = []
    for noises, density, npart in zip(noiseVals, densityVals, npartVals):
        for noise in noises:
            inputs.append(runParams(noise,norm, dt, density, npart, duration, maxruntime=minutes_per_sim))
    
    
    inputs = inputs
    p = Pool(args.proc)
    seconds = 0
    while True: # never stop!
        print printtime() + "starting {0} processes".format(len(inputs))
        ret = p.map(run_single, inputs)
        break
        
    
    elapsed_time = datetime.datetime.now() - start_time
    rest, seconds = divmod(elapsed_time.total_seconds(), 60)
    hours, minutes = divmod(rest,60)
    print printtime() + "Finished everything in {0} hours {1} minutes {2} seconds".format(int(hours),int(minutes), int(seconds))



#os.system(matlabCmd)
