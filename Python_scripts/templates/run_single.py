# this file runs a single simulation, and takes care of all the parameters except "noise"
# which is given.
import fileinput

import os
import time
import sys
import math
import datetime
import shutil
from utils import *
import subprocess

def run_single(params):
    noise = params.noise    # noise 
    norm = params.norm      # Vicksek or Flying XY (Neighbor normalization)
    dt = params.dt
    # parameters
    #precision   = 0.01	    # 1/ALIGN, NOISE*dt, dt*(V0*RSTAR) will all be no larger than "precision"
    duration    = params.duration		# ordinary time units
    density     = params.density       # how many particles per unit area
    #L           = 40.0       # box_X, sort of
    npart       = params.npart #  number of particles
    if norm == 'true':
        norm_str = 'vicsek'
    else:
        norm_str = 'flyingXY'
    TIME_BETWEEN_SAMPLES = 1.0
    maxruntime = params.maxruntime
      
    input_file_string = """\
// file names
LOG_FILE    = log$.dat		// the $ will be replaced by a time slice number
OUTPUT_DIR  = {6}        // the directory for all the *.dat files

// output format (options "csv", bin")
OUTPUT_FORMAT = bin

// physical parameters
NPART       = {0}	// number of particles
NITER       = 0	    // number of iterations
SIMTIME     = {1}   // replaces NITER
NOISE       = {2}	// noise strength
BOX_X       = {3}	// box X length
BOX_Y       = {4}	// box Y length
V0	        = 1.0	// particle speed
DT	        = {5}	// default time step
ALIGN       = 1.0	// alignment strength
DEFLECT     = 0.0	// deflection strength
RDEFLECTION = 0.0   // deflection radius
RSTAR       = 1.0   // interaction radius
SOFT_WALLS = false  // soft walls on Y
SOFT_WALLS_STRENGTH = 0.0   // strength factor for wall repulsion
OBSTACLEDEFLECTION = 0.0   // strength of deflection from obstacles
OBSTACLENUM = 0             // number of obstacles
PERIODIC_X = true
PERIODIC_Y = true
NORMALIZE_NEIGHBORS = {8}
NOISE_MODE = angular
TIME_BETWEEN_SAMPLES = {7}
HYSTER = false
MAX_RUNTIME_MINUTES = {9}


"""
    # shell command definitions
    simCmd = "MoreRollers.exe"

    # filenames
    suffix = norm_str + '_noise_' + str(noise) + '_dens_' + str(density) + '_npart_' + str(npart)
    direc = "output" + suffix
    inputFile = "input" + suffix
    error_file_name = 'errorlog' + suffix
    
    # set up derived physical quantities
    area = float(npart)/density
    box_x = math.sqrt(area)
    box_y = box_x
    
    print printtime() + "Starting simulation: D = {0} ; density = {1}".format(noise,density)
    with open(inputFile, "w") as text_file:
        text_file.write(input_file_string.format(npart,duration,noise,box_x,box_y,dt,direc,TIME_BETWEEN_SAMPLES, norm, maxruntime))
    
    # let's get this output working!
    proc = subprocess.Popen((simCmd, inputFile, '--resume'), stdout=subprocess.PIPE)
    # copy inputFile, why not
    f1=open(error_file_name, 'w')
    # take care of output
    while True:
      line = proc.stdout.readline()
      if line != '':
        prompt = '{0}_D={1}'.format(norm_str,noise).ljust(14)
        print "{0}: {1}".format(prompt, line.rstrip())
        print >>f1, line.rstrip()
        f1.flush()
      else:
        break
        
    f1.close()
    i = 0
    print direc
    for f in os.listdir(direc):
        if 'error' in f:
            print f
            i += 1
    print i
    if not direc.endswith('\\'):
        direc += '\\'
    shutil.move(inputFile, direc + inputFile)
    shutil.move(error_file_name, direc + error_file_name + '_{0:03d}'.format(i))
    print printtime() + "Finished single simulation, D = {0}".format(noise)
    return r"c'est tout"
