import os
import glob
from plots import *
import time
import sys
from utils import *
import subprocess
from multiprocessing import Process
from multiprocessing import Pool
import multiprocessing
import itertools
import argparse
import subprocess
from os.path import isfile,join

prefixOutput = 'output'

def make_single_movie(direc):
    movie_name = 'compass_' + direc.split('\\')[-1][len(prefixOutput):]
    if isfile(movie_name):
        print movie_name, ' already exists'
        return
    start_time = last_time_in_direc(direc) / 2.0
    make_movie(direc, start_time=0.0,frames=500, shortTime=True, # filename=None, 
        maxpart=None, nthframe=1, show_compass=True, filename=movie_name, show_clusters = True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make an .mp4 movie from a simulation output directory')
    parser.add_argument('directory', metavar='dir', type=str, nargs='?',default='.',
                       help='directory with simulation data')                   
    parser.add_argument('-p','--proc', metavar='processors', type=int,# action='store_const', nargs=1,
                   default=multiprocessing.cpu_count(), help='number of CPUs to use')
                       
    args = parser.parse_args()
    simDir = args.directory
    print printtime() + "Reading from directory: {0}".format(simDir)
    # directories
    output_dirs = []
    for subdir, dirs, files in os.walk(simDir):
        for dir in dirs:
            if dir.startswith(prefixOutput):
                output_dirs.append(simDir + '\\' + dir)
    #print output_dirs
    print printtime() + "Making {0} movies".format(len(output_dirs))
    # Make movies!
    p = Pool(args.proc)
    
    ret = p.map(make_single_movie, output_dirs[:5])
    #ret = map(make_single_movie, output_dirs)
    print printtime() + "Finished everything"



#os.system(matlabCmd)
