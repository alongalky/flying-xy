import sys
sys.path.insert(0, r'F:\Alon\Code\Scripts\Utility')
import numpy as np
from utils import *
from plots import *
# stopped working on this in the middle, might be useful one day


import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import sys
import re
import random
import argparse

def drawsnap(output_dir):
    try:
        input = read_input_from_dir(output_dir)
        #
        #if input.npart < 8000:
        #    return
        #
        snap = snapshot(output_dir, compass = False)
    except:
        return
    plot_snap(snap, show=False, save=True, arrows=True, color='black', 
        backgroundcolor='white', clusters = clusters)
    plt.close("all")
    # plot_snap(snap, show=True, save=True, arrows=True, color='black', 
        # backgroundcolor='white', clusters = clusters)
    return True

parser = argparse.ArgumentParser(description='Draw final snapshots')
parser.add_argument('directory', metavar='dir', type=str, nargs='?',default='.',
               help='directory with simulation data')                   
parser.add_argument('-p','--proc', metavar='processors', type=int,# action='store_const', nargs=1,
               default=multiprocessing.cpu_count(), help='number of CPUs to use')
parser.add_argument('-c','--clusters', action='store_true', 
            help='Draw clusters in the snaps')


args = parser.parse_args()
resultsDir = args.directory
clusters = args.clusters
procs = args.proc
    
if __name__ == '__main__':

    dirs = []
    for f in os.listdir(resultsDir):
        if 'output' in f and not f.endswith('.png'):
            dirs.append(join(resultsDir,f))
    #print dirs
    p = Pool(procs)

    map(drawsnap, dirs[:])
    #p.map(drawsnap, ['outputnoise_0.01_dens_0.05'])
 