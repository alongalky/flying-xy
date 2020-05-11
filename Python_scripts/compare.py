from scipy import fftpack
from scipy.spatial import Voronoi, voronoi_plot_2d
from Tkinter import *
import threading
import psutil
import gc
import numpy as np
import sys
sys.path.insert(0, r'F:\Alon\Code\Scripts\Utility')
from plots import *
from utils import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing

import os.path
from os.path import isfile,join
import random
import argparse
from scipy.optimize import curve_fit
from collections import defaultdict
import pickle
from mpl_toolkits.mplot3d import Axes3D
import struct
import subprocess

buffer_in_mb = 100

## parsing parameters
parser = argparse.ArgumentParser(description='Calculate polarization for simulation data')
parser.add_argument('directory', metavar='dir', type=str, nargs='?',default=r'.',
               help='directory with simulation data')                   
parser.add_argument('-p','--proc', metavar='processors', type=int,# action='store_const', nargs=1,
               default=multiprocessing.cpu_count(), help='number of CPUs to use')
parser.add_argument('-r','--read', action='store_true', 
               help='force rereading sim data and rewriting .npy files')
parser.add_argument('-c','--calc', action='store_true', 
               help='force recalculating means')
parser.add_argument('-D','--density', action='store_true', 
               help='plot density variance')
parser.add_argument('-E','--energy', action='store_true',
               help='Show energy data')
parser.add_argument('-t','--time', action='store_true',
               help='Plot properties vs time rather than averages')
parser.add_argument('-B','--binder', action='store_true',
               help='Plot Binder cumulant')
parser.add_argument('-y','--ylog', action='store_true',
               help='Logarithmic y axis')
parser.add_argument('-x','--xlog', action='store_true',
               help='Logarithmic x axis')
parser.add_argument('-P','--polar', action='store_true',
               help='Show global polarization data')
args = parser.parse_args()

if __name__ == '__main__':

    # if you reread you better recalculate means
    if args.read or args.time:
        args.calc = True
    start_time = datetime.datetime.now()
    resultsDir = args.directory
    procs = args.proc
    
    output_dirs = []
    for output_dir in os.listdir(resultsDir):
        if 'output' in output_dir and os.path.isdir(join(resultsDir ,output_dir)):
            ###
            ###
            output_dirs.append(join(resultsDir, output_dir))
    p = Pool(procs)
    #
    #output_dirs = output_dirs[:4]
    #
    input_list = []
    Pols = []
    
    if args.energy:
        if args.calc or not isfile('mEnergies.pickle'):
            Energies = p.map(calc_Energy, output_dirs)
            #Energies = map(calc_Energy, output_dirs)
            print printtime() + 'Calculating Energies'
            mEnergies = p.map(calc_mEnergies, Energies)
            pickle_save('mEnergies.pickle', mEnergies)
        else:
            mEnergies = pickle_load('mEnergies.pickle')
        if not args.time:
            plot_mProperty(mEnergies, title='Energy')
        else:
            plot_property_vs_time(Energies, title='Energy', ylabel='Energy', show=True, smooth=False)
    mPols = []
    if args.polar:
        if args.calc or not isfile('mPols.pickle'):
            for direc in output_dirs:
                mPols.append(calc_pol_var(direc))
                print printtime() + 'Finished {} / {} polarization calculations'.format(len(mPols), len(output_dirs))
            #Pols = map(calc_Pol, output_dirs)
            #Pols = p.map(calc_Pol, output_dirs)
            #print printtime() + 'Calculating Polarization'
            #mPols = map(calc_mPols, Pols)
            #mPols = p.map(calc_mPols, Pols)
            pickle_save('mPols.pickle', mPols)
        else:
            mPols = pickle_load('mPols.pickle')
    if args.density:
        if False:
            for direc in output_dirs[3:]:
                #calc_pdf_fourier(direc)
                calc_pdf_windows(direc)
                sys.exit()
        #pdfs = []
        vars = []
        #for direc in output_dirs[:8]:
        for direc in output_dirs[:3]:
            #pdfs.append(calc_pdf(direc))
            #print 'Finished {} / {} pdf calculations'.format(len(pdfs), len(output_dirs))
            vars.append(calc_pdf_vars(output_dir=direc))
        #plot_single_pdf(pdfs[0])
        plot_pdfs(vars=vars)
    hourminsec = hours_minutes_seconds(datetime.datetime.now() - start_time)
    print printtime() + "Total time processing data: {}:{}:{}".format(hourminsec.hours, hourminsec.minutes, hourminsec.seconds)
    
    if args.polar:
        plot_mean_pol([x for x in mPols if x.input.npart>0], show=False)
        #plot_mean_pol([x for x in mPols if x.input.npart>9000], show=False)
        if args.time:
            if len(Pols) == 0:
                # todo: make parallel
                for direc in output_dirs:
                    try:
                        Pols.append(calc_Pol(direc))
                    except:
                        continue
            plot_property_vs_time(Pols, title='Pol. vs time', ylabel=r'<\Pi>', show=True)
            #plot_pol_time(Pols, show=False)
    if args.binder:
        plot_binder(Binders, show=False)
    plt.show()