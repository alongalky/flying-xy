from matplotlib.patches import Ellipse
import sys
sys.path.insert(0, r'F:\Alon\Code\Scripts\Utility')
from matplotlib.text import Text
from matplotlib.widgets import RadioButtons,CheckButtons
from numpy.random import randn
from matplotlib import cm
from scipy import fftpack
from scipy.spatial import Voronoi, voronoi_plot_2d
from Tkinter import *
import threading
import psutil
import gc
import numpy as np
from plots import *
from utils import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from multiprocessing import Pool
import multiprocessing
import os.path
from os.path import isfile,join,isdir
import random
import argparse
from scipy.optimize import curve_fit
from collections import defaultdict
import pickle
from mpl_toolkits.mplot3d import Axes3D
import struct
import itertools
import subprocess

read = False

def plot_transition_curves(ax, filename):
    if not isfile(filename):
        print 'Warning: cannot draw transition curves, no data'
        return
    points = pickle_load(filename)
    methods = set([point.method for point in points])
    for method in methods:
        data = [point for point in points if point.method==method]
        data = sorted(data, key=lambda x: x.align)
        x_dat = [x.align for x in data]
        y_dat = [x.density for x in data]
        ax.plot(x_dat, y_dat, label = transitionMethods.names[method])
        #ax.plot(x_dat, 1/np.array(y_dat), label = transitionMethods.names[method] + ' 1/y' )
    ax.legend(loc='upper left')
    
def get_trans():
    filez=[r'align_density_trans.pickle',
        r'noise_density_trans.pickle']
    trans = []
    for file in filez:
        if not isfile(file):
            return [0.739369762117]
            continue
        data = pickle_load(file)
        data = [x for x in data if x.align == align and x.noise == noise]
        trans += [x.density for x in data]
    return trans
parser = argparse.ArgumentParser(description='Calculate polarization for simulation data')
parser.add_argument('directory', metavar='dir', type=str, nargs='?',default=r'data',
               help='directory with simulation data')                   
parser.add_argument('-p','--proc', metavar='processors', type=int,# action='store_const', nargs=1,
               default=multiprocessing.cpu_count(), help='number of CPUs to use')
parser.add_argument('-n','--noise', metavar='alignment', type=float,# action='store_const', nargs=1,
               default=2.0, help='which noise value to use (default D=2.0)')
parser.add_argument('-a','--align', metavar='alignment', type=float,# action='store_const', nargs=1,
               default=1.0, help='which A value to use (default A=1.0)')
parser.add_argument('-r','--read', action='store_true', 
               help='force rereading sim data and rewriting .npy files')
parser.add_argument('-c','--calc', action='store_true', 
               help='force recalculating means')
parser.add_argument('-N','--npart', metavar='noise value', type=float,# action='store_const', nargs=1,
               default=None, help='how many particles?')
               
parser.add_argument('-y','--ylog', action='store_true',
               help='Logarithmic y axis')
parser.add_argument('-x','--xlog', action='store_true',
               help='Logarithmic x axis')
parser.add_argument('-m','--minimal', action='store_true',
               help='Show only simple features')
args = parser.parse_args()
resultsDir = args.directory
procs = args.proc
noise = args.noise
align = args.align
read = args.read
def calc_mGyration(output_dir):
    return calc_mProp(output_dir, 'mGyration', 'mean_cluster_gyration_radius',calc_cluster_gyration_radius)
def calc_mClusterSize(output_dir):
    return calc_mProp(output_dir, 'mClusterSize', 'mean_cluster_size',calc_mean_cluster_size)
if __name__ == '__main__':
    
    ## parsing parameters
    # if you reread you better recalculate means
    #if args.read or args.time:
    #    args.calc = True
    start_time = datetime.datetime.now()
    
    
    output_dirs = []
    for output_dir in os.listdir(resultsDir)[:]:
        full_dir = join(resultsDir,output_dir)
        if not isdir(full_dir):
            continue
        try:
            input = read_input_from_dir(full_dir)
        except:
            continue
        if input.noise != noise or input.align != align:    
            continue
        if args.npart != None and input.npart != args.npart:
            continue
        output_dirs.append(full_dir)
    print 'starting work on {} directories'.format(len(output_dirs))
    p = Pool(procs)
    output_dirs = sorted(output_dirs, key=lambda x: getDensityInput(read_input_from_dir(x)))
    #Binders = [calc_binder(output_dir,read=args.read) for output_dir in output_dirs]
    #plot_mProperty_vs_dens(Binders, title=False, ylabel='$G$', vlines=get_trans(), showtitle=False)
    for att in get_order_parameters():
        #if att.pickle_name != 'abs_momentum':
            #    continue
        props = []
        for output_dir in output_dirs:
            props.append(calc_mProp_orderparam(output_dir, att, read=args.read, calc=args.calc))
        plot_mProperty_vs_dens(props, title=False, ylabel=att.ylabel, vlines=get_trans(), showtitle=False)
    # plot Binders
    means = [calc_mProp_orderparam(output_dir, orderParamter('Polarization', 'Pols', Pol_data_calc, mean_pickle_name='mPol'), read=args.read, calc=args.calc)
                for output_dir in output_dirs]
    plot_mProperty_vs_dens(means, title='Polarization Fluctuations',
        ylabel = '$<(\delta \Pi)^2>$',
        order_param = orderParamter('Polarization', 'Pols', Pol_data_calc, mean_pickle_name='mPol'),
        vlines=get_trans(), variance = True, clusters=False, show=True, compass=False)
