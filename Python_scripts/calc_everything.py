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
from os.path import isfile,join
import random
import argparse
from scipy.optimize import curve_fit
from collections import defaultdict
import pickle
from mpl_toolkits.mplot3d import Axes3D
import struct
import subprocess


# function definitions
class mean_propInfo():
  def __init__(self, mValue, simlen, input, output_dir):
     self.mValue = mValue
     self.simlen = simlen
     self.input = input
     self.output_dir = output_dir
class calcMissionInfo():
  def __init__(self, output_dir,att, read=False, resume=True):
     self.output_dir = output_dir
     self.att = att
     self.read = read
     self.resume = resume
def run_calc(mission):
    #print 'calculating ', mission.att.pickle_name, ' for ', mission.output_dir
    return calc_property(mission.output_dir,mission.att.pickle_name,mission.att.calc_data_fun, read=mission.read)
class orderParamter():
  def __init__(self, name, pickle_name, calc_data_fun):
     self.name = name
     self.pickle_name = pickle_name
     self.calc_data_fun = calc_data_fun

#current_att = 6
if __name__ == '__main__':
    
    global args
    ## parsing parameters
    parser = argparse.ArgumentParser(description='Calculate polarization for simulation data')
    parser.add_argument('directory', metavar='dir', type=str, nargs='?',default=r'.',
                   help='directory with simulation data')                   
    parser.add_argument('-p','--proc', metavar='processors', type=int,# action='store_const', nargs=1,
                   default=multiprocessing.cpu_count(), help='number of CPUs to use')
    parser.add_argument('-r','--read', action='store_true', 
                   help='force rereading sim data and rewriting .npy files')
    parser.add_argument('-P','--pols', action='store_true', 
                   help='calculate polarization')
    parser.add_argument('-D','--density', action='store_true', 
                   help='calculate density variance')
    parser.add_argument('-W','--momentum', action='store_true', 
                   help='calculate mean momentum')
    parser.add_argument('-E','--energy', action='store_true',
                   help='calculate energy data')                   
    parser.add_argument('-C','--cluster', action='store_true',
                   help='calculate clustering data')                   
    parser.add_argument('-n','--noise', metavar='alignment', type=float,# action='store_const', nargs=1,
               default=None, help='which noise value to use (default D=2.0)')
    parser.add_argument('-a','--align', metavar='alignment', type=float,# action='store_const', nargs=1,
               default=None, help='which A value to use (default A=1.0)')
    parser.add_argument('-m','--minpart', metavar='minimumpart', type=float,# action='store_const', nargs=1,
               default=None, help='minimum number of particles to take seriously')
    parser.add_argument('-M','--maxpart', metavar='maximumpart', type=float,# action='store_const', nargs=1,
               default=None, help='maximum number of particles to take seriously')
    args = parser.parse_args()

    start_time = datetime.datetime.now()
    resultsDir = args.directory
    procs = args.proc
    
    attributes = get_order_parameters()
    if args.pols:
        attributes = [orderParamter('Polarization', 'Pols', Pol_data_calc)]
    if args.density:
        attributes = [orderParamter('Density variance', 'densvar', calc_vars)]
    if args.cluster:
        attributes = [orderParamter('Clustering Stats', 'cluster_stats', calc_cluster_stats),
        orderParamter('Cluster Fraction (neighbors)', 'neighbor_pdf', calc_neighbors_pdf)]
    #print attributes
    """
    if args.energy:
        attributes.append(orderParamter('Energy', 'Energy', Energy_data_calc))
    if args.momentum:
        attributes.append(orderParamter('Mean local momentum (W)', 'momentum', calc_W))
    if args.clus:
        attributes.append(orderParamter('Mean local momentum (W)', 'momentum', calc_W))"""
    missions = [] 
    
    output_dirs = []
    inputs = []
    lastPols = []
    
    for output_dir in os.listdir(resultsDir):
        if 'output' in output_dir and os.path.isdir(join(resultsDir,output_dir)):
            try:
                input = read_input_from_dir(join(resultsDir,output_dir))
                if args.align != None and input.align != args.align:
                    continue
                if args.noise != None and input.noise != args.noise:
                    continue
                if args.minpart != None:
                    if args.minpart > input.npart: 
                        continue
                if args.maxpart != None:
                    if args.maxpart < input.npart:
                        continue
            except:
                continue
            output_dirs.append(join(resultsDir,output_dir))
    
    for att in attributes:
        for output_dir in output_dirs[:]:
            missions.append(calcMissionInfo(output_dir, att, read=args.read))
    print 'starting {} calculations on {} dirs'.format(len(missions), len(output_dirs))
    p = Pool(procs)
    p.map(run_calc,missions)
    #map(run_calc,missions[:])
