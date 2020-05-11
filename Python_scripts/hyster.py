from scipy import fftpack
from scipy.spatial import Voronoi, voronoi_plot_2d
from Tkinter import *
import threading
import psutil
import gc
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import sys
sys.path.insert(0, r'F:\Alon\Code\Scripts\Utility')
from plots import *
from utils import *
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
parser.add_argument('-P','--polar', action='store_true',
               help='Show global polarization data')
parser.add_argument('-D','--density', action='store_true',
               help='Show local density analysis')
parser.add_argument('-s','--smooth', action='store_true',
               help='Smooth results')
parser.add_argument('-y','--ylog', action='store_true',
               help='Logarithmic y axis')
parser.add_argument('-x','--xlog', action='store_true',
               help='Logarithmic x axis')
args = parser.parse_args()
# if you reread you better recalculate means
if (args.read):
    args.calc = True

# function definitions
class mean_propInfo():
  def __init__(self, mValue, simlen, input, output_dir):
     self.mValue = mValue
     self.simlen = simlen
     self.input = input
     self.output_dir = output_dir
class propInfo():
  def __init__(self, values, times, input, prop_name, output_dir):
     self.values = values
     self.times = times
     self.input = input
     self.prop_name = prop_name
     self.output_dir = output_dir
def get_value_box(title='Give Value', prompt='val'):
    root = Tk()
    root.geometry("300x280+300+300")
    
    w = Label(root, text=prompt)
    w.pack()
    
    v = StringVar()
    e = Entry(root, textvariable=v)
    e.pack()

    e.delete(0, END)
    v.set('0.0')
    b = Button(root, text="Plot", width=10, command=lambda : root.destroy())
    b.pack()
    
    e.bind("<Return>", lambda e: root.destroy())
    e.focus_set()
    e.selection_range(0, END)
    root.mainloop() 

    return v.get()

    box_x = input.box_x
    res = 50
    resolution = res * 1j
    X, Y = np.mgrid[0:box_x:resolution, 0:box_x:resolution]
    grid = np.vstack([X.ravel(), Y.ravel()])
    values = []
    for i,frame in enumerate(data):
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        thetas = frame[:,2]
        positions = np.vstack([x_pos, y_pos]).T
        test = []
        sigma = 5.0
        for point in grid.T:
            test.append(eval_pol_density(positions, thetas, point, input.box_x, sigma = sigma))
        
        values.append(np.mean(test))
        print printtime() + 'D={:.2f} Finished {} / {}'.format(input.noise, i+1, len(data))
    return values    
def plot_densvar_vs_noise_hyster(data, show=True):
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    model = 'Vicsek' if data[0].input.neighbor_norm == True else 'Flying XY'
    txt = 'Density Variance'
    title = '{}, {} Model'.format(txt, model)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('noise', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'$<\Delta \rho^2>$', fontsize=14, fontweight='bold')
    for sim in data:
        input = sim.input
        noises =  noiseFunc(sim.input, sim.times)
        y = savitzky_golay(sim.values, 1001, 4)
        ax.plot(noises, y, linewidth=0.5, linestyle='-',# markersize=4, marker='o',
            label='{} particles, t={:.1f}'.format(input.npart, input.simtime), picker=True)
    #last_time = np.min([sim.times[-1] for sim in data])
    #ax.set_xlim([0,last_time])
    plt.legend(loc='center right',prop={'size':11})
    filename = 'dens_variations_vs_noise_hyster.png'
    fig.savefig(filename, dpi=fig.dpi)    
    def onpick4(event):
        ind = event.ind
        thisline = event.artist
        line_index = ax.lines.index(thisline)
        sim = data[line_index]
        time = np.take(sim.times, ind)[0]
        snap = snapshot(sim.output_dir, time=time)
        plot_snap(snap, save=True)

    fig.canvas.mpl_connect('pick_event', onpick4)  
    if show:
        plt.show() 
      
def plot_pdfpol_vs_noise_hyster(data, show=True):
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    model = 'Vicsek' if data[0].input.neighbor_norm == True else 'Flying XY'
    txt = 'Local Polarization'
    title = '{}, {} Model'.format(txt, model)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('noise', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'$<|\rho (\vec{r}, t) \pi (\vec{r}, t)|>_{\vec{r}, t}$', fontsize=14, fontweight='bold')
    for sim in data:
        input = sim.input
        noises =  noiseFunc(sim.input, sim.times)
        y = savitzky_golay(sim.values, 1001, 4)
        ax.plot(noises, y, linewidth=0.5, linestyle='-',# markersize=4, marker='o',
            label='{} particles, t={:.1f}'.format(input.npart, input.simtime), picker=True)
    #last_time = np.min([sim.times[-1] for sim in data])
    #ax.set_xlim([0,last_time])
    plt.legend(loc='center right',prop={'size':11})
    filename = 'mean_local_pol_vs_noise_hyster.png'
    fig.savefig(filename, dpi=fig.dpi)    
    def onpick4(event):
        ind = event.ind
        thisline = event.artist
        line_index = ax.lines.index(thisline)
        sim = data[line_index]
        time = np.take(sim.times, ind)[0]
        snap = snapshot(sim.output_dir, time=time)
        plot_snap(snap, save=True)

    fig.canvas.mpl_connect('pick_event', onpick4)    
    if show:
        plt.show()

if __name__ == '__main__':
    start_time = datetime.datetime.now()
    resultsDir = args.directory
    procs = args.proc
    
    output_dirs = []
    for output_dir in os.listdir(resultsDir):
        if 'output' in output_dir and os.path.isdir(join(resultsDir,output_dir)):
            output_dirs.append(join(resultsDir,output_dir))
            
            
    p = Pool(procs)
    
    ## censurship
    #output_dirs = output_dirs[:2]
    #print output_dirs
    ##

    input_list = []
    Pols = []
    hourminsec = hours_minutes_seconds(datetime.datetime.now() - start_time)
    print printtime() + "Total time processing data: {}:{}:{}".format(hourminsec.hours, hourminsec.minutes, hourminsec.seconds)
    
    
    
    if args.density:
        Dens = []
        for direc in output_dirs[:]:
            input = read_input_from_dir(direc)
            #if input.simtime < 15000 or input.simtime > 25000:
            #    continue
            if input.simtime < 9000:
                continue
            Dens.append(calc_property(direc,'densvar',calc_vars, read=args.read))
        Dens = sorted(Dens, key=lambda x: x.input.simtime)
        for dens in Dens:
            dens.input = read_input_from_dir(dens.output_dir)
        #Dens = [x for x in Dens if x.input.simtime > 15000 and x.input.simtime < 25000]
        plot_hyster(Dens, show=False, smooth=args.smooth, ylabel=r'$<(\delta \rho)^2>$')
    elif args.polar:
        if len(Pols) == 0:
            # todo: make parallel
            for direc in output_dirs:
                Pols.append(calc_property(direc,'Pols',Pol_data_calc, read=args.read))
        Pols = [pol for pol in Pols if pol.input.simtime > 9000]
        Pols = [pol for pol in Pols if pol.input.hyster]
        Pols = sorted(Pols, key=lambda x: x.input.simtime)
        if Pols[0].input.hyster:
            plot_hyster(Pols, show=False, smooth=args.smooth)
        else:
            print 'ERROR! meant for hysterisis'
         
    plt.show()