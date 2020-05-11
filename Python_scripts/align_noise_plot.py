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


#
def calc_last_pol(output_dir):
    sim = read_data_from_dir(output_dir, start_time = last_time_in_direc(output_dir)-2.0)
    if len(sim.data) == 0:
        print 'No data found'
        return 0
    frame = sim.data[-1]
    thetas = frame[:,2]
    pol = np.abs(np.mean(np.exp(1j * thetas)))
    return pol

#current_att = 6
if __name__ == '__main__':
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
    parser.add_argument('-y','--ylog', action='store_true',
                   help='Logarithmic y axis')
    parser.add_argument('-x','--xlog', action='store_true',
                   help='Logarithmic x axis')
    args = parser.parse_args()
    # if you reread you better recalculate means
    #if args.read or args.time:
    #    args.calc = True
    start_time = datetime.datetime.now()
    resultsDir = args.directory
    resultsDir = '.\\'
    procs = args.proc
    
    attributes = get_order_parameters()
    current_atts = []
    
    
    output_dirs = []
    inputs = []
    lastPols = []
    for output_dir in os.listdir(resultsDir):
        if 'output' in output_dir and os.path.isdir(join(resultsDir,output_dir)):
            try:
                input = read_input_from_dir(join(resultsDir,output_dir))
            except:
                continue
            #if getDensityInput(input) > 10.1:
            #    continue
            output_dirs.append(join(resultsDir,output_dir))
            inputs.append(input)
            lastPols.append(calc_last_pol(output_dir))
            
    x = np.array([input.noise for input in inputs])
    y = np.array([input.align for input in inputs])
    input = inputs[0]
    density = getDensityInput(input)
    print density
    print input.npart / np.square(input.box_x)
    txt = 'Phase diagram'
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    fig.canvas.set_window_title("{} Phase diagram".format(model))
    ax1 = fig.add_subplot(111)
    title = r'{}, {} Model, $\rho={:.1f}$'.format(txt, model, density)
    ax1.set_title(title, fontsize=16)#, fontweight='bold')
    ax1.set_xlabel('Noise', fontsize=14, fontweight='bold')
    ax1.set_ylabel('A', fontsize=14, fontweight='bold')
    
    #x = [en.input.noise for en in mProp]
    #y = [en.mValue for en in mProp]
    cmap = plt.get_cmap('Blues')
    colors = cmap(lastPols)
    scat = ax1.scatter(x,y, 150, facecolors=colors, edgecolors='k', picker = True, linewidths=1, alpha=1)
    ax1.set_xlim([0,np.max(x)*1.1])
    ax1.set_ylim([0,np.max(y)*1.1])
    
    bounds = np.linspace(0,1,100)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax2 = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=[0,0.5,1], format='%.1f')

    #cb.ax.set_yticklabels(['< -1', '0', '> 1'])# vertically oriented colorbar    
    
    def onpick_polar(event):
        ind = event.ind[0]
        time = None
        input = inputs[ind]
        output_dir = output_dirs[ind]
        last_time = last_time_in_direc(output_dir)
        if last_time == 0:
            print 'No data for this point'
            return
        if event.mouseevent.button==3:
            time = get_value_box(title='Drawing snapshot for D={}, maximum time = {:.2f}'.format(input.noise, last_time), 
                prompt='Choose time:')
            try:
                time = float(time)
                if time == 0:
                    return
            except:
                time = None
        ### test
        if len(current_atts) == 2:
            props = [calc_property(output_dir,att.pickle_name,att.calc_data_fun, read=args.read) for att in current_atts]
            fig = plt.figure()
            axes = [fig.add_subplot(211),fig.add_subplot(212)]
            for prop, att, ax in zip(props, current_atts, axes):
                plot_property_vs_time([prop], title=att.name, ylabel=att.name, save=True, 
                    legend=True, show=False, smooth=att.smooth, ax=ax, figure=fig)        
        else:
            for att in current_atts:
                prop = calc_property(output_dir,att.pickle_name,att.calc_data_fun, read=args.read)
                plot_property_vs_time([prop], title=att.name, ylabel=att.name, save=True, 
                    legend=True, show=False, smooth=att.smooth)
        ###
        if len(current_atts) == 0:
            snap = snapshot(output_dir, time=time)
            plot_snap(snap, save=True, show=False)
        plt.show()
    if (args.xlog):
        ax1.set_xscale('log')
    if (args.ylog):
        ax1.set_yscale('log')
    fig.canvas.mpl_connect('pick_event', onpick_polar)
    filename = '{}_vs_noise.png'.format(txt).replace(" ","")
    fig.savefig(filename, dpi=fig.dpi)
    
    # radio!
    axcolor = 'lightgoldenrodyellow'
    rax = plt.axes([0.02, 0.12, 0.3, 0.15], axisbg=axcolor)
    #labels = [Text(text=att.name, color='red') for att in attributes]
    labels = [att.name for att in attributes]
    radio = CheckButtons(rax, labels, [False]*len(labels))
    for label in radio.labels:
        label.set_fontsize(11)
    def hzfunc(label):
        global current_atts
        current_att = next(att for att in attributes if att.name == label)
        if current_att in current_atts:
            current_atts.remove(current_att)
        else:
            current_atts.append(current_att)
    radio.on_clicked(hzfunc)

    
    plt.show()