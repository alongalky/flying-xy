from matplotlib.patches import Ellipse
import scipy
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
from os.path import isfile,join, isdir
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
parser.add_argument('-a','--align', metavar='alignment', type=float,# action='store_const', nargs=1,
               default=1.0, help='which A value to use (default A=1.0)')
parser.add_argument('-r','--read', action='store_true', 
               help='force rereading sim data and rewriting .npy files')
parser.add_argument('-c','--calc', action='store_true', 
               help='force recalculating means')
parser.add_argument('-y','--ylog', action='store_true',
               help='Logarithmic y axis')
parser.add_argument('-x','--xlog', action='store_true',
               help='Logarithmic x axis')
parser.add_argument('-v','--vicsek', action='store_true',
               help='Vicsek!')
parser.add_argument('-m','--minimal', action='store_true',
               help='Show only simple features')
args = parser.parse_args()

log_x = args.xlog
log_y = args.ylog


def plot_transition_curves(filename, align=None, txt='', cmap=plt.get_cmap('brg'), vicsek=False):
    points = []
    align_points = []
    noise_points = []
    #for filename in files:
    if isinstance(filename, list):
        for f in filename:
            if not isfile(f):
                print 'Warning: cannot draw transition curves, no data in ', filename
            else:
                if 'align' in f:
                    align_points += pickle_load(f)
                else:
                    noise_points += pickle_load(f)
                points += pickle_load(f)
    else:
        if not isfile(filename):
            print 'Warning: cannot draw transition curves, no data in ', filename
            return
        points += pickle_load(filename)

    txt = 'Transition points'
    model = 'Vicsek' if vicsek == True else 'Flying XY'
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    fig.canvas.set_window_title("{} {}".format(model, txt))
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel(r'$D/A$', fontsize=22, fontweight='bold', labelpad=-3)
    ax1.set_ylabel(r'$\rho$', fontsize=22, fontweight='bold', labelpad=0)
    
    #methods = set([point.method for point in points])
    colors = [cmap(i) for i in np.linspace(0,1,3)]
    lines = []
    all_x = []
    all_y = []
    
    # inset!
    #inset = vicsek
    inset = False
    if vicsek and inset:
        #axins = zoomed_inset_axes(ax1, 1.2, loc=3, bbox_to_anchor=(0.1, 0.3), 
        #             bbox_transform=ax1.figure.transFigure)
        axins = inset_axes(ax1, width=3, height=3, loc=3, bbox_to_anchor=(0.3, 0.3), 
                     bbox_transform=ax1.figure.transFigure)
        axins.axis([0.009, 0.4, 0.0061, 2.0])
        axins.get_xaxis().set_visible(False)
        axins.get_yaxis().set_visible(False)
        mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        #axins.xaxis.tick_top()
    
    labelz = [r'$A-\rho$', r'$D-\rho$']
    all_points = [align_points, noise_points]
    for i, points in enumerate(all_points):
        #data = [point for point in points if point.method==method]
        print len(points), labelz[i]
        points = sorted(points, key=lambda x: x.noise)
        points = sorted(points, key=lambda x: x.align)
        #data = sorted(data, key=lambda x: x.align)
        y_dat = [x.density for x in points]
        if vicsek:
            x_dat =np.array([x.noise for x in points]) / np.array([x.align for x in points])
            #x_dat=np.square(x_dat)
            #if inset:
            #    axins.plot(np.sqrt(x_dat), y_dat, label = transitionMethods.names[method]+txt, color = 'g', picker=True,
            #        linewidth=0, marker='o')
        else:
            x_dat = np.array([x.noise for x in points]) / np.array([x.align for x in points])
        all_x += list(x_dat)
        all_y += list(y_dat)
        ax1.plot(x_dat, y_dat, label = labelz[i], color = colors[i], picker=True,
            linewidth=0, marker='o')
        if inset:
            axins.plot(x_dat, y_dat, label = transitionMethods.names[method]+txt, color = colors[i], picker=True,
                linewidth=0, marker='o')
        lines.append(points)
    ax1.legend(loc='best')
    if vicsek:
        all_x = np.array(all_x)
        all_y = np.array(all_y)
        max_dens = 1.0
        #print all_y
        kosher_indices = np.logical_and( all_y <= max_dens, all_x >= 0.0)
        #print kosher_indices
        x_dat = all_x[kosher_indices]
        y_dat = all_y[kosher_indices]
        
        def parabola(x, a,b):
            return a*(x**2)
        def powerlaw(x,a,b,c):
            return a*(x**b)
        func = powerlaw
        popt, pcov = curve_fit(func, x_dat, y_dat)
        #print 'vals  = ',
        #print popt
        x = np.linspace(np.min(x_dat), 1.0,40)#np.max(x_dat), 30)
        yfit = func(x, *popt)
        #ax1.plot(x, yfit, color = 'r')
        #ax1.plot(x, 80*x**2+0.2, color = 'r')
        #ax1.plot(x, 1000*(x**6)+0.02, color = 'g')
        if inset:
            #axins.set_xscale('log')
            #axins.set_yscale('log')

            axins.plot(x, pol(x), color = 'r')
    
    #print 'regression data for line is:'
    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(all_x, all_y)
    fit, cov = np.polyfit(all_x, all_y, 1, cov=True)
    pol = np.poly1d(fit)
   # print 'rho = {} * (D/A) + {}'.format(fit[0], fit[1])
   # print 'errors = ', np.sqrt(np.diag(cov))
   # print 'R^2 = ', r_value**2
    x = np.linspace(np.min(x_dat), np.max(x_dat), 20)
    if not vicsek:
        ax1.plot(x, fit[0]*x + fit[1])
    def onpick(event):
        ind = event.ind
        thisline = event.artist
        line_index = ax1.lines.index(thisline)
        #x = np.take(thisline.get_xdata(), ind)[0]
        point = lines[line_index][ind]
        print point.noise, point.align, point.density
        method = point.method
        dirlist = []
        for output_dir in os.listdir(resultsDir):
            full_dir = join(resultsDir,output_dir)
            if not isdir(full_dir) or not 'output' in output_dir:
                continue
            try:
                input = read_input_from_dir(full_dir)
            except:
                print 'error in ', output_dir
                sys.exit()
            #if point.const_align == True:
            if point.align == input.align and point.noise == input.noise:
                dirlist.append(full_dir)
        plot_transition(point, dirlist)
        #print x,ind, thisline
    if (log_x):
        ax1.set_xscale('log')
    if (log_y):
        ax1.set_yscale('log')
    
    filename = '{}_surface_conjecture.png'.format(model.replace(' ',''))
    fig.savefig(filename, dpi=fig.dpi)    
    
    fig.canvas.mpl_connect('pick_event', onpick)


    return lines
def plot_transition(point, dirlist, show=True):
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    
    pols = map(calc_pol_var, dirlist)
    if len(pols) == 0:
        print 'Warning: empty pols in dirs: ' , dirlist
        return
    noise = pols[0].input.noise
    align = pols[0].input.align
    
    ax1.set_title(r'A = {:.2f}, D={:.2f}'.format(align, noise), fontsize=14, fontweight='bold')
    ax1.set_xlabel('Density', fontsize=14, fontweight='bold')
    
    if point.method == transitionMethods.polInflection:
        ax1.set_ylabel(r'$<\Pi>$', fontsize=14, fontweight='bold')
        pols = sorted(pols, key=lambda x: getDensityInput(x.input))
        x_dat = np.array([getDensityInput(x.input) for x in pols])
        y_dat = np.array([x.mValue for x in pols])

        ax1.set_xlim([np.min(x_dat),np.max(x_dat)])
        ax1.set_ylim([np.min(y_dat),np.max(y_dat)])
        line = ax1.plot(x_dat,y_dat, marker='o')
        #c = CubicFit(ax1, x_dat, y_dat)
    else:
        ax1.set_ylabel(r'$<(\delta\Pi)^2>$', fontsize=14, fontweight='bold')
        pols = sorted(pols, key=lambda x: getDensityInput(x.input))
        x_dat = np.array([getDensityInput(x.input) for x in pols])
        y_dat = np.array([x.var for x in pols])

        ax1.set_xlim([np.min(x_dat),np.max(x_dat)])
        ax1.set_ylim([np.min(y_dat),np.max(y_dat)])
        line = ax1.plot(x_dat,y_dat, marker='o')
    v_line = ax1.vlines(point.density, 0, np.max(y_dat), linestyle='--')
    if show:
        plt.show()
if __name__ == '__main__':
    # if you reread you better recalculate means
    #if args.read or args.time:
    #    args.calc = True
    start_time = datetime.datetime.now()
    resultsDir = args.directory
    procs = args.proc
    
    attributes = get_order_parameters()
    current_atts = []
    
    lines = plot_transition_curves(
        ['noise_density_trans.pickle'  , 'density_noise_trans.pickle' ,'align_density_trans.pickle', 'density_align_trans.pickle'], vicsek=args.vicsek)
    #lines = plot_transition_curves(
    #    ['noise_density_trans.pickle'  , 'density_noise_trans.pickle'], vicsek=args.vicsek)
    ## plot transition curve!
    plt.show()