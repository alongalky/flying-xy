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
import subprocess


def calc_last_pol(output_dir):
    sim = read_data_from_dir(output_dir, start_time = last_time_in_direc(output_dir)-2.0)
    if len(sim.data) == 0:
        print 'No data found'
        return 0
    frame = sim.data[-1]
    thetas = frame[:,2]
    pol = np.abs(np.mean(np.exp(1j * thetas)))
    return pol
    
class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
        self.ttime = -1.0

    def __call__(self, event):
        #print 'click', event
        if event.inaxes!=self.line.axes: return
        #print 'time = ', event.xdata
        self.ttime = event.xdata
        plt.close()
    
    def get_transient_time(self):
        return self.ttime

        
def choose_transient_time(prop):
    input = prop.input
    txt = prop.prop_name
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    title = '{}, {} Model'.format(txt, model)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('time', fontsize=14, fontweight='bold')
    ax.set_ylabel(txt, fontsize=14, fontweight='bold')
    line, = ax.plot(prop.times, prop.values, label='D = {:.2f}'.format(input.noise), picker = True)
    plt.legend(loc='best')
    #line, = ax.plot([0], [0])  # empty line
    linebuilder = LineBuilder(line)

    plt.show()
    return linebuilder.get_transient_time()

def find_transition_density_by_pol(pols):
    pols = sorted(pols, key=lambda x: getDensityInput(x.input))
    x_dat = np.array([getDensityInput(x.input) for x in pols])
    y_dat = np.array([x.var for x in pols])
    avg_dens = np.dot(x_dat,y_dat)/np.sum(y_dat)
    return avg_dens
def calc_transition_curve(pols):
    noises = sorted(set([pol.input.noise for pol in pols]))
    x = []
    y = []
    for noise in noises:
        data = [pol for pol in pols if pol.input.noise == noise]
        if len(data) < 5:
            continue
        x.append(noise)
        y.append(find_transition_density_by_pol(data))
    return (x,y)
def plot_transition_curves(ax, filename, sortbyrho=False):
    if not isfile(filename):
        print 'Warning: cannot draw transition curves, no data'
        return
    points = pickle_load(filename)
    methods = set([point.method for point in points])
    for method in methods:
        data = [point for point in points if point.method==method]
        data = sorted(data, key=lambda x: x.align)
        if sortbyrho:
            data = sorted(data, key=lambda x: x.density)
        x_dat = [x.align for x in data]
        y_dat = [x.density for x in data]
        ax.plot(x_dat, y_dat, label = transitionMethods.names[method], linewidth=5, color='r')
        #ax.plot(x_dat, 1/np.array(y_dat), label = transitionMethods.names[method] + ' 1/y' )
    #ax.legend(loc='upper right')
if __name__ == '__main__':
    ## parsing parameters
    parser = argparse.ArgumentParser(description='Calculate polarization for simulation data')
    parser.add_argument('directory', metavar='dir', type=str, nargs='?',default=r'.',
                   help='directory with simulation data')                   
    parser.add_argument('-p','--proc', metavar='processors', type=int,# action='store_const', nargs=1,
                   default=multiprocessing.cpu_count(), help='number of CPUs to use')
    parser.add_argument('-n','--noise', metavar='noise value', type=float,# action='store_const', nargs=1,
                   default=2.0, help='which noise to plot')
    parser.add_argument('-r','--read', action='store_true', 
                   help='force rereading sim data and rewriting .npy files')
    parser.add_argument('-c','--calc', action='store_true', 
                   help='force recalculating means')
    parser.add_argument('-y','--ylog', action='store_true',
                   help='Logarithmic y axis')
    parser.add_argument('-x','--xlog', action='store_true',
                   help='Logarithmic x axis')
    parser.add_argument('-m','--minimal', action='store_true',
                   help='Show only simple features')
    args = parser.parse_args()
    # if you reread you better recalculate means
    #if args.read or args.time:
    #    args.calc = True
    start_time = datetime.datetime.now()
    resultsDir = args.directory
    procs = args.proc
    noise = args.noise
    
    attributes = get_order_parameters()
    current_atts = []
    
    
    output_dirs = []
    inputs = []
    lastPols = []

    for output_dir in os.listdir(resultsDir):
        full_dir = join(resultsDir,output_dir)
        if not isdir(full_dir):
            continue
        try:
            input = read_input_from_dir(full_dir)
            if input.noise != noise:    
                continue
        except:
            continue
        shade = calc_pol_var(full_dir, calc=args.calc, read=args.read)
        if shade == None:
            continue
        output_dirs.append(full_dir)
        inputs.append(input)
        lastPols.append(shade)
    
    x = np.array([input.align for input in inputs])
    y = np.array([getDensityInput(input) for input in inputs])
    input = inputs[0]
    txt = 'Phase diagram'
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    fig.canvas.set_window_title("{} Phase diagram".format(model))
    ax1 = fig.add_subplot(111)
    title = '{} {}, D={}'.format(txt, model, input.noise)
    if not args.minimal:
        ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel(r'$A$', fontsize=22, fontweight='bold', labelpad=-2)
    ax1.set_ylabel(r'$\rho$', fontsize=22, fontweight='bold')
    

    cmap = plt.get_cmap('Blues')
    bounds = np.linspace(0,1,100)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax2 = fig.add_axes([0.91, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=[0,0.5,1], format='%.1f')


    ##
    normalization_factor = 1#np.max([sim.mValue for sim in lastPols])
    colors = cmap(np.array([sim.mValue for sim in lastPols])/normalization_factor)
    #normalization_factor = np.max([sim.var for sim in lastPols])
    #colors = cmap(np.array([sim.var for sim in lastPols])/normalization_factor)

    scat = ax1.scatter(x,y, 150, facecolors=colors, edgecolors='k', picker = True, linewidths=1, alpha=1)
    ax1.set_xlim([0,np.max(x)*1.1])
    ax1.set_ylim([0,np.max(y)*1.1])
    

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
    #el = Ellipse((2, -1), 0.5, 0.5)
    #ax1.add_patch(el)
    def onpick_hover(event):
        collisionFound = False
        if event.xdata != None and event.ydata != None: # mouse is inside the axes
            for kid in ax1.get_children():
                if isinstance(kid, mpl.text.Annotation):
                    kid.remove()
            for i in xrange(len(x)):
                radius = 0.035
                xmin, xmax = ax1.get_xlim()
                ymin, ymax = ax1.get_ylim()
                radius_x = radius*(xmax-xmin)
                radius_y = radius*(ymax-ymin)
                if abs(event.xdata - x[i]) < radius_x and abs(event.ydata - y[i]) < radius_y:
                    txt1 = r'Sim length = {}'.format(lastPols[i].simlen)
                    txt2 = r'  $<\Pi>=' + '{:.5f}$'.format(lastPols[i].mValue)
                    #txt2 = r'$<\Delta\Pi^2>=' + '{:.5f}$'.format(lastPols[i].var)
                    txt = txt1 + '\n' + txt2
                    #print lastPols[i].input.noise
                    ax1.annotate(txt, xy=(x[i], y[i]),  xycoords='data',
                        xytext=(-130, 40), textcoords='offset points',
                        size=20,
                        bbox=dict(boxstyle="round", fc="0.8"),
                        arrowprops=dict(arrowstyle="wedge,tail_width=0.7",
                                        fc="0.6", ec="none",
                                        #patchB=el,
                                        connectionstyle="arc3,rad=-0.3"),
                        )

                    plt.draw()
                    collisionFound = True
                    break
        if not collisionFound:
            plt.draw()
    def mousefunc(event):
        if event.button == 2:
            print event
            print scat
        
    if (args.xlog):
        ax1.set_xscale('log')
    if (args.ylog):
        ax1.set_yscale('log')
    fig.canvas.mpl_connect('pick_event', onpick_polar)
    fig.canvas.mpl_connect('button_press_event', mousefunc)
    fig.canvas.mpl_connect('motion_notify_event', onpick_hover)
    filename = '{}_vs_noise.png'.format(txt).replace(" ","")
    fig.savefig(filename, dpi=fig.dpi)
    
    # radio!
    if not args.minimal:
        axcolor = 'lightgoldenrodyellow'
        rax = plt.axes([0.7, 0.7, 0.3, 0.15], axisbg=axcolor)
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

    ## plot transition curve!
    plot_transition_curves(ax1, 'align_density_trans.pickle')
    if isfile('density_align_trans.pickle'):
        print 'yay'
        plot_transition_curves(ax1, 'density_align_trans.pickle', sortbyrho=True)
    plt.show()