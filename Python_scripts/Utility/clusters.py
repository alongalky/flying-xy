#rom matplotlib.patches import Circle, Wedge, Polygon
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
from multiprocessing import Pool
import multiprocessing
import sys
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
import networkx as nx

buffer_in_mb = 100


def distanceSq(point, rest, box_x):
    delta = np.abs(point - rest)
    delta = np.where(delta > 0.5 * box_x, box_x - delta, delta)
    return (delta ** 2).sum(axis = 1)

def neighbors(x,y,part1,part2, input):
    coord1 = np.array([x[part1], y[part1]])
    coord2 = np.array([x[part2], y[part2]])
    return distanceSq(coord1, coord2, np.array([input.box_x, input.box_x])) < 1.0
def make_graph(x,y, input):
    G=nx.Graph()
    G.add_nodes_from(range(0,input.npart))
    edges = []
    coords = np.vstack((x,y)).T
    for i in range(0,input.npart-1):
        rest = np.arange(i+1,input.npart)
        neighbors = rest[distanceSq(coords[i+1:], coords[i], np.array([input.box_x, input.box_x])) < 1.0]
        edges += [(i, j) for j in neighbors]
    G.add_edges_from(edges)
    return G
def make_edges(x,y, input):
    coords = np.vstack((x,y)).T
    edge_list = [[]] * input.npart
    for i in range(0,input.npart-1):
        rest = np.arange(i+1,input.npart)
        neighbors = rest[distanceSq(coords[i+1:], coords[i], np.array([input.box_x, input.box_x])) < 1.0]
        edge_list[i] = neighbors
    return edge_list
def cluster_analysis(x,y, input):
    start_time = datetime.datetime.now()
    G = make_graph(x,y, input)
    clusters = nx.connected_components(G)
    return clusters
def cluster_diameter(cluster, x, y, input):
    coords = np.vstack((x,y)).T
    size = len(cluster)
    dists = []
    for part1 in range(0,size-1):
        rest = cluster[part1+1:size]
        dist = distanceSq(coords[cluster[part1]], coords[rest], np.array([input.box_x, input.box_x]))
        dists.append(np.max(dist))
    return np.sqrt(np.max(dists))
def center_of_mass(cluster, x, y, input):
    coords = np.vstack((x,y)).T
    coords = coords[cluster]
    angle = coords * (2*np.pi / input.box_x)
    xi = np.cos(angle)
    ni = np.sin(angle)
    xi_bar = np.mean(xi, axis=0)
    ni_bar = np.mean(ni, axis=0)
    com_angle = np.pi + np.arctan2(-ni_bar, -xi_bar)
    com_coords =  com_angle / (2*np.pi / input.box_x)
    return (com_coords[0], com_coords[1])

def draw_cluster(output_dir,input=None):
    input = read_input_from_dir(output_dir)
    box_x = input.box_x
    sim = read_data_from_dir(output_dir, start=0, max_mem_size_mb=buffer_in_mb)
    cmap=plt.get_cmap('RdYlGn')
    for frame in sim.data[:10]:
        fig = plt.figure(figsize=(10, 10), dpi=100, facecolor='w', edgecolor='k')
        ax1 = fig.add_subplot(1,1,1)
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        theta = frame[:,2]
        ax1.quiver(x_pos/box_x, y_pos/box_x, np.cos(theta)*0.1, np.sin(theta)*0.1, color='k')
        clusters = cluster_analysis(x_pos,y_pos, input)
        big_clusters = [cluster for cluster in clusters if len(cluster) > 10]
        for i,cluster in enumerate(big_clusters):
            if len(cluster) > 10:
                #print cluster
                diam = cluster_diameter(cluster, x_pos, y_pos, input)
                #ax1.scatter(x_pos[cluster]/box_x, y_pos[cluster]/box_x, color=cmap(float(i)/len(big_clusters)))
                (x,y) = center_of_mass(cluster, x_pos, y_pos, input)
                
                #ax1.scatter(x/box_x,y/box_x,diam*10, alpha = 0.5, color=cmap(i/len(big_clusters)))
                fig.gca().add_artist(plt.Circle((x/box_x,y/box_x),0.6*diam/box_x,alpha = 0.5,color='blue'))
                
        #print 'max cluster size',np.max([len(cluster) for cluster in clusters])
        plt.show()
def make_cluster_movie(output_dir, filename=None, show=True):
    #print output_dir
    input = read_input_from_dir(output_dir)
    box_x = input.box_x
    fig = plt.figure(figsize=(10, 10), dpi=100, facecolor='w', edgecolor='k')
    #fig = plt.figure()
    ax1 = fig.add_subplot(2,2,2)
    ax2 = fig.add_subplot(2,2,1)
    ax3 = fig.add_subplot(2,2,4)
    ax4 = fig.add_subplot(2,2,3)
    #ax1.set_ylabel(txt, fontsize=14, fontweight='bold')
     
    #fig.suptitle('{} Model, t = {:4.0f}'.format(model, snap.time), fontsize=14, fontweight='bold', color=color)#, y=0.9)
    def init():
        rect = fig.patch
        rect.set_facecolor('w')
        line, = ax1.plot([], [])
        line, = ax2.plot([], [])
        line, = ax3.plot([], [])
        line, = ax4.plot([], [])
        return line,
    sim = read_data_from_dir(output_dir, start=0, max_mem_size_mb=buffer_in_mb)
    
    sim.data = sim.data[:]
    frame_count = len(sim.data)
    
    #frame_count = 
    last_time = sim.times[-1]
    times = []
    vars = []
    def animate(i):
        ax1.cla()
        ax2.cla()
        ax3.cla()
        ax4.cla()
        ax1.set_title('Gaussian base PDF', fontsize=14, fontweight='bold')
        ax2.set_title('PDF histogram', fontsize=14, fontweight='bold')
        ax3.set_title('Snapshot', fontsize=14, fontweight='bold')
        ax4.set_title('Variance', fontsize=14, fontweight='bold')
        ax2.set_xlabel('density', fontsize=10)#, fontweight='bold')
        ax4.set_xlabel('time', fontsize=14, fontweight='bold')
        ax4.set_ylabel('PDF variance', fontsize=14, fontweight='bold')
        res = 60
        resolution = res * 1j
        X, Y = np.mgrid[0:box_x:resolution, 0:box_x:resolution]
        grid = np.vstack([X.ravel(), Y.ravel()])
        frame = sim.data[i]
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        theta = frame[:,2]
        positions = np.vstack([x_pos, y_pos]).T
        test = []
        sigma = 5.0
        for point in grid.T:
            test.append(eval_density(positions, point, input.box_x, sigma = sigma))
        PDF = np.array(test)
        particle_area = np.pi * np.sqrt(2 * sigma)
        PDF /= particle_area
        
        bin_area = np.square(input.box_x) / PDF.size
        pic = PDF.reshape((res, res))
        ax1.imshow(np.rot90(pic))
        hist, bin_edges = np.histogram(PDF, bins = res)
        #density = float(input.npart) / np.square(input.box_x)
        ax2.plot(np.arange(len(hist)), hist)
        ax2.set_yscale('log')

        snap = snapShot(x_pos,y_pos, theta, [], [], sim.times[i], input)
        ax3.quiver(x_pos/box_x, y_pos/box_x, np.cos(theta)*0.1, np.sin(theta)*0.1, color='k')
        
        times.append(sim.times[i])
        vars.append(np.var(PDF.flatten()))
        ax4.plot(times, vars)
        ax4.set_xlim([0,last_time])
        if True and len(vars)>2: # Fit!
            fitfunc = lambda p, x: p[2]*np.power(x, p[0])+p[1]    # Target function
            #fitfunc = lambda p, x: x*p[0]    # Target function
            errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
            t = np.array(times)
            v = np.array(vars)
            p1 = [1.0,0.1, 1.0] # Initial guess for the parameters
            p1, success = optimize.leastsq(errfunc, p1, args=(t, v))
            ax4.plot(t, fitfunc(p1, t)) # Plot of the fit
            ax4.text(0.04, 0.1, r'$\beta+\gamma t^\alpha$, $\alpha={:0.3f}$, $\beta={:0.3f}$, $\gamma={:0.3f}$'.format(p1[0], p1[1], p1[2]), fontsize=12,
                verticalalignment='center', transform=ax4.transAxes)
        
        
        plt.draw()
        if (i % 2 == 0):
            print printtime() + 'Finished {} / {} frames for D = {:.2f}'.format(i, frame_count, input.noise)
        line, = ax3.plot([], [])
        return line,

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=frame_count, interval=1, blit=True)
    if filename != None:
        start_time = datetime.datetime.now()
        print printtime() + 'Saving movie...'
        kwargs = {'facecolor':fig.get_facecolor(), 'edgecolor':fig.get_facecolor()}
        anim.save(filename, fps=15, dpi=200, savefig_kwargs=kwargs)
        now = datetime.datetime.now()
        print 'Finished saving movie! It took {}:{:02d}'.format((now - start_time).seconds // 60,(now - start_time).seconds % 60)
    if show==True:
        plt.show()
    return
def move_pool_wrapper(output_dir):
    input = read_input_from_dir(output_dir)
    make_pdf_movie(output_dir, filename='pdf_noise{:.2f}.mp4'.format(input.noise), show=False)
    #make_pdf_movie(output_dir, filename=None, show=True)
"""
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
    parser.add_argument('-R','--rhosq', action='store_true',
                   help='Show Rho^2 data')
    parser.add_argument('-P','--polar', action='store_true',
                   help='Show polarization data')
    parser.add_argument('-E','--energy', action='store_true',
                   help='Show energy data')
    parser.add_argument('-t','--time', action='store_true',
                   help='Plot properties vs time')
    parser.add_argument('-B','--binder', action='store_true',
                   help='Plot Binder cumulant')
    parser.add_argument('-D','--pdf', action='store_true',
                   help='Plot PDF quantities')
    parser.add_argument('-y','--ylog', action='store_true',
                   help='Logarithmic y axis')
    parser.add_argument('-x','--xlog', action='store_true',
                   help='Logarithmic x axis')
    args = parser.parse_args()
    # if you reread you better recalculate means
    if (args.read):
        args.calc = True
    resultsDir = args.directory
    procs = args.proc
    
    output_dirs = []
    for output_dir in os.listdir(resultsDir):
        if 'output' in output_dir and os.path.isdir(resultsDir + '\\' + output_dir):
            output_dirs.append(resultsDir + '\\' + output_dir)
    p = Pool(procs)
    #p.map(move_pool_wrapper, output_dirs[:5])
    map(draw_cluster, output_dirs[:1])
"""
"""