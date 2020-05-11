from numpy.fft import fft2, fftshift, rfft2
from scipy import ndimage
import networkx as nx
import matplotlib.animation as animation
import os
import sys
import shutil
import time
import datetime
import glob
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from os.path import isfile,join,isdir
###
#from plots import plot_snap
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def noiseFunc(input, times):
    halfsim = input.simtime / 2.0
    times = np.array(times)
    low_times = times[times < halfsim]
    low_noises = input.noise + (low_times / halfsim)*(input.hyster_noise - input.noise)

    high_times = times[times >= halfsim] - halfsim
    high_noises = input.hyster_noise + (high_times / halfsim)*(input.noise - input.hyster_noise)
    return np.concatenate([low_noises, high_noises])

class CubicFit:
    def __init__(self, ax, x_dat, y_dat):
        self.ax = ax
        self.x_dat = x_dat
        self.y_dat = y_dat
        self.line, = ax.plot([0,1],[0,1])
        self.v_line = self.ax.vlines(0, 0, np.max(self.y_dat))
        fig = ax.get_figure()
        self.canvas = fig.canvas
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.state=0
        self.left = x_dat[0]
        self.right = x_dat[-1]
        self.x_trans = 0
        self.clicked = False
        self.update()
        plt.show()
    def update(self):
        mask = np.logical_and(self.x_dat >= self.left, self.x_dat <= self.right)
        x = self.x_dat[mask]
        y = self.y_dat[mask]
        p = np.polyfit(x, y, 3)
        xfit = np.linspace(self.x_dat[0], self.x_dat[-1], 100)
        yfit = np.polyval(p, xfit)
        self.line.set_data(xfit,yfit)
        self.v_line.remove()
        self.x_trans  = -p[1]/(3*p[0])
        self.v_line = self.ax.vlines(self.x_trans, 0, np.max(self.y_dat), linestyle='--', linewidth=1.5)
        plt.draw()
    def get_trans(self):
        if self.clicked == True:
            return self.x_trans
        else:
            return None
    def button_press_callback(self, event):
        if event.button != 3: return
        """
        if event.xdata < self.x_trans:
            self.left = event.xdata
            print 'Left = ', self.left
        else:
            self.right = event.xdata
            print 'Right = ', self.right
        """
        if self.state==0:
            self.clicked = True
            self.state = 1-self.state
            self.left = event.xdata
            print 'Left = ', self.left
        elif self.state==1:
            self.clicked = True
            self.state = 1-self.state
            self.right = event.xdata
            print 'Right = ', self.right
        self.update()
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
    title = '{}, {} Model, {} particles'.format(txt, model, input.npart)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('time', fontsize=14, fontweight='bold')
    ax.set_ylabel(txt, fontsize=14, fontweight='bold')
    line, = ax.plot(prop.times, prop.values, label=r'$A$={}, $D$={}, $\rho$={:.3f}'.format(input.align, input.noise, getDensityInput(input)), picker = True)
    plt.legend(loc='best')
    #line, = ax.plot([0], [0])  # empty line
    linebuilder = LineBuilder(line)

    plt.show()
    return linebuilder.get_transient_time()
    
def get_order_parameters():
    return [orderParamter('Polarization', 'Pols', Pol_data_calc, mean_pickle_name='mPol', ylabel=r'$<\Pi>$'),
                  orderParamter('Cluster Fraction (neighbors)', 'neighbor_pdf', calc_neighbors_pdf, mean_pickle_name='mNeighbors', ylabel=r'Cluster fraction'),
                  #orderParamter('Cluster radius of gyration', 'mean_cluster_gyration_radius', calc_cluster_gyration_radius, mean_pickle_name='mGyration', ylabel=r'$<{R_g}^2>$'),
                  #orderParamter('Cluster size', 'mean_cluster_size', calc_mean_cluster_size, mean_pickle_name='mClusterSize', ylabel=r'$<N_A>$'),
                  #orderParamter('Cluster density', 'mean_cluster_density', calc_mean_cluster_density, mean_pickle_name='mClusterDensity'),
                  #orderParamter('Cluster size PDF', 'mean_cluster_size_pdf', calc_cluster_size_pdf, mean_pickle_name='mClusterSizePDF'),
                  orderParamter('Density variance', 'densvar', calc_vars, mean_pickle_name='mDensVar', ylabel=r'$<(\delta \rho)^2>$'),
                  #orderParamter('Binder cumulant', 'binder', calc_binder, mean_pickle_name='mBinder', ylabel=r'$G$'),
                  #orderParamter('Energy', 'Energy', Energy_data_calc, mean_pickle_name='mEnergy', ylabel=r'$E_{\rm align}$'),
                  #orderParamter('Mean local momentum', 'momentum', calc_W, mean_pickle_name='mMomentum', ylabel=r'$|<W_l(\vec{r},t)>|$'),
                  #orderParamter('Mean local momentum magnitude', 'abs_momentum', calc_pol_absmean, mean_pickle_name='mAbsMomentum', ylabel=r'$<|W_l(\vec{r},t)|>$'),
                  
                  #orderParamter('Axis inhomogeneity','BandDens',calc_band_dens, smooth=True),
                  #orderParamter('Axis inhomogeneity along polarization','BandDensPol',calc_band_dens_pol)
                  ]

class transitionMethods:
    polInflection, polFluctuation = range(2)    
    names = [r'$<\Pi>$ inflection point', r'$<(\delta\Pi)^2>$ maximum']
class transitionPoint:
    def __init__(self, noise, density, align=1.0, pols=None, method=None, const_noise=False, const_align=False):
        self.noise = noise
        self.density = density
        self.align = align
        self.pols = pols
        self.method = method
        self.const_noise = const_noise
        self.const_align = const_align
    def isSame(self, other):
        return (self.noise == other.noise and    
                    self.density == other.density and    
                    self.align == other.align and    
                    self.method == other.method)
class orderParamter():
  def __init__(self, name, pickle_name, calc_data_fun, mean_pickle_name=None,smooth=False, ylabel=None):
     self.name = name
     self.pickle_name = pickle_name
     self.calc_data_fun = calc_data_fun
     self.smooth = smooth
     self.mean_pickle_name = mean_pickle_name
     self.ylabel = ylabel
class mean_propInfo():
  def __init__(self, mValue, simlen, input, output_dir, startcalctime=None, 
    name=None, var=None):
     self.mValue = mValue
     self.simlen = simlen
     self.input = input
     self.output_dir = output_dir
     self.startcalctime = startcalctime
     self.name = name
     self.var = var
class propInfo():
  def __init__(self, values, times, input, prop_name, output_dir):
     self.values = values
     self.times = times
     self.input = input
     self.prop_name = prop_name
     self.output_dir = output_dir
class clusterStats():
  def __init__(self, lengths, gyration_radii_sq, pols=None):
     self.lengths = lengths
     self.gyration_radii_sq = gyration_radii_sq
     self.pols = pols

def distanceSq(point, rest, box_x):
    delta = np.abs(point - rest)
    delta = np.where(delta > 0.5 * box_x, box_x - delta, delta)
    return (delta ** 2).sum(axis = 1)

def eval_density(positions, point, box_x, sigma=5):
    dimensions = np.array([box_x, box_x])
    delta = np.abs(positions - point)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    dists_sq = (delta ** 2).sum(axis=-1)
    gaussian_func_val = np.exp(-dists_sq / np.sqrt(2*sigma))
    # for normalization purposes, calculate single particle area (a Gaussian integral)
    particle_area = np.pi * np.sqrt(2 * sigma)
    return gaussian_func_val.sum() / particle_area

def eval_pol_density(positions, thetas, point, box_x, sigma = 5):
    # TODO: think :(
    #return None
    dimensions = np.array([box_x, box_x])
    delta = np.abs(positions - point)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    dists_sq = (delta ** 2).sum(axis=-1)
    gaussian_func_val = np.exp(-dists_sq / np.sqrt(2*sigma))
    pol = np.sum(np.exp(1j * thetas) * gaussian_func_val)
    particle_area = np.pi * np.sqrt(2 * sigma)
    ret = pol / particle_area # / gaussian_func_val.sum()
    return ret     
def calc_property(output_dir,prop_name, prop_data_func, read=False, resume=True,
        buffer_in_mb=100):
    filename = join(output_dir,prop_name + '.pickle')
     
    times = []
    Props = []
    last_calc_time = 0
    if isfile(filename) and read == False:
        try:
            prop = pickle_load(filename)
        except:
            print 'error in loading ', filename
            return calc_property(output_dir,prop_name, prop_data_func, read=True, resume=resume,
                    buffer_in_mb=buffer_in_mb)
        if not resume :
            return prop
        if len(prop.times) > 0:
            # making sure the file isn't empty
            last_calc_time = prop.times[-1]
            last_sim_time = last_time_in_direc(output_dir)
            if last_sim_time <= last_calc_time:
                print printtime(),output_dir,' property {} pre-calculated. '.format(prop_name), last_sim_time, last_calc_time
                return prop
            times = prop.times
            Props = prop.values
   
    input = read_input_from_dir(output_dir)
    print 'starting calc_property {}, output dir = {}'.format(prop_name,output_dir)
    npart = input.npart


    start_time = datetime.datetime.now()
    frames_read = 0
    files_written = 0
    time = 0
    while True:
        sim = read_data_from_dir(output_dir, start=frames_read, max_mem_size_mb=buffer_in_mb)
        if len(sim.data) ==  0:        # finished reading data dir
            break
        frames_read += len(sim.data)
        if sim.times[-1] < last_calc_time:
            continue
        first_index = 0
        if sim.times[0] < last_calc_time:
            first_index = next(tup[0] for tup in enumerate(sim.times) if tup[1] > last_calc_time)
        # debug
        if len(sim.times[first_index:]) == 0:
            print 'Error! No data in ', output_dir
            return None
        #
        times += sim.times[first_index:]
        Props.extend(prop_data_func(sim.data[first_index:], input=input))

        last_time = datetime.datetime.now()
        sec_per_buf = (last_time - start_time).microseconds / 1000000.0
        print '{} MB calculation time: {}:{:03d}. Last timestamp: {:.2f} Dir: {}'.format(buffer_in_mb,(last_time - start_time).seconds, 
            (last_time - start_time).microseconds // 1000, times[-1] if times else 0, output_dir)
        start_time = datetime.datetime.now()
    
    if len(Props) != len(times) or len(times) == 0:
        print 'Warning: directory {} is probably empty'.format(output_dir)
        return None
    ret = propInfo(Props, times, input, prop_name, output_dir)
    pickle_save(filename, ret)
    return ret
def calc_vars(data,input):
    box_x = input.box_x
    res = 50
    resolution = res * 1j
    X, Y = np.mgrid[0:box_x:resolution, 0:box_x:resolution]
    grid = np.vstack([X.ravel(), Y.ravel()])
    vars = []
    for i,frame in enumerate(data):
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        positions = np.vstack([x_pos, y_pos]).T
        test = []
        sigma = 5.0
        for point in grid.T:
            test.append(eval_density(positions, point, input.box_x, sigma = sigma))
        vars.append(np.var(test))
        if (i+1) % 20 == 0:
            print data_calc_prompt(i+1, len(data), input)
    return vars
def calc_pol_absmean(data,input):
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
        
        values.append(np.mean(np.abs(test)))
        if (i+1) % 20 == 0:
            print printtime() + 'D={:.2f} Finished {} / {}'.format(input.noise, i+1, len(data))
    return values
def calc_W(data,input):
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
        if (i+1) % 20 == 0:
            print printtime() + 'D={:.2f} dens={:.2f} Finished {} / {}'.format(input.noise, getDensityInput(input),i+1, len(data))
    return values    
def Pol_data_calc(data, input = None):
    Pols = []
    for frame in data:
        thetas = frame[:,2]
        pol = np.abs(np.mean(np.exp(1j * thetas)))
        Pols.append(pol)
        if len(Pols) % 50 == 0:
            print data_calc_prompt(len(Pols), len(data), input)
    Pols = np.array(Pols)
    return Pols
def calc_pol_var(output_dir, read=False, calc=False):
    prop_name = 'mPol'
    filename = join(output_dir,prop_name + '.pickle')
    
    last_sim_time = last_time_in_direc(output_dir)
    if last_sim_time < 20:
        print output_dir, ' is too short! not calculating pol! last_time is ', last_sim_time
        return None
    # else:
        # print last
    trans_time = None
    #print isfile(filename), read, calc
    if isfile(filename) and read == False and calc==False:
        mProp = pickle_load(filename)
        last_calc_time = mProp.simlen
        if last_sim_time <= last_calc_time:
            return mProp
        trans_time = mProp.startcalctime    
    input = read_input_from_dir(output_dir)
    pols = calc_property(output_dir,'Pols',Pol_data_calc, read=read)
    if pols == None:    # empty directory or some evil silliness
        print 'Error with calc_pol in ', output_dir
        return None
    if trans_time == None or calc==True:
        trans_time = choose_transient_time(pols)
    if trans_time >= last_sim_time:
        print '2) Error with calc_pol in ', output_dir
        return None
    ### do math:
    first_ind = next(x[0] for x in enumerate(pols.times) if x[1] > trans_time)
    
    vals = pols.values[first_ind:]
    mean_pol = np.mean(vals)
    ans = np.var(vals)
    ret = mean_propInfo(mean_pol, pols.times[-1], input, output_dir, startcalctime=trans_time,
            name = prop_name, var = ans)
    pickle_save(filename, ret)
    return ret
def calc_mProp_orderparam(output_dir, order_param,read=False, calc=False):
    if not isinstance(order_param, orderParamter):
        print 'Error! use with order parameter'
        return None
    return calc_mProp(output_dir, order_param.mean_pickle_name, order_param.pickle_name, order_param.calc_data_fun, read, calc)
def calc_mProp(output_dir, prop_name, pickleFilename,calc_data_fun,read=False, calc=False):
    filename = join(output_dir,prop_name + '.pickle')
    print filename
    last_sim_time = last_time_in_direc(output_dir)
    if last_sim_time < 20:
        return None
    # else:
        # print last
    trans_time = None
    if isfile(filename) and read == False and calc==False:
        mProp = pickle_load(filename)
        last_calc_time = mProp.simlen
        if last_sim_time <= last_calc_time:
            return mProp
        trans_time = mProp.startcalctime    
    input = read_input_from_dir(output_dir)
    pols = calc_property(output_dir,pickleFilename,calc_data_fun, read=read)
    if pols == None:    # empty directory or some evil silliness
        return None
    if trans_time == None or calc==True:
        trans_time = choose_transient_time(pols)
    if trans_time >= last_sim_time:
        return None
    ### do math:
    first_ind = next(x[0] for x in enumerate(pols.times) if x[1] > trans_time)
    
    vals = pols.values[first_ind:]
    if prop_name == 'mMomentum':
        print 'calculating mean abs momentum'
        mean_pol = np.mean(np.abs(vals))
    else:
        mean_pol = np.mean(vals)    
    ans = np.var(vals)
    ret = mean_propInfo(mean_pol, pols.times[-1], input, output_dir, startcalctime=trans_time,
            name = prop_name, var = ans)
    pickle_save(filename, ret)
    return ret
def calc_mClusterStats(output_dir,read=False, calc=False):
    prop_name = 'mClusterStats'
    filename = join(output_dir,prop_name + '.pickle')
    filename_meanclustersizenumber = join(output_dir, 'mean_cluster_size_number.pickle')
    print filename
    last_sim_time = last_time_in_direc(output_dir)
    if last_sim_time < 20:
        return None
    # else:
        # print last
    trans_time = None
    #print isfile(filename), isfile(join(output_dir,'mGyrationRadiusSq.pickle')) , read, calc
    if isfile(filename) and isfile(join(output_dir,'mGyrationRadiusSq.pickle')) and isfile(filename_meanclustersizenumber) and read == False and calc==False:
        mProp = pickle_load(filename)
        last_calc_time = mProp.simlen
        if last_sim_time <= last_calc_time:
            return mProp
        trans_time = mProp.startcalctime    
    input = read_input_from_dir(output_dir)
    #pols = calc_property(output_dir,'Pols',Pol_data_calc, read=False)
    mPol = calc_pol_var(output_dir)
    trans_time = mPol.startcalctime
    if mPol == None:    # empty directory or some evil silliness
        return None
    if trans_time >= last_sim_time:
        return None
    print printtime() + 'loading cluster stats'
    stats = calc_property(output_dir,'cluster_stats', calc_cluster_stats, read=read)
    print printtime() + 'Finished loading cluster stats'
    
    ### do math:
    #print len(stats.times), len(stats.values)
    first_ind = next(x[0] for x in enumerate(stats.times) if x[1] > trans_time)
    
    big_bin_count = np.zeros(input.npart)
    vals = stats.values[first_ind:]
    i = 0
    mean_size_nums = []
    mean_radii_sq = []
    for val in vals:  
        lengths = val.lengths
        count = np.bincount(lengths)
        count.resize(big_bin_count.shape)  # pad with zeros
        #print count.shape, big_bin_count.shape
        big_bin_count += count
        mean_size_nums.append(np.mean(lengths))
        mean_radii_sq.append(np.mean(val.gyration_radii_sq))
        i += 1
        #print i, np.mean(val.gyration_radii_sq)
    #print 'debug: calculating data, outputting'
    ans = big_bin_count / len(vals)
    ret = mean_propInfo(ans, stats.times[-1], input, output_dir, startcalctime=trans_time,
            name = prop_name)
    pickle_save(filename, ret)
    ### save mean cluster_size_number as well!
    mean_cluster_size_number = np.mean(mean_size_nums)
    ret3 = mean_propInfo(mean_cluster_size_number, stats.times[-1], input, output_dir, startcalctime=trans_time,
            name = prop_name)
    pickle_save(filename_meanclustersizenumber, ret3)
    ###
    ### save mean gyration_radius as well!
    mean_gyration_radius_sq = np.mean(mean_radii_sq)
    if np.isnan(mean_gyration_radius_sq):
        raise Exception('bad rsq value')
    ret2 = mean_propInfo(mean_gyration_radius_sq, stats.times[-1], input, output_dir, startcalctime=trans_time,
            name = prop_name)
    print 'saving gyration radius pickle'
    pickle_save(join(output_dir,'mGyrationRadiusSq.pickle'), ret2)
    ###
    try:
        # error if pols unimplemented
        print vals[0].pols[:5]
    except:
        print 'Warning: no cluster pol info!'
        return ret
    cluster_pol_pdf_unweighted = np.mean(mean_radii_sq)
    mean_pols_mProp = calc_property(output_dir,'Pols',Pol_data_calc)
    mean_pols = mean_pols_mProp.values[first_ind:]
    if len (vals) != len(mean_pols):
        raise Exception('integrity error')
    res_angle = 360
    bins = np.linspace(-np.pi, np.pi, res_angle+1)
    total_hist_unweighted = np.zeros(res_angle)
    total_hist_weighted = np.zeros(res_angle)
    total_hist_big = np.zeros(res_angle)
    i = 0
    ####
    #for time1, time2 in zip(mean_pols_mProp.times[first_ind:], stats.times[first_ind:]):
    #    print time1, time2, ' Same' if time1==time2 else 'Different'
    #return
    ###
    for pol,val in zip(mean_pols,vals): 
        #if mean
        lengths = val.lengths
        pols = val.pols
        #hist,bin_edges = np.histogram(angles, bins=bins)
        
        # normalize
        pol = pol / np.absolute(pol)
        pols = pols / np.absolute(pols)
        # rotate
        pols /= pol
        
        angles = np.angle(pols)
        hist_unweighted,bin_edges = np.histogram(angles, bins=bins)
        hist_weighted,bin_edges = np.histogram(angles, bins=bins, weights=lengths)
        
        total_hist_unweighted += hist_unweighted
        total_hist_weighted += hist_weighted      
        
        ## filter small clusters
        angles = angles[val.lengths >= 20]
        hist_big,bin_edges = np.histogram(angles, bins=bins)
        total_hist_big  += hist_big
        
    mean_hist_unweighted = total_hist_unweighted / float(len(vals))
    mean_hist_weighted = total_hist_weighted / float(len(vals))
    mean_hist_big = total_hist_big / float(len(vals))
    #print mean_hist_unweighted.shape
    ret3 = mean_propInfo(mean_hist_unweighted, stats.times[-1], input, output_dir, startcalctime=trans_time,
            name = 'unweighted_cluster_pol_pdf')
    print 'saving unweighted_cluster_pol_pdf pickle'
    pickle_save(join(output_dir,'unweighted_cluster_pol_pdf.pickle'), ret3)
    ret4 = mean_propInfo(mean_hist_weighted, stats.times[-1], input, output_dir, startcalctime=trans_time,
            name = 'weighted_cluster_pol_pdf')
    print 'saving weighted_cluster_pol_pdf pickle'
    pickle_save(join(output_dir,'weighted_cluster_pol_pdf.pickle'), ret4)
    ret5 = mean_propInfo(mean_hist_big, stats.times[-1], input, output_dir, startcalctime=trans_time,
            name = 'big_cluster_pol_pdf')
    print 'saving big_cluster_pol_pdf pickle'
    pickle_save(join(output_dir,'big_cluster_pol_pdf.pickle'), ret5)
        
    ###
    return ret
def calc_mNeighbors(output_dir, read=False, calc=False):
    prop_name = 'mNeighbors'
    filename = join(output_dir,prop_name + '.pickle')
    last_sim_time = last_time_in_direc(output_dir)
    if last_sim_time < 20:
        return None
    # else:
        # print last
    trans_time = None
    if isfile(filename) and read == False and calc==False:
        mProp = pickle_load(filename)
        last_calc_time = mProp.simlen
        if last_sim_time <= last_calc_time:
            return mProp
        trans_time = mProp.startcalctime    
    input = read_input_from_dir(output_dir)
    pols = calc_property(output_dir,'neighbor_pdf',calc_neighbors_pdf, read=read)
    if pols == None:    # empty directory or some evil silliness
        return None
    print 'debug3: mneighbors ', output_dir
    if trans_time == None or calc==True:
        trans_time = choose_transient_time(pols)
    if trans_time >= last_sim_time:
        return None
    print 'debug2: mneighbors ', output_dir
    ### do math:
    first_ind = next(x[0] for x in enumerate(pols.times) if x[1] > trans_time)
    
    vals = pols.values[first_ind:]
    mean_pol = np.mean(vals)
    ans = np.var(vals)
    ret = mean_propInfo(mean_pol, pols.times[-1], input, output_dir, startcalctime=trans_time,
            name = prop_name, var = ans)
    pickle_save(filename, ret)
    return ret
def calc_neighbors_pdf(data,input):
    box_x = input.box_x
    rstar = input.rstar
    rstarsq = np.square(input.rstar)
    npart = input.npart
    density = getDensityInput(input)
    threshold_neighbors = np.pi*rstarsq*density+1
    
    cluster_fracs = []
    print 'debug: threshold_neighbors = ', threshold_neighbors
    #for frame in (data[len(data)//2:])[:1]:
    for frame in data:
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        #theta = frame[:,2]
        coords = np.vstack((x_pos,y_pos)).T
        #neighbors = np.zeros(npart)
        cluster_frac = 0
        for i in range(input.npart):
            dists = distanceSq(coords[i], coords, box_x)
            #neighbor_thetas = theta[dists < rstarsq]
            neighbors = np.sum(dists < rstarsq) #including yourself!
            if neighbors > threshold_neighbors:
                cluster_frac += 1
        cluster_frac = 1.0*cluster_frac / npart
        cluster_fracs.append(cluster_frac)
        """
        print 'average neighborcount is ', np.mean(neighbors)
        max = np.max(neighbors)
        min = np.min(neighbors)
        extent = int(max-min)
        print 'max neighborcount is ', max
        print 'min neighborcount is ', min
        log_neighbor = np.pi*rstarsq*density
        print 'logical neighborcount is ', log_neighbor
        plt.hist(neighbors, bins=extent)
        colors = ['r' if x > log_neighbor else 'b' for x in neighbors]
        plot_snap_fake(snapShot(x_pos,y_pos, theta, [], [], 6, input), show=False, colors = colors)
        plt.show()
        """
        if len(cluster_fracs) % 50 == 0:
            print printtime() + 'Finished {} / {} in D={:.2f}'.format(len(cluster_fracs), len(data), input.noise)
    """
    x = range(len(cluster_fracs))
    plt.plot(x,cluster_fracs)
    plt.show()
    """
    return cluster_fracs   
def calc_mean_cluster_density(data, input):
    box_x = input.box_x
    rstar = input.rstar
    rstarsq = np.square(input.rstar)
    npart = input.npart
    density = getDensityInput(input)
    threshold_neighbors = np.pi*rstarsq*density+1
    
    cluster_densities = []
    #for frame in (data[len(data)//2:])[:1]:
    for frame in data:
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        #theta = frame[:,2]
        clusters = cluster_analysis(x_pos, y_pos, input)
        radii_squared = np.array([cluster_gyration_radius(x, x_pos, y_pos, input) for x in clusters])
        areas = radii_squared * np.pi
        particle_nums = np.array([len(x) for x in clusters])
        
        densities = particle_nums / areas
        avg_dens = np.mean(densities)
        
        cluster_densities.append(avg_dens)
        if len(cluster_densities) % 50 == 0:
            print data_calc_prompt(len(cluster_densities), len(data), input)
    return cluster_densities   
def calc_cluster_stats(data,input):
    box_x = input.box_x
    rstar = input.rstar
    rstarsq = np.square(input.rstar)
    npart = input.npart
    density = getDensityInput(input)
    threshold_neighbors = np.pi*rstarsq*density+1
    
    cluster_stats = []
    #for frame in (data[len(data)//2:])[:1]:
    for frame in data:
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        theta = frame[:,2]
        clusters = list(cluster_analysis(x_pos, y_pos, input))
        #print len(list(clusters))
        lengths = [len(x) for x in clusters]
        gyration_radii_sq = [cluster_gyration_radius(x, x_pos, y_pos, input) for x in clusters]
        if np.isnan(np.min(gyration_radii_sq)):
            print gyration_radii_sq
            raise Exception ('getting rsq = inf!')
        pols = [cluster_polarization(x, theta) for x in clusters]
        if len(gyration_radii_sq) != len(lengths) or len(pols) != len(lengths):
            raise Exception('mmm')
        cluster_stats.append(clusterStats(lengths, gyration_radii_sq, pols=pols))
        if len(cluster_stats) % 50 == 0:
            print data_calc_prompt(len(cluster_stats), len(data), input)
    return cluster_stats   
def calc_mean_cluster_size(data,input):
    box_x = input.box_x
    rstar = input.rstar
    rstarsq = np.square(input.rstar)
    npart = input.npart
    density = getDensityInput(input)
    threshold_neighbors = np.pi*rstarsq*density+1
    
    cluster_sizes = []
    #for frame in (data[len(data)//2:])[:1]:
    for frame in data:
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        #theta = frame[:,2]
        clusters = cluster_analysis(x_pos, y_pos, input)
        #print len(list(clusters))
        lengths = [len(x) for x in clusters]
        avg_size = np.mean(lengths)
        #print avg_size
        cluster_sizes.append(avg_size)
        #sys.exit()
        if len(cluster_sizes) % 50 == 0:
            print printtime() + 'Finished {} / {} in D={:.2f}'.format(len(cluster_sizes), len(data), input.noise)
    return cluster_sizes   
def calc_binder(output_dir, read=False):
    prop_name = 'Binder'
    filename = join(output_dir,prop_name + '.pickle')
    last_sim_time = last_time_in_direc(output_dir)
    if isfile(filename) and read == False:
        mProp = pickle_load(filename)
        last_calc_time = mProp.simlen
        if last_sim_time <= last_calc_time:
            return mProp
        trans_time = mProp.startcalctime    
    print 'calculating binder for ', output_dir
    mPol = calc_pol_var(output_dir)
    trans_time = mPol.startcalctime
    if trans_time == None:
        raise Exception("no transition time!")
    Pols = calc_property(output_dir,'Pols',Pol_data_calc)
    pols = Pols.values
    first_ind = next(x[0] for x in enumerate(Pols.times) if x[1] > trans_time)
    pols = pols[first_ind:]
    pol2 = np.square(pols)
    pol4 = np.square(pol2)
    mPol2 = np.mean(pol2)
    mPol4 = np.mean(pol4)
    value = 1 - (mPol4/(3*np.square(mPol2)))
    ret = mean_propInfo(value, Pols.times[-1], Pols.input, output_dir)
    pickle_save(filename, ret)
    return ret
def calc_mBandProfile(output_dir, read=False):
    prop_name = 'mBandProfile'
    filename = join(output_dir,prop_name + '.pickle')
    last_sim_time = last_time_in_direc(output_dir)
    if isfile(filename) and read == False:
        mProp = pickle_load(filename)
        last_calc_time = mProp.simlen
        if last_sim_time <= last_calc_time:
            return mProp
        trans_time = mProp.startcalctime   
    print 'calculating mBandProfile for ', output_dir
    mPol = calc_pol_var(output_dir)
    trans_time = mPol.startcalctime
    if mPol.input.npart < 5000:
        trans_time = 0.75*last_sim_time

    if trans_time == None:
        raise Exception("no transition time!")
    Profiles = calc_property(output_dir,'BandProfile',calc_band_profile)
    data = Profiles.values
    profiles = np.empty((len(data), len(data[0])))
    for i,prof in enumerate(data):
        profiles[i] = prof
    
    profile = np.mean(profiles,0)        
    
    ret = mean_propInfo(profile, Profiles.times[-1], Profiles.input, output_dir)
    pickle_save(filename, ret)
    return ret
def calc_band_profile(data,input):
    box_x = input.box_x
    res_length = 100
    res_width = 30
    sigma = 5.0
    profiles = []
    print len(data)
    for i,frame in enumerate(data[:]):
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        thetas = frame[:,2]
        pol = np.mean(np.exp(1j * thetas))
        reverse = False
            
        if np.abs(np.real(pol)) > np.abs(np.imag(pol)):
            positions = np.vstack([x_pos,y_pos]).T
            # going left/right
            if np.real(pol) >= 0:
                reverse = True
        else:
            positions = np.vstack([y_pos,x_pos]).T
            # going up/down
            if np.imag(pol) >= 0:
                reverse = True
        
        profile = np.empty(res_length)
        for x in range(res_length):
            y_dens = np.empty(res_width)
            for y in range(res_width):
                x_coord = x*(box_x/res_length)
                y_coord = y*(box_x/res_width)
                y_dens[y] = eval_density(positions, np.array((x_coord,y_coord)), input.box_x, sigma = sigma)
            profile[x] = np.mean(y_dens)
        if reverse:
            profile = profile[::-1]
        max_ind = np.argmax(profile)
        profile = np.roll(profile, -max_ind + res_length//2)
        profiles.append(profile)
        if (i+1) % 20 == 0:
            print data_calc_prompt(len(profiles), len(data), input)
    return profiles

def data_calc_prompt(done, out_of, input):
    return printtime() + 'Finished {} / {} in A={},D={},rho={:.3f},npart={}'.format(done, out_of, input.align, input.noise, getDensityInput(input), input.npart)
def calc_cluster_gyration_radius(data,input):
    box_x = input.box_x
    rstar = input.rstar
    rstarsq = np.square(input.rstar)
    npart = input.npart
    density = getDensityInput(input)
    threshold_neighbors = np.pi*rstarsq*density+1
    
    cluster_diameters = []
    #for frame in (data[len(data)//2:])[:1]:
    for frame in data:
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        #theta = frame[:,2]
        clusters = cluster_analysis(x_pos, y_pos, input)
        diameters = [cluster_gyration_radius(x, x_pos, y_pos, input) for x in clusters]
        avg_diam = np.mean(diameters)
        #print avg_diam
        cluster_diameters.append(avg_diam)
        if len(cluster_diameters) % 50 == 0:
            print data_calc_prompt(len(cluster_diameters), len(data), input)
    
    return cluster_diameters   
def calc_gyration_radius(data,input):
    ### TODO: untested
    box_x = input.box_x
    rstar = input.rstar
    rstarsq = np.square(input.rstar)
    npart = input.npart
    density = getDensityInput(input)
    threshold_neighbors = np.pi*rstarsq*density+1
    
    cluster_diameters = []
    #for frame in (data[len(data)//2:])[:1]:
    for frame in data:
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        #theta = frame[:,2]
        clusters = cluster_analysis(x_pos, y_pos, input)
        diameters = [cluster_diameter(x, x_pos, y_pos, input) for x in clusters]
        avg_diam = np.mean(diameters)
        #print avg_diam
        cluster_diameters.append(avg_diam)
        if len(cluster_diameters) % 50 == 0:
            print printtime() + 'Finished {} / {} in A={},D={},rho={:.3f}'.format(len(cluster_diameters), len(data), input.align, input.noise, getDensityInput(input))
    return cluster_diameters   
def energy_per_particle_pdf(data,input):
    # NEVER TESTED
    box_x = input.box_x
    rstar = input.rstar
    rstarsq = np.square(input.rstar)
    pdfs = []
    for frame in data:
        energies = []
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        theta = frame[:,2]
        coords = np.vstack((x_pos,y_pos)).T
        for i in range(input.npart-1):
            rest_coords = coords[i+1:]
            rest_thetas = theta[i+1:]
            dists = distanceSq(coords[i], rest_coords, box_x)
            neighbor_thetas = rest_thetas[dists < rstarsq]
            neighbors = len(neighbor_thetas)
            theta_diff = theta[i] - neighbor_thetas
            energy = np.cos(theta_diff).sum()
            energy_pp = (1.0+energy) / (1.0 + neighbors)
        energy /= input.npart
        energies.append(energy)
        if len(energies) % 50 == 0:
            print printtime() + 'Finished {} / {} in D={:.2f}'.format(len(energies), len(data), input.noise)
    return energies   
def calc_band_dens(data,input):
    box_x = input.box_x
    res = 20
    resolution = res * 1j
    X, Y = np.mgrid[0:box_x:resolution, 0:box_x:resolution]
    grid = np.vstack([X.ravel(), Y.ravel()])
    vars = []
    print len(data)
    for i,frame in enumerate(data[:]):
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        positions = np.vstack([x_pos, y_pos]).T
        dens = []
        sigma = 5.0
        dens = np.empty([res, res])
        for m in range(res):
            for n in range(res):
                point = (X[m][n], Y[m][n])
                dens[m][n] = eval_density(positions, point, input.box_x, sigma = sigma)
        dens = np.rot90(dens)
        #avg_dens = getDensityInput(input)
        x = np.linspace(0,input.box_x,res)
        x_var = np.var(np.mean(dens, axis=0))
        y_var = np.var(np.mean(dens, axis=1))
        ratio = x_var / y_var
        ratio = np.max([ratio, 1/ratio])
        vars.append(ratio)
        if (i+1) % 20 == 0:
            print printtime() + 'D={:5.2f}: Finished {} / {}'.format(input.noise, i+1, len(data))
    return vars
def calc_band_dens_pol(data,input):
    box_x = input.box_x
    res = 20
    resolution = res * 1j
    X, Y = np.mgrid[0:box_x:resolution, 0:box_x:resolution]
    grid = np.vstack([X.ravel(), Y.ravel()])
    vars = []
    for i,frame in enumerate(data[:]):
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        thetas = frame[:,2]
        positions = np.vstack([x_pos, y_pos]).T
        dens = []
        sigma = 5.0
        dens = np.empty([res, res])
        for m in range(res):
            for n in range(res):
                point = (X[m][n], Y[m][n])
                dens[m][n] = eval_density(positions, point, input.box_x, sigma = sigma)
        dens = np.rot90(dens)
        #avg_dens = getDensityInput(input)
        x = np.linspace(0,input.box_x,res)
        x_var = np.var(np.mean(dens, axis=0))
        y_var = np.var(np.mean(dens, axis=1))
        pol = np.mean(np.exp(1j*thetas))
        pol_x = np.abs(np.real(pol))
        pol_y = np.abs(np.imag(pol))
        ratio = x_var / y_var if pol_x > pol_y else y_var / x_var
        #print pol_x, pol_y
        ratio *= np.max([pol_x, pol_y])
        vars.append(ratio)
        if (i+1) % 20 == 0:
            print printtime() + 'D={:5.2f}: Finished {} / {}'.format(input.noise, i+1, len(data))
    return vars 
def calc_band_fourier(data, input):
    box_x = input.box_x
    res = 50
    resolution = res * 1j
    X, Y = np.mgrid[0:box_x:resolution, 0:box_x:resolution]
    grid = np.vstack([X.ravel(), Y.ravel()])
    vars = []
    print len(data)
    for i,frame in enumerate(data[406:]):
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        positions = np.vstack([x_pos, y_pos]).T
        dens = []
        sigma = 5.0
        dens = np.empty([res, res])
        for m in range(res):
            for n in range(res):
                point = (X[m][n], Y[m][n])
                dens[m][n] = eval_density(positions, point, input.box_x, sigma = sigma)
        dens = np.rot90(dens)
        
        #dens = dens.T
        
        fig = plt.figure(1)
        ax1 = fig.add_subplot(221)
        ax1.imshow( dens)
        
        avg_density = input.npart / np.square(input.box_x)
        dens -= avg_density
        F1 = fft2(dens)
        
        #F2 = fftshift( F1 )
        F2 = F1
        # the 2D power spectrum is:
        psf2D = np.abs( F2 )**2
        ax2 = fig.add_subplot(222)
        ax2.imshow( psf2D)
        ax3 = fig.add_subplot(223)
        snap = snapShot(x_pos, y_pos, frame[:,2], [],[], 6, input)
        plot_snap(snap, save=False, ax=ax3)        
        ax = fig.add_subplot(224,projection='3d')
        surf = ax.plot_surface(np.arange(res), np.arange(res), psf2D, rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)
        plt.show()
        sys.exit()
        #vars.append(np.var(test))
        if (i+1) % 20 == 0:
            print printtime() + 'D={:5.2f}: Finished {} / {}'.format(input.noise, i+1, len(data))
    return vars
def Energy_data_calc(data, input):
    box_x = input.box_x
    rstar = input.rstar
    rstarsq = np.square(input.rstar)
    energies = []
    for frame in data:
        energy = 0
        x_pos = frame[:,0]
        y_pos = frame[:,1]
        theta = frame[:,2]
        coords = np.vstack((x_pos,y_pos)).T
        for i in range(input.npart-1):
            rest_coords = coords[i+1:]
            rest_thetas = theta[i+1:]
            dists = distanceSq(coords[i], rest_coords, box_x)
            theta_diff = theta[i] - rest_thetas[dists < rstarsq]
            energy += np.cos(theta_diff).sum()
        energy /= input.npart
        energies.append(energy)
        if len(energies) % 50 == 0:
            print printtime() + 'Finished {} / {} in D={:.2f}'.format(len(energies), len(data), input.noise)
    return energies   

def pickle_save(filename, data):
    f = open(filename, 'wb')
    pickle.dump(data, f)
    f.close()
    
def pickle_load(filename):
    f = open(filename, 'rb')
    ret = pickle.load(f)
    f.close()
    return ret

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
def printtime():
   return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S> ')

class hours_minutes_secondsInfo():
  def __init__(self, hours, minutes, seconds):
     self.hours = hours
     self.minutes = minutes
     self.seconds = seconds

def hours_minutes_seconds(elapsed_time):
    rest, seconds = divmod(elapsed_time.total_seconds(), 60)
    hours, minutes = divmod(rest,60)
    return hours_minutes_secondsInfo(int(hours), int(minutes), int(seconds))
def remove_inputs_and_errors():  
    try:
        for fl in glob.glob('errorlog*') + glob.glob('input*'):
            os.remove(fl)
    except OSError:
        pass

def remove_output_dirs(prefix):
    if not (prefix == 'output'):
        print 'just what are you trying to do?'
        return
    if len([f for f in os.listdir('.') if 'output' in f]) > 0:
        if not query_yes_no('Output data will be deleted! Are you sure about this?', default="no"):
            sys.exit()
    for subdirs, dirs, files in os.walk('.'):
        for dir in dirs:
            if (dir.startswith(prefix)):
                print printtime() + "Deleting directory {0}".format(dir)
                shutil.rmtree(dir)
                
def read_obstacle_file(output_dir):
    obFiles = [name for name in os.listdir(output_dir) if 'obstacle' in name]
    if len(obFiles) == 0:
        output_dir = output_dir + '\\' + 'longtime'
        obFiles = [name for name in os.listdir(output_dir) if 'obstacle' in name]
    if len(obFiles) == 0:
        raise FileNotFoundError('Where is the obstacle file?')
    for filename in obFiles:
        fullname = output_dir + '\\' + filename
        with open(fullname, "rb") as f:
            #print 'reading from {0}'.format(filename)
            dt = np.dtype([('x','float64'),('y','float64')])
            obstacles = np.fromfile(f, dtype=dt) 
            #print 'reading {0} lines from {1}'.format(len(frame), fullname)
            xOb, yOb = obstacles['x'],obstacles['y']
    return [xOb, yOb]

class readDirInfo():
  def __init__(self,data, xOb, yOb, times):
     self.data = data
     self.xOb = xOb
     self.yOb = yOb
     self.times = times
def last_time_in_file(f, input):
    dt_time = np.dtype([('time','float64')])
    iter_size = 8*(3*input.npart + 1)
    f.seek(0, os.SEEK_END)
    filesize = f.tell()
    #print f, filesize, iter_size
    f.seek(filesize - iter_size)
    time_list = np.fromfile(f, count = 1, dtype=dt_time)
    time = time_list[0][0]
    f.seek(0)
    return time
def last_time_in_direc(output_dir, input=None):
    if input==None:
        print 'reading input from: ',output_dir
        input = read_input_from_dir(output_dir)
    datafiles = [join(output_dir,'longtime',name) for name in os.listdir(join(output_dir,'longtime')) if name.startswith('log')]
    last_times = [last_time_in_file(open(file,'rb'), input) for file in datafiles]
    return np.max(last_times) if len(last_times) > 0 else 0
def read_data_from_dir(output_dir, start=0, frames=None, max_mem_size_mb=None, start_time=None, end_time=None,**kwargs):
    start_time_jack = datetime.datetime.now()
    #for key, value in kwargs.iteritems():
    #    print "%s = %s" % (key, value)
    input_file = kwargs.get('input')
    if input_file == None:
        input = read_input_from_dir(output_dir)
    else:
        input = read_input(input_file)
    # find simulation results files (from 'longtime')
    datafiles = [name for name in os.listdir(output_dir) if name.startswith('log')]
    if len(datafiles) == 0:
        output_dir = join(output_dir ,'longtime')
        datafiles = [name for name in os.listdir(output_dir) if name.startswith('log')]
        
    npart = input.npart
    iter_size = 8*(3*npart + 1)
    file_sizes = [os.path.getsize(join(output_dir , f)) for f in datafiles]
    niter = np.sum(file_sizes) / iter_size
    # maybe old vicsek style without timestamps
    for f, size in zip (datafiles,file_sizes):
        if size % iter_size != 0:
            print 'error in directory: ', output_dir
            print 'filename: ', f
            print '{} % {} = {}'.format(size,iter_size,size%iter_size)
            raise Exception('bad file format')
            
    
    if np.sum(file_sizes) % iter_size != 0:
        old_iter_size = 8*3*npart
        if np.sum(file_sizes) % old_iter_size == 0:
            print 'yuck! old vicsek model'
        else:
            print 'hmm. exotic', output_dir
            print np.sum(file_sizes)
            print iter_size
            print np.sum(file_sizes) % iter_size
        raise Exception('bad file format')
    # first option: "start" is given as a fraction
    if isinstance(start, float):
        first_iter = int(start*niter)
    # second: specific frame to start on
    elif isinstance(start, int):
        first_iter = start
    else:
        raise ValueError("start must be an int or a float")
   
    total_iter = niter - first_iter
    # limit number of frames if "frames" given as parameter
    if frames != None:
        total_iter = np.min([total_iter, frames])
    # read less frames if we're crossing the maximum memory usage
    if (max_mem_size_mb != None):
        max_iter = int(float(max_mem_size_mb)*(1024*1024) / iter_size)
        total_iter = np.min([total_iter, max_iter])
    #print printtime() + 'Reading {} frames out of {}, starting {}'.format(total_iter, niter, first_iter)
    i = 0
    sum = 0
    for s in file_sizes:
        sum += s
        if sum > first_iter * iter_size:
            sum -= s
            offset = (first_iter * iter_size) - sum
            offset += 0     # because of the <elapsed_time> in the beginning
            break
        i += 1
    #print 'starting with {}'.format(datafiles[i])
    dt_cell = np.dtype([('x','float64'),('y','float64'),('theta','float64')])
    dt_time = np.dtype([('time','float64')])
    iterations_read = 0
    data = []
    times = []
    for filename in datafiles[i:]:

        fullname = output_dir + '\\' + filename
        with open(fullname, "rb") as f:
            # if start_time given, find the right file first
            if start_time != None:
                last_time = last_time_in_file(f,input)
                if last_time < start_time:
                    continue
            elif datafiles.index(filename) == i:
                f.seek(offset, os.SEEK_SET)  # jump ahead in first file
                #print 'reading {} as first file. offset = {} filesize = {}'.format(filename,offset, file_sizes[i])
            #hourminsec = hours_minutes_seconds(datetime.datetime.now() - start_time_jack)
            #print printtime() + "Found correct file after: {}:{}:{}".format(hourminsec.hours, hourminsec.minutes, hourminsec.seconds)
            #print printtime() + 'reading from {0}'.format(fullname)
            while True:
                hourminsec = hours_minutes_seconds(datetime.datetime.now() - start_time_jack)
                #print printtime() + "Started reading frame after: {}:{}:{}".format(hourminsec.hours, hourminsec.minutes, hourminsec.seconds)
                time_list = np.fromfile(f, count = 1, dtype=dt_time)
                hourminsec = hours_minutes_seconds(datetime.datetime.now() - start_time_jack)
                #print printtime() + "Finished reading frame after: {}:{}:{}".format(hourminsec.hours, hourminsec.minutes, hourminsec.seconds)
                if len(time_list) == 0 or iterations_read == total_iter:
                    break
                time = time_list[0][0]
                ###
                frame = np.fromfile(f, count = npart*3, dtype=np.float64)
                frame= frame.reshape((npart,3))
                ##
                if start_time != None and time < start_time:
                    continue
                if end_time != None and time > end_time:
                    break
                times.append(time)
                data.append(frame)
                iterations_read += 1
    xOb, yOb = read_obstacle_file(output_dir)
    return readDirInfo(data, xOb, yOb, times)

def read_dir_wrapped(output_dir, npart, start, finish, nthframe):
    dt_cell = np.dtype([('x','float64'),('y','float64'),('theta','float64')])
    dt_time = np.dtype([('time','float64')])
    iter_size = 8*(3*npart + 1)
    xOb,yOb = [],[]
    # read simulation results
    data = []
    times = []
    datafiles = [name for name in os.listdir(output_dir) if name.startswith('log')]
    if len(datafiles) == 0:
        output_dir = output_dir + '\\' + 'longtime'
        datafiles = [name for name in os.listdir(output_dir) if name.startswith('log')]
    file_sizes = [os.path.getsize(output_dir + '\\' +f) for f in datafiles]
    niter = np.sum(file_sizes) / iter_size
    first_iter = int(start*niter)
    total_iter = int(finish*niter) - first_iter
    print printtime() + 'Reading [{3}:{4}] out of {2} files from {0}, npart = {1}'.format(output_dir, npart, len(datafiles), start, finish)
    i = 0
    sum = 0
    for s in file_sizes:
        sum += s
        if sum > first_iter * iter_size:
            sum -= s
            offset = (first_iter * iter_size) - sum
            offset += 0     # because of the <elapsed_time> in the beginning
            break
        i += 1
    #print 'starting with {}'.format(datafiles[i])
    iterations_read = 0
    for filename in datafiles[i:]:
        fullname = output_dir + '\\' + filename
        with open(fullname, "rb") as f:
            if datafiles.index(filename) == i:
                f.seek(offset, os.SEEK_SET)  # jump ahead in first file
            #print printtime() + 'reading from {0}'.format(fullname)
            while True:
                time_list = np.fromfile(f, count = 1, dtype=dt_time)
                if len(time_list) == 0:
                    break
                time = time_list[0][0]
                times.append(time)
                frame = np.fromfile(f, count = npart, dtype=dt_cell)
                data.append(frame)
                iterations_read += 1
    xOb, yOb = read_obstacle_file(output_dir)
    return readDirInfo(data, xOb, yOb, times)

def read_dir(output_dir):
    return read_dir_wrapped(output_dir, read_input_from_dir(output_dir).npart, 0.0, 1.0, 1)

class InputInfo():
  def __init__(self,npart, box_x, noise, neighbor_norm, dt, rstar, align, simtime, hyster=False, hyster_noise=None):
     self.npart = npart
     self.box_x = box_x
     self.noise = noise 
     self.neighbor_norm = neighbor_norm
     self.dt = dt
     self.rstar = rstar
     self.align = align
     self.simtime = simtime
     self.hyster = hyster
     self.hyster_noise = hyster_noise

def read_input(input_file):
    hyster = False
    hyster_noise = None
    #print input_file
    with open(input_file, "rb") as f:
        content = f.readlines()
        for line in content:
            fields = line.strip().split()
            if len(fields) == 0:
                continue
            if 'NPART' == fields[0]:
                npart = int(fields[fields.index('=')+1])
            elif 'BOX_X' == fields[0]:
                box_x = float(fields[fields.index('=')+1])
            elif 'NOISE' == fields[0]:
                noise = float(fields[fields.index('=')+1])
            elif 'NORMALIZE_NEIGHBORS' == fields[0]:
                neighbor_norm = fields[fields.index('=')+1].strip().lower() == 'true'
            elif 'DT' == fields[0]:
                dt = float(fields[fields.index('=')+1])
            elif 'RSTAR' == fields[0]:
                rstar = float(fields[fields.index('=')+1])
            elif 'ALIGN' == fields[0]:
                align = float(fields[fields.index('=')+1])
            elif 'SIMTIME' == fields[0]:
                simtime = float(fields[fields.index('=')+1])
            elif 'HYSTER' == fields[0]:
                hyster = fields[fields.index('=')+1].strip().lower() == 'true'
            elif 'HYSTER_NOISE' == fields[0]:
                hyster_noise = float(fields[fields.index('=')+1])
    return InputInfo(npart, box_x, noise, neighbor_norm, dt, rstar, align, simtime,
        hyster=hyster, hyster_noise=hyster_noise)

def find_input_file(dir):
    input_file = 'ERROR'
    for f in os.listdir(dir):
        if 'input' in f:
            input_file = dir + '\\' + f
            break
    if input_file == 'ERROR':
        if dir.endswith('longtime'):
            dir = dir[:-len(r'\longtime')]
            for f in os.listdir(dir):
                if 'input' in f:
                    input_file = join(dir,f)
                    return read_input(input_file)
        raise FileNotFoundError('Error: No input file found!')
    return input_file
def read_input_from_dir(dir):
    input_file = 'ERROR'
    if not isdir(dir):
        # try to look for it in the current directory
        print 'warning: ', dir ,' does not exist'
        outdirname = dir.split(os.sep)[-1]
        print outdirname, isdir(outdirname)
        if isdir(outdirname):
            dir = outdirname
        else:
            raise FileNotFoundError('Error: directory not found!')
    for f in os.listdir(dir):
        if 'input' in f:
            input_file = join(dir , f)
            break
    if input_file == 'ERROR':
        if dir.endswith('longtime'):
            dir = dir[:-len(r'\longtime')]
            for f in os.listdir(dir):
                if 'input' in f:
                    input_file = join(dir,f)
                    return read_input(input_file)
        raise FileNotFoundError('Error: No input file found!')
    return read_input(input_file)

class snapShot():
  def __init__(self,x,y, theta, xOb, yOb, time, input, pol=None):
     self.x = x
     self.y = y
     self.theta = theta
     self.xOb = xOb
     self.yOb = yOb
     self.time = time
     self.input = input
     self.pol = pol

def snapshot(output_dir, time=None, compass = True):
    if not isdir(output_dir):
        # try to look for it in the current directory
        outdirname = output_dir.split(os.sep)[-1]
        if isdir(outdirname):
            output_dir = outdirname
        else:
            raise FileNotFoundError('Error: directory not found!')
    input = read_input_from_dir(output_dir)
    if time != None:
        sim = read_data_from_dir(output_dir, frames=1, start_time=time)
        if len(sim.data) == 0:
            print "error!"
            return
    else:
        start_percent = 0.999
        sim = []
        while True:
            sim = read_dir_wrapped(output_dir, input.npart, start_percent, 1.0, 1)
            if len(sim.data) > 0:
                break
            start_percent /= 10.0
            if start_percent < 0.1:
                raise IOError('cannot locate file in {}'.format(output_dir))
            print start_percent
    last_frame = sim.data[-1]
    x = np.array([val[0] for val in last_frame])
    y = np.array([val[1] for val in last_frame])
    theta = np.array([val[2] for val in last_frame])
    xOb, yOb = read_obstacle_file(output_dir)
    time = sim.times[-1]
    print 'real time is ', time
    if compass == True:
        pol = np.mean(np.exp(1j * theta))
    else:
        pol = None
    return snapShot(x,y,theta,xOb,yOb,time, input, pol = pol)

def getDensity(sim):
    return getDensityInput(sim.input)
def getDensityInput(input):
    return input.npart / np.square(input.box_x)
    
# A smoothing function I copied from http://wiki.scipy.org/Cookbook/SavitzkyGolay
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
    
"""
Useful plotting snippet
        
        fig = plt.figure(1)
        ax1 = fig.add_subplot(221)
        ax1.imshow( dens)
        
        ax2 = fig.add_subplot(222)
        x = np.linspace(0,input.box_x,res)
        y = np.mean(dens, axis=0)
        ax2.plot(x,y)
        ax3 = fig.add_subplot(223)
        snap = snapShot(x_pos, y_pos, frame[:,2], [],[], 6, input)
        plot_snap(snap, save=False, ax=ax3)        
        ax4 = fig.add_subplot(224,projection='3d')
        surf = ax4.plot_surface(X, Y, psf2D, rstride=1, cstride=1, cmap=cm.coolwarm,
            linewidth=0, antialiased=False)

        plt.show()
        sys.exit()
        """
def plot_snap_fake(snap, show=True, save=False, draw=False, ax=None, arrows=True, color='black', 
        figure=None, backgroundcolor='white', clusters=False, compass_color='r', cluster_color = 'r', colors=None):
    input = snap.input
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    density = snap.input.npart / np.square(snap.input.box_x)
    filename = 'snapshot_{}_dens{:.2f}_noise{}_npart{}_time{}'.format(model, density,input.noise,input.npart,int(snap.time)).replace(' ','').replace('.','_') + '.png'
    xedge = 0.01

    if ax == None:
        fig = plt.figure(num=None, figsize=(10, 10), dpi=100, facecolor=backgroundcolor, edgecolor='k')
        ax1 = fig.add_subplot(111)
        fig.suptitle('{} Model, t = {:4.0f}'.format(model, snap.time), fontsize=14, fontweight='bold', color=color)#, y=0.9)
        #fig.suptitle('L = {}'.format(input.box_x), fontsize=14, fontweight='bold')
        ax1.set_title(r'D = {}, {} particles, L = {:.0f}, $\rho$ = {:g}'.format(input.noise, input.npart, input.box_x, density), fontsize=14, color=color)
        plt.subplots_adjust(bottom=0.01, left=xedge, right=1-xedge, top=0.91)
    else:
        ax1 = ax
        fig = figure
        plt.subplots_adjust(bottom=0.01, left=xedge, right=1-xedge, top=0.99)
    #fig.patch.set_facecolor('k')
    box_x = input.box_x
    x,y,theta,xOb,yOb = snap.x, snap.y, snap.theta, snap.xOb, snap.yOb
    # draw box - vertical lines
    border_color = color
    plt.plot([0, 0], [0, 1], color=border_color, linestyle='--', linewidth=2)
    plt.plot([1, 1], [0, 1], color=border_color, linestyle='--', linewidth=2)
    # draw box - horizontal lines
    plt.plot([0, 1], [0, 0], color=border_color, linestyle='--', linewidth=4)
    plt.plot([0, 1], [1, 1], color=border_color, linestyle='--', linewidth=4)
    #Q = plt.quiver(x, y, np.cos(theta),np.sin(theta), width=0.1, scale=7,units='inches', 
    #    color = 'black', headlength=5, headwidth =3, facecolor='black',  edgecolors=('black'))
    factor = 0.1
    """
    for x_val, y_val, theta_val,  in zip(x/box_x,y/box_x,theta):
        ax1.add_artist(plt.Circle((x_val,y_val), 1.0/box_x,ec="k", color='r', alpha = 0.5)) """
    if colors != None:
        color = colors
    if (arrows):
        ax1.quiver(x/box_x, y/box_x, np.cos(theta)*factor, np.sin(theta)*factor, color=color)
    else:
        ax1.scatter(x/box_x, y/box_x, color=color, marker="o", s=0.3, alpha=0.8)
    #ax.scatter(xOb, yOb, s=40, c='m', marker="o", units='x')
    for X, Y in zip(xOb, yOb):
        ax1.add_artist(plt.Circle((X/box_x,Y/box_x), 1/100,ec="none", color='m')) 
    
    if clusters:
        x = np.array(x)
        y = np.array(y)
        clusters = list(cluster_analysis(x,y, input))
        big_clusters = [cluster for cluster in clusters if len(cluster) > 10]
        for cluster in big_clusters:
            diam = cluster_diameter(cluster, x, y, input)
            (x_com,y_com) = center_of_mass(cluster, x, y, input)
            ax1.add_artist(plt.Circle((x_com/box_x,y_com/box_x),0.6*diam/box_x,alpha = 0.5,color=cluster_color))
        
    
    if snap.pol != None:
        off = 0.75
        size = 0.15
        factor = 0.25
        pol = snap.pol
        im = plt.imread(r'G:\Alon\Workspace\01042015\compass\compass.png')
        implot = plt.imshow(im, extent=(off, off+size, off, off+size))
        norm = np.abs(pol)
        pol = pol / norm
        pol *= factor
        
        rect_xy = (0.91, off)
        rect_width = 0.02
        rect_height = size
        ax1.add_artist(plt.Rectangle(rect_xy, rect_width, rect_height, ec="none", color=compass_color, alpha=0.1)) 
        ax1.add_artist(plt.Rectangle(rect_xy, rect_width, rect_height*norm, ec="none", color=compass_color)) 
        plt.arrow(off+size/2.0,off+size/2.0, np.real(pol)*factor, np.imag(pol)*factor, 
            width=0.003, length_includes_head=True, head_width=0.02, color=compass_color, alpha=0.9)
        
    ax1.set_autoscaley_on(False)
    ax1.set_autoscalex_on(False)
    plt.axis('off')
    plt.axis('equal')
    if show and ax==None:
        def onclick(event):
            kids = ax1.get_children()
            for kid in kids:
                if isinstance(kid, Circle):
                    kid.set(Visible=not kid.get_visible())
            plt.draw()
        fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
    if save:
        fig.savefig(filename, dpi=fig.dpi, facecolor = fig.get_facecolor(), edgecolor = fig.get_facecolor())
    if draw:
        plt.draw()
    #plt.close(fig)
    return True        
    
### cluster functions
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
    size = len(cluster)
    if size <2:
        return 0.0
    coords = np.vstack((x,y)).T
    dists = []
    for part1 in range(0,size-1):
        rest = cluster[part1+1:size]
        dist = distanceSq(coords[cluster[part1]], coords[rest], np.array([input.box_x, input.box_x]))
        dists.append(np.max(dist))
    return np.sqrt(np.max(dists))
def cluster_gyration_radius(cluster, x, y, input):
    size = len(cluster)
    if size <2:
        return 0.0
    reduced_x, reduced_y = x[cluster], y[cluster]
    coords = np.vstack((reduced_x,reduced_y)).T
    c_o_m_coords = np.array(center_of_mass(cluster, x, y, input)).T
    dists = distanceSq(c_o_m_coords, coords, np.array([input.box_x, input.box_x]))
    gyration_radius = np.sum(dists) / size
    #print 'debug R^2 = ', gyration_radius
    return gyration_radius
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
def cluster_polarization(cluster, theta):
    angles = theta[cluster]
    pol = np.mean(np.exp(1j * angles))
    return pol
