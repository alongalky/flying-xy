import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import poisson
import random
from Tkinter import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from collections import defaultdict
from utils import *
#from clusters import *

def plot_hyster(data, show=True, smooth=True, ylabel=r'$<\Pi>$'):
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    model = 'Vicsek' if data[0].input.neighbor_norm == True else 'Flying XY'
    txt = 'Polarization'
    title = '{}, {} Model'.format(txt, model)
    #ax.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel('$D$', fontsize=22, fontweight='bold', labelpad=-2)
    ax1.set_ylabel(ylabel, fontsize=22, fontweight='bold')
    for pol in data:
        x = noiseFunc(pol.input, pol.times)
        y = pol.values
        #def savitzky_golay(y, window_size, order, deriv=0, rate=1):
        cmap = plt.get_cmap('jet')
        col = [cmap(i) for i in np.linspace(0,1,len(x))]
        if smooth:
            y = savitzky_golay(y, 201, 4)
        ax1.plot(x,y,label='t = {}'.format(int(pol.input.simtime)), #marker='o', markersize = 0,
            picker=True)#color=col,picker = False, edgecolor=None)
    #plt.legend(loc='best')
    filename = '{}_npart{}_global_pol_vs_noise_hyster.png'.format(model.replace(' ',''), data[0].input.npart)
    fig.savefig(filename, dpi=fig.dpi)    
    def onpick4(event):
        ind = event.ind
        thisline = event.artist
        line_index = ax1.lines.index(thisline)
        sim = data[line_index]
        time = np.take(sim.times, ind)[0]
        snap = snapshot(sim.output_dir, time=time)
        plot_snap(snap, save=True)

    fig.canvas.mpl_connect('pick_event', onpick4)        
    if show:
        plt.show()


def seperate_by_npart(energies_flattened):
    dict = defaultdict(list)
    for energies in energies_flattened:
        npart = energies.input.npart
        dict[npart].append(energies)
    energies_list = []
    for npart, sims in dict.iteritems():
        sims = sorted(sims, key=lambda sim: sim.input.noise)
        energies_list.append((npart, sims))
    # sort by density
    energies_list = sorted(energies_list, key=lambda val: val[0])
    ret = [val[1] for val in energies_list]
    return ret
    txt = 'Polarization'
    title = '{} vs time'.format(txt)
    print printtime() + 'Plotting {}'.format(title)
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax2 = fig.add_subplot(111)
    ax2.set_title(title, fontsize=14, fontweight='bold')
    ax2.set_xlabel('Time', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Polarization', fontsize=14, fontweight='bold')
        
    #selected_Pols = random.sample(Pols, np.min([4,len(Pols)]))
    #selected_Pols = sorted(selected_Pols, key=lambda pol: pol.input.noise)
    selected_Pols = sorted(Pols, key=lambda pol: pol.input.npart / np.square(pol.input.box_x))
    #selected_inputs = [input for input in flat_inputs if input.npart == 2048 and input.noise==0.335]
    #
    plot_output_dirs2 = []
    cmap = plt.get_cmap('hsv')
    colors = [cmap(i) for i in np.linspace(0,1,len(selected_Pols))]
    for color, pol in zip(colors,selected_Pols):
        input = pol.input
        npart = input.npart
        density = npart / np.square(input.box_x)
        model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
        max_set_size = 500
        nthframe = np.max([int(len(pol.times) / max_set_size), 1])
        x = pol.times[::nthframe]
        y = pol.values[::nthframe]
        #color=colors[len(plot_output_dirs2)]
        plot_output_dirs2.append(pol.output_dir)
        ax2.plot(x,y, '-', linewidth=2, label=r'dens={:.1f},D={:.1f}'.format(density,input.noise),# marker='o', markersize = 4,
            color=color, picker = True)
    def onpick4(event):
        ind = event.ind
        thisline = event.artist
        line_index = ax2.lines.index(thisline)
        time = np.take(thisline.get_xdata(), ind)[0]
        if event.mouseevent.button == 3:
            print 'select second point for movie making (doesn\'t really work)'
        else:
            snap = snapshot(plot_output_dirs2[line_index], time=time)
            plot_snap(snap, save=True)

    fig.canvas.mpl_connect('pick_event', onpick4)
    filename = 'Polarization_vs_time.png'
    ax2.legend(prop={'size':12},loc='best')
    fig.savefig(filename, dpi=fig.dpi)
    if show:
        plt.show()
def plot_property_vs_time(Props, title='Property Name', ylabel=None, show=True, smooth=False, cmap=plt.cm.hsv,
        log_x = False, log_y = False, legend_loc='best', save=True, linewidth=1, legend=True, fitpairs=False,
        ax=None, figure=None, clusters=False, compass=True):
    txt = title
    tit = '{} vs time'.format(txt)
    print printtime() + 'Plotting {}'.format(tit)
    if figure == None or ax==None:
        fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
        ax1 = fig.add_subplot(111)
    else:
        fig = figure
        ax1 = ax
    ax1.set_title(tit, fontsize=14, fontweight='bold')
    ax1.set_xlabel('Time', fontsize=14, fontweight='bold')
    if ylabel != None:
        ax1.set_ylabel(ylabel, fontsize=14, fontweight='bold')
    selected_Pols = Props
    selected_Pols = sorted(Props, key=lambda pol: pol.input.npart / np.square(pol.input.box_x))
    plot_output_dirs2 = []
    #cmap = plt.get_cmap('Dark2')
    if fitpairs:
        colors = [cmap(i) for i in np.linspace(0,1,len(selected_Pols)/2)]
        colors = [x for pair in zip(colors,colors) for x in pair]
    else:
        colors = [cmap(i) for i in np.linspace(0,1,len(selected_Pols))]
    if len(selected_Pols) < 4:
        colors = ['r','b','g','m']
    for color, pol in zip(colors,selected_Pols):
        input = pol.input
        npart = input.npart
        density = npart / np.square(input.box_x)
        model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
        # max_set_size = 500
        # nthframe = np.max([int(len(pol.times) / max_set_size), 1])
        nthframe =1
        x = pol.times[::nthframe]
        y = pol.values[::nthframe]
        print pol.prop_name
        if pol.prop_name == 'momentum':
            print 'mMomemntum'
            y = np.abs(y)
        plot_output_dirs2.append(pol.output_dir)
        
        #ax1.plot(x,y, '-', linewidth=0, label=r'D={}'.format(input.noise), marker='o', markersize = 4,
        #    color=color, picker = True)
        if smooth and len(x) > 500:
            window_size = np.max([1001, int(len(x)/4)*2 + 1])
            y = savitzky_golay(y, window_size, 4)

        ax1.plot(x,y, '-', linewidth=linewidth, label=r'$A={},D={},\rho={:.2f},N={}$'.format(input.align,input.noise, density,input.npart),# marker='o', markersize = 4,
            color=color, picker = True)
    def onpick4(event):
        ind = event.ind
        thisline = event.artist
        line_index = ax1.lines.index(thisline)
        time = np.take(thisline.get_xdata(), ind)[0]
        if event.mouseevent.button == 3:
            print 'select second point for movie making (doesn\'t really work)'
        else:
            snap = snapshot(plot_output_dirs2[line_index], time=time, compass=compass)
            plot_snap(snap, save=True, clusters=clusters)
    
    if (log_x):
        ax1.set_xscale('log')
    if (log_y):
        ax1.set_yscale('log')

    fig.canvas.mpl_connect('pick_event', onpick4)
    if legend:
        if fitpairs:
            print 'Fit pairs debug'
            handles, labels = ax1.get_legend_handles_labels()
            ax1.legend(handles[::2], labels[::2], loc=legend_loc)
        else:
            ax1.legend(prop={'size':12},loc=legend_loc)
    if save:
        input = Props[0].input
        filename = '{}_align_{:.2f}_noise_{:.3f}_dens_{:.3f}_vs_time.png'.format(txt,
            input.align, input.noise, getDensityInput(input)).replace(" ","").replace('.','_').replace(',','_')
        fig.savefig(filename, dpi=fig.dpi)
    if show:
        plt.show()    
def plot_mProperty_vs_dens(mProp, title=None, log_x = False, log_y = False, 
        show=True, variance=False, ylabel=None, order_param=None, vlines=[], clusters=False, compass=True, showtitle=True, linewidth=0):
    input = mProp[0].input
    if title==None:
        title = mProp[0].name
    txt = title
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    title = '{}'.format(txt)#, model,input.noise)
    if showtitle:
        ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel(r'$\rho$', fontsize=22, labelpad=-2)#fontweight='bold')
    ax1.set_ylabel(txt if ylabel==None else ylabel, fontsize=22)#, fontweight='bold')
    
    x = [getDensityInput(en.input) for en in mProp]
    if variance:
        y = [en.var for en in mProp]
        y_dat = np.array(y)
        x_dat = np.array(x)
        def gaussian(x, a,b,c):
            return a*np.exp(-(b*(x**2)) / np.sqrt(2*c))
        def pois(x, a,b):
            return a*poisson.pmf(x,b)
        def cauch(x, a,b):
            return a/(1+x_dat*b)
        func = gaussian
        popt, pcov = curve_fit(func, x_dat, y_dat)
        yfit = func(x_dat, *popt)
        #ax1.plot(x_dat,yfit, '--', linewidth=1, label=r'fit', marker='o', markersize = 4,
        #    color='r', picker = True)
            
        ###
        avg_x = np.dot(x_dat,y_dat)/np.sum(y_dat)
        #plt.vlines(avg_x, 0, np.max(y_dat), colors=u'k', linestyles=u'solid', label=u'')
        #ax1.plot((avg_x, 0), (avg_x, ax1.get_ylim()[1]), 'k-')
    else:
        #print mProp[0].name
        if mProp[0].name == 'mMomentum':
            ax1.set_ylabel(r'$|<W_l(\vec{r},t)>|/\rho$')
            y = [en.mValue/getDensityInput(en.input) for en in mProp]
            #y = [en.mValue for en in mProp]
        elif mProp[0].name == 'mDensVar':
            y = [en.mValue/getDensityInput(en.input) for en in mProp]
            ax1.set_ylabel(r'$<(\delta \rho)^2>/\rho$')
        elif mProp[0].name == 'mEnergy':
            y = [en.mValue/getDensityInput(en.input) for en in mProp]
            ax1.set_ylabel(r'$E_{\rm align}/\rho$')
        elif mProp[0].name == 'mAbsMomentum':
            y = [en.mValue/getDensityInput(en.input) for en in mProp]
            ax1.set_ylabel(r'$<|W_l(\vec{r},t)|>/\rho$')
        else:
            y = [en.mValue for en in mProp]
    ax1.plot(x,y, '--', linewidth=linewidth, label=r'npart = {}'.format(mProp[0].input.npart), marker='o', markersize = 8,
        color='b', picker = True)
    def onpick_polar(event):
        ind = event.ind
        if order_param == None:
            time = None
            if event.mouseevent.button==3:
                input = mProp[ind].input
                time = get_value_box(title='Drawing snapshot for D={}, maximum time = {:.2f}'.format(input.noise, mProp[ind].simlen), 
                    prompt='Choose time:')
                try:
                    time = float(time)
                    if time == 0:
                        return
                except:
                    time = None
            snap = snapshot(mProp[ind].output_dir, time=time, compass=False)
            plot_snap(snap, save=True, clusters=clusters)
        else:
            output_dir = mProp[ind].output_dir
            prop = calc_property(output_dir,order_param.pickle_name,order_param.calc_data_fun, read=False)
            plot_property_vs_time([prop], title=order_param.name, ylabel=order_param.name,save=True, 
                    legend=True, show=True, smooth=order_param.smooth, clusters=clusters, compass=compass)
    for point in vlines:
        plt.vlines(point, 0, ax1.get_ylim()[1], linestyles='--', label=u'')

    if (log_x):
        ax1.set_xscale('log')
    if (log_y):
        ax1.set_yscale('log')
    fig.canvas.mpl_connect('pick_event', onpick_polar)
    filename = '{}_vs_dens.png'.format(txt).replace(" ","").replace('.','_').replace(',','_')
    fig.savefig(filename, dpi=fig.dpi)
    if show:
        plt.show()
    return
def plot_mProperty_vs_noise(mProp, title='Property Name', log_x = False, log_y = False, 
        show=True, variance=False, ylabel=None):
    input = mProp[0].input
    txt = title
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    title = '{}, {}, $\\rho$={}'.format(txt, model,getDensityInput(input))
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel('Noise', fontsize=14, fontweight='bold')
    ax1.set_ylabel(txt if ylabel==None else ylabel, fontsize=14, fontweight='bold')
    
    x = [en.input.noise for en in mProp]
    if variance:
        y = [en.var for en in mProp]
    else:
        y = [en.mValue for en in mProp]
    ax1.plot(x,y, '--', linewidth=1, label=r'npart = {}'.format(mProp[0].input.npart), marker='o', markersize = 4,
        color='b', picker = True)
    def onpick_polar(event):
        ind = event.ind
        time = None
        if event.mouseevent.button==3:
            input = mProp[ind].input
            time = get_value_box(title='Drawing snapshot for D={}, maximum time = {:.2f}'.format(input.noise, mProp[ind].simlen), 
                prompt='Choose time:')
            try:
                time = float(time)
                if time == 0:
                    return
            except:
                time = None
        snap = snapshot(mProp[ind].output_dir, time=time)
        plot_snap(snap, save=True)

    if (log_x):
        ax1.set_xscale('log')
    if (log_y):
        ax1.set_yscale('log')
    fig.canvas.mpl_connect('pick_event', onpick_polar)
    filename = '{}_vs_noise.png'.format(txt).replace(" ","")
    fig.savefig(filename, dpi=fig.dpi)
    if show:
        plt.show()
    return
def plot_mProperty(mProp, title='Property Name', log_x = False, log_y = False, show=True):
    input = mProp[0].input
    txt = title
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    title = '{}, {} Model, {} particles'.format(txt, model, input.npart)
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel('Noise', fontsize=14, fontweight='bold')
    ax1.set_ylabel(txt, fontsize=14, fontweight='bold')
    
    x = [en.input.noise for en in mProp]
    y = [en.mValue for en in mProp]
    ax1.plot(x,y, '--', linewidth=1, label=r'npart = {}'.format(mProp[0].input.npart), marker='o', markersize = 4,
        color='b', picker = True)
    def onpick_polar(event):
        ind = event.ind
        time = None
        if event.mouseevent.button==3:
            input = mProp[ind].input
            time = get_value_box(title='Drawing snapshot for D={}, maximum time = {:.2f}'.format(input.noise, mProp[ind].simlen), 
                prompt='Choose time:')
            try:
                time = float(time)
                if time == 0:
                    return
            except:
                time = None
        snap = snapshot(mProp[ind].output_dir, time=time)
        plot_snap(snap, save=True)

    if (log_x):
        ax1.set_xscale('log')
    if (log_y):
        ax1.set_yscale('log')
    fig.canvas.mpl_connect('pick_event', onpick_polar)
    filename = '{}_vs_noise.png'.format(txt).replace(" ","")
    fig.savefig(filename, dpi=fig.dpi)
    if show:
        plt.show()
    return
def plot_mean_pol(mPols, log_x = False, log_y = False, show=True):
    txt = 'Polarization'
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    model = 'Vicsek' if mPols[0].input.neighbor_norm == True else 'Flying XY'
    title = '{}, {} Model'.format(txt, model)
    #ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel('$D$', fontsize=22, fontweight='bold', labelpad=-2)
    ax1.set_ylabel(r'$<\Pi>$', fontsize=22, fontweight='bold')
    
    mPols_list = seperate_by_npart(mPols)
    #ax1.set_xlim([0,0.6])
    cmap = plt.get_cmap('jet')
    colors = [cmap(i) for i in np.linspace(0,0.9,len(mPols_list))]
    #
    #inset = True
    inset = False
    if inset:
        colors = ['b']*10
        axins = zoomed_inset_axes(ax1, 7, loc=1, bbox_to_anchor=(0.65, 0.9), 
                     bbox_transform=ax1.figure.transFigure)
        axins.axis([0.23, 0.5, 0, 0.1])
        axins.get_xaxis().set_visible(False)
        axins.get_yaxis().set_visible(False)
        mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    linestyles = ['-', '--','-.',':']
    for n, mPols in enumerate(mPols_list):
        dirs = [val.output_dir for val in mPols] # dirs
        x = [val.input.noise for val in mPols] # noise
        y = [val.mValue for val in mPols] # mPol
        if inset:
            ax1.plot(x,y, linestyle=linestyles[n%len(linestyles)], linewidth=2.5, label=r'N = {}'.format(mPols[0].input.npart), 
                color=colors[n], picker = True)
        else:
            ax1.plot(x,y, linestyle=linestyles[n%len(linestyles)], linewidth=0, label=r'N = {}'.format(mPols[0].input.npart), 
                marker='o', markersize = 8,
                color=colors[n], picker = True)
        if inset:
            axins.plot(x,y, linestyle=linestyles[n%len(linestyles)], linewidth=3, label=r'N = {}'.format(mPols[0].input.npart), #marker='o', markersize = 4,
                color=colors[n],picker = False)
    def onpick_polar(event):
        ind = event.ind
        thisline = event.artist
        line_index = ax1.lines.index(thisline)
        #print 'Plotting ', plot_output_dirs[line_index][ind]
        time = None
        direc = mPols_list[line_index][ind].output_dir
        if not isdir (direc):
            direc = direc.split(os.sep)[-1]
        if not isdir (direc):
            raise FileNotFoundError('Error: directory not found!')
        if event.mouseevent.button==2:
            input = read_input_from_dir(direc)
            pol = calc_property(direc,'Pols',Pol_data_calc, read=False)
            plot_property_vs_time([pol])
            return True
        if event.mouseevent.button==3:
            input = mPols_list[line_index][ind].input
            time = get_value_box(title='Drawing snapshot for D={}, maximum time = {:.2f}'.format(input.noise, mPols_list[line_index][ind].simlen), 
                prompt='Choose time:')
            try:
                time = float(time)
                if time == 0:
                    return
            except:
                time = None
        snap = snapshot(mPols_list[line_index][ind].output_dir, time=time)
        plot_snap(snap, save=True)#, color = 'w', backgroundcolor='k', clusters=True)

    if (log_x):
        ax1.set_xscale('log')
    if (log_y):
        ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.set_ylim([0,1])
    if inset:
        axins.xaxis.tick_top()
    
    fig.canvas.mpl_connect('pick_event', onpick_polar)
    filename = 'Polarization_vs_noise{}.png'.format(model)
    fig.savefig(filename, dpi=fig.dpi)
    if show:
        plt.show()
    return
def plot_binder(binders, log_x = False, log_y = False, xlim=None, show=True): 
    txt = 'Binder Cumulant'
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    model = 'Vicsek' if binders[0].input.neighbor_norm == True else 'Flying XY'
    title = '{}, {} Model'.format(txt, model)
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel('Noise', fontsize=14, fontweight='bold')
    ax1.set_ylabel(txt, fontsize=14, fontweight='bold')
    
    binders_list = seperate_by_npart(binders)
    #ax1.set_xlim([0,0.6])
    cmap = plt.get_cmap('Dark2')
    colors = [cmap(i) for i in np.linspace(0,1,len(binders_list))]
    for n, mPols in enumerate(binders_list):
        dirs = [val.output_dir for val in mPols] # dirs
        x = [val.input.noise for val in mPols] # noise
        y = [val.mValue for val in mPols] # mPol
        ax1.plot(x,y, '--', linewidth=1, label=r'npart = {}'.format(mPols[0].input.npart), marker='o', markersize = 4,
            color=colors[n], picker = True)
    def onpick_polar(event):
        ind = event.ind
        thisline = event.artist
        line_index = ax1.lines.index(thisline)
        #print 'Plotting ', plot_output_dirs[line_index][ind]
        time = None
        if event.mouseevent.button==3:
            input = binders_list[line_index][ind].input
            time = get_value_box(title='Drawing snapshot for D={}, npart={}'.format(input.noise, input.npart), prompt='Choose time:')
            try:
                time = float(time)
                print 'good! time = {}'.format(time)
                if time == 0:
                    return
            except:
                time = None
        snap = snapshot(binders_list[line_index][ind].output_dir, time=time)
        print 'Plotting t = ', snap.time
        plot_snap(snap, save=True, color = 'w', backgroundcolor='k')

    if (log_x):
        ax1.set_xscale('log')
    if (log_y):
        ax1.set_yscale('log')
    
    if xlim != None:
        ax1.set_xlim(xlim)

    fig.canvas.mpl_connect('pick_event', onpick_polar)
    filename = 'Binder_vs_noise.png'
    fig.savefig(filename, dpi=fig.dpi)
    
    if show:
        plt.show()
    return
"""
def plot_pdfs(vars=vars, pdfs=None, show=True):
    txt = 'PDF Variance'
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    model = 'Vicsek' if vars[0].input.neighbor_norm == True else 'Flying XY'
    title = '{}, {} Model'.format(txt, model)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('time', fontsize=14, fontweight='bold')
    ax.set_ylabel(txt, fontsize=14, fontweight='bold')
    plot_output_dirs2 = []
    for var in vars[:3]:
        print 'plotting ' , var.input.noise
        #plot_pdf_variance(var, ax=ax, show=False)
        ax.plot(var.times, var.values, label='D = {:.2f}'.format(var.input.noise), picker = True)
        if True:
            # power law fitting:
            fitfunc = lambda p, x: np.power(x, p[0])    # Target function
            errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
            p0 = [1.0] # Initial guess for the parameters
            p1, success = optimize.leastsq(errfunc, p0, args=(var.times, var.values))
            #time = linspace(Tx.min(), Tx.max(), 100)
            p1 = [0.8]
            print var.input.noise, ' alpha = ', p1[0]
            ax.plot(var.times, fitfunc(p1, var.times), "-", label='D={:.2f} fit'.format(var.input.noise)) # Plot of the fit
        plot_output_dirs2.append(var.output_dir)
    def onpick2(event):
        select_mode = False
        ind = event.ind
        thisline = event.artist
        line_index = ax.lines.index(thisline)
        #density = densities[line_index]
        time = np.take(thisline.get_xdata(), ind)[0]
        if event.mouseevent.button == 3:
            print 'select second point for movie making (doesn\'t really work)'
        else:
            #print 'trying to print line_index = {}, time = {}'.format(line_index, time)
            snap = snapshot(plot_output_dirs2[line_index], time=time)
            plot_snap(snap, save=True)
    plt.xlim([0,1500])
    fig.canvas.mpl_connect('pick_event', onpick2)
    plt.legend()
    fig.savefig('pdf_variance.png', dpi=fig.dpi)
    if show:
        plt.show()
   """ 
def plot_pdf_variance(var, show=True, ax=None):
    if ax == None:
        input = pdf.input
        txt = 'PDF Variance'
        fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
        model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
        title = '{}, {} Model'.format(txt, model)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('time', fontsize=14, fontweight='bold')
        ax.set_ylabel(txt, fontsize=14, fontweight='bold')
    ax.plot(var.times, var.values, label='D = {:.2f}'.format(var.input.noise), picker = True)
            
            
def plot_single_pdf(pdf, show=True):
    var = []
    print 'Plotting!'
    for i, hist in enumerate(pdf.values[:10]):
        # Fourier
        if False:
            test = np.fft.fft2(hist)
            print test
        # Time evolving variance
        if False:
            print 'Plotting' ,i
            input = pdf.input
            txt = 'PDF distribution'
            fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
            ax1 = fig.add_subplot(111)
            model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
            title = '{}, {} Model'.format(txt, model)
            ax1.set_title(title, fontsize=14, fontweight='bold')
            ax1.set_xlabel('particle density', fontsize=14, fontweight='bold')
            ax1.set_ylabel(txt, fontsize=14, fontweight='bold')
            
            # actual calculation
            data = hist.reshape((-1,))  # flatten histogram
            hist1d = np.bincount(data.astype(int))
            edges = np.arange(float(len(hist1d)))
            
            hist1d = np.divide(hist1d, float(hist.size)) # normalize
            edges *= hist.size
            edges /= np.square(input.box_x)
            ax1.plot(edges, hist1d)
            #ax1.plot(edges, hist1d*np.power(edges,2))
            fig.savefig('pdf_noise{:.2f}_{:06d}.png'.format(pdf.input.noise,i), dpi=fig.dpi)
            plt.close()
        if False:
            x = np.linspace(0, pdf.input.box_x, hist.shape[0])
            plt.pcolormesh(x, x, hist)
            plt.colorbar()
            plt.xlim([0,pdf.input.box_x])
            plt.ylim([0,pdf.input.box_x])
            fig.savefig('pdf_noise{:.2f}_{:06d}.png'.format(pdf.input.noise,i), dpi=fig.dpi)
            plt.close()


            

def get_value_box(title='Give Value', prompt='val'):
    root = Tk()
    root.geometry("300x280+300+300")

    w = Label(root, text=title)
    w.pack()
    
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
def plot_snap_torus(snap, show=True, save=False, draw=False, ax=None, arrows=True, color='black', 
        figure=None, backgroundcolor='white', clusters=False, compass_color='r', cluster_color = 'r'):
    input = snap.input
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    density = snap.input.npart / np.square(snap.input.box_x)
    filename = 'snapshot_{}_dens{:.2f}_noise{}_align{}_npart{}_time{}'.format(model, density,input.noise,input.align,input.npart,int(snap.time)).replace(' ','').replace('.','_') + '.png'
    xedge = 0.01

    if ax == None:
        #fig = plt.figure(num=None, figsize=(10, 10), dpi=100, facecolor=backgroundcolor, edgecolor='k')
        fig = plt.figure(figsize=(8, 8), dpi=150)
        ax1 = fig.gca(projection = '3d')
        ax1.set_xlim3d(-1, 1)
        ax1.set_ylim3d(-1, 1)
        ax1.set_zlim3d(-1, 1)
        ax1.set_axis_off()
        ax1.view_init(elev=30., azim=0)
    else:
        ax1 = ax
        fig = figure
        ax.set_axis_off()
    box_x = input.box_x
    x,y,theta,xOb,yOb = snap.x, snap.y, snap.theta, snap.xOb, snap.yOb
    # draw box - vertical lines
    theta = (snap.y / box_x) * 2*np.pi 
    phi = (snap.x / box_x) * 2*np.pi 
    r, R = .40, 1.
    X = (R + r * np.cos(phi)) * np.cos(theta)
    Y = (R + r * np.cos(phi)) * np.sin(theta)
    Z = r * np.sin(phi)
    if (arrows):
        tx = -np.sin(theta)
        ty = np.cos(theta)
        tz = 0
        sx = np.cos(theta)*(-np.sin(phi))
        sy = np.sin(theta)*(-np.sin(phi))
        sz = np.cos(phi)
        b,a = np.cos(snap.theta), np.sin(snap.theta)
        U,V,W = a*tx+b*sx, a*ty+b*sy,a*tz+b*sz
        ax1.quiver(X,Y,Z,U,V,W, length=0.1, color=color)
    else:
        ax1.set_axis_off()
        ax1.scatter(X,Y,Z, color=color, marker="o", s=0.3, alpha=0.8)
        ax1.set_axis_off()
    ax1.axis('off')
    ax1.set_axis_off()
    if save:
        fig.savefig(filename, dpi=fig.dpi, facecolor = fig.get_facecolor(), edgecolor = fig.get_facecolor())
    if draw:
        ax1.set_axis_off()
        plt.draw()
    #plt.close(fig)
    return True
def plot_snap(snap, show=True, save=False, draw=False, ax=None, arrows=True, color='black', 
        figure=None, backgroundcolor='white', clusters=False, compass_color='r', cluster_color = 'r'):
    input = snap.input
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    density = snap.input.npart / np.square(snap.input.box_x)
    filename = 'snapshot_{}_dens{:.2f}_noise{}_align{}_npart{}_time{}'.format(model, density,input.noise,input.align,input.npart,int(snap.time)).replace(' ','').replace('.','_') + '.png'
    xedge = 0.01

    if ax == None:
        fig = plt.figure(num=None, figsize=(10, 10), dpi=100, facecolor=backgroundcolor, edgecolor='k')
        ax1 = fig.add_subplot(111)
        fig.suptitle('{} Model, t = {:4.0f}'.format(model, snap.time), fontsize=14, fontweight='bold', color=color)#, y=0.9)
        #fig.suptitle('L = {}'.format(input.box_x), fontsize=14, fontweight='bold')
        ax1.set_title(r'A={}, D={}, $\rho$={:g}, {} particles, L = {:.0f}'.format(input.align,input.noise, density, input.npart, input.box_x), fontsize=14, color=color)
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
            diam = np.sqrt(cluster_gyration_radius(cluster, x, y, input))
            (x_com,y_com) = center_of_mass(cluster, x, y, input)
            ax1.add_artist(plt.Circle((x_com/box_x,y_com/box_x),diam/box_x,alpha = 0.5,color=cluster_color))
        
    
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
        alpha = 0.1
        if color=='black':
            alpha = 0.8
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
def make_movie(output_dir, start_time=0.0, end_time=None,frames=None, shortTime=False, input_file=None, 
    input=None, filename=None, maxpart=None, nthframe=1, show_compass=True, show_clusters = False, fps=25, black=False,torus=False):
    if filename != None:
        if filename.endswith('.mp4'):
            movie_filename = filename
        else:
            movie_filename = filename + '.mp4'
    if input_file != None:
        input = read_input(input_file)
    else:
        input = read_input_from_dir(output_dir)

    x,y,theta = [],[],[]
    start = start_time
        
    print printtime() + 'Making movie! {}'.format(filename)
    print printtime() + 'Data directory: {0}'.format(output_dir)
    if maxpart != None:
        maxpart = np.min([input.npart, maxpart])
    else:
        maxpart = input.npart

    start_time_jack = datetime.datetime.now()
    ###
    #frames = 10
    if end_time != None:
        ret = read_data_from_dir(output_dir, start_time=start, end_time=end_time, shortTime=True)
        if frames != None:
            total_frames = len(ret.data)
            nthframe = np.max([total_frames // frames, 1])
            print 'nthframe = ', nthframe, ' total_frames = ', total_frames
    elif frames != None:
        ret = read_data_from_dir(output_dir, start_time=start, shortTime=True, frames=frames*nthframe)
        
    times = ret.times
    data = ret.data
    xOb = ret.xOb
    yOb = ret.yOb
    print printtime() +'make_movie: iterations read = {}'.format(len(data))
    print printtime() +'start = {} end = {}'.format(start_time, end_time)
    box_x = input.box_x
    print printtime() + 'box_x = {0} ; npart = {1}'.format(box_x, input.npart)
    now = datetime.datetime.now()
    print printtime() + 'Finished reading data! It took {}:{:02d}'.format((now - start_time_jack).seconds // 60,(now - start_time_jack).seconds % 60)
    if torus:
        fig = plt.figure(figsize=(10, 10), dpi=100)
        #ax = fig.gca(projection = '3d', animated=True)
        ax = fig.add_subplot(1,1,1,projection = '3d')#, animated=True
        ax.set_axis_off()
        ax.set_xlim3d(-1, 1)
        ax.set_ylim3d(-1, 1)
        ax.set_zlim3d(-1, 1)
        ax.view_init(elev=30., azim=0)
        ax.set_axis_off()
    else:
        if black:
            fig = plt.figure(figsize=(10, 10), dpi=100, facecolor='k', edgecolor='k')
        else:
            fig = plt.figure(figsize=(10, 10), dpi=100, facecolor='w')
        ax = fig.add_subplot(1,1,1, animated=True)
        plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.05)
        ax = plt.axes(xlim=(0, box_x), ylim=(0, box_x))

    line, = ax.plot([], [])    
    #ax.set_axis_off()
    def init():
        #line, = ax.plot([], [])
        if torus:
            ax.set_axis_off()
            #rect = fig.patch
            #rect.set_facecolor('w')
            ax.set_axis_off()

        else:
            rect = fig.patch
            if black:
                rect.set_facecolor('k')
            else:
                rect.set_facecolor('w')
        
        return line,

    frame_list = range(len(data))[::nthframe]
    print 'going to write {} frames'.format(len(frame_list))
    data = [frame[:maxpart] for frame in data]
    def animate(i):
        #line, = ax.plot([], [])
        x,y,theta= [],[],[]
        frame = frame_list[i]
        if i % 20 == 0 or torus:
            print printtime() + 'writing frame {0} / {1}'.format(i, len(frame_list))
        
        x = data[frame][:,0]
        y = data[frame][:,1]
        theta = data[frame][:,2]
        # TODO: precalc Pols
        pol = None
        if (show_compass==True):
            pol = np.mean(np.exp(1j * theta))
        snap = snapShot(x,y, theta, xOb, yOb, times[frame], input, pol=pol)
        
        if torus:   
            #fig.clf()
            ax = fig.add_subplot(1,1,1,projection = '3d')
            ax.set_xlim3d(-1, 1)
            ax.set_ylim3d(-1, 1)
            ax.set_zlim3d(-1, 1)
            ax.view_init(elev=30., azim=0)
            ax.set_axis_off()
            plot_snap_torus(snap, show=False, save=False, draw=True, ax=ax, arrows=False,figure=fig)
        elif black:
            ax.cla()
            plot_snap(snap, show=False, save=False, draw=True, ax=ax, arrows=True,figure=fig, color='white', clusters = show_clusters)
        else:
            ax.cla()
            plot_snap(snap, show=False, save=False, draw=True, ax=ax, arrows=True,figure=fig,clusters = show_clusters)
        return line,

    print 'Making movie...'
    ax.set_axis_off()
    anim = animation.FuncAnimation(fig, animate, #init_func=init,
    #                               frames=len(data), interval=20, blit=True)
                                   frames=len(frame_list), interval=1, blit=True)
    
    if not filename == None:
        print printtime() + 'Saving movie...'
        kwargs = {'facecolor':fig.get_facecolor(), 'edgecolor':fig.get_facecolor()}
        anim.save(movie_filename, fps=fps, dpi=200, savefig_kwargs=kwargs)
        now = datetime.datetime.now()
        print 'Finished saving movie! It took {}:{:02d}'.format((now - start_time_jack).seconds // 60,(now - start_time_jack).seconds % 60)
    else:
        print printtime() + 'Showing movie...'
        plt.show()    