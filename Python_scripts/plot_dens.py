import sys
sys.path.insert(0, r'F:\Alon\Code\Scripts\Utility')
import numpy as np
from plots import *
from utils import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
from os.path import isfile,join,isdir
import argparse
def plot_mProperty_vs_noise_groupbynpart(mProp, title='Property Name', log_x = False, log_y = False, 
        show=True, variance=False, ylabel=None, legend=True):
    input = mProp[0].input
    txt = title
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    title = '{}, {}, D={}'.format(txt, model,input.noise)
    #ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel(r'$D$', fontsize=22, fontweight='bold', labelpad=-2)
    ax1.set_ylabel(r'$<(\delta \rho)^2>$', fontsize=22)#, fontweight='bold')

    inset = False
    if inset:
        axins = zoomed_inset_axes(ax1, 7, loc=1, bbox_to_anchor=(0.65, 0.7), 
                     bbox_transform=ax1.figure.transFigure)
        axins.axis([0.09, 0.57, 0, 0.08])
        axins.get_xaxis().set_visible(False)
        axins.get_yaxis().set_visible(False)
        mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        axins.xaxis.tick_top()

    linestyles = ['-', '--','-.',':']
    colors = ['b']*30
    pols_bynpart = {}
    for mPol in mProp:
        npart = mPol.input.npart
        if npart < 500:
            npart = 250
        if not npart in pols_bynpart:
            pols_bynpart[npart] = [mPol]
        else:
            pols_bynpart[npart].append(mPol)
    
    #pols_bynpart = sorted(pols_bynpart)
    ###
    #print pols_bynpart
    #nparts = [500,2000,4000,8000]
    #pols_bynpart = {k:v for (k,v) in pols_bynpart.iteritems() if k in nparts} # filter only align values with enough data points (> 5)

    ###
    mProps = []
    for n,npart in enumerate( sorted(pols_bynpart)):
        pols_bynpart[npart] = sorted(pols_bynpart[npart], key=lambda x: x.input.noise)
        mProp = pols_bynpart[npart]
        mProps.append(mProp)
        x = [en.input.noise for en in mProp]
        y = [en.mValue for en in mProp]
        if inset:
            ax1.plot(x,y, linestyle=linestyles[n%len(linestyles)], linewidth=2.5, label=r'N = {}'.format(npart), 
                color=colors[n], picker = True)
            axins.plot(x,y, linestyle=linestyles[n%len(linestyles)], linewidth=3, label=r'N = {}'.format(npart), #marker='o', markersize = 4,
                color=colors[n],picker = False)
        else:
            ax1.plot(x,y, '--', linewidth=0, label=r'npart = {}'.format(mProp[0].input.npart), marker='o', markersize = 8,
                picker = True)
    def onpick_polar(event):
        ind = event.ind
        time = None
        thisline = event.artist
        line_index = ax1.lines.index(thisline)
        #npart = pols_bynpart.iteritems()[line_index]
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
        """
        print 'debug:'
        output_dir = mProps[line_index][ind].output_dir
        input = read_input_from_dir(output_dir)
        print 'output_dir = ',output_dir,', density = ', getDensityInput(input)
        """
        snap = snapshot(mProps[line_index][ind].output_dir, time=time)
        plot_snap(snap, save=True)

    if (log_x):
        ax1.set_xscale('log')
    if (log_y):
        ax1.set_yscale('log')
    fig.canvas.mpl_connect('pick_event', onpick_polar)
    if legend and len(pols_bynpart)>1:
        ax1.legend(loc='right')
    filename = '{}_vs_dens.png'.format(txt).replace(" ","")
    fig.savefig(filename, dpi=fig.dpi)
    if show:
        plt.show()
    return
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='draw transition curve')
    parser.add_argument('directory', metavar='dir', type=str, nargs='?',default=r'data',
                   help='directory with simulation data')                   
    parser.add_argument('-a','--align', metavar='alignment', type=float,# action='store_const', nargs=1,
                   default=1.0, help='which A value to use (default A=1.0)')
    parser.add_argument('-n','--noise', metavar='noise value', type=float,# action='store_const', nargs=1,
                   default=2.0, help='which noise to plot')
    parser.add_argument('-r','--read', action='store_true', 
                   help='force rereading sim data and rewriting .npy files')
    parser.add_argument('-c','--calc', action='store_true', 
                   help='force recalculating means')
    args = parser.parse_args()
    resultsDir = args.directory
    #print resultsDir
    pickleFilename = 'mDensVars.pickle'
    if args.read or not isfile(pickleFilename):
        mDensVars = []
        output_dirs = [join(resultsDir, x) for x in os.listdir(resultsDir)]
        print output_dirs
        for full_dir in output_dirs:
            if 'output' in full_dir and isdir(full_dir):
                try:
                    input = read_input_from_dir(full_dir)
                    #if input.align != args.align or input.noise != args.noise:
                    #    continue
                    #if input.npart != 8000:
                    #    continue
                    print full_dir
                    pol = calc_mProp_orderparam(full_dir, orderParamter('Density variance', 'densvar', calc_vars, mean_pickle_name='mDensVar', ylabel=r'$<(\delta \rho)^2>$'),
                                read=False, calc=args.calc)
                    mDensVars.append(pol)
                except:
                    continue
        pickle_save(pickleFilename, mDensVars)
    else:
        mDensVars = pickle_load(pickleFilename)
    plot_mProperty_vs_noise_groupbynpart(mDensVars, title='Density Variations', log_x = False, log_y = False, 
        show=True, variance=False, ylabel=None)
