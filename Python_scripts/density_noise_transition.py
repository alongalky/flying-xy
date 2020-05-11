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

        
   
def find_transition_noise_by_pol_inflexion(pols):
    if len(pols) < 5:
        return None
    
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    density = getDensityInput(pols[0].input)
    align = pols[0].input.align
    print align
    ax1.set_title(r'A = {:.2f}, $\rho$={:.3f}'.format(align, density), fontsize=14, fontweight='bold')
    ax1.set_xlabel('$D$', fontsize=14, fontweight='bold')
    ax1.set_ylabel(r'$<\Pi>$', fontsize=14, fontweight='bold')
    pols = sorted(pols, key=lambda x: x.input.noise)
    x_dat = np.array([x.input.noise for x in pols])
    y_dat = np.array([x.mValue for x in pols])

    ax1.set_xlim([np.min(x_dat),np.max(x_dat)])
    ax1.set_ylim([np.min(y_dat),np.max(y_dat)])
    line = ax1.plot(x_dat,y_dat, marker='o')
    c = CubicFit(ax1, x_dat, y_dat)
    trans = c.get_trans()
    if trans == None:
        return None
    return transitionPoint(noise=c.get_trans(), density=density, align=align, pols=pols, method=transitionMethods.polInflection,
        const_align=True)
    
def find_transition_noise_by_pol(pols):
    if len(pols) < 5:
        return None
    pols = sorted(pols, key=lambda x: x.input.noise)
    x_dat = np.array([x.input.noise for x in pols])
    y_dat = np.array([x.var for x in pols])
    avg_noise = np.dot(x_dat,y_dat)/np.sum(y_dat)
    ###
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)
    #plt.plot(x_dat, y_dat)
    #plt.vlines(avg_noise, 0, ax1.get_ylim()[1], linestyles='--')
    #plt.show()
    # TODO: possibly use method to select range for averaging
    ###
    density = getDensityInput(pols[0].input)
    align = pols[0].input.align
    return transitionPoint(noise=avg_noise, density=density, align=align, pols=pols, method=transitionMethods.polFluctuation,
        const_align=True)

def updateTransitions(pickleFilename, points):
    if not isfile(pickleFilename):
        pickle_save(pickleFilename, points)
    else:
        old_points = pickle_load(pickleFilename)
        ans = list(old_points) # duplicate
        for new_point in points:
            for old_point in old_points:
                if old_point.isSame(new_point):
                    ans.remove(old_point)
                    ans.append(new_point)
        pickle_save(pickleFilename, points)
                    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='draw transition curve')
    parser.add_argument('directory', metavar='dir', type=str, nargs='?',default=r'.',
                   help='directory with simulation data')                   
    parser.add_argument('-a','--align', metavar='alignment', type=float,# action='store_const', nargs=1,
                   default=1.0, help='which A value to use (default A=1.0)')
    parser.add_argument('-W','--momentum', action='store_true', 
                   help='calculate W (may be slow)')
    parser.add_argument('-r','--read', action='store_true', 
                   help='force rereading sim data and rewriting .npy files')
    parser.add_argument('-c','--calc', action='store_true', 
                   help='force recalculating means')
    args = parser.parse_args()
    resultsDir = args.directory
    print resultsDir
    pickleFilename = 'density_noise_trans.pickle'
    if args.read or not isfile(pickleFilename):
        pols = {}
        Ws = {}
        for output_dir in os.listdir(resultsDir):
            full_dir = join(resultsDir,output_dir)
            print full_dir
            if 'output' in output_dir and isdir(full_dir):
                try:
                    input = read_input_from_dir(full_dir)
                    if input.align != args.align:
                        continue
                except:
                    continue
                pol = calc_pol_var(full_dir, calc=args.calc)
                if pol == None:
                    continue
                dens = getDensityInput(input)
                if dens in pols:
                    pols[dens].append(pol)
                else:
                    pols[dens] = [pol]
                
        pols = {k:v for (k,v) in pols.iteritems() if len(v)>5} # filter only align values with enough data points (> 5)
        densities = sorted(pols)
        print densities
        points = []
        for density in densities[:]:
            points.append(find_transition_noise_by_pol_inflexion(pols[density]))
            #points.append(find_transition_noise_by_pol(pols[density]))
            #points.append(find_transition_density_by_W(pols[noise]))
        points = [point for point in points if point != None]
        updateTransitions(pickleFilename, points)
    else:
        points = pickle_load(pickleFilename)
    txt = 'Transition'
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    input = points[0].pols[0].input
    model = 'Vicsek' if input.neighbor_norm == True else 'Flying XY'
    title = txt#'{}, {}, D={}'.format(txt, model,input.noise)
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlabel('Noise', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Density', fontsize=14, fontweight='bold')
    
    #points1 = [point for point in points if point.method == transitionMethods.polFluctuation]
    #x1 = [point.noise for point in points1]
    #y1 = [point.density for point in points1]
    #ax1.plot(x1,y1,label = transitionMethods.names[transitionMethods.polFluctuation])
    points2 = [point for point in points if point.method == transitionMethods.polInflection]
    x2 = [point.noise for point in points2]
    y2 = [point.density for point in points2]
    ax1.plot(x2,y2,label = transitionMethods.names[transitionMethods.polInflection])
        
    plt.legend()
    plt.show()