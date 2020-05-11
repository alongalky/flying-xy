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

        
   

def find_transition_density_by_pol_inflexion(pols):
    if len(pols) < 5:
        return None
    
    fig = plt.figure(num=None, figsize=(8, 5), dpi=150, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(111)
    noise = pols[0].input.noise
    align = pols[0].input.align
    print align
    ax1.set_title(r'A = {:.2f}, D={:.2f}'.format(align, noise), fontsize=22, fontweight='bold')
    ax1.set_xlabel(r'$\rho$', fontsize=22, fontweight='bold', labelpad=-2)
    ax1.set_ylabel(r'$<\Pi>$', fontsize=22, fontweight='bold')
    pols = sorted(pols, key=lambda x: getDensityInput(x.input))
    x_dat = np.array([getDensityInput(x.input) for x in pols])
    y_dat = np.array([x.mValue for x in pols])

    ax1.set_xlim([np.min(x_dat),np.max(x_dat)])
    ax1.set_ylim([np.min(y_dat),np.max(y_dat)])
    line = ax1.plot(x_dat,y_dat, marker='o')
    c = CubicFit(ax1, x_dat, y_dat)
    trans = c.get_trans()
    if trans == None:
        return None
    return transitionPoint(noise=noise, density=c.get_trans(), align=align, pols=pols, method=transitionMethods.polInflection,
        const_noise=True)
    
def find_transition_density_by_pol(pols):
    if len(pols) < 5:
        return None
    pols = sorted(pols, key=lambda x: getDensityInput(x.input))
    x_dat = np.array([getDensityInput(x.input) for x in pols])
    y_dat = np.array([x.var for x in pols])
    avg_dens = np.dot(x_dat,y_dat)/np.sum(y_dat)
    
    noise = pols[0].input.noise
    align = pols[0].input.align
    return transitionPoint(noise=noise, density=avg_dens, align=align, pols=pols, method=transitionMethods.polFluctuation,
        const_noise=True)

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
    parser.add_argument('-n','--noise', metavar='noise value', type=float,# action='store_const', nargs=1,
                   default=2.0, help='which noise to plot')
    parser.add_argument('-W','--momentum', action='store_true', 
                   help='calculate W (may be slow)')
    parser.add_argument('-r','--read', action='store_true', 
                   help='force rereading sim data and rewriting .npy files')
    parser.add_argument('-c','--calc', action='store_true', 
                   help='force recalculating means')
    args = parser.parse_args()
    resultsDir = args.directory
    print resultsDir
    pickleFilename = 'align_density_trans.pickle'
    if args.read or not isfile(pickleFilename):
        pols = {}
        Ws = {}
        for output_dir in os.listdir(resultsDir):
            full_dir = join(resultsDir,output_dir)
            print full_dir
            if 'output' in output_dir and isdir(full_dir):
                try:
                    input = read_input_from_dir(full_dir)
                    if input.noise != args.noise:
                        continue
                except:
                    print 'horrible error in ', full_dir
                    continue
                pol = calc_pol_var(full_dir, calc=args.calc)
                if pol == None:
                    continue
                align = input.align
                if align in pols:
                    pols[align].append(pol)
                else:
                    pols[align] = [pol]
                
        pols = {k:v for (k,v) in pols.iteritems() if len(v)>5} # filter only align values with enough data points (> 5)
        aligns = sorted(pols)
        print aligns
        points = []
        for align in aligns[:]:
            points.append(find_transition_density_by_pol_inflexion(pols[align]))
            #points.append(find_transition_density_by_pol(pols[align]))
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
    ax1.set_xlabel('Align', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Density', fontsize=14, fontweight='bold')
    
    points1 = [point for point in points if point.method == transitionMethods.polFluctuation]
    #x1 = [point.align for point in points1]
    #y1 = [point.density for point in points1]
    points2 = [point for point in points if point.method == transitionMethods.polInflection]
    x2 = [point.align for point in points2]
    y2 = [point.density for point in points2]
    #ax1.plot(x1,y1,label = transitionMethods.names[transitionMethods.polFluctuation])
    ax1.plot(x2,y2,label = transitionMethods.names[transitionMethods.polInflection])
        
    plt.legend()
    plt.show()
