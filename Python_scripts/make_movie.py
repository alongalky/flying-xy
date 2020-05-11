import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
sys.path.insert(0, r'F:\Alon\Code\Scripts\Utility')
import struct
import os
from utils import *
from plots import *
import argparse
import datetime

parser = argparse.ArgumentParser(description='Make an .mp4 movie from a simulation output directory')
parser.add_argument('directory', metavar='dir', type=str, nargs=1,
                   help='directory with simulation data')                   
parser.add_argument('-s','--save', metavar='out', #type=str, action='store_const', nargs=1,
                   default=None, help='movie file name (will not output if none given')
parser.add_argument('-i','--input', metavar='input', #type=str, action='store_const', nargs=1,
                   default=None, help='input file name (optional)')
parser.add_argument('-p','--play', action='store_true', 
                   default=None, help='play movie')
parser.add_argument('-T','--torus', action='store_true', 
                   default=None, help='bagel!')
parser.add_argument('-t','--time', metavar='start_time', type=float,# action='store_const', nargs=1,
               default=0.0, help='time to start')
parser.add_argument('-f','--frames', metavar='length_inframes', type=int,# action='store_const', nargs=1,
               default=100, help='Movie length in frames')
parser.add_argument('-e','--endtime', metavar='end_time', type=int,# action='store_const', nargs=1,
               default=None, help='time to end')
parser.add_argument('-n','--npart', metavar='particle number', type=int,# action='store_const', nargs=1,
               default=None, help='maximum number of particles on frame')
parser.add_argument('-v','--nthframe', metavar='every nth frame / movie speed', type=int,# action='store_const', nargs=1,
               default=1, help='skip every nth frame to get faseter movie')
args = parser.parse_args()
output_dir      = args.directory[0]
play            = args.play
nthframe        = args.nthframe
movie_filename = args.save
if not args.save == None:
    name = args.save
    if name.endswith('.mp4'):
        movie_filename = name
    else:
        movie_filename = name + '.mp4'

make_movie(output_dir, start_time=args.time, frames=args.frames, filename=args.save, torus=args.torus)
sys.exit()        
######################
if not args.input == None:
    input_file      = args.input
else:
    # look for an input file in output_dir
    print 'debug: {0}'.format(output_dir)
    for f in os.listdir(output_dir):
        if 'input' in f:
            input_file = output_dir + '\\' + f
            break
    if input_file == None:
        print printtime() + 'Error: No input file found!'
        sys.exit()

x,y,theta = [],[],[]
start = args.time
    
print printtime() + 'Making movie!'
print printtime() + 'Data directory: {0}'.format(output_dir)
print printtime() + 'input_file: {0}'.format(input_file)
input = read_input(input_file)
#def read_data_from_dir(output_dir, start=0, frames=None, max_mem_size_mb=None, start_time=None, shortTime=False, **kwargs):
if args.npart != None:
    maxpart = np.min([args.npart, input.npart])
else:
    maxpart = input.npart

start_time = datetime.datetime.now()
ret = read_data_from_dir(output_dir, start_time=start, end_time=args.endtime, shortTime=True, frames=args.frames*nthframe)
times = ret.times
data = ret.data
xOb = ret.xOb
yOb = ret.yOb
print 'make_movie: iterations read = {}'.format(len(data))
box_x = input.box_x
print printtime() + 'box_x = {0} ; npart = {1}'.format(box_x, input.npart)
now = datetime.datetime.now()
print 'Finished reading data! It took {}:{:02d}'.format((now - start_time).seconds // 60,(now - start_time).seconds % 60)
fig = plt.figure(figsize=(10, 10), dpi=100, facecolor='k', edgecolor='k')
ax = fig.add_subplot(1,1,1, animated=True, axisbg='red')
plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.05)
ax = plt.axes(xlim=(0, box_x), ylim=(0, box_x))

def init():
    rect = fig.patch
    rect.set_facecolor('k')
    line, = ax.plot([], [])
    return line,

frame_list = range(len(data))[::nthframe]
print 'going to write {} frames'.format(len(frame_list))
print np.shape(data)
print type(data[:])
data = [frame[:maxpart] for frame in data]
print np.shape(data)
def animate(i):
    x,y,theta= [],[],[]
    frame = frame_list[i]
    if i % 20 == 0:
        print printtime() + 'writing frame {0} / {1}'.format(i, len(frame_list))
    
    x = data[frame][:,0]
    y = data[frame][:,1]
    theta = data[frame][:,2]
    snap = snapShot(x,y, theta, xOb, yOb, times[frame], input)
    ax.cla()
    plot_snap(snap, show=False, save=False, draw=True, ax=ax, arrows=True,figure=fig, color='white')
    line, = ax.plot([], [])
    return line,

print 'Making movie...'
anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=len(data), interval=20, blit=True)
                               frames=len(frame_list), interval=1, blit=True)
if not movie_filename == None:
    print printtime() + 'Saving movie...'
    kwargs = {'facecolor':fig.get_facecolor(), 'edgecolor':fig.get_facecolor()}
    anim.save(movie_filename, fps=25, dpi=200, savefig_kwargs=kwargs)
    now = datetime.datetime.now()
    print 'Finished saving movie! It took {}:{:02d}'.format((now - start_time).seconds // 60,(now - start_time).seconds % 60)

if play :
    print printtime() + 'Showing movie...'
    plt.show()

