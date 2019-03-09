#
# simple plotting
#

import re
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# ------------------------------------------------------------------------------

# check args and open param file
if len(sys.argv)<2:
    print "Error: please provide a output directory."
    exit(1)

def get_data(frame):
    fn = '{0}/frame{1}.dat'.format(outdir ,str(frame*ninfo))
    dat = open(fn, mode='r').read()
    arr = np.array(struct.unpack('d'*(len(dat)//8), dat))
    return arr.reshape((-1, 3))

outdir = sys.argv[1]
dat = open(outdir + '/parameters').read()

# load params of the form 'name = value'
def load(name, dat):
    m = re.search(r'\s+{}\s+=\s+(.*)'.format(name), dat)
    return m.group(1)

npart  = int(load('npart', dat))
nsteps = int(load('nsteps', dat))
ninfo  = int(load('ninfo', dat))

# ------------------------------------------------------------------------------
# plot

def plot_frame(frame):
    fig.clf()
    particles = get_data(frame)
    for i in range(npart):
        print particles[i]
        plt.plot(particles[i,0], particles[i,1], 'or')
        plt.axes().set_aspect('equal', adjustable='box')
        plt.xlim([-5, 5])
        plt.ylim([-5, 5])

fig = plt.figure(figsize=(12,6))
an = ani.FuncAnimation(fig, plot_frame,
                       frames = np.arange(0, nsteps//ninfo),
                       interval = 100, blit = False)
plt.show(); exit(0)

writer = ani.writers['ffmpeg'](fps=10, bitrate=1800)
an.save('movie.mp4', writer=writer)
