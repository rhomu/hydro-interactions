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

outdir = sys.argv[1]
dat = open(outdir + '/parameters').read()

# load params of the form 'name = value'
def load(name, dat):
    m = re.search(r'\s+{}\s+=\s+(.*)'.format(name), dat)
    return m.group(1)

npart  = int(load('npart', dat))
nsteps = int(load('nsteps', dat))
ninfo  = int(load('ninfo', dat))
nbins  = int(load('nbins', dat))
hmax   = float(load('hmax', dat))

# ------------------------------------------------------------------------------
# plot

dat = open(outdir+'/histogram.dat', mode='r').read()
arr = np.array(struct.unpack('d'*(len(dat)//8), dat))

fig = plt.figure(figsize=(12,6))
x = np.linspace(0, hmax*(1+1./nbins), num = nbins+1) + hmax/nbins/2.
plt.plot(x, arr)
plt.show(); exit(0)

##writer = ani.writers['ffmpeg'](fps=10, bitrate=1800)
#an.save('movie.mp4', writer=writer)
