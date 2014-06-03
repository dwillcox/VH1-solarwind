from yt.mods import *
import Scientific.IO.NetCDF
import numpy as np
from Scientific.IO.NetCDF import *
from yt.utilities.physical_constants import cm_per_kpc
import h5py
import pylab as P
import os

filePrefix = 'solarwind'
fileSuffix = ['1000','2000','4000','6000','8000','9999']
fileNames = [filePrefix+fS for fS in fileSuffix]

# Lists to hold profiles, labels, and plot specifications.
profiles = []
labels = []

#pf = load(fileNames)
# Loop over each dataset in the time-series.
# Create a data container to hold the whole dataset.
#ad = pf.h.sphere([0.0,0.0,0.0],(2.5e8,"cm"))
# Create a 1d profile of density vs. radius.
#profiles.append(create_profile(ad, ["Radius"],
#                                   fields=["dens"],n=512,
#                                   weight_field=None,
#                                   accumulation=False))
# Add labels
#labels.append("t = %.2f" % pf.current_time)

# Create the profile plot from the list of profiles.
#plot = ProfilePlot.from_profiles(profiles, labels=labels)

# Save the image.
#plot.save()

bbox = np.array([[7.0,50.0],[0.0,1.0*3.14],[0.0,2.0*3.14]])
scale = 100.*cm_per_kpc
field = 'Xvelocity'

fig,pax = P.subplots(2)
for fN in fileNames:
	pfile = NetCDFFile(fN)
	data = {}
	for k in pfile.variables.keys():
		data[k] = pfile.variables[k].getValue()
	pf = load_uniform_grid(data,data[field].shape,scale,bbox=bbox,nprocs=1,periodicity=(False,False,False))
	axcut = 0 # take a line cut along the x axis
	loccut = (0.0,0.0) # place cut ray at origin
	ray = pf.h.ortho_ray(axcut, loccut)
	pax[0].plot(ray['x'], ray['X'],
		label='t = %.2f' % pf.current_time)
	pax[1].semilogy(ray['x'], ray['Xvelocity'],
                label='t = %.2f' % pf.current_time)
pax[0].set_ylabel('Radial Velocity')
pax[0].set_xlabel('r')
#pax[0].set_xlim([0,0.5e8])
#pax[0].set_ylim([1e9,4e9])
#pax[0].legend(loc = 'upper right')
pax[1].set_ylabel('Radial Velocity')
pax[1].set_xlabel('r')
pax[1].set_xlim([0,1e9])
pax[1].legend(loc = 'upper right')
titobj = pax[0].set_title('Radial Velocity Lineout along radius dimension')
P.tight_layout()
titobj.set_y(1.09)
fig.subplots_adjust(top=0.86) 
P.savefig("solarwind_rsweep_ts_Rvel.png")

