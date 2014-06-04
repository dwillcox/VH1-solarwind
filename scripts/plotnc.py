from yt.mods import *
import Scientific.IO.NetCDF
from numpy import *
from Scientific.IO.NetCDF import *
from yt.utilities.physical_constants import cm_per_kpc
import h5py

file = NetCDFFile('sod-mpi1032.nc','r')
data = {}
dkeys = file.variables.keys()
for k in dkeys:
	data[k] = file.variables[k].getValue()

# To plot 2D or 1D data, set the size of the unused dims to 1 in dshape
# Otherwise, treat it like 3D data

bbox = array([[0.0,1.0],[0.0,1.0],[0.0,0.0]])
dshape2d = data['Density'].shape
dshape = (dshape2d[0],dshape2d[1],1) # 2D geometry
scalemultiplier = 1.0

pf = load_uniform_grid(data,dshape,scalemultiplier,bbox=bbox,nprocs=1,periodicity=(False,False,False))

for k in dkeys:
	slc = SlicePlot(pf,"z",[k])
	slc.save()

